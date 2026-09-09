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


#include "SimdThreeCenterElectronRepulsionVrrRecSOI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;

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

    const auto *sni0_0 = buffer.data(sni0 + 0);
    const auto *sni0_3 = buffer.data(sni0 + 3);
    const auto *sni0_5 = buffer.data(sni0 + 5);
    const auto *sni0_6 = buffer.data(sni0 + 6);
    const auto *sni0_9 = buffer.data(sni0 + 9);
    const auto *sni0_10 = buffer.data(sni0 + 10);
    const auto *sni0_12 = buffer.data(sni0 + 12);
    const auto *sni0_14 = buffer.data(sni0 + 14);
    const auto *sni0_21 = buffer.data(sni0 + 21);
    const auto *sni0_27 = buffer.data(sni0 + 27);
    const auto *sni0_31 = buffer.data(sni0 + 31);
    const auto *sni0_34 = buffer.data(sni0 + 34);
    const auto *sni0_38 = buffer.data(sni0 + 38);
    const auto *sni0_56 = buffer.data(sni0 + 56);
    const auto *sni0_61 = buffer.data(sni0 + 61);
    const auto *sni0_65 = buffer.data(sni0 + 65);

    const auto *snh_0 = buffer.data(snh + 0);
    const auto *snh_1 = buffer.data(snh + 1);
    const auto *snh_2 = buffer.data(snh + 2);
    const auto *snh_3 = buffer.data(snh + 3);
    const auto *snh_5 = buffer.data(snh + 5);
    const auto *snh_6 = buffer.data(snh + 6);
    const auto *snh_7 = buffer.data(snh + 7);
    const auto *snh_8 = buffer.data(snh + 8);
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
    const auto *snh_47 = buffer.data(snh + 47);
    const auto *snh_57 = buffer.data(snh + 57);
    const auto *snh_58 = buffer.data(snh + 58);
    const auto *snh_59 = buffer.data(snh + 59);
    const auto *snh_60 = buffer.data(snh + 60);
    const auto *snh_61 = buffer.data(snh + 61);
    const auto *snh_62 = buffer.data(snh + 62);
    const auto *snh_63 = buffer.data(snh + 63);
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

    const auto *sni1_0 = buffer.data(sni1 + 0);
    const auto *sni1_3 = buffer.data(sni1 + 3);
    const auto *sni1_5 = buffer.data(sni1 + 5);
    const auto *sni1_6 = buffer.data(sni1 + 6);
    const auto *sni1_9 = buffer.data(sni1 + 9);
    const auto *sni1_10 = buffer.data(sni1 + 10);
    const auto *sni1_12 = buffer.data(sni1 + 12);
    const auto *sni1_14 = buffer.data(sni1 + 14);
    const auto *sni1_21 = buffer.data(sni1 + 21);
    const auto *sni1_27 = buffer.data(sni1 + 27);
    const auto *sni1_31 = buffer.data(sni1 + 31);
    const auto *sni1_34 = buffer.data(sni1 + 34);
    const auto *sni1_38 = buffer.data(sni1 + 38);
    const auto *sni1_56 = buffer.data(sni1 + 56);
    const auto *sni1_61 = buffer.data(sni1 + 61);
    const auto *sni1_65 = buffer.data(sni1 + 65);

    const auto *sog0_0 = buffer.data(sog0 + 0);
    const auto *sog0_3 = buffer.data(sog0 + 3);
    const auto *sog0_5 = buffer.data(sog0 + 5);
    const auto *sog0_6 = buffer.data(sog0 + 6);
    const auto *sog0_9 = buffer.data(sog0 + 9);
    const auto *sog0_10 = buffer.data(sog0 + 10);
    const auto *sog0_12 = buffer.data(sog0 + 12);
    const auto *sog0_13 = buffer.data(sog0 + 13);
    const auto *sog0_14 = buffer.data(sog0 + 14);
    const auto *sog0_25 = buffer.data(sog0 + 25);
    const auto *sog0_27 = buffer.data(sog0 + 27);
    const auto *sog0_28 = buffer.data(sog0 + 28);
    const auto *sog0_29 = buffer.data(sog0 + 29);
    const auto *sog0_42 = buffer.data(sog0 + 42);
    const auto *sog0_43 = buffer.data(sog0 + 43);
    const auto *sog0_44 = buffer.data(sog0 + 44);
    const auto *sog0_45 = buffer.data(sog0 + 45);
    const auto *sog0_48 = buffer.data(sog0 + 48);
    const auto *sog0_50 = buffer.data(sog0 + 50);
    const auto *sog0_51 = buffer.data(sog0 + 51);
    const auto *sog0_54 = buffer.data(sog0 + 54);
    const auto *sog0_55 = buffer.data(sog0 + 55);
    const auto *sog0_57 = buffer.data(sog0 + 57);
    const auto *sog0_58 = buffer.data(sog0 + 58);
    const auto *sog0_59 = buffer.data(sog0 + 59);

    const auto *sog1_0 = buffer.data(sog1 + 0);
    const auto *sog1_3 = buffer.data(sog1 + 3);
    const auto *sog1_5 = buffer.data(sog1 + 5);
    const auto *sog1_6 = buffer.data(sog1 + 6);
    const auto *sog1_9 = buffer.data(sog1 + 9);
    const auto *sog1_10 = buffer.data(sog1 + 10);
    const auto *sog1_12 = buffer.data(sog1 + 12);
    const auto *sog1_13 = buffer.data(sog1 + 13);
    const auto *sog1_14 = buffer.data(sog1 + 14);
    const auto *sog1_25 = buffer.data(sog1 + 25);
    const auto *sog1_27 = buffer.data(sog1 + 27);
    const auto *sog1_28 = buffer.data(sog1 + 28);
    const auto *sog1_29 = buffer.data(sog1 + 29);
    const auto *sog1_42 = buffer.data(sog1 + 42);
    const auto *sog1_43 = buffer.data(sog1 + 43);
    const auto *sog1_44 = buffer.data(sog1 + 44);
    const auto *sog1_45 = buffer.data(sog1 + 45);
    const auto *sog1_48 = buffer.data(sog1 + 48);
    const auto *sog1_50 = buffer.data(sog1 + 50);
    const auto *sog1_51 = buffer.data(sog1 + 51);
    const auto *sog1_54 = buffer.data(sog1 + 54);
    const auto *sog1_55 = buffer.data(sog1 + 55);
    const auto *sog1_57 = buffer.data(sog1 + 57);
    const auto *sog1_58 = buffer.data(sog1 + 58);
    const auto *sog1_59 = buffer.data(sog1 + 59);

    const auto *soh_0 = buffer.data(soh + 0);
    const auto *soh_2 = buffer.data(soh + 2);
    const auto *soh_3 = buffer.data(soh + 3);
    const auto *soh_5 = buffer.data(soh + 5);
    const auto *soh_6 = buffer.data(soh + 6);
    const auto *soh_9 = buffer.data(soh + 9);
    const auto *soh_10 = buffer.data(soh + 10);
    const auto *soh_12 = buffer.data(soh + 12);
    const auto *soh_14 = buffer.data(soh + 14);
    const auto *soh_15 = buffer.data(soh + 15);
    const auto *soh_16 = buffer.data(soh + 16);
    const auto *soh_17 = buffer.data(soh + 17);
    const auto *soh_18 = buffer.data(soh + 18);
    const auto *soh_19 = buffer.data(soh + 19);
    const auto *soh_20 = buffer.data(soh + 20);
    const auto *soh_21 = buffer.data(soh + 21);
    const auto *soh_23 = buffer.data(soh + 23);
    const auto *soh_24 = buffer.data(soh + 24);
    const auto *soh_26 = buffer.data(soh + 26);
    const auto *soh_27 = buffer.data(soh + 27);
    const auto *soh_30 = buffer.data(soh + 30);
    const auto *soh_36 = buffer.data(soh + 36);
    const auto *soh_37 = buffer.data(soh + 37);
    const auto *soh_38 = buffer.data(soh + 38);
    const auto *soh_39 = buffer.data(soh + 39);
    const auto *soh_40 = buffer.data(soh + 40);
    const auto *soh_41 = buffer.data(soh + 41);
    const auto *soh_42 = buffer.data(soh + 42);
    const auto *soh_44 = buffer.data(soh + 44);
    const auto *soh_45 = buffer.data(soh + 45);
    const auto *soh_47 = buffer.data(soh + 47);
    const auto *soh_48 = buffer.data(soh + 48);
    const auto *soh_51 = buffer.data(soh + 51);
    const auto *soh_57 = buffer.data(soh + 57);
    const auto *soh_58 = buffer.data(soh + 58);
    const auto *soh_59 = buffer.data(soh + 59);
    const auto *soh_60 = buffer.data(soh + 60);
    const auto *soh_61 = buffer.data(soh + 61);
    const auto *soh_62 = buffer.data(soh + 62);
    const auto *soh_63 = buffer.data(soh + 63);
    const auto *soh_65 = buffer.data(soh + 65);
    const auto *soh_66 = buffer.data(soh + 66);
    const auto *soh_68 = buffer.data(soh + 68);
    const auto *soh_69 = buffer.data(soh + 69);
    const auto *soh_72 = buffer.data(soh + 72);
    const auto *soh_73 = buffer.data(soh + 73);
    const auto *soh_75 = buffer.data(soh + 75);
    const auto *soh_77 = buffer.data(soh + 77);
    const auto *soh_78 = buffer.data(soh + 78);
    const auto *soh_79 = buffer.data(soh + 79);
    const auto *soh_80 = buffer.data(soh + 80);
    const auto *soh_81 = buffer.data(soh + 81);
    const auto *soh_82 = buffer.data(soh + 82);
    const auto *soh_83 = buffer.data(soh + 83);
    const auto *soh_84 = buffer.data(soh + 84);
    const auto *soh_86 = buffer.data(soh + 86);
    const auto *soh_87 = buffer.data(soh + 87);
    const auto *soh_89 = buffer.data(soh + 89);
    const auto *soh_90 = buffer.data(soh + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, snh_0, snh_3, sog0_0, sog0_3, \
                         sog1_0, sog1_3, soh_0, soh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * snh_0[k]
                 + f_1 * sog0_0[k]
                 - f_2 * sog1_0[k]
                 + f_3 * pc_x[k] * soh_0[k];

        t_1[k] = f_3 * pc_y[k] * soh_0[k];

        t_2[k] = f_3 * pc_z[k] * soh_0[k];

        t_3[k] = f_0 * snh_3[k]
                 + f_4 * sog0_3[k]
                 - f_5 * sog1_3[k]
                 + f_3 * pc_x[k] * soh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, snh_5, snh_6, sog0_5, sog0_6, sog1_5, \
                         sog1_6, soh_2, soh_5, soh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * soh_2[k];

        t_5[k] = f_0 * snh_5[k]
                 + f_4 * sog0_5[k]
                 - f_5 * sog1_5[k]
                 + f_3 * pc_x[k] * soh_5[k];

        t_6[k] = f_0 * snh_6[k]
                 + f_6 * sog0_6[k]
                 - f_7 * sog1_6[k]
                 + f_3 * pc_x[k] * soh_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, snh_9, sog0_9, sog1_9, soh_3, soh_5, \
                         soh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * soh_3[k];

        t_8[k] = f_3 * pc_y[k] * soh_5[k];

        t_9[k] = f_0 * snh_9[k]
                 + f_6 * sog0_9[k]
                 - f_7 * sog1_9[k]
                 + f_3 * pc_x[k] * soh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, snh_10, snh_12, sog0_10, sog0_12, \
                         sog1_10, sog1_12, soh_6, soh_10, soh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * snh_10[k]
                  + f_8 * sog0_10[k]
                  - f_9 * sog1_10[k]
                  + f_3 * pc_x[k] * soh_10[k];

        t_11[k] = f_3 * pc_z[k] * soh_6[k];

        t_12[k] = f_0 * snh_12[k]
                  + f_8 * sog0_12[k]
                  - f_9 * sog1_12[k]
                  + f_3 * pc_x[k] * soh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, snh_14, snh_15, snh_16, sog0_14, \
                         sog1_14, soh_9, soh_14, soh_15, soh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * soh_9[k];

        t_14[k] = f_0 * snh_14[k]
                  + f_8 * sog0_14[k]
                  - f_9 * sog1_14[k]
                  + f_3 * pc_x[k] * soh_14[k];

        t_15[k] = f_0 * snh_15[k]
                  + f_3 * pc_x[k] * soh_15[k];

        t_16[k] = f_0 * snh_16[k]
                  + f_3 * pc_x[k] * soh_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, snh_17, snh_18, snh_19, snh_20, soh_17, \
                         soh_18, soh_19, soh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * snh_17[k]
                  + f_3 * pc_x[k] * soh_17[k];

        t_18[k] = f_0 * snh_18[k]
                  + f_3 * pc_x[k] * soh_18[k];

        t_19[k] = f_0 * snh_19[k]
                  + f_3 * pc_x[k] * soh_19[k];

        t_20[k] = f_0 * snh_20[k]
                  + f_3 * pc_x[k] * soh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, sog0_10, sog0_12, sog0_13, \
                         sog1_10, sog1_12, sog1_13, soh_15, soh_17, \
                         soh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * sog0_10[k]
                  - f_2 * sog1_10[k]
                  + f_3 * pc_y[k] * soh_15[k];

        t_22[k] = f_3 * pc_z[k] * soh_15[k];

        t_23[k] = f_4 * sog0_12[k]
                  - f_5 * sog1_12[k]
                  + f_3 * pc_y[k] * soh_17[k];

        t_24[k] = f_6 * sog0_13[k]
                  - f_7 * sog1_13[k]
                  + f_3 * pc_y[k] * soh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_y, pc_y, pc_z, sni0_0, snh_0, \
                         sni1_0, sog0_14, sog1_14, soh_19, soh_20, \
                         soh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_8 * sog0_14[k]
                  - f_9 * sog1_14[k]
                  + f_3 * pc_y[k] * soh_19[k];

        t_26[k] = f_3 * pc_y[k] * soh_20[k];

        t_27[k] = f_1 * sog0_14[k]
                  - f_2 * sog1_14[k]
                  + f_3 * pc_z[k] * soh_20[k];

        t_28[k] = pb_y[k] * sni0_0[k]
                  - f_10 * pc_y[k] * sni1_0[k];

        t_29[k] = f_11 * snh_0[k]
                  + f_3 * pc_y[k] * soh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_y, pc_z, sni0_3, sni0_5, snh_1, \
                         snh_2, sni1_3, sni1_5, soh_21, soh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * soh_21[k];

        t_31[k] = pb_y[k] * sni0_3[k]
                  + f_12 * snh_1[k]
                  - f_10 * pc_y[k] * sni1_3[k];

        t_32[k] = f_11 * snh_2[k]
                  + f_3 * pc_y[k] * soh_23[k];

        t_33[k] = pb_y[k] * sni0_5[k]
                  - f_10 * pc_y[k] * sni1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_y, pc_y, pc_z, sni0_6, sni0_9, snh_3, \
                         snh_5, sni1_6, sni1_9, soh_24, soh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * sni0_6[k]
                  + f_13 * snh_3[k]
                  - f_10 * pc_y[k] * sni1_6[k];

        t_35[k] = f_3 * pc_z[k] * soh_24[k];

        t_36[k] = f_11 * snh_5[k]
                  + f_3 * pc_y[k] * soh_26[k];

        t_37[k] = pb_y[k] * sni0_9[k]
                  - f_10 * pc_y[k] * sni1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sni0_10, sni0_12, snh_6, \
                         snh_8, snh_9, sni1_10, sni1_12, soh_27, \
                         soh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * sni0_10[k]
                  + f_14 * snh_6[k]
                  - f_10 * pc_y[k] * sni1_10[k];

        t_39[k] = f_3 * pc_z[k] * soh_27[k];

        t_40[k] = pb_y[k] * sni0_12[k]
                  + f_12 * snh_8[k]
                  - f_10 * pc_y[k] * sni1_12[k];

        t_41[k] = f_11 * snh_9[k]
                  + f_3 * pc_y[k] * soh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, sni0_14, snh_36, snh_37, \
                         snh_38, sni1_14, soh_36, soh_37, soh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * sni0_14[k]
                  - f_10 * pc_y[k] * sni1_14[k];

        t_43[k] = f_15 * snh_36[k]
                  + f_3 * pc_x[k] * soh_36[k];

        t_44[k] = f_15 * snh_37[k]
                  + f_3 * pc_x[k] * soh_37[k];

        t_45[k] = f_15 * snh_38[k]
                  + f_3 * pc_x[k] * soh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, snh_15, snh_39, snh_40, snh_41, \
                         sog0_25, sog1_25, soh_36, soh_39, soh_40, \
                         soh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * snh_39[k]
                  + f_3 * pc_x[k] * soh_39[k];

        t_47[k] = f_15 * snh_40[k]
                  + f_3 * pc_x[k] * soh_40[k];

        t_48[k] = f_15 * snh_41[k]
                  + f_3 * pc_x[k] * soh_41[k];

        t_49[k] = f_11 * snh_15[k]
                  + f_1 * sog0_25[k]
                  - f_2 * sog1_25[k]
                  + f_3 * pc_y[k] * soh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pc_y, pc_z, snh_17, snh_18, sog0_27, sog0_28, \
                         sog1_27, sog1_28, soh_36, soh_38, soh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * soh_36[k];

        t_51[k] = f_11 * snh_17[k]
                  + f_4 * sog0_27[k]
                  - f_5 * sog1_27[k]
                  + f_3 * pc_y[k] * soh_38[k];

        t_52[k] = f_11 * snh_18[k]
                  + f_6 * sog0_28[k]
                  - f_7 * sog1_28[k]
                  + f_3 * pc_y[k] * soh_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_y, sni0_27, snh_19, snh_20, sni1_27, \
                         sog0_29, sog1_29, soh_40, soh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * snh_19[k]
                  + f_8 * sog0_29[k]
                  - f_9 * sog1_29[k]
                  + f_3 * pc_y[k] * soh_40[k];

        t_54[k] = f_11 * snh_20[k]
                  + f_3 * pc_y[k] * soh_41[k];

        t_55[k] = pb_y[k] * sni0_27[k]
                  - f_10 * pc_y[k] * sni1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, sni0_0, sni0_3, \
                         snh_0, sni1_0, sni1_3, soh_42, soh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * sni0_0[k]
                  - f_10 * pc_z[k] * sni1_0[k];

        t_57[k] = f_3 * pc_y[k] * soh_42[k];

        t_58[k] = f_11 * snh_0[k]
                  + f_3 * pc_z[k] * soh_42[k];

        t_59[k] = pb_z[k] * sni0_3[k]
                  - f_10 * pc_z[k] * sni1_3[k];

        t_60[k] = f_3 * pc_y[k] * soh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sni0_5, sni0_6, snh_2, \
                         snh_3, sni1_5, sni1_6, soh_45, soh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_z[k] * sni0_5[k]
                  + f_12 * snh_2[k]
                  - f_10 * pc_z[k] * sni1_5[k];

        t_62[k] = pb_z[k] * sni0_6[k]
                  - f_10 * pc_z[k] * sni1_6[k];

        t_63[k] = f_11 * snh_3[k]
                  + f_3 * pc_z[k] * soh_45[k];

        t_64[k] = f_3 * pc_y[k] * soh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_z, pc_z, sni0_9, sni0_10, sni0_12, snh_5, \
                         snh_6, snh_7, sni1_9, sni1_10, sni1_12, \
                         soh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_z[k] * sni0_9[k]
                  + f_13 * snh_5[k]
                  - f_10 * pc_z[k] * sni1_9[k];

        t_66[k] = pb_z[k] * sni0_10[k]
                  - f_10 * pc_z[k] * sni1_10[k];

        t_67[k] = f_11 * snh_6[k]
                  + f_3 * pc_z[k] * soh_48[k];

        t_68[k] = pb_z[k] * sni0_12[k]
                  + f_12 * snh_7[k]
                  - f_10 * pc_z[k] * sni1_12[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_z, pc_x, pc_y, pc_z, sni0_14, snh_9, \
                         snh_57, snh_58, sni1_14, soh_51, soh_57, \
                         soh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_3 * pc_y[k] * soh_51[k];

        t_70[k] = pb_z[k] * sni0_14[k]
                  + f_14 * snh_9[k]
                  - f_10 * pc_z[k] * sni1_14[k];

        t_71[k] = f_15 * snh_57[k]
                  + f_3 * pc_x[k] * soh_57[k];

        t_72[k] = f_15 * snh_58[k]
                  + f_3 * pc_x[k] * soh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, snh_59, snh_60, snh_61, snh_62, soh_59, \
                         soh_60, soh_61, soh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * snh_59[k]
                  + f_3 * pc_x[k] * soh_59[k];

        t_74[k] = f_15 * snh_60[k]
                  + f_3 * pc_x[k] * soh_60[k];

        t_75[k] = f_15 * snh_61[k]
                  + f_3 * pc_x[k] * soh_61[k];

        t_76[k] = f_15 * snh_62[k]
                  + f_3 * pc_x[k] * soh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, sni0_21, snh_15, sni1_21, \
                         sog0_42, sog1_42, soh_57, soh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * sni0_21[k]
                  - f_10 * pc_z[k] * sni1_21[k];

        t_78[k] = f_11 * snh_15[k]
                  + f_3 * pc_z[k] * soh_57[k];

        t_79[k] = f_4 * sog0_42[k]
                  - f_5 * sog1_42[k]
                  + f_3 * pc_y[k] * soh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, snh_20, sog0_43, sog0_44, \
                         sog1_43, sog1_44, soh_60, soh_61, soh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * sog0_43[k]
                  - f_7 * sog1_43[k]
                  + f_3 * pc_y[k] * soh_60[k];

        t_81[k] = f_8 * sog0_44[k]
                  - f_9 * sog1_44[k]
                  + f_3 * pc_y[k] * soh_61[k];

        t_82[k] = f_3 * pc_y[k] * soh_62[k];

        t_83[k] = f_11 * snh_20[k]
                  + f_1 * sog0_44[k]
                  - f_2 * sog1_44[k]
                  + f_3 * pc_z[k] * soh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, snh_21, snh_63, snh_66, \
                         sog0_45, sog0_48, sog1_45, sog1_48, soh_63, \
                         soh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_16 * snh_63[k]
                  + f_1 * sog0_45[k]
                  - f_2 * sog1_45[k]
                  + f_3 * pc_x[k] * soh_63[k];

        t_85[k] = f_12 * snh_21[k]
                  + f_3 * pc_y[k] * soh_63[k];

        t_86[k] = f_3 * pc_z[k] * soh_63[k];

        t_87[k] = f_16 * snh_66[k]
                  + f_4 * sog0_48[k]
                  - f_5 * sog1_48[k]
                  + f_3 * pc_x[k] * soh_66[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pc_x, pc_y, snh_23, snh_68, snh_69, sog0_50, \
                         sog0_51, sog1_50, sog1_51, soh_65, soh_68, \
                         soh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * snh_23[k]
                  + f_3 * pc_y[k] * soh_65[k];

        t_89[k] = f_16 * snh_68[k]
                  + f_4 * sog0_50[k]
                  - f_5 * sog1_50[k]
                  + f_3 * pc_x[k] * soh_68[k];

        t_90[k] = f_16 * snh_69[k]
                  + f_6 * sog0_51[k]
                  - f_7 * sog1_51[k]
                  + f_3 * pc_x[k] * soh_69[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pc_x, pc_y, pc_z, snh_26, snh_72, sog0_54, sog1_54, \
                         soh_66, soh_68, soh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_z[k] * soh_66[k];

        t_92[k] = f_12 * snh_26[k]
                  + f_3 * pc_y[k] * soh_68[k];

        t_93[k] = f_16 * snh_72[k]
                  + f_6 * sog0_54[k]
                  - f_7 * sog1_54[k]
                  + f_3 * pc_x[k] * soh_72[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pc_x, pc_z, snh_73, snh_75, sog0_55, sog0_57, \
                         sog1_55, sog1_57, soh_69, soh_73, soh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_16 * snh_73[k]
                  + f_8 * sog0_55[k]
                  - f_9 * sog1_55[k]
                  + f_3 * pc_x[k] * soh_73[k];

        t_95[k] = f_3 * pc_z[k] * soh_69[k];

        t_96[k] = f_16 * snh_75[k]
                  + f_8 * sog0_57[k]
                  - f_9 * sog1_57[k]
                  + f_3 * pc_x[k] * soh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, snh_30, snh_77, snh_78, snh_79, \
                         sog0_59, sog1_59, soh_72, soh_77, soh_78, \
                         soh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_12 * snh_30[k]
                  + f_3 * pc_y[k] * soh_72[k];

        t_98[k] = f_16 * snh_77[k]
                  + f_8 * sog0_59[k]
                  - f_9 * sog1_59[k]
                  + f_3 * pc_x[k] * soh_77[k];

        t_99[k] = f_16 * snh_78[k]
                  + f_3 * pc_x[k] * soh_78[k];

        t_100[k] = f_16 * snh_79[k]
                   + f_3 * pc_x[k] * soh_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, snh_80, snh_81, snh_82, snh_83, \
                         soh_80, soh_81, soh_82, soh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_16 * snh_80[k]
                   + f_3 * pc_x[k] * soh_80[k];

        t_102[k] = f_16 * snh_81[k]
                   + f_3 * pc_x[k] * soh_81[k];

        t_103[k] = f_16 * snh_82[k]
                   + f_3 * pc_x[k] * soh_82[k];

        t_104[k] = f_16 * snh_83[k]
                   + f_3 * pc_x[k] * soh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, snh_36, snh_38, sog0_55, sog0_57, \
                         sog1_55, sog1_57, soh_78, soh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * snh_36[k]
                   + f_1 * sog0_55[k]
                   - f_2 * sog1_55[k]
                   + f_3 * pc_y[k] * soh_78[k];

        t_106[k] = f_3 * pc_z[k] * soh_78[k];

        t_107[k] = f_12 * snh_38[k]
                   + f_4 * sog0_57[k]
                   - f_5 * sog1_57[k]
                   + f_3 * pc_y[k] * soh_80[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_y, pc_z, snh_39, snh_40, snh_41, \
                         sog0_58, sog0_59, sog1_58, sog1_59, soh_81, soh_82, \
                         soh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_12 * snh_39[k]
                   + f_6 * sog0_58[k]
                   - f_7 * sog1_58[k]
                   + f_3 * pc_y[k] * soh_81[k];

        t_109[k] = f_12 * snh_40[k]
                   + f_8 * sog0_59[k]
                   - f_9 * sog1_59[k]
                   + f_3 * pc_y[k] * soh_82[k];

        t_110[k] = f_12 * snh_41[k]
                   + f_3 * pc_y[k] * soh_83[k];

        t_111[k] = f_1 * sog0_59[k]
                   - f_2 * sog1_59[k]
                   + f_3 * pc_z[k] * soh_83[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pb_y, pb_z, pc_y, pc_z, sni0_31, sni0_56, \
                         snh_21, snh_42, sni1_31, sni1_56, soh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pb_y[k] * sni0_56[k]
                   - f_10 * pc_y[k] * sni1_56[k];

        t_113[k] = f_11 * snh_42[k]
                   + f_3 * pc_y[k] * soh_84[k];

        t_114[k] = f_11 * snh_21[k]
                   + f_3 * pc_z[k] * soh_84[k];

        t_115[k] = pb_z[k] * sni0_31[k]
                   - f_10 * pc_z[k] * sni1_31[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pb_z, pc_y, pc_z, sni0_34, sni0_61, \
                         snh_24, snh_44, sni1_34, sni1_61, soh_86, \
                         soh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_11 * snh_44[k]
                   + f_3 * pc_y[k] * soh_86[k];

        t_117[k] = pb_y[k] * sni0_61[k]
                   - f_10 * pc_y[k] * sni1_61[k];

        t_118[k] = pb_z[k] * sni0_34[k]
                   - f_10 * pc_z[k] * sni1_34[k];

        t_119[k] = f_11 * snh_24[k]
                   + f_3 * pc_z[k] * soh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, pc_y, pc_z, sni0_38, sni0_65, \
                         snh_27, snh_47, sni1_38, sni1_65, soh_89, \
                         soh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_11 * snh_47[k]
                   + f_3 * pc_y[k] * soh_89[k];

        t_121[k] = pb_y[k] * sni0_65[k]
                   - f_10 * pc_y[k] * sni1_65[k];

        t_122[k] = pb_z[k] * sni0_38[k]
                   - f_10 * pc_z[k] * sni1_38[k];

        t_123[k] = f_11 * snh_27[k]
                   + f_3 * pc_z[k] * soh_90[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;

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

    const auto *sni0_49 = buffer.data(sni0 + 49);
    const auto *sni0_68 = buffer.data(sni0 + 68);
    const auto *sni0_70 = buffer.data(sni0 + 70);
    const auto *sni0_83 = buffer.data(sni0 + 83);
    const auto *sni0_84 = buffer.data(sni0 + 84);
    const auto *sni0_87 = buffer.data(sni0 + 87);
    const auto *sni0_90 = buffer.data(sni0 + 90);
    const auto *sni0_94 = buffer.data(sni0 + 94);
    const auto *sni0_96 = buffer.data(sni0 + 96);
    const auto *sni0_105 = buffer.data(sni0 + 105);
    const auto *sni0_140 = buffer.data(sni0 + 140);
    const auto *sni0_143 = buffer.data(sni0 + 143);
    const auto *sni0_145 = buffer.data(sni0 + 145);
    const auto *sni0_146 = buffer.data(sni0 + 146);
    const auto *sni0_149 = buffer.data(sni0 + 149);
    const auto *sni0_150 = buffer.data(sni0 + 150);
    const auto *sni0_152 = buffer.data(sni0 + 152);
    const auto *sni0_154 = buffer.data(sni0 + 154);

    const auto *snh_36 = buffer.data(snh + 36);
    const auto *snh_42 = buffer.data(snh + 42);
    const auto *snh_45 = buffer.data(snh + 45);
    const auto *snh_48 = buffer.data(snh + 48);
    const auto *snh_50 = buffer.data(snh + 50);
    const auto *snh_51 = buffer.data(snh + 51);
    const auto *snh_57 = buffer.data(snh + 57);
    const auto *snh_59 = buffer.data(snh + 59);
    const auto *snh_60 = buffer.data(snh + 60);
    const auto *snh_61 = buffer.data(snh + 61);
    const auto *snh_62 = buffer.data(snh + 62);
    const auto *snh_63 = buffer.data(snh + 63);
    const auto *snh_65 = buffer.data(snh + 65);
    const auto *snh_66 = buffer.data(snh + 66);
    const auto *snh_68 = buffer.data(snh + 68);
    const auto *snh_69 = buffer.data(snh + 69);
    const auto *snh_70 = buffer.data(snh + 70);
    const auto *snh_72 = buffer.data(snh + 72);
    const auto *snh_78 = buffer.data(snh + 78);
    const auto *snh_80 = buffer.data(snh + 80);
    const auto *snh_81 = buffer.data(snh + 81);
    const auto *snh_82 = buffer.data(snh + 82);
    const auto *snh_83 = buffer.data(snh + 83);
    const auto *snh_84 = buffer.data(snh + 84);
    const auto *snh_86 = buffer.data(snh + 86);
    const auto *snh_87 = buffer.data(snh + 87);
    const auto *snh_89 = buffer.data(snh + 89);
    const auto *snh_90 = buffer.data(snh + 90);
    const auto *snh_93 = buffer.data(snh + 93);
    const auto *snh_99 = buffer.data(snh + 99);
    const auto *snh_100 = buffer.data(snh + 100);
    const auto *snh_101 = buffer.data(snh + 101);
    const auto *snh_102 = buffer.data(snh + 102);
    const auto *snh_103 = buffer.data(snh + 103);
    const auto *snh_104 = buffer.data(snh + 104);
    const auto *snh_105 = buffer.data(snh + 105);
    const auto *snh_106 = buffer.data(snh + 106);
    const auto *snh_107 = buffer.data(snh + 107);
    const auto *snh_108 = buffer.data(snh + 108);
    const auto *snh_110 = buffer.data(snh + 110);
    const auto *snh_111 = buffer.data(snh + 111);
    const auto *snh_113 = buffer.data(snh + 113);
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
    const auto *snh_152 = buffer.data(snh + 152);
    const auto *snh_156 = buffer.data(snh + 156);
    const auto *snh_161 = buffer.data(snh + 161);
    const auto *snh_162 = buffer.data(snh + 162);
    const auto *snh_163 = buffer.data(snh + 163);
    const auto *snh_164 = buffer.data(snh + 164);
    const auto *snh_165 = buffer.data(snh + 165);
    const auto *snh_166 = buffer.data(snh + 166);
    const auto *snh_167 = buffer.data(snh + 167);
    const auto *snh_183 = buffer.data(snh + 183);

    const auto *sni1_49 = buffer.data(sni1 + 49);
    const auto *sni1_68 = buffer.data(sni1 + 68);
    const auto *sni1_70 = buffer.data(sni1 + 70);
    const auto *sni1_83 = buffer.data(sni1 + 83);
    const auto *sni1_84 = buffer.data(sni1 + 84);
    const auto *sni1_87 = buffer.data(sni1 + 87);
    const auto *sni1_90 = buffer.data(sni1 + 90);
    const auto *sni1_94 = buffer.data(sni1 + 94);
    const auto *sni1_96 = buffer.data(sni1 + 96);
    const auto *sni1_105 = buffer.data(sni1 + 105);
    const auto *sni1_140 = buffer.data(sni1 + 140);
    const auto *sni1_143 = buffer.data(sni1 + 143);
    const auto *sni1_145 = buffer.data(sni1 + 145);
    const auto *sni1_146 = buffer.data(sni1 + 146);
    const auto *sni1_149 = buffer.data(sni1 + 149);
    const auto *sni1_150 = buffer.data(sni1 + 150);
    const auto *sni1_152 = buffer.data(sni1 + 152);
    const auto *sni1_154 = buffer.data(sni1 + 154);

    const auto *sog0_72 = buffer.data(sog0 + 72);
    const auto *sog0_73 = buffer.data(sog0 + 73);
    const auto *sog0_74 = buffer.data(sog0 + 74);
    const auto *sog0_75 = buffer.data(sog0 + 75);
    const auto *sog0_78 = buffer.data(sog0 + 78);
    const auto *sog0_80 = buffer.data(sog0 + 80);
    const auto *sog0_81 = buffer.data(sog0 + 81);
    const auto *sog0_84 = buffer.data(sog0 + 84);
    const auto *sog0_85 = buffer.data(sog0 + 85);
    const auto *sog0_87 = buffer.data(sog0 + 87);
    const auto *sog0_88 = buffer.data(sog0 + 88);
    const auto *sog0_89 = buffer.data(sog0 + 89);
    const auto *sog0_90 = buffer.data(sog0 + 90);
    const auto *sog0_93 = buffer.data(sog0 + 93);
    const auto *sog0_95 = buffer.data(sog0 + 95);
    const auto *sog0_96 = buffer.data(sog0 + 96);
    const auto *sog0_99 = buffer.data(sog0 + 99);
    const auto *sog0_100 = buffer.data(sog0 + 100);
    const auto *sog0_102 = buffer.data(sog0 + 102);
    const auto *sog0_103 = buffer.data(sog0 + 103);
    const auto *sog0_104 = buffer.data(sog0 + 104);
    const auto *sog0_110 = buffer.data(sog0 + 110);
    const auto *sog0_114 = buffer.data(sog0 + 114);
    const auto *sog0_117 = buffer.data(sog0 + 117);
    const auto *sog0_118 = buffer.data(sog0 + 118);
    const auto *sog0_119 = buffer.data(sog0 + 119);

    const auto *sog1_72 = buffer.data(sog1 + 72);
    const auto *sog1_73 = buffer.data(sog1 + 73);
    const auto *sog1_74 = buffer.data(sog1 + 74);
    const auto *sog1_75 = buffer.data(sog1 + 75);
    const auto *sog1_78 = buffer.data(sog1 + 78);
    const auto *sog1_80 = buffer.data(sog1 + 80);
    const auto *sog1_81 = buffer.data(sog1 + 81);
    const auto *sog1_84 = buffer.data(sog1 + 84);
    const auto *sog1_85 = buffer.data(sog1 + 85);
    const auto *sog1_87 = buffer.data(sog1 + 87);
    const auto *sog1_88 = buffer.data(sog1 + 88);
    const auto *sog1_89 = buffer.data(sog1 + 89);
    const auto *sog1_90 = buffer.data(sog1 + 90);
    const auto *sog1_93 = buffer.data(sog1 + 93);
    const auto *sog1_95 = buffer.data(sog1 + 95);
    const auto *sog1_96 = buffer.data(sog1 + 96);
    const auto *sog1_99 = buffer.data(sog1 + 99);
    const auto *sog1_100 = buffer.data(sog1 + 100);
    const auto *sog1_102 = buffer.data(sog1 + 102);
    const auto *sog1_103 = buffer.data(sog1 + 103);
    const auto *sog1_104 = buffer.data(sog1 + 104);
    const auto *sog1_110 = buffer.data(sog1 + 110);
    const auto *sog1_114 = buffer.data(sog1 + 114);
    const auto *sog1_117 = buffer.data(sog1 + 117);
    const auto *sog1_118 = buffer.data(sog1 + 118);
    const auto *sog1_119 = buffer.data(sog1 + 119);

    const auto *soh_93 = buffer.data(soh + 93);
    const auto *soh_99 = buffer.data(soh + 99);
    const auto *soh_100 = buffer.data(soh + 100);
    const auto *soh_101 = buffer.data(soh + 101);
    const auto *soh_102 = buffer.data(soh + 102);
    const auto *soh_103 = buffer.data(soh + 103);
    const auto *soh_104 = buffer.data(soh + 104);
    const auto *soh_105 = buffer.data(soh + 105);
    const auto *soh_107 = buffer.data(soh + 107);
    const auto *soh_108 = buffer.data(soh + 108);
    const auto *soh_110 = buffer.data(soh + 110);
    const auto *soh_111 = buffer.data(soh + 111);
    const auto *soh_114 = buffer.data(soh + 114);
    const auto *soh_115 = buffer.data(soh + 115);
    const auto *soh_117 = buffer.data(soh + 117);
    const auto *soh_119 = buffer.data(soh + 119);
    const auto *soh_120 = buffer.data(soh + 120);
    const auto *soh_121 = buffer.data(soh + 121);
    const auto *soh_122 = buffer.data(soh + 122);
    const auto *soh_123 = buffer.data(soh + 123);
    const auto *soh_124 = buffer.data(soh + 124);
    const auto *soh_125 = buffer.data(soh + 125);
    const auto *soh_126 = buffer.data(soh + 126);
    const auto *soh_128 = buffer.data(soh + 128);
    const auto *soh_129 = buffer.data(soh + 129);
    const auto *soh_131 = buffer.data(soh + 131);
    const auto *soh_132 = buffer.data(soh + 132);
    const auto *soh_135 = buffer.data(soh + 135);
    const auto *soh_136 = buffer.data(soh + 136);
    const auto *soh_138 = buffer.data(soh + 138);
    const auto *soh_140 = buffer.data(soh + 140);
    const auto *soh_141 = buffer.data(soh + 141);
    const auto *soh_142 = buffer.data(soh + 142);
    const auto *soh_143 = buffer.data(soh + 143);
    const auto *soh_144 = buffer.data(soh + 144);
    const auto *soh_145 = buffer.data(soh + 145);
    const auto *soh_146 = buffer.data(soh + 146);
    const auto *soh_147 = buffer.data(soh + 147);
    const auto *soh_149 = buffer.data(soh + 149);
    const auto *soh_150 = buffer.data(soh + 150);
    const auto *soh_152 = buffer.data(soh + 152);
    const auto *soh_153 = buffer.data(soh + 153);
    const auto *soh_156 = buffer.data(soh + 156);
    const auto *soh_161 = buffer.data(soh + 161);
    const auto *soh_162 = buffer.data(soh + 162);
    const auto *soh_163 = buffer.data(soh + 163);
    const auto *soh_164 = buffer.data(soh + 164);
    const auto *soh_165 = buffer.data(soh + 165);
    const auto *soh_166 = buffer.data(soh + 166);
    const auto *soh_167 = buffer.data(soh + 167);
    const auto *soh_168 = buffer.data(soh + 168);
    const auto *soh_170 = buffer.data(soh + 170);
    const auto *soh_171 = buffer.data(soh + 171);
    const auto *soh_173 = buffer.data(soh + 173);
    const auto *soh_174 = buffer.data(soh + 174);
    const auto *soh_177 = buffer.data(soh + 177);
    const auto *soh_183 = buffer.data(soh + 183);

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_y, pc_x, pc_y, sni0_68, sni0_70, \
                         snh_50, snh_51, snh_99, sni1_68, sni1_70, soh_93, \
                         soh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * sni0_68[k]
                   + f_12 * snh_50[k]
                   - f_10 * pc_y[k] * sni1_68[k];

        t_125[k] = f_11 * snh_51[k]
                   + f_3 * pc_y[k] * soh_93[k];

        t_126[k] = pb_y[k] * sni0_70[k]
                   - f_10 * pc_y[k] * sni1_70[k];

        t_127[k] = f_16 * snh_99[k]
                   + f_3 * pc_x[k] * soh_99[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, snh_100, snh_101, snh_102, \
                         snh_103, snh_104, soh_100, soh_101, soh_102, soh_103, \
                         soh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_16 * snh_100[k]
                   + f_3 * pc_x[k] * soh_100[k];

        t_129[k] = f_16 * snh_101[k]
                   + f_3 * pc_x[k] * soh_101[k];

        t_130[k] = f_16 * snh_102[k]
                   + f_3 * pc_x[k] * soh_102[k];

        t_131[k] = f_16 * snh_103[k]
                   + f_3 * pc_x[k] * soh_103[k];

        t_132[k] = f_16 * snh_104[k]
                   + f_3 * pc_x[k] * soh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, pc_y, pc_z, sni0_49, snh_36, snh_59, \
                         sni1_49, sog0_72, sog1_72, soh_99, soh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_z[k] * sni0_49[k]
                   - f_10 * pc_z[k] * sni1_49[k];

        t_134[k] = f_11 * snh_36[k]
                   + f_3 * pc_z[k] * soh_99[k];

        t_135[k] = f_11 * snh_59[k]
                   + f_4 * sog0_72[k]
                   - f_5 * sog1_72[k]
                   + f_3 * pc_y[k] * soh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, snh_60, snh_61, snh_62, sog0_73, sog0_74, \
                         sog1_73, sog1_74, soh_102, soh_103, soh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * snh_60[k]
                   + f_6 * sog0_73[k]
                   - f_7 * sog1_73[k]
                   + f_3 * pc_y[k] * soh_102[k];

        t_137[k] = f_11 * snh_61[k]
                   + f_8 * sog0_74[k]
                   - f_9 * sog1_74[k]
                   + f_3 * pc_y[k] * soh_103[k];

        t_138[k] = f_11 * snh_62[k]
                   + f_3 * pc_y[k] * soh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_y, pc_x, pc_y, pc_z, sni0_83, snh_42, \
                         snh_105, sni1_83, sog0_75, sog1_75, soh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * sni0_83[k]
                   - f_10 * pc_y[k] * sni1_83[k];

        t_140[k] = f_16 * snh_105[k]
                   + f_1 * sog0_75[k]
                   - f_2 * sog1_75[k]
                   + f_3 * pc_x[k] * soh_105[k];

        t_141[k] = f_3 * pc_y[k] * soh_105[k];

        t_142[k] = f_12 * snh_42[k]
                   + f_3 * pc_z[k] * soh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, snh_108, snh_110, sog0_78, sog0_80, \
                         sog1_78, sog1_80, soh_107, soh_108, soh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_16 * snh_108[k]
                   + f_4 * sog0_78[k]
                   - f_5 * sog1_78[k]
                   + f_3 * pc_x[k] * soh_108[k];

        t_144[k] = f_3 * pc_y[k] * soh_107[k];

        t_145[k] = f_16 * snh_110[k]
                   + f_4 * sog0_80[k]
                   - f_5 * sog1_80[k]
                   + f_3 * pc_x[k] * soh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_x, pc_y, pc_z, snh_45, snh_111, sog0_81, \
                         sog1_81, soh_108, soh_110, soh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_16 * snh_111[k]
                   + f_6 * sog0_81[k]
                   - f_7 * sog1_81[k]
                   + f_3 * pc_x[k] * soh_111[k];

        t_147[k] = f_12 * snh_45[k]
                   + f_3 * pc_z[k] * soh_108[k];

        t_148[k] = f_3 * pc_y[k] * soh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_z, snh_48, snh_114, snh_115, sog0_84, \
                         sog0_85, sog1_84, sog1_85, soh_111, soh_114, \
                         soh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_16 * snh_114[k]
                   + f_6 * sog0_84[k]
                   - f_7 * sog1_84[k]
                   + f_3 * pc_x[k] * soh_114[k];

        t_150[k] = f_16 * snh_115[k]
                   + f_8 * sog0_85[k]
                   - f_9 * sog1_85[k]
                   + f_3 * pc_x[k] * soh_115[k];

        t_151[k] = f_12 * snh_48[k]
                   + f_3 * pc_z[k] * soh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pc_x, pc_y, snh_117, snh_119, sog0_87, sog0_89, \
                         sog1_87, sog1_89, soh_114, soh_117, soh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_16 * snh_117[k]
                   + f_8 * sog0_87[k]
                   - f_9 * sog1_87[k]
                   + f_3 * pc_x[k] * soh_117[k];

        t_153[k] = f_3 * pc_y[k] * soh_114[k];

        t_154[k] = f_16 * snh_119[k]
                   + f_8 * sog0_89[k]
                   - f_9 * sog1_89[k]
                   + f_3 * pc_x[k] * soh_119[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, snh_120, snh_121, snh_122, \
                         snh_123, snh_124, soh_120, soh_121, soh_122, soh_123, \
                         soh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_16 * snh_120[k]
                   + f_3 * pc_x[k] * soh_120[k];

        t_156[k] = f_16 * snh_121[k]
                   + f_3 * pc_x[k] * soh_121[k];

        t_157[k] = f_16 * snh_122[k]
                   + f_3 * pc_x[k] * soh_122[k];

        t_158[k] = f_16 * snh_123[k]
                   + f_3 * pc_x[k] * soh_123[k];

        t_159[k] = f_16 * snh_124[k]
                   + f_3 * pc_x[k] * soh_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, snh_57, snh_125, \
                         sog0_85, sog0_87, sog1_85, sog1_87, soh_120, soh_122, \
                         soh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * snh_125[k]
                   + f_3 * pc_x[k] * soh_125[k];

        t_161[k] = f_1 * sog0_85[k]
                   - f_2 * sog1_85[k]
                   + f_3 * pc_y[k] * soh_120[k];

        t_162[k] = f_12 * snh_57[k]
                   + f_3 * pc_z[k] * soh_120[k];

        t_163[k] = f_4 * sog0_87[k]
                   - f_5 * sog1_87[k]
                   + f_3 * pc_y[k] * soh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, snh_62, sog0_88, sog0_89, \
                         sog1_88, sog1_89, soh_123, soh_124, soh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * sog0_88[k]
                   - f_7 * sog1_88[k]
                   + f_3 * pc_y[k] * soh_123[k];

        t_165[k] = f_8 * sog0_89[k]
                   - f_9 * sog1_89[k]
                   + f_3 * pc_y[k] * soh_124[k];

        t_166[k] = f_3 * pc_y[k] * soh_125[k];

        t_167[k] = f_12 * snh_62[k]
                   + f_1 * sog0_89[k]
                   - f_2 * sog1_89[k]
                   + f_3 * pc_z[k] * soh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, snh_63, snh_126, \
                         snh_129, sog0_90, sog0_93, sog1_90, sog1_93, soh_126, \
                         soh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_17 * snh_126[k]
                   + f_1 * sog0_90[k]
                   - f_2 * sog1_90[k]
                   + f_3 * pc_x[k] * soh_126[k];

        t_169[k] = f_13 * snh_63[k]
                   + f_3 * pc_y[k] * soh_126[k];

        t_170[k] = f_3 * pc_z[k] * soh_126[k];

        t_171[k] = f_17 * snh_129[k]
                   + f_4 * sog0_93[k]
                   - f_5 * sog1_93[k]
                   + f_3 * pc_x[k] * soh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pc_x, pc_y, snh_65, snh_131, snh_132, sog0_95, \
                         sog0_96, sog1_95, sog1_96, soh_128, soh_131, \
                         soh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_13 * snh_65[k]
                   + f_3 * pc_y[k] * soh_128[k];

        t_173[k] = f_17 * snh_131[k]
                   + f_4 * sog0_95[k]
                   - f_5 * sog1_95[k]
                   + f_3 * pc_x[k] * soh_131[k];

        t_174[k] = f_17 * snh_132[k]
                   + f_6 * sog0_96[k]
                   - f_7 * sog1_96[k]
                   + f_3 * pc_x[k] * soh_132[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pc_x, pc_y, pc_z, snh_68, snh_135, sog0_99, \
                         sog1_99, soh_129, soh_131, soh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_3 * pc_z[k] * soh_129[k];

        t_176[k] = f_13 * snh_68[k]
                   + f_3 * pc_y[k] * soh_131[k];

        t_177[k] = f_17 * snh_135[k]
                   + f_6 * sog0_99[k]
                   - f_7 * sog1_99[k]
                   + f_3 * pc_x[k] * soh_135[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_z, snh_136, snh_138, sog0_100, \
                         sog0_102, sog1_100, sog1_102, soh_132, soh_136, \
                         soh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_17 * snh_136[k]
                   + f_8 * sog0_100[k]
                   - f_9 * sog1_100[k]
                   + f_3 * pc_x[k] * soh_136[k];

        t_179[k] = f_3 * pc_z[k] * soh_132[k];

        t_180[k] = f_17 * snh_138[k]
                   + f_8 * sog0_102[k]
                   - f_9 * sog1_102[k]
                   + f_3 * pc_x[k] * soh_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, snh_72, snh_140, snh_141, \
                         snh_142, sog0_104, sog1_104, soh_135, soh_140, soh_141, \
                         soh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_13 * snh_72[k]
                   + f_3 * pc_y[k] * soh_135[k];

        t_182[k] = f_17 * snh_140[k]
                   + f_8 * sog0_104[k]
                   - f_9 * sog1_104[k]
                   + f_3 * pc_x[k] * soh_140[k];

        t_183[k] = f_17 * snh_141[k]
                   + f_3 * pc_x[k] * soh_141[k];

        t_184[k] = f_17 * snh_142[k]
                   + f_3 * pc_x[k] * soh_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, snh_143, snh_144, snh_145, snh_146, \
                         soh_143, soh_144, soh_145, soh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_17 * snh_143[k]
                   + f_3 * pc_x[k] * soh_143[k];

        t_186[k] = f_17 * snh_144[k]
                   + f_3 * pc_x[k] * soh_144[k];

        t_187[k] = f_17 * snh_145[k]
                   + f_3 * pc_x[k] * soh_145[k];

        t_188[k] = f_17 * snh_146[k]
                   + f_3 * pc_x[k] * soh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, pc_z, snh_78, snh_80, sog0_100, sog0_102, \
                         sog1_100, sog1_102, soh_141, soh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * snh_78[k]
                   + f_1 * sog0_100[k]
                   - f_2 * sog1_100[k]
                   + f_3 * pc_y[k] * soh_141[k];

        t_190[k] = f_3 * pc_z[k] * soh_141[k];

        t_191[k] = f_13 * snh_80[k]
                   + f_4 * sog0_102[k]
                   - f_5 * sog1_102[k]
                   + f_3 * pc_y[k] * soh_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, snh_81, snh_82, snh_83, \
                         sog0_103, sog0_104, sog1_103, sog1_104, soh_144, soh_145, \
                         soh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_13 * snh_81[k]
                   + f_6 * sog0_103[k]
                   - f_7 * sog1_103[k]
                   + f_3 * pc_y[k] * soh_144[k];

        t_193[k] = f_13 * snh_82[k]
                   + f_8 * sog0_104[k]
                   - f_9 * sog1_104[k]
                   + f_3 * pc_y[k] * soh_145[k];

        t_194[k] = f_13 * snh_83[k]
                   + f_3 * pc_y[k] * soh_146[k];

        t_195[k] = f_1 * sog0_104[k]
                   - f_2 * sog1_104[k]
                   + f_3 * pc_z[k] * soh_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_z, pc_y, pc_z, sni0_84, sni0_87, \
                         snh_63, snh_84, sni1_84, sni1_87, soh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_z[k] * sni0_84[k]
                   - f_10 * pc_z[k] * sni1_84[k];

        t_197[k] = f_12 * snh_84[k]
                   + f_3 * pc_y[k] * soh_147[k];

        t_198[k] = f_11 * snh_63[k]
                   + f_3 * pc_z[k] * soh_147[k];

        t_199[k] = pb_z[k] * sni0_87[k]
                   - f_10 * pc_z[k] * sni1_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_z, pc_x, pc_y, pc_z, sni0_90, snh_86, \
                         snh_152, sni1_90, sog0_110, sog1_110, soh_149, \
                         soh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * snh_86[k]
                   + f_3 * pc_y[k] * soh_149[k];

        t_201[k] = f_17 * snh_152[k]
                   + f_4 * sog0_110[k]
                   - f_5 * sog1_110[k]
                   + f_3 * pc_x[k] * soh_152[k];

        t_202[k] = pb_z[k] * sni0_90[k]
                   - f_10 * pc_z[k] * sni1_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pc_x, pc_y, pc_z, snh_66, snh_89, snh_156, \
                         sog0_114, sog1_114, soh_150, soh_152, \
                         soh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_11 * snh_66[k]
                   + f_3 * pc_z[k] * soh_150[k];

        t_204[k] = f_12 * snh_89[k]
                   + f_3 * pc_y[k] * soh_152[k];

        t_205[k] = f_17 * snh_156[k]
                   + f_6 * sog0_114[k]
                   - f_7 * sog1_114[k]
                   + f_3 * pc_x[k] * soh_156[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_z, pc_y, pc_z, sni0_94, sni0_96, \
                         snh_69, snh_70, snh_93, sni1_94, sni1_96, soh_153, \
                         soh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pb_z[k] * sni0_94[k]
                   - f_10 * pc_z[k] * sni1_94[k];

        t_207[k] = f_11 * snh_69[k]
                   + f_3 * pc_z[k] * soh_153[k];

        t_208[k] = pb_z[k] * sni0_96[k]
                   + f_12 * snh_70[k]
                   - f_10 * pc_z[k] * sni1_96[k];

        t_209[k] = f_12 * snh_93[k]
                   + f_3 * pc_y[k] * soh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, snh_161, snh_162, snh_163, snh_164, \
                         sog0_119, sog1_119, soh_161, soh_162, soh_163, \
                         soh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * snh_161[k]
                   + f_8 * sog0_119[k]
                   - f_9 * sog1_119[k]
                   + f_3 * pc_x[k] * soh_161[k];

        t_211[k] = f_17 * snh_162[k]
                   + f_3 * pc_x[k] * soh_162[k];

        t_212[k] = f_17 * snh_163[k]
                   + f_3 * pc_x[k] * soh_163[k];

        t_213[k] = f_17 * snh_164[k]
                   + f_3 * pc_x[k] * soh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, sni0_105, snh_165, \
                         snh_166, snh_167, sni1_105, soh_165, soh_166, \
                         soh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_17 * snh_165[k]
                   + f_3 * pc_x[k] * soh_165[k];

        t_215[k] = f_17 * snh_166[k]
                   + f_3 * pc_x[k] * soh_166[k];

        t_216[k] = f_17 * snh_167[k]
                   + f_3 * pc_x[k] * soh_167[k];

        t_217[k] = pb_z[k] * sni0_105[k]
                   - f_10 * pc_z[k] * sni1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, snh_78, snh_101, snh_102, sog0_117, \
                         sog0_118, sog1_117, sog1_118, soh_162, soh_164, \
                         soh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * snh_78[k]
                   + f_3 * pc_z[k] * soh_162[k];

        t_219[k] = f_12 * snh_101[k]
                   + f_4 * sog0_117[k]
                   - f_5 * sog1_117[k]
                   + f_3 * pc_y[k] * soh_164[k];

        t_220[k] = f_12 * snh_102[k]
                   + f_6 * sog0_118[k]
                   - f_7 * sog1_118[k]
                   + f_3 * pc_y[k] * soh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, pc_z, sni0_140, snh_83, \
                         snh_103, snh_104, sni1_140, sog0_119, sog1_119, soh_166, \
                         soh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * snh_103[k]
                   + f_8 * sog0_119[k]
                   - f_9 * sog1_119[k]
                   + f_3 * pc_y[k] * soh_166[k];

        t_222[k] = f_12 * snh_104[k]
                   + f_3 * pc_y[k] * soh_167[k];

        t_223[k] = f_11 * snh_83[k]
                   + f_1 * sog0_119[k]
                   - f_2 * sog1_119[k]
                   + f_3 * pc_z[k] * soh_167[k];

        t_224[k] = pb_y[k] * sni0_140[k]
                   - f_10 * pc_y[k] * sni1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_y, pc_y, pc_z, sni0_143, snh_84, \
                         snh_105, snh_106, snh_107, sni1_143, soh_168, \
                         soh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * snh_105[k]
                   + f_3 * pc_y[k] * soh_168[k];

        t_226[k] = f_12 * snh_84[k]
                   + f_3 * pc_z[k] * soh_168[k];

        t_227[k] = pb_y[k] * sni0_143[k]
                   + f_12 * snh_106[k]
                   - f_10 * pc_y[k] * sni1_143[k];

        t_228[k] = f_11 * snh_107[k]
                   + f_3 * pc_y[k] * soh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_y, pc_y, pc_z, sni0_145, sni0_146, \
                         snh_87, snh_108, snh_110, sni1_145, sni1_146, soh_171, \
                         soh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * sni0_145[k]
                   - f_10 * pc_y[k] * sni1_145[k];

        t_230[k] = pb_y[k] * sni0_146[k]
                   + f_13 * snh_108[k]
                   - f_10 * pc_y[k] * sni1_146[k];

        t_231[k] = f_12 * snh_87[k]
                   + f_3 * pc_z[k] * soh_171[k];

        t_232[k] = f_11 * snh_110[k]
                   + f_3 * pc_y[k] * soh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pc_y, pc_z, sni0_149, sni0_150, snh_90, \
                         snh_111, sni1_149, sni1_150, soh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * sni0_149[k]
                   - f_10 * pc_y[k] * sni1_149[k];

        t_234[k] = pb_y[k] * sni0_150[k]
                   + f_14 * snh_111[k]
                   - f_10 * pc_y[k] * sni1_150[k];

        t_235[k] = f_12 * snh_90[k]
                   + f_3 * pc_z[k] * soh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_y, pc_x, pc_y, sni0_152, sni0_154, \
                         snh_113, snh_114, snh_183, sni1_152, sni1_154, soh_177, \
                         soh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * sni0_152[k]
                   + f_12 * snh_113[k]
                   - f_10 * pc_y[k] * sni1_152[k];

        t_237[k] = f_11 * snh_114[k]
                   + f_3 * pc_y[k] * soh_177[k];

        t_238[k] = pb_y[k] * sni0_154[k]
                   - f_10 * pc_y[k] * sni1_154[k];

        t_239[k] = f_17 * snh_183[k]
                   + f_3 * pc_x[k] * soh_183[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;

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

    const auto *sni0_167 = buffer.data(sni0 + 167);
    const auto *sni0_168 = buffer.data(sni0 + 168);
    const auto *sni0_171 = buffer.data(sni0 + 171);
    const auto *sni0_174 = buffer.data(sni0 + 174);
    const auto *sni0_178 = buffer.data(sni0 + 178);
    const auto *sni0_180 = buffer.data(sni0 + 180);
    const auto *sni0_189 = buffer.data(sni0 + 189);

    const auto *snh_99 = buffer.data(snh + 99);
    const auto *snh_105 = buffer.data(snh + 105);
    const auto *snh_108 = buffer.data(snh + 108);
    const auto *snh_111 = buffer.data(snh + 111);
    const auto *snh_120 = buffer.data(snh + 120);
    const auto *snh_122 = buffer.data(snh + 122);
    const auto *snh_123 = buffer.data(snh + 123);
    const auto *snh_124 = buffer.data(snh + 124);
    const auto *snh_125 = buffer.data(snh + 125);
    const auto *snh_126 = buffer.data(snh + 126);
    const auto *snh_128 = buffer.data(snh + 128);
    const auto *snh_129 = buffer.data(snh + 129);
    const auto *snh_131 = buffer.data(snh + 131);
    const auto *snh_132 = buffer.data(snh + 132);
    const auto *snh_133 = buffer.data(snh + 133);
    const auto *snh_135 = buffer.data(snh + 135);
    const auto *snh_141 = buffer.data(snh + 141);
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
    const auto *snh_164 = buffer.data(snh + 164);
    const auto *snh_165 = buffer.data(snh + 165);
    const auto *snh_166 = buffer.data(snh + 166);
    const auto *snh_167 = buffer.data(snh + 167);
    const auto *snh_168 = buffer.data(snh + 168);
    const auto *snh_170 = buffer.data(snh + 170);
    const auto *snh_173 = buffer.data(snh + 173);
    const auto *snh_177 = buffer.data(snh + 177);
    const auto *snh_184 = buffer.data(snh + 184);
    const auto *snh_185 = buffer.data(snh + 185);
    const auto *snh_186 = buffer.data(snh + 186);
    const auto *snh_187 = buffer.data(snh + 187);
    const auto *snh_188 = buffer.data(snh + 188);
    const auto *snh_189 = buffer.data(snh + 189);
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
    const auto *snh_236 = buffer.data(snh + 236);
    const auto *snh_240 = buffer.data(snh + 240);
    const auto *snh_245 = buffer.data(snh + 245);
    const auto *snh_246 = buffer.data(snh + 246);
    const auto *snh_247 = buffer.data(snh + 247);
    const auto *snh_248 = buffer.data(snh + 248);
    const auto *snh_249 = buffer.data(snh + 249);
    const auto *snh_250 = buffer.data(snh + 250);
    const auto *snh_251 = buffer.data(snh + 251);
    const auto *snh_252 = buffer.data(snh + 252);
    const auto *snh_255 = buffer.data(snh + 255);
    const auto *snh_257 = buffer.data(snh + 257);
    const auto *snh_258 = buffer.data(snh + 258);
    const auto *snh_261 = buffer.data(snh + 261);
    const auto *snh_262 = buffer.data(snh + 262);
    const auto *snh_264 = buffer.data(snh + 264);
    const auto *snh_266 = buffer.data(snh + 266);

    const auto *sni1_167 = buffer.data(sni1 + 167);
    const auto *sni1_168 = buffer.data(sni1 + 168);
    const auto *sni1_171 = buffer.data(sni1 + 171);
    const auto *sni1_174 = buffer.data(sni1 + 174);
    const auto *sni1_178 = buffer.data(sni1 + 178);
    const auto *sni1_180 = buffer.data(sni1 + 180);
    const auto *sni1_189 = buffer.data(sni1 + 189);

    const auto *sog0_130 = buffer.data(sog0 + 130);
    const auto *sog0_132 = buffer.data(sog0 + 132);
    const auto *sog0_133 = buffer.data(sog0 + 133);
    const auto *sog0_134 = buffer.data(sog0 + 134);
    const auto *sog0_135 = buffer.data(sog0 + 135);
    const auto *sog0_138 = buffer.data(sog0 + 138);
    const auto *sog0_140 = buffer.data(sog0 + 140);
    const auto *sog0_141 = buffer.data(sog0 + 141);
    const auto *sog0_144 = buffer.data(sog0 + 144);
    const auto *sog0_145 = buffer.data(sog0 + 145);
    const auto *sog0_147 = buffer.data(sog0 + 147);
    const auto *sog0_148 = buffer.data(sog0 + 148);
    const auto *sog0_149 = buffer.data(sog0 + 149);
    const auto *sog0_150 = buffer.data(sog0 + 150);
    const auto *sog0_153 = buffer.data(sog0 + 153);
    const auto *sog0_155 = buffer.data(sog0 + 155);
    const auto *sog0_156 = buffer.data(sog0 + 156);
    const auto *sog0_159 = buffer.data(sog0 + 159);
    const auto *sog0_160 = buffer.data(sog0 + 160);
    const auto *sog0_162 = buffer.data(sog0 + 162);
    const auto *sog0_163 = buffer.data(sog0 + 163);
    const auto *sog0_164 = buffer.data(sog0 + 164);
    const auto *sog0_170 = buffer.data(sog0 + 170);
    const auto *sog0_174 = buffer.data(sog0 + 174);
    const auto *sog0_177 = buffer.data(sog0 + 177);
    const auto *sog0_178 = buffer.data(sog0 + 178);
    const auto *sog0_179 = buffer.data(sog0 + 179);
    const auto *sog0_180 = buffer.data(sog0 + 180);
    const auto *sog0_183 = buffer.data(sog0 + 183);
    const auto *sog0_185 = buffer.data(sog0 + 185);
    const auto *sog0_186 = buffer.data(sog0 + 186);
    const auto *sog0_189 = buffer.data(sog0 + 189);
    const auto *sog0_190 = buffer.data(sog0 + 190);
    const auto *sog0_192 = buffer.data(sog0 + 192);
    const auto *sog0_194 = buffer.data(sog0 + 194);

    const auto *sog1_130 = buffer.data(sog1 + 130);
    const auto *sog1_132 = buffer.data(sog1 + 132);
    const auto *sog1_133 = buffer.data(sog1 + 133);
    const auto *sog1_134 = buffer.data(sog1 + 134);
    const auto *sog1_135 = buffer.data(sog1 + 135);
    const auto *sog1_138 = buffer.data(sog1 + 138);
    const auto *sog1_140 = buffer.data(sog1 + 140);
    const auto *sog1_141 = buffer.data(sog1 + 141);
    const auto *sog1_144 = buffer.data(sog1 + 144);
    const auto *sog1_145 = buffer.data(sog1 + 145);
    const auto *sog1_147 = buffer.data(sog1 + 147);
    const auto *sog1_148 = buffer.data(sog1 + 148);
    const auto *sog1_149 = buffer.data(sog1 + 149);
    const auto *sog1_150 = buffer.data(sog1 + 150);
    const auto *sog1_153 = buffer.data(sog1 + 153);
    const auto *sog1_155 = buffer.data(sog1 + 155);
    const auto *sog1_156 = buffer.data(sog1 + 156);
    const auto *sog1_159 = buffer.data(sog1 + 159);
    const auto *sog1_160 = buffer.data(sog1 + 160);
    const auto *sog1_162 = buffer.data(sog1 + 162);
    const auto *sog1_163 = buffer.data(sog1 + 163);
    const auto *sog1_164 = buffer.data(sog1 + 164);
    const auto *sog1_170 = buffer.data(sog1 + 170);
    const auto *sog1_174 = buffer.data(sog1 + 174);
    const auto *sog1_177 = buffer.data(sog1 + 177);
    const auto *sog1_178 = buffer.data(sog1 + 178);
    const auto *sog1_179 = buffer.data(sog1 + 179);
    const auto *sog1_180 = buffer.data(sog1 + 180);
    const auto *sog1_183 = buffer.data(sog1 + 183);
    const auto *sog1_185 = buffer.data(sog1 + 185);
    const auto *sog1_186 = buffer.data(sog1 + 186);
    const auto *sog1_189 = buffer.data(sog1 + 189);
    const auto *sog1_190 = buffer.data(sog1 + 190);
    const auto *sog1_192 = buffer.data(sog1 + 192);
    const auto *sog1_194 = buffer.data(sog1 + 194);

    const auto *soh_183 = buffer.data(soh + 183);
    const auto *soh_184 = buffer.data(soh + 184);
    const auto *soh_185 = buffer.data(soh + 185);
    const auto *soh_186 = buffer.data(soh + 186);
    const auto *soh_187 = buffer.data(soh + 187);
    const auto *soh_188 = buffer.data(soh + 188);
    const auto *soh_189 = buffer.data(soh + 189);
    const auto *soh_191 = buffer.data(soh + 191);
    const auto *soh_192 = buffer.data(soh + 192);
    const auto *soh_194 = buffer.data(soh + 194);
    const auto *soh_195 = buffer.data(soh + 195);
    const auto *soh_198 = buffer.data(soh + 198);
    const auto *soh_199 = buffer.data(soh + 199);
    const auto *soh_201 = buffer.data(soh + 201);
    const auto *soh_203 = buffer.data(soh + 203);
    const auto *soh_204 = buffer.data(soh + 204);
    const auto *soh_205 = buffer.data(soh + 205);
    const auto *soh_206 = buffer.data(soh + 206);
    const auto *soh_207 = buffer.data(soh + 207);
    const auto *soh_208 = buffer.data(soh + 208);
    const auto *soh_209 = buffer.data(soh + 209);
    const auto *soh_210 = buffer.data(soh + 210);
    const auto *soh_212 = buffer.data(soh + 212);
    const auto *soh_213 = buffer.data(soh + 213);
    const auto *soh_215 = buffer.data(soh + 215);
    const auto *soh_216 = buffer.data(soh + 216);
    const auto *soh_219 = buffer.data(soh + 219);
    const auto *soh_220 = buffer.data(soh + 220);
    const auto *soh_222 = buffer.data(soh + 222);
    const auto *soh_224 = buffer.data(soh + 224);
    const auto *soh_225 = buffer.data(soh + 225);
    const auto *soh_226 = buffer.data(soh + 226);
    const auto *soh_227 = buffer.data(soh + 227);
    const auto *soh_228 = buffer.data(soh + 228);
    const auto *soh_229 = buffer.data(soh + 229);
    const auto *soh_230 = buffer.data(soh + 230);
    const auto *soh_231 = buffer.data(soh + 231);
    const auto *soh_233 = buffer.data(soh + 233);
    const auto *soh_234 = buffer.data(soh + 234);
    const auto *soh_236 = buffer.data(soh + 236);
    const auto *soh_237 = buffer.data(soh + 237);
    const auto *soh_240 = buffer.data(soh + 240);
    const auto *soh_245 = buffer.data(soh + 245);
    const auto *soh_246 = buffer.data(soh + 246);
    const auto *soh_247 = buffer.data(soh + 247);
    const auto *soh_248 = buffer.data(soh + 248);
    const auto *soh_249 = buffer.data(soh + 249);
    const auto *soh_250 = buffer.data(soh + 250);
    const auto *soh_251 = buffer.data(soh + 251);
    const auto *soh_252 = buffer.data(soh + 252);
    const auto *soh_254 = buffer.data(soh + 254);
    const auto *soh_255 = buffer.data(soh + 255);
    const auto *soh_257 = buffer.data(soh + 257);
    const auto *soh_258 = buffer.data(soh + 258);
    const auto *soh_261 = buffer.data(soh + 261);
    const auto *soh_262 = buffer.data(soh + 262);
    const auto *soh_264 = buffer.data(soh + 264);
    const auto *soh_266 = buffer.data(soh + 266);

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, snh_184, snh_185, snh_186, \
                         snh_187, snh_188, soh_184, soh_185, soh_186, soh_187, \
                         soh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * snh_184[k]
                   + f_3 * pc_x[k] * soh_184[k];

        t_241[k] = f_17 * snh_185[k]
                   + f_3 * pc_x[k] * soh_185[k];

        t_242[k] = f_17 * snh_186[k]
                   + f_3 * pc_x[k] * soh_186[k];

        t_243[k] = f_17 * snh_187[k]
                   + f_3 * pc_x[k] * soh_187[k];

        t_244[k] = f_17 * snh_188[k]
                   + f_3 * pc_x[k] * soh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, snh_99, snh_120, snh_122, sog0_130, \
                         sog0_132, sog1_130, sog1_132, soh_183, \
                         soh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * snh_120[k]
                   + f_1 * sog0_130[k]
                   - f_2 * sog1_130[k]
                   + f_3 * pc_y[k] * soh_183[k];

        t_246[k] = f_12 * snh_99[k]
                   + f_3 * pc_z[k] * soh_183[k];

        t_247[k] = f_11 * snh_122[k]
                   + f_4 * sog0_132[k]
                   - f_5 * sog1_132[k]
                   + f_3 * pc_y[k] * soh_185[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, snh_123, snh_124, snh_125, sog0_133, \
                         sog0_134, sog1_133, sog1_134, soh_186, soh_187, \
                         soh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * snh_123[k]
                   + f_6 * sog0_133[k]
                   - f_7 * sog1_133[k]
                   + f_3 * pc_y[k] * soh_186[k];

        t_249[k] = f_11 * snh_124[k]
                   + f_8 * sog0_134[k]
                   - f_9 * sog1_134[k]
                   + f_3 * pc_y[k] * soh_187[k];

        t_250[k] = f_11 * snh_125[k]
                   + f_3 * pc_y[k] * soh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_x, pc_y, pc_z, sni0_167, \
                         snh_105, snh_189, sni1_167, sog0_135, sog1_135, \
                         soh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * sni0_167[k]
                   - f_10 * pc_y[k] * sni1_167[k];

        t_252[k] = f_17 * snh_189[k]
                   + f_1 * sog0_135[k]
                   - f_2 * sog1_135[k]
                   + f_3 * pc_x[k] * soh_189[k];

        t_253[k] = f_3 * pc_y[k] * soh_189[k];

        t_254[k] = f_13 * snh_105[k]
                   + f_3 * pc_z[k] * soh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, snh_192, snh_194, sog0_138, \
                         sog0_140, sog1_138, sog1_140, soh_191, soh_192, \
                         soh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_17 * snh_192[k]
                   + f_4 * sog0_138[k]
                   - f_5 * sog1_138[k]
                   + f_3 * pc_x[k] * soh_192[k];

        t_256[k] = f_3 * pc_y[k] * soh_191[k];

        t_257[k] = f_17 * snh_194[k]
                   + f_4 * sog0_140[k]
                   - f_5 * sog1_140[k]
                   + f_3 * pc_x[k] * soh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, snh_108, snh_195, sog0_141, \
                         sog1_141, soh_192, soh_194, soh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_17 * snh_195[k]
                   + f_6 * sog0_141[k]
                   - f_7 * sog1_141[k]
                   + f_3 * pc_x[k] * soh_195[k];

        t_259[k] = f_13 * snh_108[k]
                   + f_3 * pc_z[k] * soh_192[k];

        t_260[k] = f_3 * pc_y[k] * soh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, snh_111, snh_198, snh_199, sog0_144, \
                         sog0_145, sog1_144, sog1_145, soh_195, soh_198, \
                         soh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_17 * snh_198[k]
                   + f_6 * sog0_144[k]
                   - f_7 * sog1_144[k]
                   + f_3 * pc_x[k] * soh_198[k];

        t_262[k] = f_17 * snh_199[k]
                   + f_8 * sog0_145[k]
                   - f_9 * sog1_145[k]
                   + f_3 * pc_x[k] * soh_199[k];

        t_263[k] = f_13 * snh_111[k]
                   + f_3 * pc_z[k] * soh_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, snh_201, snh_203, sog0_147, \
                         sog0_149, sog1_147, sog1_149, soh_198, soh_201, \
                         soh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * snh_201[k]
                   + f_8 * sog0_147[k]
                   - f_9 * sog1_147[k]
                   + f_3 * pc_x[k] * soh_201[k];

        t_265[k] = f_3 * pc_y[k] * soh_198[k];

        t_266[k] = f_17 * snh_203[k]
                   + f_8 * sog0_149[k]
                   - f_9 * sog1_149[k]
                   + f_3 * pc_x[k] * soh_203[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pc_x, snh_204, snh_205, snh_206, \
                         snh_207, snh_208, soh_204, soh_205, soh_206, soh_207, \
                         soh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_17 * snh_204[k]
                   + f_3 * pc_x[k] * soh_204[k];

        t_268[k] = f_17 * snh_205[k]
                   + f_3 * pc_x[k] * soh_205[k];

        t_269[k] = f_17 * snh_206[k]
                   + f_3 * pc_x[k] * soh_206[k];

        t_270[k] = f_17 * snh_207[k]
                   + f_3 * pc_x[k] * soh_207[k];

        t_271[k] = f_17 * snh_208[k]
                   + f_3 * pc_x[k] * soh_208[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, snh_120, snh_209, \
                         sog0_145, sog0_147, sog1_145, sog1_147, soh_204, soh_206, \
                         soh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_17 * snh_209[k]
                   + f_3 * pc_x[k] * soh_209[k];

        t_273[k] = f_1 * sog0_145[k]
                   - f_2 * sog1_145[k]
                   + f_3 * pc_y[k] * soh_204[k];

        t_274[k] = f_13 * snh_120[k]
                   + f_3 * pc_z[k] * soh_204[k];

        t_275[k] = f_4 * sog0_147[k]
                   - f_5 * sog1_147[k]
                   + f_3 * pc_y[k] * soh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, snh_125, sog0_148, sog0_149, \
                         sog1_148, sog1_149, soh_207, soh_208, \
                         soh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * sog0_148[k]
                   - f_7 * sog1_148[k]
                   + f_3 * pc_y[k] * soh_207[k];

        t_277[k] = f_8 * sog0_149[k]
                   - f_9 * sog1_149[k]
                   + f_3 * pc_y[k] * soh_208[k];

        t_278[k] = f_3 * pc_y[k] * soh_209[k];

        t_279[k] = f_13 * snh_125[k]
                   + f_1 * sog0_149[k]
                   - f_2 * sog1_149[k]
                   + f_3 * pc_z[k] * soh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, snh_126, snh_210, \
                         snh_213, sog0_150, sog0_153, sog1_150, sog1_153, soh_210, \
                         soh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_18 * snh_210[k]
                   + f_1 * sog0_150[k]
                   - f_2 * sog1_150[k]
                   + f_3 * pc_x[k] * soh_210[k];

        t_281[k] = f_14 * snh_126[k]
                   + f_3 * pc_y[k] * soh_210[k];

        t_282[k] = f_3 * pc_z[k] * soh_210[k];

        t_283[k] = f_18 * snh_213[k]
                   + f_4 * sog0_153[k]
                   - f_5 * sog1_153[k]
                   + f_3 * pc_x[k] * soh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_y, snh_128, snh_215, snh_216, sog0_155, \
                         sog0_156, sog1_155, sog1_156, soh_212, soh_215, \
                         soh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_14 * snh_128[k]
                   + f_3 * pc_y[k] * soh_212[k];

        t_285[k] = f_18 * snh_215[k]
                   + f_4 * sog0_155[k]
                   - f_5 * sog1_155[k]
                   + f_3 * pc_x[k] * soh_215[k];

        t_286[k] = f_18 * snh_216[k]
                   + f_6 * sog0_156[k]
                   - f_7 * sog1_156[k]
                   + f_3 * pc_x[k] * soh_216[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, pc_z, snh_131, snh_219, sog0_159, \
                         sog1_159, soh_213, soh_215, soh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * soh_213[k];

        t_288[k] = f_14 * snh_131[k]
                   + f_3 * pc_y[k] * soh_215[k];

        t_289[k] = f_18 * snh_219[k]
                   + f_6 * sog0_159[k]
                   - f_7 * sog1_159[k]
                   + f_3 * pc_x[k] * soh_219[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, snh_220, snh_222, sog0_160, \
                         sog0_162, sog1_160, sog1_162, soh_216, soh_220, \
                         soh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_18 * snh_220[k]
                   + f_8 * sog0_160[k]
                   - f_9 * sog1_160[k]
                   + f_3 * pc_x[k] * soh_220[k];

        t_291[k] = f_3 * pc_z[k] * soh_216[k];

        t_292[k] = f_18 * snh_222[k]
                   + f_8 * sog0_162[k]
                   - f_9 * sog1_162[k]
                   + f_3 * pc_x[k] * soh_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_y, snh_135, snh_224, snh_225, \
                         snh_226, sog0_164, sog1_164, soh_219, soh_224, soh_225, \
                         soh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * snh_135[k]
                   + f_3 * pc_y[k] * soh_219[k];

        t_294[k] = f_18 * snh_224[k]
                   + f_8 * sog0_164[k]
                   - f_9 * sog1_164[k]
                   + f_3 * pc_x[k] * soh_224[k];

        t_295[k] = f_18 * snh_225[k]
                   + f_3 * pc_x[k] * soh_225[k];

        t_296[k] = f_18 * snh_226[k]
                   + f_3 * pc_x[k] * soh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, snh_227, snh_228, snh_229, snh_230, \
                         soh_227, soh_228, soh_229, soh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_18 * snh_227[k]
                   + f_3 * pc_x[k] * soh_227[k];

        t_298[k] = f_18 * snh_228[k]
                   + f_3 * pc_x[k] * soh_228[k];

        t_299[k] = f_18 * snh_229[k]
                   + f_3 * pc_x[k] * soh_229[k];

        t_300[k] = f_18 * snh_230[k]
                   + f_3 * pc_x[k] * soh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pc_y, pc_z, snh_141, snh_143, sog0_160, \
                         sog0_162, sog1_160, sog1_162, soh_225, \
                         soh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * snh_141[k]
                   + f_1 * sog0_160[k]
                   - f_2 * sog1_160[k]
                   + f_3 * pc_y[k] * soh_225[k];

        t_302[k] = f_3 * pc_z[k] * soh_225[k];

        t_303[k] = f_14 * snh_143[k]
                   + f_4 * sog0_162[k]
                   - f_5 * sog1_162[k]
                   + f_3 * pc_y[k] * soh_227[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_y, pc_z, snh_144, snh_145, snh_146, \
                         sog0_163, sog0_164, sog1_163, sog1_164, soh_228, soh_229, \
                         soh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * snh_144[k]
                   + f_6 * sog0_163[k]
                   - f_7 * sog1_163[k]
                   + f_3 * pc_y[k] * soh_228[k];

        t_305[k] = f_14 * snh_145[k]
                   + f_8 * sog0_164[k]
                   - f_9 * sog1_164[k]
                   + f_3 * pc_y[k] * soh_229[k];

        t_306[k] = f_14 * snh_146[k]
                   + f_3 * pc_y[k] * soh_230[k];

        t_307[k] = f_1 * sog0_164[k]
                   - f_2 * sog1_164[k]
                   + f_3 * pc_z[k] * soh_230[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, sni0_168, sni0_171, \
                         snh_126, snh_147, sni1_168, sni1_171, \
                         soh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * sni0_168[k]
                   - f_10 * pc_z[k] * sni1_168[k];

        t_309[k] = f_13 * snh_147[k]
                   + f_3 * pc_y[k] * soh_231[k];

        t_310[k] = f_11 * snh_126[k]
                   + f_3 * pc_z[k] * soh_231[k];

        t_311[k] = pb_z[k] * sni0_171[k]
                   - f_10 * pc_z[k] * sni1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_z, pc_x, pc_y, pc_z, sni0_174, snh_149, \
                         snh_236, sni1_174, sog0_170, sog1_170, soh_233, \
                         soh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * snh_149[k]
                   + f_3 * pc_y[k] * soh_233[k];

        t_313[k] = f_18 * snh_236[k]
                   + f_4 * sog0_170[k]
                   - f_5 * sog1_170[k]
                   + f_3 * pc_x[k] * soh_236[k];

        t_314[k] = pb_z[k] * sni0_174[k]
                   - f_10 * pc_z[k] * sni1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pc_x, pc_y, pc_z, snh_129, snh_152, snh_240, \
                         sog0_174, sog1_174, soh_234, soh_236, \
                         soh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * snh_129[k]
                   + f_3 * pc_z[k] * soh_234[k];

        t_316[k] = f_13 * snh_152[k]
                   + f_3 * pc_y[k] * soh_236[k];

        t_317[k] = f_18 * snh_240[k]
                   + f_6 * sog0_174[k]
                   - f_7 * sog1_174[k]
                   + f_3 * pc_x[k] * soh_240[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pb_z, pc_y, pc_z, sni0_178, sni0_180, \
                         snh_132, snh_133, snh_156, sni1_178, sni1_180, soh_237, \
                         soh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * sni0_178[k]
                   - f_10 * pc_z[k] * sni1_178[k];

        t_319[k] = f_11 * snh_132[k]
                   + f_3 * pc_z[k] * soh_237[k];

        t_320[k] = pb_z[k] * sni0_180[k]
                   + f_12 * snh_133[k]
                   - f_10 * pc_z[k] * sni1_180[k];

        t_321[k] = f_13 * snh_156[k]
                   + f_3 * pc_y[k] * soh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, snh_245, snh_246, snh_247, snh_248, \
                         sog0_179, sog1_179, soh_245, soh_246, soh_247, \
                         soh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_18 * snh_245[k]
                   + f_8 * sog0_179[k]
                   - f_9 * sog1_179[k]
                   + f_3 * pc_x[k] * soh_245[k];

        t_323[k] = f_18 * snh_246[k]
                   + f_3 * pc_x[k] * soh_246[k];

        t_324[k] = f_18 * snh_247[k]
                   + f_3 * pc_x[k] * soh_247[k];

        t_325[k] = f_18 * snh_248[k]
                   + f_3 * pc_x[k] * soh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, sni0_189, snh_249, \
                         snh_250, snh_251, sni1_189, soh_249, soh_250, \
                         soh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_18 * snh_249[k]
                   + f_3 * pc_x[k] * soh_249[k];

        t_327[k] = f_18 * snh_250[k]
                   + f_3 * pc_x[k] * soh_250[k];

        t_328[k] = f_18 * snh_251[k]
                   + f_3 * pc_x[k] * soh_251[k];

        t_329[k] = pb_z[k] * sni0_189[k]
                   - f_10 * pc_z[k] * sni1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, snh_141, snh_164, snh_165, sog0_177, \
                         sog0_178, sog1_177, sog1_178, soh_246, soh_248, \
                         soh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * snh_141[k]
                   + f_3 * pc_z[k] * soh_246[k];

        t_331[k] = f_13 * snh_164[k]
                   + f_4 * sog0_177[k]
                   - f_5 * sog1_177[k]
                   + f_3 * pc_y[k] * soh_248[k];

        t_332[k] = f_13 * snh_165[k]
                   + f_6 * sog0_178[k]
                   - f_7 * sog1_178[k]
                   + f_3 * pc_y[k] * soh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, snh_146, snh_166, snh_167, sog0_179, \
                         sog1_179, soh_250, soh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * snh_166[k]
                   + f_8 * sog0_179[k]
                   - f_9 * sog1_179[k]
                   + f_3 * pc_y[k] * soh_250[k];

        t_334[k] = f_13 * snh_167[k]
                   + f_3 * pc_y[k] * soh_251[k];

        t_335[k] = f_11 * snh_146[k]
                   + f_1 * sog0_179[k]
                   - f_2 * sog1_179[k]
                   + f_3 * pc_z[k] * soh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, snh_147, snh_168, snh_252, \
                         sog0_180, sog1_180, soh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_18 * snh_252[k]
                   + f_1 * sog0_180[k]
                   - f_2 * sog1_180[k]
                   + f_3 * pc_x[k] * soh_252[k];

        t_337[k] = f_12 * snh_168[k]
                   + f_3 * pc_y[k] * soh_252[k];

        t_338[k] = f_12 * snh_147[k]
                   + f_3 * pc_z[k] * soh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, snh_170, snh_255, snh_257, sog0_183, \
                         sog0_185, sog1_183, sog1_185, soh_254, soh_255, \
                         soh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_18 * snh_255[k]
                   + f_4 * sog0_183[k]
                   - f_5 * sog1_183[k]
                   + f_3 * pc_x[k] * soh_255[k];

        t_340[k] = f_12 * snh_170[k]
                   + f_3 * pc_y[k] * soh_254[k];

        t_341[k] = f_18 * snh_257[k]
                   + f_4 * sog0_185[k]
                   - f_5 * sog1_185[k]
                   + f_3 * pc_x[k] * soh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, snh_150, snh_173, snh_258, \
                         sog0_186, sog1_186, soh_255, soh_257, \
                         soh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_18 * snh_258[k]
                   + f_6 * sog0_186[k]
                   - f_7 * sog1_186[k]
                   + f_3 * pc_x[k] * soh_258[k];

        t_343[k] = f_12 * snh_150[k]
                   + f_3 * pc_z[k] * soh_255[k];

        t_344[k] = f_12 * snh_173[k]
                   + f_3 * pc_y[k] * soh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, snh_153, snh_261, snh_262, sog0_189, \
                         sog0_190, sog1_189, sog1_190, soh_258, soh_261, \
                         soh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * snh_261[k]
                   + f_6 * sog0_189[k]
                   - f_7 * sog1_189[k]
                   + f_3 * pc_x[k] * soh_261[k];

        t_346[k] = f_18 * snh_262[k]
                   + f_8 * sog0_190[k]
                   - f_9 * sog1_190[k]
                   + f_3 * pc_x[k] * soh_262[k];

        t_347[k] = f_12 * snh_153[k]
                   + f_3 * pc_z[k] * soh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, snh_177, snh_264, snh_266, sog0_192, \
                         sog0_194, sog1_192, sog1_194, soh_261, soh_264, \
                         soh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_18 * snh_264[k]
                   + f_8 * sog0_192[k]
                   - f_9 * sog1_192[k]
                   + f_3 * pc_x[k] * soh_264[k];

        t_349[k] = f_12 * snh_177[k]
                   + f_3 * pc_y[k] * soh_261[k];

        t_350[k] = f_18 * snh_266[k]
                   + f_8 * sog0_194[k]
                   - f_9 * sog1_194[k]
                   + f_3 * pc_x[k] * soh_266[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_252 = buffer.data(sni0 + 252);
    const auto *sni0_255 = buffer.data(sni0 + 255);
    const auto *sni0_257 = buffer.data(sni0 + 257);
    const auto *sni0_258 = buffer.data(sni0 + 258);
    const auto *sni0_261 = buffer.data(sni0 + 261);
    const auto *sni0_262 = buffer.data(sni0 + 262);
    const auto *sni0_264 = buffer.data(sni0 + 264);
    const auto *sni0_266 = buffer.data(sni0 + 266);
    const auto *sni0_279 = buffer.data(sni0 + 279);
    const auto *sni0_280 = buffer.data(sni0 + 280);
    const auto *sni0_283 = buffer.data(sni0 + 283);
    const auto *sni0_286 = buffer.data(sni0 + 286);
    const auto *sni0_290 = buffer.data(sni0 + 290);
    const auto *sni0_292 = buffer.data(sni0 + 292);

    const auto *snh_162 = buffer.data(snh + 162);
    const auto *snh_167 = buffer.data(snh + 167);
    const auto *snh_168 = buffer.data(snh + 168);
    const auto *snh_171 = buffer.data(snh + 171);
    const auto *snh_174 = buffer.data(snh + 174);
    const auto *snh_183 = buffer.data(snh + 183);
    const auto *snh_185 = buffer.data(snh + 185);
    const auto *snh_186 = buffer.data(snh + 186);
    const auto *snh_187 = buffer.data(snh + 187);
    const auto *snh_188 = buffer.data(snh + 188);
    const auto *snh_189 = buffer.data(snh + 189);
    const auto *snh_190 = buffer.data(snh + 190);
    const auto *snh_191 = buffer.data(snh + 191);
    const auto *snh_192 = buffer.data(snh + 192);
    const auto *snh_194 = buffer.data(snh + 194);
    const auto *snh_195 = buffer.data(snh + 195);
    const auto *snh_197 = buffer.data(snh + 197);
    const auto *snh_198 = buffer.data(snh + 198);
    const auto *snh_204 = buffer.data(snh + 204);
    const auto *snh_206 = buffer.data(snh + 206);
    const auto *snh_207 = buffer.data(snh + 207);
    const auto *snh_208 = buffer.data(snh + 208);
    const auto *snh_209 = buffer.data(snh + 209);
    const auto *snh_210 = buffer.data(snh + 210);
    const auto *snh_212 = buffer.data(snh + 212);
    const auto *snh_213 = buffer.data(snh + 213);
    const auto *snh_215 = buffer.data(snh + 215);
    const auto *snh_216 = buffer.data(snh + 216);
    const auto *snh_217 = buffer.data(snh + 217);
    const auto *snh_219 = buffer.data(snh + 219);
    const auto *snh_225 = buffer.data(snh + 225);
    const auto *snh_227 = buffer.data(snh + 227);
    const auto *snh_228 = buffer.data(snh + 228);
    const auto *snh_229 = buffer.data(snh + 229);
    const auto *snh_230 = buffer.data(snh + 230);
    const auto *snh_231 = buffer.data(snh + 231);
    const auto *snh_233 = buffer.data(snh + 233);
    const auto *snh_236 = buffer.data(snh + 236);
    const auto *snh_240 = buffer.data(snh + 240);
    const auto *snh_267 = buffer.data(snh + 267);
    const auto *snh_268 = buffer.data(snh + 268);
    const auto *snh_269 = buffer.data(snh + 269);
    const auto *snh_270 = buffer.data(snh + 270);
    const auto *snh_271 = buffer.data(snh + 271);
    const auto *snh_272 = buffer.data(snh + 272);
    const auto *snh_288 = buffer.data(snh + 288);
    const auto *snh_289 = buffer.data(snh + 289);
    const auto *snh_290 = buffer.data(snh + 290);
    const auto *snh_291 = buffer.data(snh + 291);
    const auto *snh_292 = buffer.data(snh + 292);
    const auto *snh_293 = buffer.data(snh + 293);
    const auto *snh_294 = buffer.data(snh + 294);
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
    const auto *snh_341 = buffer.data(snh + 341);
    const auto *snh_345 = buffer.data(snh + 345);
    const auto *snh_350 = buffer.data(snh + 350);
    const auto *snh_351 = buffer.data(snh + 351);
    const auto *snh_352 = buffer.data(snh + 352);
    const auto *snh_353 = buffer.data(snh + 353);

    const auto *sni1_252 = buffer.data(sni1 + 252);
    const auto *sni1_255 = buffer.data(sni1 + 255);
    const auto *sni1_257 = buffer.data(sni1 + 257);
    const auto *sni1_258 = buffer.data(sni1 + 258);
    const auto *sni1_261 = buffer.data(sni1 + 261);
    const auto *sni1_262 = buffer.data(sni1 + 262);
    const auto *sni1_264 = buffer.data(sni1 + 264);
    const auto *sni1_266 = buffer.data(sni1 + 266);
    const auto *sni1_279 = buffer.data(sni1 + 279);
    const auto *sni1_280 = buffer.data(sni1 + 280);
    const auto *sni1_283 = buffer.data(sni1 + 283);
    const auto *sni1_286 = buffer.data(sni1 + 286);
    const auto *sni1_290 = buffer.data(sni1 + 290);
    const auto *sni1_292 = buffer.data(sni1 + 292);

    const auto *sog0_190 = buffer.data(sog0 + 190);
    const auto *sog0_192 = buffer.data(sog0 + 192);
    const auto *sog0_193 = buffer.data(sog0 + 193);
    const auto *sog0_194 = buffer.data(sog0 + 194);
    const auto *sog0_205 = buffer.data(sog0 + 205);
    const auto *sog0_207 = buffer.data(sog0 + 207);
    const auto *sog0_208 = buffer.data(sog0 + 208);
    const auto *sog0_209 = buffer.data(sog0 + 209);
    const auto *sog0_210 = buffer.data(sog0 + 210);
    const auto *sog0_213 = buffer.data(sog0 + 213);
    const auto *sog0_215 = buffer.data(sog0 + 215);
    const auto *sog0_216 = buffer.data(sog0 + 216);
    const auto *sog0_219 = buffer.data(sog0 + 219);
    const auto *sog0_220 = buffer.data(sog0 + 220);
    const auto *sog0_222 = buffer.data(sog0 + 222);
    const auto *sog0_223 = buffer.data(sog0 + 223);
    const auto *sog0_224 = buffer.data(sog0 + 224);
    const auto *sog0_225 = buffer.data(sog0 + 225);
    const auto *sog0_228 = buffer.data(sog0 + 228);
    const auto *sog0_230 = buffer.data(sog0 + 230);
    const auto *sog0_231 = buffer.data(sog0 + 231);
    const auto *sog0_234 = buffer.data(sog0 + 234);
    const auto *sog0_235 = buffer.data(sog0 + 235);
    const auto *sog0_237 = buffer.data(sog0 + 237);
    const auto *sog0_238 = buffer.data(sog0 + 238);
    const auto *sog0_239 = buffer.data(sog0 + 239);
    const auto *sog0_245 = buffer.data(sog0 + 245);
    const auto *sog0_249 = buffer.data(sog0 + 249);
    const auto *sog0_254 = buffer.data(sog0 + 254);

    const auto *sog1_190 = buffer.data(sog1 + 190);
    const auto *sog1_192 = buffer.data(sog1 + 192);
    const auto *sog1_193 = buffer.data(sog1 + 193);
    const auto *sog1_194 = buffer.data(sog1 + 194);
    const auto *sog1_205 = buffer.data(sog1 + 205);
    const auto *sog1_207 = buffer.data(sog1 + 207);
    const auto *sog1_208 = buffer.data(sog1 + 208);
    const auto *sog1_209 = buffer.data(sog1 + 209);
    const auto *sog1_210 = buffer.data(sog1 + 210);
    const auto *sog1_213 = buffer.data(sog1 + 213);
    const auto *sog1_215 = buffer.data(sog1 + 215);
    const auto *sog1_216 = buffer.data(sog1 + 216);
    const auto *sog1_219 = buffer.data(sog1 + 219);
    const auto *sog1_220 = buffer.data(sog1 + 220);
    const auto *sog1_222 = buffer.data(sog1 + 222);
    const auto *sog1_223 = buffer.data(sog1 + 223);
    const auto *sog1_224 = buffer.data(sog1 + 224);
    const auto *sog1_225 = buffer.data(sog1 + 225);
    const auto *sog1_228 = buffer.data(sog1 + 228);
    const auto *sog1_230 = buffer.data(sog1 + 230);
    const auto *sog1_231 = buffer.data(sog1 + 231);
    const auto *sog1_234 = buffer.data(sog1 + 234);
    const auto *sog1_235 = buffer.data(sog1 + 235);
    const auto *sog1_237 = buffer.data(sog1 + 237);
    const auto *sog1_238 = buffer.data(sog1 + 238);
    const auto *sog1_239 = buffer.data(sog1 + 239);
    const auto *sog1_245 = buffer.data(sog1 + 245);
    const auto *sog1_249 = buffer.data(sog1 + 249);
    const auto *sog1_254 = buffer.data(sog1 + 254);

    const auto *soh_267 = buffer.data(soh + 267);
    const auto *soh_268 = buffer.data(soh + 268);
    const auto *soh_269 = buffer.data(soh + 269);
    const auto *soh_270 = buffer.data(soh + 270);
    const auto *soh_271 = buffer.data(soh + 271);
    const auto *soh_272 = buffer.data(soh + 272);
    const auto *soh_273 = buffer.data(soh + 273);
    const auto *soh_275 = buffer.data(soh + 275);
    const auto *soh_276 = buffer.data(soh + 276);
    const auto *soh_278 = buffer.data(soh + 278);
    const auto *soh_279 = buffer.data(soh + 279);
    const auto *soh_282 = buffer.data(soh + 282);
    const auto *soh_288 = buffer.data(soh + 288);
    const auto *soh_289 = buffer.data(soh + 289);
    const auto *soh_290 = buffer.data(soh + 290);
    const auto *soh_291 = buffer.data(soh + 291);
    const auto *soh_292 = buffer.data(soh + 292);
    const auto *soh_293 = buffer.data(soh + 293);
    const auto *soh_294 = buffer.data(soh + 294);
    const auto *soh_296 = buffer.data(soh + 296);
    const auto *soh_297 = buffer.data(soh + 297);
    const auto *soh_299 = buffer.data(soh + 299);
    const auto *soh_300 = buffer.data(soh + 300);
    const auto *soh_303 = buffer.data(soh + 303);
    const auto *soh_304 = buffer.data(soh + 304);
    const auto *soh_306 = buffer.data(soh + 306);
    const auto *soh_308 = buffer.data(soh + 308);
    const auto *soh_309 = buffer.data(soh + 309);
    const auto *soh_310 = buffer.data(soh + 310);
    const auto *soh_311 = buffer.data(soh + 311);
    const auto *soh_312 = buffer.data(soh + 312);
    const auto *soh_313 = buffer.data(soh + 313);
    const auto *soh_314 = buffer.data(soh + 314);
    const auto *soh_315 = buffer.data(soh + 315);
    const auto *soh_317 = buffer.data(soh + 317);
    const auto *soh_318 = buffer.data(soh + 318);
    const auto *soh_320 = buffer.data(soh + 320);
    const auto *soh_321 = buffer.data(soh + 321);
    const auto *soh_324 = buffer.data(soh + 324);
    const auto *soh_325 = buffer.data(soh + 325);
    const auto *soh_327 = buffer.data(soh + 327);
    const auto *soh_329 = buffer.data(soh + 329);
    const auto *soh_330 = buffer.data(soh + 330);
    const auto *soh_331 = buffer.data(soh + 331);
    const auto *soh_332 = buffer.data(soh + 332);
    const auto *soh_333 = buffer.data(soh + 333);
    const auto *soh_334 = buffer.data(soh + 334);
    const auto *soh_335 = buffer.data(soh + 335);
    const auto *soh_336 = buffer.data(soh + 336);
    const auto *soh_338 = buffer.data(soh + 338);
    const auto *soh_339 = buffer.data(soh + 339);
    const auto *soh_341 = buffer.data(soh + 341);
    const auto *soh_342 = buffer.data(soh + 342);
    const auto *soh_345 = buffer.data(soh + 345);
    const auto *soh_350 = buffer.data(soh + 350);
    const auto *soh_351 = buffer.data(soh + 351);
    const auto *soh_352 = buffer.data(soh + 352);
    const auto *soh_353 = buffer.data(soh + 353);

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, snh_267, snh_268, snh_269, \
                         snh_270, snh_271, soh_267, soh_268, soh_269, soh_270, \
                         soh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_18 * snh_267[k]
                   + f_3 * pc_x[k] * soh_267[k];

        t_352[k] = f_18 * snh_268[k]
                   + f_3 * pc_x[k] * soh_268[k];

        t_353[k] = f_18 * snh_269[k]
                   + f_3 * pc_x[k] * soh_269[k];

        t_354[k] = f_18 * snh_270[k]
                   + f_3 * pc_x[k] * soh_270[k];

        t_355[k] = f_18 * snh_271[k]
                   + f_3 * pc_x[k] * soh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, snh_162, snh_183, snh_272, \
                         sog0_190, sog1_190, soh_267, soh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_18 * snh_272[k]
                   + f_3 * pc_x[k] * soh_272[k];

        t_357[k] = f_12 * snh_183[k]
                   + f_1 * sog0_190[k]
                   - f_2 * sog1_190[k]
                   + f_3 * pc_y[k] * soh_267[k];

        t_358[k] = f_12 * snh_162[k]
                   + f_3 * pc_z[k] * soh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, snh_185, snh_186, snh_187, sog0_192, \
                         sog0_193, sog0_194, sog1_192, sog1_193, sog1_194, soh_269, soh_270, \
                         soh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * snh_185[k]
                   + f_4 * sog0_192[k]
                   - f_5 * sog1_192[k]
                   + f_3 * pc_y[k] * soh_269[k];

        t_360[k] = f_12 * snh_186[k]
                   + f_6 * sog0_193[k]
                   - f_7 * sog1_193[k]
                   + f_3 * pc_y[k] * soh_270[k];

        t_361[k] = f_12 * snh_187[k]
                   + f_8 * sog0_194[k]
                   - f_9 * sog1_194[k]
                   + f_3 * pc_y[k] * soh_271[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_y, pc_y, pc_z, sni0_252, snh_167, \
                         snh_188, snh_189, sni1_252, sog0_194, sog1_194, soh_272, \
                         soh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * snh_188[k]
                   + f_3 * pc_y[k] * soh_272[k];

        t_363[k] = f_12 * snh_167[k]
                   + f_1 * sog0_194[k]
                   - f_2 * sog1_194[k]
                   + f_3 * pc_z[k] * soh_272[k];

        t_364[k] = pb_y[k] * sni0_252[k]
                   - f_10 * pc_y[k] * sni1_252[k];

        t_365[k] = f_11 * snh_189[k]
                   + f_3 * pc_y[k] * soh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_y, pc_y, pc_z, sni0_255, sni0_257, \
                         snh_168, snh_190, snh_191, sni1_255, sni1_257, soh_273, \
                         soh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * snh_168[k]
                   + f_3 * pc_z[k] * soh_273[k];

        t_367[k] = pb_y[k] * sni0_255[k]
                   + f_12 * snh_190[k]
                   - f_10 * pc_y[k] * sni1_255[k];

        t_368[k] = f_11 * snh_191[k]
                   + f_3 * pc_y[k] * soh_275[k];

        t_369[k] = pb_y[k] * sni0_257[k]
                   - f_10 * pc_y[k] * sni1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_y, pc_y, pc_z, sni0_258, sni0_261, \
                         snh_171, snh_192, snh_194, sni1_258, sni1_261, soh_276, \
                         soh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_y[k] * sni0_258[k]
                   + f_13 * snh_192[k]
                   - f_10 * pc_y[k] * sni1_258[k];

        t_371[k] = f_13 * snh_171[k]
                   + f_3 * pc_z[k] * soh_276[k];

        t_372[k] = f_11 * snh_194[k]
                   + f_3 * pc_y[k] * soh_278[k];

        t_373[k] = pb_y[k] * sni0_261[k]
                   - f_10 * pc_y[k] * sni1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pb_y, pc_y, pc_z, sni0_262, sni0_264, snh_174, \
                         snh_195, snh_197, sni1_262, sni1_264, \
                         soh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_y[k] * sni0_262[k]
                   + f_14 * snh_195[k]
                   - f_10 * pc_y[k] * sni1_262[k];

        t_375[k] = f_13 * snh_174[k]
                   + f_3 * pc_z[k] * soh_279[k];

        t_376[k] = pb_y[k] * sni0_264[k]
                   + f_12 * snh_197[k]
                   - f_10 * pc_y[k] * sni1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pb_y, pc_x, pc_y, sni0_266, snh_198, \
                         snh_288, snh_289, sni1_266, soh_282, soh_288, \
                         soh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * snh_198[k]
                   + f_3 * pc_y[k] * soh_282[k];

        t_378[k] = pb_y[k] * sni0_266[k]
                   - f_10 * pc_y[k] * sni1_266[k];

        t_379[k] = f_18 * snh_288[k]
                   + f_3 * pc_x[k] * soh_288[k];

        t_380[k] = f_18 * snh_289[k]
                   + f_3 * pc_x[k] * soh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, snh_290, snh_291, snh_292, snh_293, \
                         soh_290, soh_291, soh_292, soh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_18 * snh_290[k]
                   + f_3 * pc_x[k] * soh_290[k];

        t_382[k] = f_18 * snh_291[k]
                   + f_3 * pc_x[k] * soh_291[k];

        t_383[k] = f_18 * snh_292[k]
                   + f_3 * pc_x[k] * soh_292[k];

        t_384[k] = f_18 * snh_293[k]
                   + f_3 * pc_x[k] * soh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, snh_183, snh_204, snh_206, sog0_205, \
                         sog0_207, sog1_205, sog1_207, soh_288, \
                         soh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * snh_204[k]
                   + f_1 * sog0_205[k]
                   - f_2 * sog1_205[k]
                   + f_3 * pc_y[k] * soh_288[k];

        t_386[k] = f_13 * snh_183[k]
                   + f_3 * pc_z[k] * soh_288[k];

        t_387[k] = f_11 * snh_206[k]
                   + f_4 * sog0_207[k]
                   - f_5 * sog1_207[k]
                   + f_3 * pc_y[k] * soh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, snh_207, snh_208, snh_209, sog0_208, \
                         sog0_209, sog1_208, sog1_209, soh_291, soh_292, \
                         soh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * snh_207[k]
                   + f_6 * sog0_208[k]
                   - f_7 * sog1_208[k]
                   + f_3 * pc_y[k] * soh_291[k];

        t_389[k] = f_11 * snh_208[k]
                   + f_8 * sog0_209[k]
                   - f_9 * sog1_209[k]
                   + f_3 * pc_y[k] * soh_292[k];

        t_390[k] = f_11 * snh_209[k]
                   + f_3 * pc_y[k] * soh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pb_y, pc_x, pc_y, pc_z, sni0_279, \
                         snh_189, snh_294, sni1_279, sog0_210, sog1_210, \
                         soh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pb_y[k] * sni0_279[k]
                   - f_10 * pc_y[k] * sni1_279[k];

        t_392[k] = f_18 * snh_294[k]
                   + f_1 * sog0_210[k]
                   - f_2 * sog1_210[k]
                   + f_3 * pc_x[k] * soh_294[k];

        t_393[k] = f_3 * pc_y[k] * soh_294[k];

        t_394[k] = f_14 * snh_189[k]
                   + f_3 * pc_z[k] * soh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, snh_297, snh_299, sog0_213, \
                         sog0_215, sog1_213, sog1_215, soh_296, soh_297, \
                         soh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_18 * snh_297[k]
                   + f_4 * sog0_213[k]
                   - f_5 * sog1_213[k]
                   + f_3 * pc_x[k] * soh_297[k];

        t_396[k] = f_3 * pc_y[k] * soh_296[k];

        t_397[k] = f_18 * snh_299[k]
                   + f_4 * sog0_215[k]
                   - f_5 * sog1_215[k]
                   + f_3 * pc_x[k] * soh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_x, pc_y, pc_z, snh_192, snh_300, sog0_216, \
                         sog1_216, soh_297, soh_299, soh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_18 * snh_300[k]
                   + f_6 * sog0_216[k]
                   - f_7 * sog1_216[k]
                   + f_3 * pc_x[k] * soh_300[k];

        t_399[k] = f_14 * snh_192[k]
                   + f_3 * pc_z[k] * soh_297[k];

        t_400[k] = f_3 * pc_y[k] * soh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, snh_195, snh_303, snh_304, sog0_219, \
                         sog0_220, sog1_219, sog1_220, soh_300, soh_303, \
                         soh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_18 * snh_303[k]
                   + f_6 * sog0_219[k]
                   - f_7 * sog1_219[k]
                   + f_3 * pc_x[k] * soh_303[k];

        t_402[k] = f_18 * snh_304[k]
                   + f_8 * sog0_220[k]
                   - f_9 * sog1_220[k]
                   + f_3 * pc_x[k] * soh_304[k];

        t_403[k] = f_14 * snh_195[k]
                   + f_3 * pc_z[k] * soh_300[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pc_x, pc_y, snh_306, snh_308, sog0_222, \
                         sog0_224, sog1_222, sog1_224, soh_303, soh_306, \
                         soh_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_18 * snh_306[k]
                   + f_8 * sog0_222[k]
                   - f_9 * sog1_222[k]
                   + f_3 * pc_x[k] * soh_306[k];

        t_405[k] = f_3 * pc_y[k] * soh_303[k];

        t_406[k] = f_18 * snh_308[k]
                   + f_8 * sog0_224[k]
                   - f_9 * sog1_224[k]
                   + f_3 * pc_x[k] * soh_308[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pc_x, snh_309, snh_310, snh_311, \
                         snh_312, snh_313, soh_309, soh_310, soh_311, soh_312, \
                         soh_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_18 * snh_309[k]
                   + f_3 * pc_x[k] * soh_309[k];

        t_408[k] = f_18 * snh_310[k]
                   + f_3 * pc_x[k] * soh_310[k];

        t_409[k] = f_18 * snh_311[k]
                   + f_3 * pc_x[k] * soh_311[k];

        t_410[k] = f_18 * snh_312[k]
                   + f_3 * pc_x[k] * soh_312[k];

        t_411[k] = f_18 * snh_313[k]
                   + f_3 * pc_x[k] * soh_313[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pc_x, pc_y, pc_z, snh_204, snh_314, \
                         sog0_220, sog0_222, sog1_220, sog1_222, soh_309, soh_311, \
                         soh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_18 * snh_314[k]
                   + f_3 * pc_x[k] * soh_314[k];

        t_413[k] = f_1 * sog0_220[k]
                   - f_2 * sog1_220[k]
                   + f_3 * pc_y[k] * soh_309[k];

        t_414[k] = f_14 * snh_204[k]
                   + f_3 * pc_z[k] * soh_309[k];

        t_415[k] = f_4 * sog0_222[k]
                   - f_5 * sog1_222[k]
                   + f_3 * pc_y[k] * soh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, snh_209, sog0_223, sog0_224, \
                         sog1_223, sog1_224, soh_312, soh_313, \
                         soh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * sog0_223[k]
                   - f_7 * sog1_223[k]
                   + f_3 * pc_y[k] * soh_312[k];

        t_417[k] = f_8 * sog0_224[k]
                   - f_9 * sog1_224[k]
                   + f_3 * pc_y[k] * soh_313[k];

        t_418[k] = f_3 * pc_y[k] * soh_314[k];

        t_419[k] = f_14 * snh_209[k]
                   + f_1 * sog0_224[k]
                   - f_2 * sog1_224[k]
                   + f_3 * pc_z[k] * soh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, snh_210, snh_315, \
                         snh_318, sog0_225, sog0_228, sog1_225, sog1_228, soh_315, \
                         soh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_19 * snh_315[k]
                   + f_1 * sog0_225[k]
                   - f_2 * sog1_225[k]
                   + f_3 * pc_x[k] * soh_315[k];

        t_421[k] = f_20 * snh_210[k]
                   + f_3 * pc_y[k] * soh_315[k];

        t_422[k] = f_3 * pc_z[k] * soh_315[k];

        t_423[k] = f_19 * snh_318[k]
                   + f_4 * sog0_228[k]
                   - f_5 * sog1_228[k]
                   + f_3 * pc_x[k] * soh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pc_x, pc_y, snh_212, snh_320, snh_321, sog0_230, \
                         sog0_231, sog1_230, sog1_231, soh_317, soh_320, \
                         soh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_20 * snh_212[k]
                   + f_3 * pc_y[k] * soh_317[k];

        t_425[k] = f_19 * snh_320[k]
                   + f_4 * sog0_230[k]
                   - f_5 * sog1_230[k]
                   + f_3 * pc_x[k] * soh_320[k];

        t_426[k] = f_19 * snh_321[k]
                   + f_6 * sog0_231[k]
                   - f_7 * sog1_231[k]
                   + f_3 * pc_x[k] * soh_321[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pc_x, pc_y, pc_z, snh_215, snh_324, sog0_234, \
                         sog1_234, soh_318, soh_320, soh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_z[k] * soh_318[k];

        t_428[k] = f_20 * snh_215[k]
                   + f_3 * pc_y[k] * soh_320[k];

        t_429[k] = f_19 * snh_324[k]
                   + f_6 * sog0_234[k]
                   - f_7 * sog1_234[k]
                   + f_3 * pc_x[k] * soh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pc_x, pc_z, snh_325, snh_327, sog0_235, \
                         sog0_237, sog1_235, sog1_237, soh_321, soh_325, \
                         soh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_19 * snh_325[k]
                   + f_8 * sog0_235[k]
                   - f_9 * sog1_235[k]
                   + f_3 * pc_x[k] * soh_325[k];

        t_431[k] = f_3 * pc_z[k] * soh_321[k];

        t_432[k] = f_19 * snh_327[k]
                   + f_8 * sog0_237[k]
                   - f_9 * sog1_237[k]
                   + f_3 * pc_x[k] * soh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, snh_219, snh_329, snh_330, \
                         snh_331, sog0_239, sog1_239, soh_324, soh_329, soh_330, \
                         soh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_20 * snh_219[k]
                   + f_3 * pc_y[k] * soh_324[k];

        t_434[k] = f_19 * snh_329[k]
                   + f_8 * sog0_239[k]
                   - f_9 * sog1_239[k]
                   + f_3 * pc_x[k] * soh_329[k];

        t_435[k] = f_19 * snh_330[k]
                   + f_3 * pc_x[k] * soh_330[k];

        t_436[k] = f_19 * snh_331[k]
                   + f_3 * pc_x[k] * soh_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, snh_332, snh_333, snh_334, snh_335, \
                         soh_332, soh_333, soh_334, soh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_19 * snh_332[k]
                   + f_3 * pc_x[k] * soh_332[k];

        t_438[k] = f_19 * snh_333[k]
                   + f_3 * pc_x[k] * soh_333[k];

        t_439[k] = f_19 * snh_334[k]
                   + f_3 * pc_x[k] * soh_334[k];

        t_440[k] = f_19 * snh_335[k]
                   + f_3 * pc_x[k] * soh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, pc_z, snh_225, snh_227, sog0_235, \
                         sog0_237, sog1_235, sog1_237, soh_330, \
                         soh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_20 * snh_225[k]
                   + f_1 * sog0_235[k]
                   - f_2 * sog1_235[k]
                   + f_3 * pc_y[k] * soh_330[k];

        t_442[k] = f_3 * pc_z[k] * soh_330[k];

        t_443[k] = f_20 * snh_227[k]
                   + f_4 * sog0_237[k]
                   - f_5 * sog1_237[k]
                   + f_3 * pc_y[k] * soh_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, snh_228, snh_229, snh_230, \
                         sog0_238, sog0_239, sog1_238, sog1_239, soh_333, soh_334, \
                         soh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_20 * snh_228[k]
                   + f_6 * sog0_238[k]
                   - f_7 * sog1_238[k]
                   + f_3 * pc_y[k] * soh_333[k];

        t_445[k] = f_20 * snh_229[k]
                   + f_8 * sog0_239[k]
                   - f_9 * sog1_239[k]
                   + f_3 * pc_y[k] * soh_334[k];

        t_446[k] = f_20 * snh_230[k]
                   + f_3 * pc_y[k] * soh_335[k];

        t_447[k] = f_1 * sog0_239[k]
                   - f_2 * sog1_239[k]
                   + f_3 * pc_z[k] * soh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pb_z, pc_y, pc_z, sni0_280, sni0_283, \
                         snh_210, snh_231, sni1_280, sni1_283, \
                         soh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_z[k] * sni0_280[k]
                   - f_10 * pc_z[k] * sni1_280[k];

        t_449[k] = f_14 * snh_231[k]
                   + f_3 * pc_y[k] * soh_336[k];

        t_450[k] = f_11 * snh_210[k]
                   + f_3 * pc_z[k] * soh_336[k];

        t_451[k] = pb_z[k] * sni0_283[k]
                   - f_10 * pc_z[k] * sni1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pb_z, pc_x, pc_y, pc_z, sni0_286, snh_233, \
                         snh_341, sni1_286, sog0_245, sog1_245, soh_338, \
                         soh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * snh_233[k]
                   + f_3 * pc_y[k] * soh_338[k];

        t_453[k] = f_19 * snh_341[k]
                   + f_4 * sog0_245[k]
                   - f_5 * sog1_245[k]
                   + f_3 * pc_x[k] * soh_341[k];

        t_454[k] = pb_z[k] * sni0_286[k]
                   - f_10 * pc_z[k] * sni1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, pc_y, pc_z, snh_213, snh_236, snh_345, \
                         sog0_249, sog1_249, soh_339, soh_341, \
                         soh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * snh_213[k]
                   + f_3 * pc_z[k] * soh_339[k];

        t_456[k] = f_14 * snh_236[k]
                   + f_3 * pc_y[k] * soh_341[k];

        t_457[k] = f_19 * snh_345[k]
                   + f_6 * sog0_249[k]
                   - f_7 * sog1_249[k]
                   + f_3 * pc_x[k] * soh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pb_z, pc_y, pc_z, sni0_290, sni0_292, \
                         snh_216, snh_217, snh_240, sni1_290, sni1_292, soh_342, \
                         soh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pb_z[k] * sni0_290[k]
                   - f_10 * pc_z[k] * sni1_290[k];

        t_459[k] = f_11 * snh_216[k]
                   + f_3 * pc_z[k] * soh_342[k];

        t_460[k] = pb_z[k] * sni0_292[k]
                   + f_12 * snh_217[k]
                   - f_10 * pc_z[k] * sni1_292[k];

        t_461[k] = f_14 * snh_240[k]
                   + f_3 * pc_y[k] * soh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, snh_350, snh_351, snh_352, snh_353, \
                         sog0_254, sog1_254, soh_350, soh_351, soh_352, \
                         soh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_19 * snh_350[k]
                   + f_8 * sog0_254[k]
                   - f_9 * sog1_254[k]
                   + f_3 * pc_x[k] * soh_350[k];

        t_463[k] = f_19 * snh_351[k]
                   + f_3 * pc_x[k] * soh_351[k];

        t_464[k] = f_19 * snh_352[k]
                   + f_3 * pc_x[k] * soh_352[k];

        t_465[k] = f_19 * snh_353[k]
                   + f_3 * pc_x[k] * soh_353[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_301 = buffer.data(sni0 + 301);
    const auto *sni0_392 = buffer.data(sni0 + 392);
    const auto *sni0_395 = buffer.data(sni0 + 395);
    const auto *sni0_397 = buffer.data(sni0 + 397);
    const auto *sni0_398 = buffer.data(sni0 + 398);
    const auto *sni0_401 = buffer.data(sni0 + 401);
    const auto *sni0_402 = buffer.data(sni0 + 402);
    const auto *sni0_404 = buffer.data(sni0 + 404);
    const auto *sni0_406 = buffer.data(sni0 + 406);
    const auto *sni0_419 = buffer.data(sni0 + 419);

    const auto *snh_225 = buffer.data(snh + 225);
    const auto *snh_230 = buffer.data(snh + 230);
    const auto *snh_231 = buffer.data(snh + 231);
    const auto *snh_234 = buffer.data(snh + 234);
    const auto *snh_237 = buffer.data(snh + 237);
    const auto *snh_246 = buffer.data(snh + 246);
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
    const auto *snh_267 = buffer.data(snh + 267);
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
    const auto *snh_290 = buffer.data(snh + 290);
    const auto *snh_291 = buffer.data(snh + 291);
    const auto *snh_292 = buffer.data(snh + 292);
    const auto *snh_293 = buffer.data(snh + 293);
    const auto *snh_294 = buffer.data(snh + 294);
    const auto *snh_295 = buffer.data(snh + 295);
    const auto *snh_296 = buffer.data(snh + 296);
    const auto *snh_297 = buffer.data(snh + 297);
    const auto *snh_299 = buffer.data(snh + 299);
    const auto *snh_300 = buffer.data(snh + 300);
    const auto *snh_302 = buffer.data(snh + 302);
    const auto *snh_303 = buffer.data(snh + 303);
    const auto *snh_309 = buffer.data(snh + 309);
    const auto *snh_311 = buffer.data(snh + 311);
    const auto *snh_312 = buffer.data(snh + 312);
    const auto *snh_313 = buffer.data(snh + 313);
    const auto *snh_314 = buffer.data(snh + 314);
    const auto *snh_354 = buffer.data(snh + 354);
    const auto *snh_355 = buffer.data(snh + 355);
    const auto *snh_356 = buffer.data(snh + 356);
    const auto *snh_357 = buffer.data(snh + 357);
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
    const auto *snh_414 = buffer.data(snh + 414);
    const auto *snh_415 = buffer.data(snh + 415);
    const auto *snh_416 = buffer.data(snh + 416);
    const auto *snh_417 = buffer.data(snh + 417);
    const auto *snh_418 = buffer.data(snh + 418);
    const auto *snh_419 = buffer.data(snh + 419);
    const auto *snh_420 = buffer.data(snh + 420);
    const auto *snh_423 = buffer.data(snh + 423);
    const auto *snh_425 = buffer.data(snh + 425);
    const auto *snh_426 = buffer.data(snh + 426);
    const auto *snh_429 = buffer.data(snh + 429);
    const auto *snh_430 = buffer.data(snh + 430);
    const auto *snh_432 = buffer.data(snh + 432);
    const auto *snh_434 = buffer.data(snh + 434);

    const auto *sni1_301 = buffer.data(sni1 + 301);
    const auto *sni1_392 = buffer.data(sni1 + 392);
    const auto *sni1_395 = buffer.data(sni1 + 395);
    const auto *sni1_397 = buffer.data(sni1 + 397);
    const auto *sni1_398 = buffer.data(sni1 + 398);
    const auto *sni1_401 = buffer.data(sni1 + 401);
    const auto *sni1_402 = buffer.data(sni1 + 402);
    const auto *sni1_404 = buffer.data(sni1 + 404);
    const auto *sni1_406 = buffer.data(sni1 + 406);
    const auto *sni1_419 = buffer.data(sni1 + 419);

    const auto *sog0_252 = buffer.data(sog0 + 252);
    const auto *sog0_253 = buffer.data(sog0 + 253);
    const auto *sog0_254 = buffer.data(sog0 + 254);
    const auto *sog0_255 = buffer.data(sog0 + 255);
    const auto *sog0_258 = buffer.data(sog0 + 258);
    const auto *sog0_260 = buffer.data(sog0 + 260);
    const auto *sog0_261 = buffer.data(sog0 + 261);
    const auto *sog0_264 = buffer.data(sog0 + 264);
    const auto *sog0_265 = buffer.data(sog0 + 265);
    const auto *sog0_267 = buffer.data(sog0 + 267);
    const auto *sog0_268 = buffer.data(sog0 + 268);
    const auto *sog0_269 = buffer.data(sog0 + 269);
    const auto *sog0_270 = buffer.data(sog0 + 270);
    const auto *sog0_273 = buffer.data(sog0 + 273);
    const auto *sog0_275 = buffer.data(sog0 + 275);
    const auto *sog0_276 = buffer.data(sog0 + 276);
    const auto *sog0_279 = buffer.data(sog0 + 279);
    const auto *sog0_280 = buffer.data(sog0 + 280);
    const auto *sog0_282 = buffer.data(sog0 + 282);
    const auto *sog0_283 = buffer.data(sog0 + 283);
    const auto *sog0_284 = buffer.data(sog0 + 284);
    const auto *sog0_295 = buffer.data(sog0 + 295);
    const auto *sog0_297 = buffer.data(sog0 + 297);
    const auto *sog0_298 = buffer.data(sog0 + 298);
    const auto *sog0_299 = buffer.data(sog0 + 299);
    const auto *sog0_300 = buffer.data(sog0 + 300);
    const auto *sog0_303 = buffer.data(sog0 + 303);
    const auto *sog0_305 = buffer.data(sog0 + 305);
    const auto *sog0_306 = buffer.data(sog0 + 306);
    const auto *sog0_309 = buffer.data(sog0 + 309);
    const auto *sog0_310 = buffer.data(sog0 + 310);
    const auto *sog0_312 = buffer.data(sog0 + 312);
    const auto *sog0_314 = buffer.data(sog0 + 314);

    const auto *sog1_252 = buffer.data(sog1 + 252);
    const auto *sog1_253 = buffer.data(sog1 + 253);
    const auto *sog1_254 = buffer.data(sog1 + 254);
    const auto *sog1_255 = buffer.data(sog1 + 255);
    const auto *sog1_258 = buffer.data(sog1 + 258);
    const auto *sog1_260 = buffer.data(sog1 + 260);
    const auto *sog1_261 = buffer.data(sog1 + 261);
    const auto *sog1_264 = buffer.data(sog1 + 264);
    const auto *sog1_265 = buffer.data(sog1 + 265);
    const auto *sog1_267 = buffer.data(sog1 + 267);
    const auto *sog1_268 = buffer.data(sog1 + 268);
    const auto *sog1_269 = buffer.data(sog1 + 269);
    const auto *sog1_270 = buffer.data(sog1 + 270);
    const auto *sog1_273 = buffer.data(sog1 + 273);
    const auto *sog1_275 = buffer.data(sog1 + 275);
    const auto *sog1_276 = buffer.data(sog1 + 276);
    const auto *sog1_279 = buffer.data(sog1 + 279);
    const auto *sog1_280 = buffer.data(sog1 + 280);
    const auto *sog1_282 = buffer.data(sog1 + 282);
    const auto *sog1_283 = buffer.data(sog1 + 283);
    const auto *sog1_284 = buffer.data(sog1 + 284);
    const auto *sog1_295 = buffer.data(sog1 + 295);
    const auto *sog1_297 = buffer.data(sog1 + 297);
    const auto *sog1_298 = buffer.data(sog1 + 298);
    const auto *sog1_299 = buffer.data(sog1 + 299);
    const auto *sog1_300 = buffer.data(sog1 + 300);
    const auto *sog1_303 = buffer.data(sog1 + 303);
    const auto *sog1_305 = buffer.data(sog1 + 305);
    const auto *sog1_306 = buffer.data(sog1 + 306);
    const auto *sog1_309 = buffer.data(sog1 + 309);
    const auto *sog1_310 = buffer.data(sog1 + 310);
    const auto *sog1_312 = buffer.data(sog1 + 312);
    const auto *sog1_314 = buffer.data(sog1 + 314);

    const auto *soh_351 = buffer.data(soh + 351);
    const auto *soh_353 = buffer.data(soh + 353);
    const auto *soh_354 = buffer.data(soh + 354);
    const auto *soh_355 = buffer.data(soh + 355);
    const auto *soh_356 = buffer.data(soh + 356);
    const auto *soh_357 = buffer.data(soh + 357);
    const auto *soh_359 = buffer.data(soh + 359);
    const auto *soh_360 = buffer.data(soh + 360);
    const auto *soh_362 = buffer.data(soh + 362);
    const auto *soh_363 = buffer.data(soh + 363);
    const auto *soh_366 = buffer.data(soh + 366);
    const auto *soh_367 = buffer.data(soh + 367);
    const auto *soh_369 = buffer.data(soh + 369);
    const auto *soh_371 = buffer.data(soh + 371);
    const auto *soh_372 = buffer.data(soh + 372);
    const auto *soh_373 = buffer.data(soh + 373);
    const auto *soh_374 = buffer.data(soh + 374);
    const auto *soh_375 = buffer.data(soh + 375);
    const auto *soh_376 = buffer.data(soh + 376);
    const auto *soh_377 = buffer.data(soh + 377);
    const auto *soh_378 = buffer.data(soh + 378);
    const auto *soh_380 = buffer.data(soh + 380);
    const auto *soh_381 = buffer.data(soh + 381);
    const auto *soh_383 = buffer.data(soh + 383);
    const auto *soh_384 = buffer.data(soh + 384);
    const auto *soh_387 = buffer.data(soh + 387);
    const auto *soh_388 = buffer.data(soh + 388);
    const auto *soh_390 = buffer.data(soh + 390);
    const auto *soh_392 = buffer.data(soh + 392);
    const auto *soh_393 = buffer.data(soh + 393);
    const auto *soh_394 = buffer.data(soh + 394);
    const auto *soh_395 = buffer.data(soh + 395);
    const auto *soh_396 = buffer.data(soh + 396);
    const auto *soh_397 = buffer.data(soh + 397);
    const auto *soh_398 = buffer.data(soh + 398);
    const auto *soh_399 = buffer.data(soh + 399);
    const auto *soh_401 = buffer.data(soh + 401);
    const auto *soh_402 = buffer.data(soh + 402);
    const auto *soh_404 = buffer.data(soh + 404);
    const auto *soh_405 = buffer.data(soh + 405);
    const auto *soh_408 = buffer.data(soh + 408);
    const auto *soh_414 = buffer.data(soh + 414);
    const auto *soh_415 = buffer.data(soh + 415);
    const auto *soh_416 = buffer.data(soh + 416);
    const auto *soh_417 = buffer.data(soh + 417);
    const auto *soh_418 = buffer.data(soh + 418);
    const auto *soh_419 = buffer.data(soh + 419);
    const auto *soh_420 = buffer.data(soh + 420);
    const auto *soh_422 = buffer.data(soh + 422);
    const auto *soh_423 = buffer.data(soh + 423);
    const auto *soh_425 = buffer.data(soh + 425);
    const auto *soh_426 = buffer.data(soh + 426);
    const auto *soh_429 = buffer.data(soh + 429);
    const auto *soh_430 = buffer.data(soh + 430);
    const auto *soh_432 = buffer.data(soh + 432);
    const auto *soh_434 = buffer.data(soh + 434);

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pb_z, pc_x, pc_z, sni0_301, snh_354, \
                         snh_355, snh_356, sni1_301, soh_354, soh_355, \
                         soh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_19 * snh_354[k]
                   + f_3 * pc_x[k] * soh_354[k];

        t_467[k] = f_19 * snh_355[k]
                   + f_3 * pc_x[k] * soh_355[k];

        t_468[k] = f_19 * snh_356[k]
                   + f_3 * pc_x[k] * soh_356[k];

        t_469[k] = pb_z[k] * sni0_301[k]
                   - f_10 * pc_z[k] * sni1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, snh_225, snh_248, snh_249, sog0_252, \
                         sog0_253, sog1_252, sog1_253, soh_351, soh_353, \
                         soh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * snh_225[k]
                   + f_3 * pc_z[k] * soh_351[k];

        t_471[k] = f_14 * snh_248[k]
                   + f_4 * sog0_252[k]
                   - f_5 * sog1_252[k]
                   + f_3 * pc_y[k] * soh_353[k];

        t_472[k] = f_14 * snh_249[k]
                   + f_6 * sog0_253[k]
                   - f_7 * sog1_253[k]
                   + f_3 * pc_y[k] * soh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, snh_230, snh_250, snh_251, sog0_254, \
                         sog1_254, soh_355, soh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * snh_250[k]
                   + f_8 * sog0_254[k]
                   - f_9 * sog1_254[k]
                   + f_3 * pc_y[k] * soh_355[k];

        t_474[k] = f_14 * snh_251[k]
                   + f_3 * pc_y[k] * soh_356[k];

        t_475[k] = f_11 * snh_230[k]
                   + f_1 * sog0_254[k]
                   - f_2 * sog1_254[k]
                   + f_3 * pc_z[k] * soh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, snh_231, snh_252, snh_357, \
                         sog0_255, sog1_255, soh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_19 * snh_357[k]
                   + f_1 * sog0_255[k]
                   - f_2 * sog1_255[k]
                   + f_3 * pc_x[k] * soh_357[k];

        t_477[k] = f_13 * snh_252[k]
                   + f_3 * pc_y[k] * soh_357[k];

        t_478[k] = f_12 * snh_231[k]
                   + f_3 * pc_z[k] * soh_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, snh_254, snh_360, snh_362, sog0_258, \
                         sog0_260, sog1_258, sog1_260, soh_359, soh_360, \
                         soh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_19 * snh_360[k]
                   + f_4 * sog0_258[k]
                   - f_5 * sog1_258[k]
                   + f_3 * pc_x[k] * soh_360[k];

        t_480[k] = f_13 * snh_254[k]
                   + f_3 * pc_y[k] * soh_359[k];

        t_481[k] = f_19 * snh_362[k]
                   + f_4 * sog0_260[k]
                   - f_5 * sog1_260[k]
                   + f_3 * pc_x[k] * soh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, snh_234, snh_257, snh_363, \
                         sog0_261, sog1_261, soh_360, soh_362, \
                         soh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_19 * snh_363[k]
                   + f_6 * sog0_261[k]
                   - f_7 * sog1_261[k]
                   + f_3 * pc_x[k] * soh_363[k];

        t_483[k] = f_12 * snh_234[k]
                   + f_3 * pc_z[k] * soh_360[k];

        t_484[k] = f_13 * snh_257[k]
                   + f_3 * pc_y[k] * soh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, snh_237, snh_366, snh_367, sog0_264, \
                         sog0_265, sog1_264, sog1_265, soh_363, soh_366, \
                         soh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_19 * snh_366[k]
                   + f_6 * sog0_264[k]
                   - f_7 * sog1_264[k]
                   + f_3 * pc_x[k] * soh_366[k];

        t_486[k] = f_19 * snh_367[k]
                   + f_8 * sog0_265[k]
                   - f_9 * sog1_265[k]
                   + f_3 * pc_x[k] * soh_367[k];

        t_487[k] = f_12 * snh_237[k]
                   + f_3 * pc_z[k] * soh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, snh_261, snh_369, snh_371, sog0_267, \
                         sog0_269, sog1_267, sog1_269, soh_366, soh_369, \
                         soh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_19 * snh_369[k]
                   + f_8 * sog0_267[k]
                   - f_9 * sog1_267[k]
                   + f_3 * pc_x[k] * soh_369[k];

        t_489[k] = f_13 * snh_261[k]
                   + f_3 * pc_y[k] * soh_366[k];

        t_490[k] = f_19 * snh_371[k]
                   + f_8 * sog0_269[k]
                   - f_9 * sog1_269[k]
                   + f_3 * pc_x[k] * soh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, snh_372, snh_373, snh_374, \
                         snh_375, snh_376, soh_372, soh_373, soh_374, soh_375, \
                         soh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_19 * snh_372[k]
                   + f_3 * pc_x[k] * soh_372[k];

        t_492[k] = f_19 * snh_373[k]
                   + f_3 * pc_x[k] * soh_373[k];

        t_493[k] = f_19 * snh_374[k]
                   + f_3 * pc_x[k] * soh_374[k];

        t_494[k] = f_19 * snh_375[k]
                   + f_3 * pc_x[k] * soh_375[k];

        t_495[k] = f_19 * snh_376[k]
                   + f_3 * pc_x[k] * soh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, snh_246, snh_267, snh_377, \
                         sog0_265, sog1_265, soh_372, soh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_19 * snh_377[k]
                   + f_3 * pc_x[k] * soh_377[k];

        t_497[k] = f_13 * snh_267[k]
                   + f_1 * sog0_265[k]
                   - f_2 * sog1_265[k]
                   + f_3 * pc_y[k] * soh_372[k];

        t_498[k] = f_12 * snh_246[k]
                   + f_3 * pc_z[k] * soh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, snh_269, snh_270, snh_271, sog0_267, \
                         sog0_268, sog0_269, sog1_267, sog1_268, sog1_269, soh_374, soh_375, \
                         soh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * snh_269[k]
                   + f_4 * sog0_267[k]
                   - f_5 * sog1_267[k]
                   + f_3 * pc_y[k] * soh_374[k];

        t_500[k] = f_13 * snh_270[k]
                   + f_6 * sog0_268[k]
                   - f_7 * sog1_268[k]
                   + f_3 * pc_y[k] * soh_375[k];

        t_501[k] = f_13 * snh_271[k]
                   + f_8 * sog0_269[k]
                   - f_9 * sog1_269[k]
                   + f_3 * pc_y[k] * soh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, snh_251, snh_272, snh_378, \
                         sog0_269, sog0_270, sog1_269, sog1_270, soh_377, \
                         soh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * snh_272[k]
                   + f_3 * pc_y[k] * soh_377[k];

        t_503[k] = f_12 * snh_251[k]
                   + f_1 * sog0_269[k]
                   - f_2 * sog1_269[k]
                   + f_3 * pc_z[k] * soh_377[k];

        t_504[k] = f_19 * snh_378[k]
                   + f_1 * sog0_270[k]
                   - f_2 * sog1_270[k]
                   + f_3 * pc_x[k] * soh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, snh_252, snh_273, \
                         snh_275, snh_381, sog0_273, sog1_273, soh_378, soh_380, \
                         soh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * snh_273[k]
                   + f_3 * pc_y[k] * soh_378[k];

        t_506[k] = f_13 * snh_252[k]
                   + f_3 * pc_z[k] * soh_378[k];

        t_507[k] = f_19 * snh_381[k]
                   + f_4 * sog0_273[k]
                   - f_5 * sog1_273[k]
                   + f_3 * pc_x[k] * soh_381[k];

        t_508[k] = f_12 * snh_275[k]
                   + f_3 * pc_y[k] * soh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, snh_255, snh_383, snh_384, sog0_275, \
                         sog0_276, sog1_275, sog1_276, soh_381, soh_383, \
                         soh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_19 * snh_383[k]
                   + f_4 * sog0_275[k]
                   - f_5 * sog1_275[k]
                   + f_3 * pc_x[k] * soh_383[k];

        t_510[k] = f_19 * snh_384[k]
                   + f_6 * sog0_276[k]
                   - f_7 * sog1_276[k]
                   + f_3 * pc_x[k] * soh_384[k];

        t_511[k] = f_13 * snh_255[k]
                   + f_3 * pc_z[k] * soh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, snh_278, snh_387, snh_388, sog0_279, \
                         sog0_280, sog1_279, sog1_280, soh_383, soh_387, \
                         soh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * snh_278[k]
                   + f_3 * pc_y[k] * soh_383[k];

        t_513[k] = f_19 * snh_387[k]
                   + f_6 * sog0_279[k]
                   - f_7 * sog1_279[k]
                   + f_3 * pc_x[k] * soh_387[k];

        t_514[k] = f_19 * snh_388[k]
                   + f_8 * sog0_280[k]
                   - f_9 * sog1_280[k]
                   + f_3 * pc_x[k] * soh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, snh_258, snh_282, snh_390, \
                         sog0_282, sog1_282, soh_384, soh_387, \
                         soh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * snh_258[k]
                   + f_3 * pc_z[k] * soh_384[k];

        t_516[k] = f_19 * snh_390[k]
                   + f_8 * sog0_282[k]
                   - f_9 * sog1_282[k]
                   + f_3 * pc_x[k] * soh_390[k];

        t_517[k] = f_12 * snh_282[k]
                   + f_3 * pc_y[k] * soh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, snh_392, snh_393, snh_394, snh_395, \
                         sog0_284, sog1_284, soh_392, soh_393, soh_394, \
                         soh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_19 * snh_392[k]
                   + f_8 * sog0_284[k]
                   - f_9 * sog1_284[k]
                   + f_3 * pc_x[k] * soh_392[k];

        t_519[k] = f_19 * snh_393[k]
                   + f_3 * pc_x[k] * soh_393[k];

        t_520[k] = f_19 * snh_394[k]
                   + f_3 * pc_x[k] * soh_394[k];

        t_521[k] = f_19 * snh_395[k]
                   + f_3 * pc_x[k] * soh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, snh_288, snh_396, snh_397, \
                         snh_398, sog0_280, sog1_280, soh_393, soh_396, soh_397, \
                         soh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_19 * snh_396[k]
                   + f_3 * pc_x[k] * soh_396[k];

        t_523[k] = f_19 * snh_397[k]
                   + f_3 * pc_x[k] * soh_397[k];

        t_524[k] = f_19 * snh_398[k]
                   + f_3 * pc_x[k] * soh_398[k];

        t_525[k] = f_12 * snh_288[k]
                   + f_1 * sog0_280[k]
                   - f_2 * sog1_280[k]
                   + f_3 * pc_y[k] * soh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, snh_267, snh_290, snh_291, sog0_282, \
                         sog0_283, sog1_282, sog1_283, soh_393, soh_395, \
                         soh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * snh_267[k]
                   + f_3 * pc_z[k] * soh_393[k];

        t_527[k] = f_12 * snh_290[k]
                   + f_4 * sog0_282[k]
                   - f_5 * sog1_282[k]
                   + f_3 * pc_y[k] * soh_395[k];

        t_528[k] = f_12 * snh_291[k]
                   + f_6 * sog0_283[k]
                   - f_7 * sog1_283[k]
                   + f_3 * pc_y[k] * soh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, sni0_392, snh_272, \
                         snh_292, snh_293, sni1_392, sog0_284, sog1_284, soh_397, \
                         soh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * snh_292[k]
                   + f_8 * sog0_284[k]
                   - f_9 * sog1_284[k]
                   + f_3 * pc_y[k] * soh_397[k];

        t_530[k] = f_12 * snh_293[k]
                   + f_3 * pc_y[k] * soh_398[k];

        t_531[k] = f_13 * snh_272[k]
                   + f_1 * sog0_284[k]
                   - f_2 * sog1_284[k]
                   + f_3 * pc_z[k] * soh_398[k];

        t_532[k] = pb_y[k] * sni0_392[k]
                   - f_10 * pc_y[k] * sni1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pb_y, pc_y, pc_z, sni0_395, snh_273, \
                         snh_294, snh_295, snh_296, sni1_395, soh_399, \
                         soh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * snh_294[k]
                   + f_3 * pc_y[k] * soh_399[k];

        t_534[k] = f_14 * snh_273[k]
                   + f_3 * pc_z[k] * soh_399[k];

        t_535[k] = pb_y[k] * sni0_395[k]
                   + f_12 * snh_295[k]
                   - f_10 * pc_y[k] * sni1_395[k];

        t_536[k] = f_11 * snh_296[k]
                   + f_3 * pc_y[k] * soh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pb_y, pc_y, pc_z, sni0_397, sni0_398, \
                         snh_276, snh_297, snh_299, sni1_397, sni1_398, soh_402, \
                         soh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pb_y[k] * sni0_397[k]
                   - f_10 * pc_y[k] * sni1_397[k];

        t_538[k] = pb_y[k] * sni0_398[k]
                   + f_13 * snh_297[k]
                   - f_10 * pc_y[k] * sni1_398[k];

        t_539[k] = f_14 * snh_276[k]
                   + f_3 * pc_z[k] * soh_402[k];

        t_540[k] = f_11 * snh_299[k]
                   + f_3 * pc_y[k] * soh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pb_y, pc_y, pc_z, sni0_401, sni0_402, snh_279, \
                         snh_300, sni1_401, sni1_402, soh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pb_y[k] * sni0_401[k]
                   - f_10 * pc_y[k] * sni1_401[k];

        t_542[k] = pb_y[k] * sni0_402[k]
                   + f_14 * snh_300[k]
                   - f_10 * pc_y[k] * sni1_402[k];

        t_543[k] = f_14 * snh_279[k]
                   + f_3 * pc_z[k] * soh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_y, pc_x, pc_y, sni0_404, sni0_406, \
                         snh_302, snh_303, snh_414, sni1_404, sni1_406, soh_408, \
                         soh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pb_y[k] * sni0_404[k]
                   + f_12 * snh_302[k]
                   - f_10 * pc_y[k] * sni1_404[k];

        t_545[k] = f_11 * snh_303[k]
                   + f_3 * pc_y[k] * soh_408[k];

        t_546[k] = pb_y[k] * sni0_406[k]
                   - f_10 * pc_y[k] * sni1_406[k];

        t_547[k] = f_19 * snh_414[k]
                   + f_3 * pc_x[k] * soh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, snh_415, snh_416, snh_417, \
                         snh_418, snh_419, soh_415, soh_416, soh_417, soh_418, \
                         soh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_19 * snh_415[k]
                   + f_3 * pc_x[k] * soh_415[k];

        t_549[k] = f_19 * snh_416[k]
                   + f_3 * pc_x[k] * soh_416[k];

        t_550[k] = f_19 * snh_417[k]
                   + f_3 * pc_x[k] * soh_417[k];

        t_551[k] = f_19 * snh_418[k]
                   + f_3 * pc_x[k] * soh_418[k];

        t_552[k] = f_19 * snh_419[k]
                   + f_3 * pc_x[k] * soh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, snh_288, snh_309, snh_311, sog0_295, \
                         sog0_297, sog1_295, sog1_297, soh_414, \
                         soh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * snh_309[k]
                   + f_1 * sog0_295[k]
                   - f_2 * sog1_295[k]
                   + f_3 * pc_y[k] * soh_414[k];

        t_554[k] = f_14 * snh_288[k]
                   + f_3 * pc_z[k] * soh_414[k];

        t_555[k] = f_11 * snh_311[k]
                   + f_4 * sog0_297[k]
                   - f_5 * sog1_297[k]
                   + f_3 * pc_y[k] * soh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, snh_312, snh_313, snh_314, sog0_298, \
                         sog0_299, sog1_298, sog1_299, soh_417, soh_418, \
                         soh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * snh_312[k]
                   + f_6 * sog0_298[k]
                   - f_7 * sog1_298[k]
                   + f_3 * pc_y[k] * soh_417[k];

        t_557[k] = f_11 * snh_313[k]
                   + f_8 * sog0_299[k]
                   - f_9 * sog1_299[k]
                   + f_3 * pc_y[k] * soh_418[k];

        t_558[k] = f_11 * snh_314[k]
                   + f_3 * pc_y[k] * soh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_y, pc_x, pc_y, pc_z, sni0_419, \
                         snh_294, snh_420, sni1_419, sog0_300, sog1_300, \
                         soh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pb_y[k] * sni0_419[k]
                   - f_10 * pc_y[k] * sni1_419[k];

        t_560[k] = f_19 * snh_420[k]
                   + f_1 * sog0_300[k]
                   - f_2 * sog1_300[k]
                   + f_3 * pc_x[k] * soh_420[k];

        t_561[k] = f_3 * pc_y[k] * soh_420[k];

        t_562[k] = f_20 * snh_294[k]
                   + f_3 * pc_z[k] * soh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, snh_423, snh_425, sog0_303, \
                         sog0_305, sog1_303, sog1_305, soh_422, soh_423, \
                         soh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_19 * snh_423[k]
                   + f_4 * sog0_303[k]
                   - f_5 * sog1_303[k]
                   + f_3 * pc_x[k] * soh_423[k];

        t_564[k] = f_3 * pc_y[k] * soh_422[k];

        t_565[k] = f_19 * snh_425[k]
                   + f_4 * sog0_305[k]
                   - f_5 * sog1_305[k]
                   + f_3 * pc_x[k] * soh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_x, pc_y, pc_z, snh_297, snh_426, sog0_306, \
                         sog1_306, soh_423, soh_425, soh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_19 * snh_426[k]
                   + f_6 * sog0_306[k]
                   - f_7 * sog1_306[k]
                   + f_3 * pc_x[k] * soh_426[k];

        t_567[k] = f_20 * snh_297[k]
                   + f_3 * pc_z[k] * soh_423[k];

        t_568[k] = f_3 * pc_y[k] * soh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, snh_300, snh_429, snh_430, sog0_309, \
                         sog0_310, sog1_309, sog1_310, soh_426, soh_429, \
                         soh_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_19 * snh_429[k]
                   + f_6 * sog0_309[k]
                   - f_7 * sog1_309[k]
                   + f_3 * pc_x[k] * soh_429[k];

        t_570[k] = f_19 * snh_430[k]
                   + f_8 * sog0_310[k]
                   - f_9 * sog1_310[k]
                   + f_3 * pc_x[k] * soh_430[k];

        t_571[k] = f_20 * snh_300[k]
                   + f_3 * pc_z[k] * soh_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_y, snh_432, snh_434, sog0_312, \
                         sog0_314, sog1_312, sog1_314, soh_429, soh_432, \
                         soh_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_19 * snh_432[k]
                   + f_8 * sog0_312[k]
                   - f_9 * sog1_312[k]
                   + f_3 * pc_x[k] * soh_432[k];

        t_573[k] = f_3 * pc_y[k] * soh_429[k];

        t_574[k] = f_19 * snh_434[k]
                   + f_8 * sog0_314[k]
                   - f_9 * sog1_314[k]
                   + f_3 * pc_x[k] * soh_434[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_420 = buffer.data(sni0 + 420);
    const auto *sni0_423 = buffer.data(sni0 + 423);
    const auto *sni0_426 = buffer.data(sni0 + 426);
    const auto *sni0_430 = buffer.data(sni0 + 430);
    const auto *sni0_432 = buffer.data(sni0 + 432);
    const auto *sni0_441 = buffer.data(sni0 + 441);

    const auto *snh_309 = buffer.data(snh + 309);
    const auto *snh_314 = buffer.data(snh + 314);
    const auto *snh_315 = buffer.data(snh + 315);
    const auto *snh_317 = buffer.data(snh + 317);
    const auto *snh_318 = buffer.data(snh + 318);
    const auto *snh_320 = buffer.data(snh + 320);
    const auto *snh_321 = buffer.data(snh + 321);
    const auto *snh_322 = buffer.data(snh + 322);
    const auto *snh_324 = buffer.data(snh + 324);
    const auto *snh_330 = buffer.data(snh + 330);
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
    const auto *snh_372 = buffer.data(snh + 372);
    const auto *snh_374 = buffer.data(snh + 374);
    const auto *snh_375 = buffer.data(snh + 375);
    const auto *snh_376 = buffer.data(snh + 376);
    const auto *snh_377 = buffer.data(snh + 377);
    const auto *snh_378 = buffer.data(snh + 378);
    const auto *snh_380 = buffer.data(snh + 380);
    const auto *snh_383 = buffer.data(snh + 383);
    const auto *snh_387 = buffer.data(snh + 387);
    const auto *snh_435 = buffer.data(snh + 435);
    const auto *snh_436 = buffer.data(snh + 436);
    const auto *snh_437 = buffer.data(snh + 437);
    const auto *snh_438 = buffer.data(snh + 438);
    const auto *snh_439 = buffer.data(snh + 439);
    const auto *snh_440 = buffer.data(snh + 440);
    const auto *snh_441 = buffer.data(snh + 441);
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
    const auto *snh_467 = buffer.data(snh + 467);
    const auto *snh_471 = buffer.data(snh + 471);
    const auto *snh_476 = buffer.data(snh + 476);
    const auto *snh_477 = buffer.data(snh + 477);
    const auto *snh_478 = buffer.data(snh + 478);
    const auto *snh_479 = buffer.data(snh + 479);
    const auto *snh_480 = buffer.data(snh + 480);
    const auto *snh_481 = buffer.data(snh + 481);
    const auto *snh_482 = buffer.data(snh + 482);
    const auto *snh_483 = buffer.data(snh + 483);
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
    const auto *snh_507 = buffer.data(snh + 507);
    const auto *snh_509 = buffer.data(snh + 509);
    const auto *snh_510 = buffer.data(snh + 510);
    const auto *snh_513 = buffer.data(snh + 513);
    const auto *snh_514 = buffer.data(snh + 514);
    const auto *snh_516 = buffer.data(snh + 516);

    const auto *sni1_420 = buffer.data(sni1 + 420);
    const auto *sni1_423 = buffer.data(sni1 + 423);
    const auto *sni1_426 = buffer.data(sni1 + 426);
    const auto *sni1_430 = buffer.data(sni1 + 430);
    const auto *sni1_432 = buffer.data(sni1 + 432);
    const auto *sni1_441 = buffer.data(sni1 + 441);

    const auto *sog0_310 = buffer.data(sog0 + 310);
    const auto *sog0_312 = buffer.data(sog0 + 312);
    const auto *sog0_313 = buffer.data(sog0 + 313);
    const auto *sog0_314 = buffer.data(sog0 + 314);
    const auto *sog0_315 = buffer.data(sog0 + 315);
    const auto *sog0_318 = buffer.data(sog0 + 318);
    const auto *sog0_320 = buffer.data(sog0 + 320);
    const auto *sog0_321 = buffer.data(sog0 + 321);
    const auto *sog0_324 = buffer.data(sog0 + 324);
    const auto *sog0_325 = buffer.data(sog0 + 325);
    const auto *sog0_327 = buffer.data(sog0 + 327);
    const auto *sog0_328 = buffer.data(sog0 + 328);
    const auto *sog0_329 = buffer.data(sog0 + 329);
    const auto *sog0_335 = buffer.data(sog0 + 335);
    const auto *sog0_339 = buffer.data(sog0 + 339);
    const auto *sog0_342 = buffer.data(sog0 + 342);
    const auto *sog0_343 = buffer.data(sog0 + 343);
    const auto *sog0_344 = buffer.data(sog0 + 344);
    const auto *sog0_345 = buffer.data(sog0 + 345);
    const auto *sog0_348 = buffer.data(sog0 + 348);
    const auto *sog0_350 = buffer.data(sog0 + 350);
    const auto *sog0_351 = buffer.data(sog0 + 351);
    const auto *sog0_354 = buffer.data(sog0 + 354);
    const auto *sog0_355 = buffer.data(sog0 + 355);
    const auto *sog0_357 = buffer.data(sog0 + 357);
    const auto *sog0_358 = buffer.data(sog0 + 358);
    const auto *sog0_359 = buffer.data(sog0 + 359);
    const auto *sog0_360 = buffer.data(sog0 + 360);
    const auto *sog0_363 = buffer.data(sog0 + 363);
    const auto *sog0_365 = buffer.data(sog0 + 365);
    const auto *sog0_366 = buffer.data(sog0 + 366);
    const auto *sog0_369 = buffer.data(sog0 + 369);
    const auto *sog0_370 = buffer.data(sog0 + 370);
    const auto *sog0_372 = buffer.data(sog0 + 372);

    const auto *sog1_310 = buffer.data(sog1 + 310);
    const auto *sog1_312 = buffer.data(sog1 + 312);
    const auto *sog1_313 = buffer.data(sog1 + 313);
    const auto *sog1_314 = buffer.data(sog1 + 314);
    const auto *sog1_315 = buffer.data(sog1 + 315);
    const auto *sog1_318 = buffer.data(sog1 + 318);
    const auto *sog1_320 = buffer.data(sog1 + 320);
    const auto *sog1_321 = buffer.data(sog1 + 321);
    const auto *sog1_324 = buffer.data(sog1 + 324);
    const auto *sog1_325 = buffer.data(sog1 + 325);
    const auto *sog1_327 = buffer.data(sog1 + 327);
    const auto *sog1_328 = buffer.data(sog1 + 328);
    const auto *sog1_329 = buffer.data(sog1 + 329);
    const auto *sog1_335 = buffer.data(sog1 + 335);
    const auto *sog1_339 = buffer.data(sog1 + 339);
    const auto *sog1_342 = buffer.data(sog1 + 342);
    const auto *sog1_343 = buffer.data(sog1 + 343);
    const auto *sog1_344 = buffer.data(sog1 + 344);
    const auto *sog1_345 = buffer.data(sog1 + 345);
    const auto *sog1_348 = buffer.data(sog1 + 348);
    const auto *sog1_350 = buffer.data(sog1 + 350);
    const auto *sog1_351 = buffer.data(sog1 + 351);
    const auto *sog1_354 = buffer.data(sog1 + 354);
    const auto *sog1_355 = buffer.data(sog1 + 355);
    const auto *sog1_357 = buffer.data(sog1 + 357);
    const auto *sog1_358 = buffer.data(sog1 + 358);
    const auto *sog1_359 = buffer.data(sog1 + 359);
    const auto *sog1_360 = buffer.data(sog1 + 360);
    const auto *sog1_363 = buffer.data(sog1 + 363);
    const auto *sog1_365 = buffer.data(sog1 + 365);
    const auto *sog1_366 = buffer.data(sog1 + 366);
    const auto *sog1_369 = buffer.data(sog1 + 369);
    const auto *sog1_370 = buffer.data(sog1 + 370);
    const auto *sog1_372 = buffer.data(sog1 + 372);

    const auto *soh_435 = buffer.data(soh + 435);
    const auto *soh_436 = buffer.data(soh + 436);
    const auto *soh_437 = buffer.data(soh + 437);
    const auto *soh_438 = buffer.data(soh + 438);
    const auto *soh_439 = buffer.data(soh + 439);
    const auto *soh_440 = buffer.data(soh + 440);
    const auto *soh_441 = buffer.data(soh + 441);
    const auto *soh_443 = buffer.data(soh + 443);
    const auto *soh_444 = buffer.data(soh + 444);
    const auto *soh_446 = buffer.data(soh + 446);
    const auto *soh_447 = buffer.data(soh + 447);
    const auto *soh_450 = buffer.data(soh + 450);
    const auto *soh_451 = buffer.data(soh + 451);
    const auto *soh_453 = buffer.data(soh + 453);
    const auto *soh_455 = buffer.data(soh + 455);
    const auto *soh_456 = buffer.data(soh + 456);
    const auto *soh_457 = buffer.data(soh + 457);
    const auto *soh_458 = buffer.data(soh + 458);
    const auto *soh_459 = buffer.data(soh + 459);
    const auto *soh_460 = buffer.data(soh + 460);
    const auto *soh_461 = buffer.data(soh + 461);
    const auto *soh_462 = buffer.data(soh + 462);
    const auto *soh_464 = buffer.data(soh + 464);
    const auto *soh_465 = buffer.data(soh + 465);
    const auto *soh_467 = buffer.data(soh + 467);
    const auto *soh_468 = buffer.data(soh + 468);
    const auto *soh_471 = buffer.data(soh + 471);
    const auto *soh_476 = buffer.data(soh + 476);
    const auto *soh_477 = buffer.data(soh + 477);
    const auto *soh_478 = buffer.data(soh + 478);
    const auto *soh_479 = buffer.data(soh + 479);
    const auto *soh_480 = buffer.data(soh + 480);
    const auto *soh_481 = buffer.data(soh + 481);
    const auto *soh_482 = buffer.data(soh + 482);
    const auto *soh_483 = buffer.data(soh + 483);
    const auto *soh_485 = buffer.data(soh + 485);
    const auto *soh_486 = buffer.data(soh + 486);
    const auto *soh_488 = buffer.data(soh + 488);
    const auto *soh_489 = buffer.data(soh + 489);
    const auto *soh_492 = buffer.data(soh + 492);
    const auto *soh_493 = buffer.data(soh + 493);
    const auto *soh_495 = buffer.data(soh + 495);
    const auto *soh_497 = buffer.data(soh + 497);
    const auto *soh_498 = buffer.data(soh + 498);
    const auto *soh_499 = buffer.data(soh + 499);
    const auto *soh_500 = buffer.data(soh + 500);
    const auto *soh_501 = buffer.data(soh + 501);
    const auto *soh_502 = buffer.data(soh + 502);
    const auto *soh_503 = buffer.data(soh + 503);
    const auto *soh_504 = buffer.data(soh + 504);
    const auto *soh_506 = buffer.data(soh + 506);
    const auto *soh_507 = buffer.data(soh + 507);
    const auto *soh_509 = buffer.data(soh + 509);
    const auto *soh_510 = buffer.data(soh + 510);
    const auto *soh_513 = buffer.data(soh + 513);
    const auto *soh_514 = buffer.data(soh + 514);
    const auto *soh_516 = buffer.data(soh + 516);

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, snh_435, snh_436, snh_437, \
                         snh_438, snh_439, soh_435, soh_436, soh_437, soh_438, \
                         soh_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_19 * snh_435[k]
                   + f_3 * pc_x[k] * soh_435[k];

        t_576[k] = f_19 * snh_436[k]
                   + f_3 * pc_x[k] * soh_436[k];

        t_577[k] = f_19 * snh_437[k]
                   + f_3 * pc_x[k] * soh_437[k];

        t_578[k] = f_19 * snh_438[k]
                   + f_3 * pc_x[k] * soh_438[k];

        t_579[k] = f_19 * snh_439[k]
                   + f_3 * pc_x[k] * soh_439[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, pc_z, snh_309, snh_440, \
                         sog0_310, sog0_312, sog1_310, sog1_312, soh_435, soh_437, \
                         soh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_19 * snh_440[k]
                   + f_3 * pc_x[k] * soh_440[k];

        t_581[k] = f_1 * sog0_310[k]
                   - f_2 * sog1_310[k]
                   + f_3 * pc_y[k] * soh_435[k];

        t_582[k] = f_20 * snh_309[k]
                   + f_3 * pc_z[k] * soh_435[k];

        t_583[k] = f_4 * sog0_312[k]
                   - f_5 * sog1_312[k]
                   + f_3 * pc_y[k] * soh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, snh_314, sog0_313, sog0_314, \
                         sog1_313, sog1_314, soh_438, soh_439, \
                         soh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * sog0_313[k]
                   - f_7 * sog1_313[k]
                   + f_3 * pc_y[k] * soh_438[k];

        t_585[k] = f_8 * sog0_314[k]
                   - f_9 * sog1_314[k]
                   + f_3 * pc_y[k] * soh_439[k];

        t_586[k] = f_3 * pc_y[k] * soh_440[k];

        t_587[k] = f_20 * snh_314[k]
                   + f_1 * sog0_314[k]
                   - f_2 * sog1_314[k]
                   + f_3 * pc_z[k] * soh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, snh_315, snh_441, \
                         snh_444, sog0_315, sog0_318, sog1_315, sog1_318, soh_441, \
                         soh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_20 * snh_441[k]
                   + f_1 * sog0_315[k]
                   - f_2 * sog1_315[k]
                   + f_3 * pc_x[k] * soh_441[k];

        t_589[k] = f_19 * snh_315[k]
                   + f_3 * pc_y[k] * soh_441[k];

        t_590[k] = f_3 * pc_z[k] * soh_441[k];

        t_591[k] = f_20 * snh_444[k]
                   + f_4 * sog0_318[k]
                   - f_5 * sog1_318[k]
                   + f_3 * pc_x[k] * soh_444[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, pc_x, pc_y, snh_317, snh_446, snh_447, sog0_320, \
                         sog0_321, sog1_320, sog1_321, soh_443, soh_446, \
                         soh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_19 * snh_317[k]
                   + f_3 * pc_y[k] * soh_443[k];

        t_593[k] = f_20 * snh_446[k]
                   + f_4 * sog0_320[k]
                   - f_5 * sog1_320[k]
                   + f_3 * pc_x[k] * soh_446[k];

        t_594[k] = f_20 * snh_447[k]
                   + f_6 * sog0_321[k]
                   - f_7 * sog1_321[k]
                   + f_3 * pc_x[k] * soh_447[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_x, pc_y, pc_z, snh_320, snh_450, sog0_324, \
                         sog1_324, soh_444, soh_446, soh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_3 * pc_z[k] * soh_444[k];

        t_596[k] = f_19 * snh_320[k]
                   + f_3 * pc_y[k] * soh_446[k];

        t_597[k] = f_20 * snh_450[k]
                   + f_6 * sog0_324[k]
                   - f_7 * sog1_324[k]
                   + f_3 * pc_x[k] * soh_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_z, snh_451, snh_453, sog0_325, \
                         sog0_327, sog1_325, sog1_327, soh_447, soh_451, \
                         soh_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_20 * snh_451[k]
                   + f_8 * sog0_325[k]
                   - f_9 * sog1_325[k]
                   + f_3 * pc_x[k] * soh_451[k];

        t_599[k] = f_3 * pc_z[k] * soh_447[k];

        t_600[k] = f_20 * snh_453[k]
                   + f_8 * sog0_327[k]
                   - f_9 * sog1_327[k]
                   + f_3 * pc_x[k] * soh_453[k];
    }

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, snh_324, snh_455, snh_456, \
                         snh_457, sog0_329, sog1_329, soh_450, soh_455, soh_456, \
                         soh_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_19 * snh_324[k]
                   + f_3 * pc_y[k] * soh_450[k];

        t_602[k] = f_20 * snh_455[k]
                   + f_8 * sog0_329[k]
                   - f_9 * sog1_329[k]
                   + f_3 * pc_x[k] * soh_455[k];

        t_603[k] = f_20 * snh_456[k]
                   + f_3 * pc_x[k] * soh_456[k];

        t_604[k] = f_20 * snh_457[k]
                   + f_3 * pc_x[k] * soh_457[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, snh_458, snh_459, snh_460, snh_461, \
                         soh_458, soh_459, soh_460, soh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_20 * snh_458[k]
                   + f_3 * pc_x[k] * soh_458[k];

        t_606[k] = f_20 * snh_459[k]
                   + f_3 * pc_x[k] * soh_459[k];

        t_607[k] = f_20 * snh_460[k]
                   + f_3 * pc_x[k] * soh_460[k];

        t_608[k] = f_20 * snh_461[k]
                   + f_3 * pc_x[k] * soh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_y, pc_z, snh_330, snh_332, sog0_325, \
                         sog0_327, sog1_325, sog1_327, soh_456, \
                         soh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_19 * snh_330[k]
                   + f_1 * sog0_325[k]
                   - f_2 * sog1_325[k]
                   + f_3 * pc_y[k] * soh_456[k];

        t_610[k] = f_3 * pc_z[k] * soh_456[k];

        t_611[k] = f_19 * snh_332[k]
                   + f_4 * sog0_327[k]
                   - f_5 * sog1_327[k]
                   + f_3 * pc_y[k] * soh_458[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_y, pc_z, snh_333, snh_334, snh_335, \
                         sog0_328, sog0_329, sog1_328, sog1_329, soh_459, soh_460, \
                         soh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_19 * snh_333[k]
                   + f_6 * sog0_328[k]
                   - f_7 * sog1_328[k]
                   + f_3 * pc_y[k] * soh_459[k];

        t_613[k] = f_19 * snh_334[k]
                   + f_8 * sog0_329[k]
                   - f_9 * sog1_329[k]
                   + f_3 * pc_y[k] * soh_460[k];

        t_614[k] = f_19 * snh_335[k]
                   + f_3 * pc_y[k] * soh_461[k];

        t_615[k] = f_1 * sog0_329[k]
                   - f_2 * sog1_329[k]
                   + f_3 * pc_z[k] * soh_461[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_z, pc_y, pc_z, sni0_420, sni0_423, \
                         snh_315, snh_336, sni1_420, sni1_423, \
                         soh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_z[k] * sni0_420[k]
                   - f_10 * pc_z[k] * sni1_420[k];

        t_617[k] = f_20 * snh_336[k]
                   + f_3 * pc_y[k] * soh_462[k];

        t_618[k] = f_11 * snh_315[k]
                   + f_3 * pc_z[k] * soh_462[k];

        t_619[k] = pb_z[k] * sni0_423[k]
                   - f_10 * pc_z[k] * sni1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pb_z, pc_x, pc_y, pc_z, sni0_426, snh_338, \
                         snh_467, sni1_426, sog0_335, sog1_335, soh_464, \
                         soh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_20 * snh_338[k]
                   + f_3 * pc_y[k] * soh_464[k];

        t_621[k] = f_20 * snh_467[k]
                   + f_4 * sog0_335[k]
                   - f_5 * sog1_335[k]
                   + f_3 * pc_x[k] * soh_467[k];

        t_622[k] = pb_z[k] * sni0_426[k]
                   - f_10 * pc_z[k] * sni1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, pc_y, pc_z, snh_318, snh_341, snh_471, \
                         sog0_339, sog1_339, soh_465, soh_467, \
                         soh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * snh_318[k]
                   + f_3 * pc_z[k] * soh_465[k];

        t_624[k] = f_20 * snh_341[k]
                   + f_3 * pc_y[k] * soh_467[k];

        t_625[k] = f_20 * snh_471[k]
                   + f_6 * sog0_339[k]
                   - f_7 * sog1_339[k]
                   + f_3 * pc_x[k] * soh_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pb_z, pc_y, pc_z, sni0_430, sni0_432, \
                         snh_321, snh_322, snh_345, sni1_430, sni1_432, soh_468, \
                         soh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pb_z[k] * sni0_430[k]
                   - f_10 * pc_z[k] * sni1_430[k];

        t_627[k] = f_11 * snh_321[k]
                   + f_3 * pc_z[k] * soh_468[k];

        t_628[k] = pb_z[k] * sni0_432[k]
                   + f_12 * snh_322[k]
                   - f_10 * pc_z[k] * sni1_432[k];

        t_629[k] = f_20 * snh_345[k]
                   + f_3 * pc_y[k] * soh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, snh_476, snh_477, snh_478, snh_479, \
                         sog0_344, sog1_344, soh_476, soh_477, soh_478, \
                         soh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_20 * snh_476[k]
                   + f_8 * sog0_344[k]
                   - f_9 * sog1_344[k]
                   + f_3 * pc_x[k] * soh_476[k];

        t_631[k] = f_20 * snh_477[k]
                   + f_3 * pc_x[k] * soh_477[k];

        t_632[k] = f_20 * snh_478[k]
                   + f_3 * pc_x[k] * soh_478[k];

        t_633[k] = f_20 * snh_479[k]
                   + f_3 * pc_x[k] * soh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_z, pc_x, pc_z, sni0_441, snh_480, \
                         snh_481, snh_482, sni1_441, soh_480, soh_481, \
                         soh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_20 * snh_480[k]
                   + f_3 * pc_x[k] * soh_480[k];

        t_635[k] = f_20 * snh_481[k]
                   + f_3 * pc_x[k] * soh_481[k];

        t_636[k] = f_20 * snh_482[k]
                   + f_3 * pc_x[k] * soh_482[k];

        t_637[k] = pb_z[k] * sni0_441[k]
                   - f_10 * pc_z[k] * sni1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, snh_330, snh_353, snh_354, sog0_342, \
                         sog0_343, sog1_342, sog1_343, soh_477, soh_479, \
                         soh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * snh_330[k]
                   + f_3 * pc_z[k] * soh_477[k];

        t_639[k] = f_20 * snh_353[k]
                   + f_4 * sog0_342[k]
                   - f_5 * sog1_342[k]
                   + f_3 * pc_y[k] * soh_479[k];

        t_640[k] = f_20 * snh_354[k]
                   + f_6 * sog0_343[k]
                   - f_7 * sog1_343[k]
                   + f_3 * pc_y[k] * soh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, snh_335, snh_355, snh_356, sog0_344, \
                         sog1_344, soh_481, soh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_20 * snh_355[k]
                   + f_8 * sog0_344[k]
                   - f_9 * sog1_344[k]
                   + f_3 * pc_y[k] * soh_481[k];

        t_642[k] = f_20 * snh_356[k]
                   + f_3 * pc_y[k] * soh_482[k];

        t_643[k] = f_11 * snh_335[k]
                   + f_1 * sog0_344[k]
                   - f_2 * sog1_344[k]
                   + f_3 * pc_z[k] * soh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, snh_336, snh_357, snh_483, \
                         sog0_345, sog1_345, soh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_20 * snh_483[k]
                   + f_1 * sog0_345[k]
                   - f_2 * sog1_345[k]
                   + f_3 * pc_x[k] * soh_483[k];

        t_645[k] = f_14 * snh_357[k]
                   + f_3 * pc_y[k] * soh_483[k];

        t_646[k] = f_12 * snh_336[k]
                   + f_3 * pc_z[k] * soh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, snh_359, snh_486, snh_488, sog0_348, \
                         sog0_350, sog1_348, sog1_350, soh_485, soh_486, \
                         soh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_20 * snh_486[k]
                   + f_4 * sog0_348[k]
                   - f_5 * sog1_348[k]
                   + f_3 * pc_x[k] * soh_486[k];

        t_648[k] = f_14 * snh_359[k]
                   + f_3 * pc_y[k] * soh_485[k];

        t_649[k] = f_20 * snh_488[k]
                   + f_4 * sog0_350[k]
                   - f_5 * sog1_350[k]
                   + f_3 * pc_x[k] * soh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, snh_339, snh_362, snh_489, \
                         sog0_351, sog1_351, soh_486, soh_488, \
                         soh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_20 * snh_489[k]
                   + f_6 * sog0_351[k]
                   - f_7 * sog1_351[k]
                   + f_3 * pc_x[k] * soh_489[k];

        t_651[k] = f_12 * snh_339[k]
                   + f_3 * pc_z[k] * soh_486[k];

        t_652[k] = f_14 * snh_362[k]
                   + f_3 * pc_y[k] * soh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, snh_342, snh_492, snh_493, sog0_354, \
                         sog0_355, sog1_354, sog1_355, soh_489, soh_492, \
                         soh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_20 * snh_492[k]
                   + f_6 * sog0_354[k]
                   - f_7 * sog1_354[k]
                   + f_3 * pc_x[k] * soh_492[k];

        t_654[k] = f_20 * snh_493[k]
                   + f_8 * sog0_355[k]
                   - f_9 * sog1_355[k]
                   + f_3 * pc_x[k] * soh_493[k];

        t_655[k] = f_12 * snh_342[k]
                   + f_3 * pc_z[k] * soh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, snh_366, snh_495, snh_497, sog0_357, \
                         sog0_359, sog1_357, sog1_359, soh_492, soh_495, \
                         soh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_20 * snh_495[k]
                   + f_8 * sog0_357[k]
                   - f_9 * sog1_357[k]
                   + f_3 * pc_x[k] * soh_495[k];

        t_657[k] = f_14 * snh_366[k]
                   + f_3 * pc_y[k] * soh_492[k];

        t_658[k] = f_20 * snh_497[k]
                   + f_8 * sog0_359[k]
                   - f_9 * sog1_359[k]
                   + f_3 * pc_x[k] * soh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, snh_498, snh_499, snh_500, \
                         snh_501, snh_502, soh_498, soh_499, soh_500, soh_501, \
                         soh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_20 * snh_498[k]
                   + f_3 * pc_x[k] * soh_498[k];

        t_660[k] = f_20 * snh_499[k]
                   + f_3 * pc_x[k] * soh_499[k];

        t_661[k] = f_20 * snh_500[k]
                   + f_3 * pc_x[k] * soh_500[k];

        t_662[k] = f_20 * snh_501[k]
                   + f_3 * pc_x[k] * soh_501[k];

        t_663[k] = f_20 * snh_502[k]
                   + f_3 * pc_x[k] * soh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, snh_351, snh_372, snh_503, \
                         sog0_355, sog1_355, soh_498, soh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_20 * snh_503[k]
                   + f_3 * pc_x[k] * soh_503[k];

        t_665[k] = f_14 * snh_372[k]
                   + f_1 * sog0_355[k]
                   - f_2 * sog1_355[k]
                   + f_3 * pc_y[k] * soh_498[k];

        t_666[k] = f_12 * snh_351[k]
                   + f_3 * pc_z[k] * soh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, snh_374, snh_375, snh_376, sog0_357, \
                         sog0_358, sog0_359, sog1_357, sog1_358, sog1_359, soh_500, soh_501, \
                         soh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * snh_374[k]
                   + f_4 * sog0_357[k]
                   - f_5 * sog1_357[k]
                   + f_3 * pc_y[k] * soh_500[k];

        t_668[k] = f_14 * snh_375[k]
                   + f_6 * sog0_358[k]
                   - f_7 * sog1_358[k]
                   + f_3 * pc_y[k] * soh_501[k];

        t_669[k] = f_14 * snh_376[k]
                   + f_8 * sog0_359[k]
                   - f_9 * sog1_359[k]
                   + f_3 * pc_y[k] * soh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, snh_356, snh_377, snh_504, \
                         sog0_359, sog0_360, sog1_359, sog1_360, soh_503, \
                         soh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * snh_377[k]
                   + f_3 * pc_y[k] * soh_503[k];

        t_671[k] = f_12 * snh_356[k]
                   + f_1 * sog0_359[k]
                   - f_2 * sog1_359[k]
                   + f_3 * pc_z[k] * soh_503[k];

        t_672[k] = f_20 * snh_504[k]
                   + f_1 * sog0_360[k]
                   - f_2 * sog1_360[k]
                   + f_3 * pc_x[k] * soh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, snh_357, snh_378, \
                         snh_380, snh_507, sog0_363, sog1_363, soh_504, soh_506, \
                         soh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * snh_378[k]
                   + f_3 * pc_y[k] * soh_504[k];

        t_674[k] = f_13 * snh_357[k]
                   + f_3 * pc_z[k] * soh_504[k];

        t_675[k] = f_20 * snh_507[k]
                   + f_4 * sog0_363[k]
                   - f_5 * sog1_363[k]
                   + f_3 * pc_x[k] * soh_507[k];

        t_676[k] = f_13 * snh_380[k]
                   + f_3 * pc_y[k] * soh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, snh_360, snh_509, snh_510, sog0_365, \
                         sog0_366, sog1_365, sog1_366, soh_507, soh_509, \
                         soh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_20 * snh_509[k]
                   + f_4 * sog0_365[k]
                   - f_5 * sog1_365[k]
                   + f_3 * pc_x[k] * soh_509[k];

        t_678[k] = f_20 * snh_510[k]
                   + f_6 * sog0_366[k]
                   - f_7 * sog1_366[k]
                   + f_3 * pc_x[k] * soh_510[k];

        t_679[k] = f_13 * snh_360[k]
                   + f_3 * pc_z[k] * soh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, snh_383, snh_513, snh_514, sog0_369, \
                         sog0_370, sog1_369, sog1_370, soh_509, soh_513, \
                         soh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * snh_383[k]
                   + f_3 * pc_y[k] * soh_509[k];

        t_681[k] = f_20 * snh_513[k]
                   + f_6 * sog0_369[k]
                   - f_7 * sog1_369[k]
                   + f_3 * pc_x[k] * soh_513[k];

        t_682[k] = f_20 * snh_514[k]
                   + f_8 * sog0_370[k]
                   - f_9 * sog1_370[k]
                   + f_3 * pc_x[k] * soh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, snh_363, snh_387, snh_516, \
                         sog0_372, sog1_372, soh_510, soh_513, \
                         soh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * snh_363[k]
                   + f_3 * pc_z[k] * soh_510[k];

        t_684[k] = f_20 * snh_516[k]
                   + f_8 * sog0_372[k]
                   - f_9 * sog1_372[k]
                   + f_3 * pc_x[k] * soh_516[k];

        t_685[k] = f_13 * snh_387[k]
                   + f_3 * pc_y[k] * soh_513[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_560 = buffer.data(sni0 + 560);
    const auto *sni0_563 = buffer.data(sni0 + 563);
    const auto *sni0_565 = buffer.data(sni0 + 565);
    const auto *sni0_566 = buffer.data(sni0 + 566);
    const auto *sni0_569 = buffer.data(sni0 + 569);
    const auto *sni0_570 = buffer.data(sni0 + 570);
    const auto *sni0_572 = buffer.data(sni0 + 572);
    const auto *sni0_574 = buffer.data(sni0 + 574);
    const auto *sni0_587 = buffer.data(sni0 + 587);

    const auto *snh_372 = buffer.data(snh + 372);
    const auto *snh_377 = buffer.data(snh + 377);
    const auto *snh_378 = buffer.data(snh + 378);
    const auto *snh_381 = buffer.data(snh + 381);
    const auto *snh_384 = buffer.data(snh + 384);
    const auto *snh_393 = buffer.data(snh + 393);
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
    const auto *snh_416 = buffer.data(snh + 416);
    const auto *snh_417 = buffer.data(snh + 417);
    const auto *snh_418 = buffer.data(snh + 418);
    const auto *snh_419 = buffer.data(snh + 419);
    const auto *snh_420 = buffer.data(snh + 420);
    const auto *snh_421 = buffer.data(snh + 421);
    const auto *snh_422 = buffer.data(snh + 422);
    const auto *snh_423 = buffer.data(snh + 423);
    const auto *snh_425 = buffer.data(snh + 425);
    const auto *snh_426 = buffer.data(snh + 426);
    const auto *snh_428 = buffer.data(snh + 428);
    const auto *snh_429 = buffer.data(snh + 429);
    const auto *snh_435 = buffer.data(snh + 435);
    const auto *snh_437 = buffer.data(snh + 437);
    const auto *snh_438 = buffer.data(snh + 438);
    const auto *snh_439 = buffer.data(snh + 439);
    const auto *snh_440 = buffer.data(snh + 440);
    const auto *snh_441 = buffer.data(snh + 441);
    const auto *snh_443 = buffer.data(snh + 443);
    const auto *snh_446 = buffer.data(snh + 446);
    const auto *snh_518 = buffer.data(snh + 518);
    const auto *snh_519 = buffer.data(snh + 519);
    const auto *snh_520 = buffer.data(snh + 520);
    const auto *snh_521 = buffer.data(snh + 521);
    const auto *snh_522 = buffer.data(snh + 522);
    const auto *snh_523 = buffer.data(snh + 523);
    const auto *snh_524 = buffer.data(snh + 524);
    const auto *snh_525 = buffer.data(snh + 525);
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
    const auto *snh_561 = buffer.data(snh + 561);
    const auto *snh_562 = buffer.data(snh + 562);
    const auto *snh_563 = buffer.data(snh + 563);
    const auto *snh_564 = buffer.data(snh + 564);
    const auto *snh_565 = buffer.data(snh + 565);
    const auto *snh_566 = buffer.data(snh + 566);
    const auto *snh_567 = buffer.data(snh + 567);
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
    const auto *snh_591 = buffer.data(snh + 591);
    const auto *snh_593 = buffer.data(snh + 593);
    const auto *snh_594 = buffer.data(snh + 594);
    const auto *snh_597 = buffer.data(snh + 597);
    const auto *snh_598 = buffer.data(snh + 598);
    const auto *snh_600 = buffer.data(snh + 600);

    const auto *sni1_560 = buffer.data(sni1 + 560);
    const auto *sni1_563 = buffer.data(sni1 + 563);
    const auto *sni1_565 = buffer.data(sni1 + 565);
    const auto *sni1_566 = buffer.data(sni1 + 566);
    const auto *sni1_569 = buffer.data(sni1 + 569);
    const auto *sni1_570 = buffer.data(sni1 + 570);
    const auto *sni1_572 = buffer.data(sni1 + 572);
    const auto *sni1_574 = buffer.data(sni1 + 574);
    const auto *sni1_587 = buffer.data(sni1 + 587);

    const auto *sog0_370 = buffer.data(sog0 + 370);
    const auto *sog0_372 = buffer.data(sog0 + 372);
    const auto *sog0_373 = buffer.data(sog0 + 373);
    const auto *sog0_374 = buffer.data(sog0 + 374);
    const auto *sog0_375 = buffer.data(sog0 + 375);
    const auto *sog0_378 = buffer.data(sog0 + 378);
    const auto *sog0_380 = buffer.data(sog0 + 380);
    const auto *sog0_381 = buffer.data(sog0 + 381);
    const auto *sog0_384 = buffer.data(sog0 + 384);
    const auto *sog0_385 = buffer.data(sog0 + 385);
    const auto *sog0_387 = buffer.data(sog0 + 387);
    const auto *sog0_388 = buffer.data(sog0 + 388);
    const auto *sog0_389 = buffer.data(sog0 + 389);
    const auto *sog0_400 = buffer.data(sog0 + 400);
    const auto *sog0_402 = buffer.data(sog0 + 402);
    const auto *sog0_403 = buffer.data(sog0 + 403);
    const auto *sog0_404 = buffer.data(sog0 + 404);
    const auto *sog0_405 = buffer.data(sog0 + 405);
    const auto *sog0_408 = buffer.data(sog0 + 408);
    const auto *sog0_410 = buffer.data(sog0 + 410);
    const auto *sog0_411 = buffer.data(sog0 + 411);
    const auto *sog0_414 = buffer.data(sog0 + 414);
    const auto *sog0_415 = buffer.data(sog0 + 415);
    const auto *sog0_417 = buffer.data(sog0 + 417);
    const auto *sog0_418 = buffer.data(sog0 + 418);
    const auto *sog0_419 = buffer.data(sog0 + 419);
    const auto *sog0_420 = buffer.data(sog0 + 420);
    const auto *sog0_423 = buffer.data(sog0 + 423);
    const auto *sog0_425 = buffer.data(sog0 + 425);
    const auto *sog0_426 = buffer.data(sog0 + 426);
    const auto *sog0_429 = buffer.data(sog0 + 429);
    const auto *sog0_430 = buffer.data(sog0 + 430);
    const auto *sog0_432 = buffer.data(sog0 + 432);

    const auto *sog1_370 = buffer.data(sog1 + 370);
    const auto *sog1_372 = buffer.data(sog1 + 372);
    const auto *sog1_373 = buffer.data(sog1 + 373);
    const auto *sog1_374 = buffer.data(sog1 + 374);
    const auto *sog1_375 = buffer.data(sog1 + 375);
    const auto *sog1_378 = buffer.data(sog1 + 378);
    const auto *sog1_380 = buffer.data(sog1 + 380);
    const auto *sog1_381 = buffer.data(sog1 + 381);
    const auto *sog1_384 = buffer.data(sog1 + 384);
    const auto *sog1_385 = buffer.data(sog1 + 385);
    const auto *sog1_387 = buffer.data(sog1 + 387);
    const auto *sog1_388 = buffer.data(sog1 + 388);
    const auto *sog1_389 = buffer.data(sog1 + 389);
    const auto *sog1_400 = buffer.data(sog1 + 400);
    const auto *sog1_402 = buffer.data(sog1 + 402);
    const auto *sog1_403 = buffer.data(sog1 + 403);
    const auto *sog1_404 = buffer.data(sog1 + 404);
    const auto *sog1_405 = buffer.data(sog1 + 405);
    const auto *sog1_408 = buffer.data(sog1 + 408);
    const auto *sog1_410 = buffer.data(sog1 + 410);
    const auto *sog1_411 = buffer.data(sog1 + 411);
    const auto *sog1_414 = buffer.data(sog1 + 414);
    const auto *sog1_415 = buffer.data(sog1 + 415);
    const auto *sog1_417 = buffer.data(sog1 + 417);
    const auto *sog1_418 = buffer.data(sog1 + 418);
    const auto *sog1_419 = buffer.data(sog1 + 419);
    const auto *sog1_420 = buffer.data(sog1 + 420);
    const auto *sog1_423 = buffer.data(sog1 + 423);
    const auto *sog1_425 = buffer.data(sog1 + 425);
    const auto *sog1_426 = buffer.data(sog1 + 426);
    const auto *sog1_429 = buffer.data(sog1 + 429);
    const auto *sog1_430 = buffer.data(sog1 + 430);
    const auto *sog1_432 = buffer.data(sog1 + 432);

    const auto *soh_518 = buffer.data(soh + 518);
    const auto *soh_519 = buffer.data(soh + 519);
    const auto *soh_520 = buffer.data(soh + 520);
    const auto *soh_521 = buffer.data(soh + 521);
    const auto *soh_522 = buffer.data(soh + 522);
    const auto *soh_523 = buffer.data(soh + 523);
    const auto *soh_524 = buffer.data(soh + 524);
    const auto *soh_525 = buffer.data(soh + 525);
    const auto *soh_527 = buffer.data(soh + 527);
    const auto *soh_528 = buffer.data(soh + 528);
    const auto *soh_530 = buffer.data(soh + 530);
    const auto *soh_531 = buffer.data(soh + 531);
    const auto *soh_534 = buffer.data(soh + 534);
    const auto *soh_535 = buffer.data(soh + 535);
    const auto *soh_537 = buffer.data(soh + 537);
    const auto *soh_539 = buffer.data(soh + 539);
    const auto *soh_540 = buffer.data(soh + 540);
    const auto *soh_541 = buffer.data(soh + 541);
    const auto *soh_542 = buffer.data(soh + 542);
    const auto *soh_543 = buffer.data(soh + 543);
    const auto *soh_544 = buffer.data(soh + 544);
    const auto *soh_545 = buffer.data(soh + 545);
    const auto *soh_546 = buffer.data(soh + 546);
    const auto *soh_548 = buffer.data(soh + 548);
    const auto *soh_549 = buffer.data(soh + 549);
    const auto *soh_551 = buffer.data(soh + 551);
    const auto *soh_552 = buffer.data(soh + 552);
    const auto *soh_555 = buffer.data(soh + 555);
    const auto *soh_561 = buffer.data(soh + 561);
    const auto *soh_562 = buffer.data(soh + 562);
    const auto *soh_563 = buffer.data(soh + 563);
    const auto *soh_564 = buffer.data(soh + 564);
    const auto *soh_565 = buffer.data(soh + 565);
    const auto *soh_566 = buffer.data(soh + 566);
    const auto *soh_567 = buffer.data(soh + 567);
    const auto *soh_569 = buffer.data(soh + 569);
    const auto *soh_570 = buffer.data(soh + 570);
    const auto *soh_572 = buffer.data(soh + 572);
    const auto *soh_573 = buffer.data(soh + 573);
    const auto *soh_576 = buffer.data(soh + 576);
    const auto *soh_577 = buffer.data(soh + 577);
    const auto *soh_579 = buffer.data(soh + 579);
    const auto *soh_581 = buffer.data(soh + 581);
    const auto *soh_582 = buffer.data(soh + 582);
    const auto *soh_583 = buffer.data(soh + 583);
    const auto *soh_584 = buffer.data(soh + 584);
    const auto *soh_585 = buffer.data(soh + 585);
    const auto *soh_586 = buffer.data(soh + 586);
    const auto *soh_587 = buffer.data(soh + 587);
    const auto *soh_588 = buffer.data(soh + 588);
    const auto *soh_590 = buffer.data(soh + 590);
    const auto *soh_591 = buffer.data(soh + 591);
    const auto *soh_593 = buffer.data(soh + 593);
    const auto *soh_594 = buffer.data(soh + 594);
    const auto *soh_597 = buffer.data(soh + 597);
    const auto *soh_598 = buffer.data(soh + 598);
    const auto *soh_600 = buffer.data(soh + 600);

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, snh_518, snh_519, snh_520, snh_521, \
                         sog0_374, sog1_374, soh_518, soh_519, soh_520, \
                         soh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_20 * snh_518[k]
                   + f_8 * sog0_374[k]
                   - f_9 * sog1_374[k]
                   + f_3 * pc_x[k] * soh_518[k];

        t_687[k] = f_20 * snh_519[k]
                   + f_3 * pc_x[k] * soh_519[k];

        t_688[k] = f_20 * snh_520[k]
                   + f_3 * pc_x[k] * soh_520[k];

        t_689[k] = f_20 * snh_521[k]
                   + f_3 * pc_x[k] * soh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, snh_393, snh_522, snh_523, \
                         snh_524, sog0_370, sog1_370, soh_519, soh_522, soh_523, \
                         soh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_20 * snh_522[k]
                   + f_3 * pc_x[k] * soh_522[k];

        t_691[k] = f_20 * snh_523[k]
                   + f_3 * pc_x[k] * soh_523[k];

        t_692[k] = f_20 * snh_524[k]
                   + f_3 * pc_x[k] * soh_524[k];

        t_693[k] = f_13 * snh_393[k]
                   + f_1 * sog0_370[k]
                   - f_2 * sog1_370[k]
                   + f_3 * pc_y[k] * soh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, snh_372, snh_395, snh_396, sog0_372, \
                         sog0_373, sog1_372, sog1_373, soh_519, soh_521, \
                         soh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * snh_372[k]
                   + f_3 * pc_z[k] * soh_519[k];

        t_695[k] = f_13 * snh_395[k]
                   + f_4 * sog0_372[k]
                   - f_5 * sog1_372[k]
                   + f_3 * pc_y[k] * soh_521[k];

        t_696[k] = f_13 * snh_396[k]
                   + f_6 * sog0_373[k]
                   - f_7 * sog1_373[k]
                   + f_3 * pc_y[k] * soh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, snh_377, snh_397, snh_398, sog0_374, \
                         sog1_374, soh_523, soh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * snh_397[k]
                   + f_8 * sog0_374[k]
                   - f_9 * sog1_374[k]
                   + f_3 * pc_y[k] * soh_523[k];

        t_698[k] = f_13 * snh_398[k]
                   + f_3 * pc_y[k] * soh_524[k];

        t_699[k] = f_13 * snh_377[k]
                   + f_1 * sog0_374[k]
                   - f_2 * sog1_374[k]
                   + f_3 * pc_z[k] * soh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, snh_378, snh_399, snh_525, \
                         sog0_375, sog1_375, soh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_20 * snh_525[k]
                   + f_1 * sog0_375[k]
                   - f_2 * sog1_375[k]
                   + f_3 * pc_x[k] * soh_525[k];

        t_701[k] = f_12 * snh_399[k]
                   + f_3 * pc_y[k] * soh_525[k];

        t_702[k] = f_14 * snh_378[k]
                   + f_3 * pc_z[k] * soh_525[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, snh_401, snh_528, snh_530, sog0_378, \
                         sog0_380, sog1_378, sog1_380, soh_527, soh_528, \
                         soh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_20 * snh_528[k]
                   + f_4 * sog0_378[k]
                   - f_5 * sog1_378[k]
                   + f_3 * pc_x[k] * soh_528[k];

        t_704[k] = f_12 * snh_401[k]
                   + f_3 * pc_y[k] * soh_527[k];

        t_705[k] = f_20 * snh_530[k]
                   + f_4 * sog0_380[k]
                   - f_5 * sog1_380[k]
                   + f_3 * pc_x[k] * soh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, snh_381, snh_404, snh_531, \
                         sog0_381, sog1_381, soh_528, soh_530, \
                         soh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_20 * snh_531[k]
                   + f_6 * sog0_381[k]
                   - f_7 * sog1_381[k]
                   + f_3 * pc_x[k] * soh_531[k];

        t_707[k] = f_14 * snh_381[k]
                   + f_3 * pc_z[k] * soh_528[k];

        t_708[k] = f_12 * snh_404[k]
                   + f_3 * pc_y[k] * soh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, snh_384, snh_534, snh_535, sog0_384, \
                         sog0_385, sog1_384, sog1_385, soh_531, soh_534, \
                         soh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_20 * snh_534[k]
                   + f_6 * sog0_384[k]
                   - f_7 * sog1_384[k]
                   + f_3 * pc_x[k] * soh_534[k];

        t_710[k] = f_20 * snh_535[k]
                   + f_8 * sog0_385[k]
                   - f_9 * sog1_385[k]
                   + f_3 * pc_x[k] * soh_535[k];

        t_711[k] = f_14 * snh_384[k]
                   + f_3 * pc_z[k] * soh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, snh_408, snh_537, snh_539, sog0_387, \
                         sog0_389, sog1_387, sog1_389, soh_534, soh_537, \
                         soh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_20 * snh_537[k]
                   + f_8 * sog0_387[k]
                   - f_9 * sog1_387[k]
                   + f_3 * pc_x[k] * soh_537[k];

        t_713[k] = f_12 * snh_408[k]
                   + f_3 * pc_y[k] * soh_534[k];

        t_714[k] = f_20 * snh_539[k]
                   + f_8 * sog0_389[k]
                   - f_9 * sog1_389[k]
                   + f_3 * pc_x[k] * soh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, snh_540, snh_541, snh_542, \
                         snh_543, snh_544, soh_540, soh_541, soh_542, soh_543, \
                         soh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_20 * snh_540[k]
                   + f_3 * pc_x[k] * soh_540[k];

        t_716[k] = f_20 * snh_541[k]
                   + f_3 * pc_x[k] * soh_541[k];

        t_717[k] = f_20 * snh_542[k]
                   + f_3 * pc_x[k] * soh_542[k];

        t_718[k] = f_20 * snh_543[k]
                   + f_3 * pc_x[k] * soh_543[k];

        t_719[k] = f_20 * snh_544[k]
                   + f_3 * pc_x[k] * soh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, snh_393, snh_414, snh_545, \
                         sog0_385, sog1_385, soh_540, soh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_20 * snh_545[k]
                   + f_3 * pc_x[k] * soh_545[k];

        t_721[k] = f_12 * snh_414[k]
                   + f_1 * sog0_385[k]
                   - f_2 * sog1_385[k]
                   + f_3 * pc_y[k] * soh_540[k];

        t_722[k] = f_14 * snh_393[k]
                   + f_3 * pc_z[k] * soh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, snh_416, snh_417, snh_418, sog0_387, \
                         sog0_388, sog0_389, sog1_387, sog1_388, sog1_389, soh_542, soh_543, \
                         soh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * snh_416[k]
                   + f_4 * sog0_387[k]
                   - f_5 * sog1_387[k]
                   + f_3 * pc_y[k] * soh_542[k];

        t_724[k] = f_12 * snh_417[k]
                   + f_6 * sog0_388[k]
                   - f_7 * sog1_388[k]
                   + f_3 * pc_y[k] * soh_543[k];

        t_725[k] = f_12 * snh_418[k]
                   + f_8 * sog0_389[k]
                   - f_9 * sog1_389[k]
                   + f_3 * pc_y[k] * soh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pb_y, pc_y, pc_z, sni0_560, snh_398, \
                         snh_419, snh_420, sni1_560, sog0_389, sog1_389, soh_545, \
                         soh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * snh_419[k]
                   + f_3 * pc_y[k] * soh_545[k];

        t_727[k] = f_14 * snh_398[k]
                   + f_1 * sog0_389[k]
                   - f_2 * sog1_389[k]
                   + f_3 * pc_z[k] * soh_545[k];

        t_728[k] = pb_y[k] * sni0_560[k]
                   - f_10 * pc_y[k] * sni1_560[k];

        t_729[k] = f_11 * snh_420[k]
                   + f_3 * pc_y[k] * soh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pb_y, pc_y, pc_z, sni0_563, sni0_565, \
                         snh_399, snh_421, snh_422, sni1_563, sni1_565, soh_546, \
                         soh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_20 * snh_399[k]
                   + f_3 * pc_z[k] * soh_546[k];

        t_731[k] = pb_y[k] * sni0_563[k]
                   + f_12 * snh_421[k]
                   - f_10 * pc_y[k] * sni1_563[k];

        t_732[k] = f_11 * snh_422[k]
                   + f_3 * pc_y[k] * soh_548[k];

        t_733[k] = pb_y[k] * sni0_565[k]
                   - f_10 * pc_y[k] * sni1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pb_y, pc_y, pc_z, sni0_566, sni0_569, \
                         snh_402, snh_423, snh_425, sni1_566, sni1_569, soh_549, \
                         soh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_y[k] * sni0_566[k]
                   + f_13 * snh_423[k]
                   - f_10 * pc_y[k] * sni1_566[k];

        t_735[k] = f_20 * snh_402[k]
                   + f_3 * pc_z[k] * soh_549[k];

        t_736[k] = f_11 * snh_425[k]
                   + f_3 * pc_y[k] * soh_551[k];

        t_737[k] = pb_y[k] * sni0_569[k]
                   - f_10 * pc_y[k] * sni1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pb_y, pc_y, pc_z, sni0_570, sni0_572, snh_405, \
                         snh_426, snh_428, sni1_570, sni1_572, \
                         soh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pb_y[k] * sni0_570[k]
                   + f_14 * snh_426[k]
                   - f_10 * pc_y[k] * sni1_570[k];

        t_739[k] = f_20 * snh_405[k]
                   + f_3 * pc_z[k] * soh_552[k];

        t_740[k] = pb_y[k] * sni0_572[k]
                   + f_12 * snh_428[k]
                   - f_10 * pc_y[k] * sni1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pb_y, pc_x, pc_y, sni0_574, snh_429, \
                         snh_561, snh_562, sni1_574, soh_555, soh_561, \
                         soh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * snh_429[k]
                   + f_3 * pc_y[k] * soh_555[k];

        t_742[k] = pb_y[k] * sni0_574[k]
                   - f_10 * pc_y[k] * sni1_574[k];

        t_743[k] = f_20 * snh_561[k]
                   + f_3 * pc_x[k] * soh_561[k];

        t_744[k] = f_20 * snh_562[k]
                   + f_3 * pc_x[k] * soh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, snh_563, snh_564, snh_565, snh_566, \
                         soh_563, soh_564, soh_565, soh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_20 * snh_563[k]
                   + f_3 * pc_x[k] * soh_563[k];

        t_746[k] = f_20 * snh_564[k]
                   + f_3 * pc_x[k] * soh_564[k];

        t_747[k] = f_20 * snh_565[k]
                   + f_3 * pc_x[k] * soh_565[k];

        t_748[k] = f_20 * snh_566[k]
                   + f_3 * pc_x[k] * soh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, snh_414, snh_435, snh_437, sog0_400, \
                         sog0_402, sog1_400, sog1_402, soh_561, \
                         soh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * snh_435[k]
                   + f_1 * sog0_400[k]
                   - f_2 * sog1_400[k]
                   + f_3 * pc_y[k] * soh_561[k];

        t_750[k] = f_20 * snh_414[k]
                   + f_3 * pc_z[k] * soh_561[k];

        t_751[k] = f_11 * snh_437[k]
                   + f_4 * sog0_402[k]
                   - f_5 * sog1_402[k]
                   + f_3 * pc_y[k] * soh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, snh_438, snh_439, snh_440, sog0_403, \
                         sog0_404, sog1_403, sog1_404, soh_564, soh_565, \
                         soh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * snh_438[k]
                   + f_6 * sog0_403[k]
                   - f_7 * sog1_403[k]
                   + f_3 * pc_y[k] * soh_564[k];

        t_753[k] = f_11 * snh_439[k]
                   + f_8 * sog0_404[k]
                   - f_9 * sog1_404[k]
                   + f_3 * pc_y[k] * soh_565[k];

        t_754[k] = f_11 * snh_440[k]
                   + f_3 * pc_y[k] * soh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pb_y, pc_x, pc_y, pc_z, sni0_587, \
                         snh_420, snh_567, sni1_587, sog0_405, sog1_405, \
                         soh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pb_y[k] * sni0_587[k]
                   - f_10 * pc_y[k] * sni1_587[k];

        t_756[k] = f_20 * snh_567[k]
                   + f_1 * sog0_405[k]
                   - f_2 * sog1_405[k]
                   + f_3 * pc_x[k] * soh_567[k];

        t_757[k] = f_3 * pc_y[k] * soh_567[k];

        t_758[k] = f_19 * snh_420[k]
                   + f_3 * pc_z[k] * soh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, snh_570, snh_572, sog0_408, \
                         sog0_410, sog1_408, sog1_410, soh_569, soh_570, \
                         soh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_20 * snh_570[k]
                   + f_4 * sog0_408[k]
                   - f_5 * sog1_408[k]
                   + f_3 * pc_x[k] * soh_570[k];

        t_760[k] = f_3 * pc_y[k] * soh_569[k];

        t_761[k] = f_20 * snh_572[k]
                   + f_4 * sog0_410[k]
                   - f_5 * sog1_410[k]
                   + f_3 * pc_x[k] * soh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_x, pc_y, pc_z, snh_423, snh_573, sog0_411, \
                         sog1_411, soh_570, soh_572, soh_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_20 * snh_573[k]
                   + f_6 * sog0_411[k]
                   - f_7 * sog1_411[k]
                   + f_3 * pc_x[k] * soh_573[k];

        t_763[k] = f_19 * snh_423[k]
                   + f_3 * pc_z[k] * soh_570[k];

        t_764[k] = f_3 * pc_y[k] * soh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_z, snh_426, snh_576, snh_577, sog0_414, \
                         sog0_415, sog1_414, sog1_415, soh_573, soh_576, \
                         soh_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_20 * snh_576[k]
                   + f_6 * sog0_414[k]
                   - f_7 * sog1_414[k]
                   + f_3 * pc_x[k] * soh_576[k];

        t_766[k] = f_20 * snh_577[k]
                   + f_8 * sog0_415[k]
                   - f_9 * sog1_415[k]
                   + f_3 * pc_x[k] * soh_577[k];

        t_767[k] = f_19 * snh_426[k]
                   + f_3 * pc_z[k] * soh_573[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, pc_x, pc_y, snh_579, snh_581, sog0_417, \
                         sog0_419, sog1_417, sog1_419, soh_576, soh_579, \
                         soh_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_20 * snh_579[k]
                   + f_8 * sog0_417[k]
                   - f_9 * sog1_417[k]
                   + f_3 * pc_x[k] * soh_579[k];

        t_769[k] = f_3 * pc_y[k] * soh_576[k];

        t_770[k] = f_20 * snh_581[k]
                   + f_8 * sog0_419[k]
                   - f_9 * sog1_419[k]
                   + f_3 * pc_x[k] * soh_581[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, t_775, pc_x, snh_582, snh_583, snh_584, \
                         snh_585, snh_586, soh_582, soh_583, soh_584, soh_585, \
                         soh_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_20 * snh_582[k]
                   + f_3 * pc_x[k] * soh_582[k];

        t_772[k] = f_20 * snh_583[k]
                   + f_3 * pc_x[k] * soh_583[k];

        t_773[k] = f_20 * snh_584[k]
                   + f_3 * pc_x[k] * soh_584[k];

        t_774[k] = f_20 * snh_585[k]
                   + f_3 * pc_x[k] * soh_585[k];

        t_775[k] = f_20 * snh_586[k]
                   + f_3 * pc_x[k] * soh_586[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pc_x, pc_y, pc_z, snh_435, snh_587, \
                         sog0_415, sog0_417, sog1_415, sog1_417, soh_582, soh_584, \
                         soh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_20 * snh_587[k]
                   + f_3 * pc_x[k] * soh_587[k];

        t_777[k] = f_1 * sog0_415[k]
                   - f_2 * sog1_415[k]
                   + f_3 * pc_y[k] * soh_582[k];

        t_778[k] = f_19 * snh_435[k]
                   + f_3 * pc_z[k] * soh_582[k];

        t_779[k] = f_4 * sog0_417[k]
                   - f_5 * sog1_417[k]
                   + f_3 * pc_y[k] * soh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, snh_440, sog0_418, sog0_419, \
                         sog1_418, sog1_419, soh_585, soh_586, \
                         soh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * sog0_418[k]
                   - f_7 * sog1_418[k]
                   + f_3 * pc_y[k] * soh_585[k];

        t_781[k] = f_8 * sog0_419[k]
                   - f_9 * sog1_419[k]
                   + f_3 * pc_y[k] * soh_586[k];

        t_782[k] = f_3 * pc_y[k] * soh_587[k];

        t_783[k] = f_19 * snh_440[k]
                   + f_1 * sog0_419[k]
                   - f_2 * sog1_419[k]
                   + f_3 * pc_z[k] * soh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, pc_y, pc_z, snh_441, snh_588, \
                         snh_591, sog0_420, sog0_423, sog1_420, sog1_423, soh_588, \
                         soh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_14 * snh_588[k]
                   + f_1 * sog0_420[k]
                   - f_2 * sog1_420[k]
                   + f_3 * pc_x[k] * soh_588[k];

        t_785[k] = f_18 * snh_441[k]
                   + f_3 * pc_y[k] * soh_588[k];

        t_786[k] = f_3 * pc_z[k] * soh_588[k];

        t_787[k] = f_14 * snh_591[k]
                   + f_4 * sog0_423[k]
                   - f_5 * sog1_423[k]
                   + f_3 * pc_x[k] * soh_591[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, pc_y, snh_443, snh_593, snh_594, sog0_425, \
                         sog0_426, sog1_425, sog1_426, soh_590, soh_593, \
                         soh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_18 * snh_443[k]
                   + f_3 * pc_y[k] * soh_590[k];

        t_789[k] = f_14 * snh_593[k]
                   + f_4 * sog0_425[k]
                   - f_5 * sog1_425[k]
                   + f_3 * pc_x[k] * soh_593[k];

        t_790[k] = f_14 * snh_594[k]
                   + f_6 * sog0_426[k]
                   - f_7 * sog1_426[k]
                   + f_3 * pc_x[k] * soh_594[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, pc_x, pc_y, pc_z, snh_446, snh_597, sog0_429, \
                         sog1_429, soh_591, soh_593, soh_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_3 * pc_z[k] * soh_591[k];

        t_792[k] = f_18 * snh_446[k]
                   + f_3 * pc_y[k] * soh_593[k];

        t_793[k] = f_14 * snh_597[k]
                   + f_6 * sog0_429[k]
                   - f_7 * sog1_429[k]
                   + f_3 * pc_x[k] * soh_597[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pc_x, pc_z, snh_598, snh_600, sog0_430, \
                         sog0_432, sog1_430, sog1_432, soh_594, soh_598, \
                         soh_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_14 * snh_598[k]
                   + f_8 * sog0_430[k]
                   - f_9 * sog1_430[k]
                   + f_3 * pc_x[k] * soh_598[k];

        t_795[k] = f_3 * pc_z[k] * soh_594[k];

        t_796[k] = f_14 * snh_600[k]
                   + f_8 * sog0_432[k]
                   - f_9 * sog1_432[k]
                   + f_3 * pc_x[k] * soh_600[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_588 = buffer.data(sni0 + 588);
    const auto *sni0_591 = buffer.data(sni0 + 591);
    const auto *sni0_594 = buffer.data(sni0 + 594);
    const auto *sni0_598 = buffer.data(sni0 + 598);
    const auto *sni0_600 = buffer.data(sni0 + 600);
    const auto *sni0_609 = buffer.data(sni0 + 609);

    const auto *snh_441 = buffer.data(snh + 441);
    const auto *snh_444 = buffer.data(snh + 444);
    const auto *snh_447 = buffer.data(snh + 447);
    const auto *snh_448 = buffer.data(snh + 448);
    const auto *snh_450 = buffer.data(snh + 450);
    const auto *snh_456 = buffer.data(snh + 456);
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
    const auto *snh_477 = buffer.data(snh + 477);
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
    const auto *snh_498 = buffer.data(snh + 498);
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
    const auto *snh_519 = buffer.data(snh + 519);
    const auto *snh_521 = buffer.data(snh + 521);
    const auto *snh_522 = buffer.data(snh + 522);
    const auto *snh_523 = buffer.data(snh + 523);
    const auto *snh_524 = buffer.data(snh + 524);
    const auto *snh_525 = buffer.data(snh + 525);
    const auto *snh_527 = buffer.data(snh + 527);
    const auto *snh_530 = buffer.data(snh + 530);
    const auto *snh_602 = buffer.data(snh + 602);
    const auto *snh_603 = buffer.data(snh + 603);
    const auto *snh_604 = buffer.data(snh + 604);
    const auto *snh_605 = buffer.data(snh + 605);
    const auto *snh_606 = buffer.data(snh + 606);
    const auto *snh_607 = buffer.data(snh + 607);
    const auto *snh_608 = buffer.data(snh + 608);
    const auto *snh_614 = buffer.data(snh + 614);
    const auto *snh_618 = buffer.data(snh + 618);
    const auto *snh_623 = buffer.data(snh + 623);
    const auto *snh_624 = buffer.data(snh + 624);
    const auto *snh_625 = buffer.data(snh + 625);
    const auto *snh_626 = buffer.data(snh + 626);
    const auto *snh_627 = buffer.data(snh + 627);
    const auto *snh_628 = buffer.data(snh + 628);
    const auto *snh_629 = buffer.data(snh + 629);
    const auto *snh_630 = buffer.data(snh + 630);
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
    const auto *snh_675 = buffer.data(snh + 675);
    const auto *snh_677 = buffer.data(snh + 677);
    const auto *snh_678 = buffer.data(snh + 678);
    const auto *snh_681 = buffer.data(snh + 681);
    const auto *snh_682 = buffer.data(snh + 682);

    const auto *sni1_588 = buffer.data(sni1 + 588);
    const auto *sni1_591 = buffer.data(sni1 + 591);
    const auto *sni1_594 = buffer.data(sni1 + 594);
    const auto *sni1_598 = buffer.data(sni1 + 598);
    const auto *sni1_600 = buffer.data(sni1 + 600);
    const auto *sni1_609 = buffer.data(sni1 + 609);

    const auto *sog0_430 = buffer.data(sog0 + 430);
    const auto *sog0_432 = buffer.data(sog0 + 432);
    const auto *sog0_433 = buffer.data(sog0 + 433);
    const auto *sog0_434 = buffer.data(sog0 + 434);
    const auto *sog0_440 = buffer.data(sog0 + 440);
    const auto *sog0_444 = buffer.data(sog0 + 444);
    const auto *sog0_447 = buffer.data(sog0 + 447);
    const auto *sog0_448 = buffer.data(sog0 + 448);
    const auto *sog0_449 = buffer.data(sog0 + 449);
    const auto *sog0_450 = buffer.data(sog0 + 450);
    const auto *sog0_453 = buffer.data(sog0 + 453);
    const auto *sog0_455 = buffer.data(sog0 + 455);
    const auto *sog0_456 = buffer.data(sog0 + 456);
    const auto *sog0_459 = buffer.data(sog0 + 459);
    const auto *sog0_460 = buffer.data(sog0 + 460);
    const auto *sog0_462 = buffer.data(sog0 + 462);
    const auto *sog0_463 = buffer.data(sog0 + 463);
    const auto *sog0_464 = buffer.data(sog0 + 464);
    const auto *sog0_465 = buffer.data(sog0 + 465);
    const auto *sog0_468 = buffer.data(sog0 + 468);
    const auto *sog0_470 = buffer.data(sog0 + 470);
    const auto *sog0_471 = buffer.data(sog0 + 471);
    const auto *sog0_474 = buffer.data(sog0 + 474);
    const auto *sog0_475 = buffer.data(sog0 + 475);
    const auto *sog0_477 = buffer.data(sog0 + 477);
    const auto *sog0_478 = buffer.data(sog0 + 478);
    const auto *sog0_479 = buffer.data(sog0 + 479);
    const auto *sog0_480 = buffer.data(sog0 + 480);
    const auto *sog0_483 = buffer.data(sog0 + 483);
    const auto *sog0_485 = buffer.data(sog0 + 485);
    const auto *sog0_486 = buffer.data(sog0 + 486);
    const auto *sog0_489 = buffer.data(sog0 + 489);
    const auto *sog0_490 = buffer.data(sog0 + 490);

    const auto *sog1_430 = buffer.data(sog1 + 430);
    const auto *sog1_432 = buffer.data(sog1 + 432);
    const auto *sog1_433 = buffer.data(sog1 + 433);
    const auto *sog1_434 = buffer.data(sog1 + 434);
    const auto *sog1_440 = buffer.data(sog1 + 440);
    const auto *sog1_444 = buffer.data(sog1 + 444);
    const auto *sog1_447 = buffer.data(sog1 + 447);
    const auto *sog1_448 = buffer.data(sog1 + 448);
    const auto *sog1_449 = buffer.data(sog1 + 449);
    const auto *sog1_450 = buffer.data(sog1 + 450);
    const auto *sog1_453 = buffer.data(sog1 + 453);
    const auto *sog1_455 = buffer.data(sog1 + 455);
    const auto *sog1_456 = buffer.data(sog1 + 456);
    const auto *sog1_459 = buffer.data(sog1 + 459);
    const auto *sog1_460 = buffer.data(sog1 + 460);
    const auto *sog1_462 = buffer.data(sog1 + 462);
    const auto *sog1_463 = buffer.data(sog1 + 463);
    const auto *sog1_464 = buffer.data(sog1 + 464);
    const auto *sog1_465 = buffer.data(sog1 + 465);
    const auto *sog1_468 = buffer.data(sog1 + 468);
    const auto *sog1_470 = buffer.data(sog1 + 470);
    const auto *sog1_471 = buffer.data(sog1 + 471);
    const auto *sog1_474 = buffer.data(sog1 + 474);
    const auto *sog1_475 = buffer.data(sog1 + 475);
    const auto *sog1_477 = buffer.data(sog1 + 477);
    const auto *sog1_478 = buffer.data(sog1 + 478);
    const auto *sog1_479 = buffer.data(sog1 + 479);
    const auto *sog1_480 = buffer.data(sog1 + 480);
    const auto *sog1_483 = buffer.data(sog1 + 483);
    const auto *sog1_485 = buffer.data(sog1 + 485);
    const auto *sog1_486 = buffer.data(sog1 + 486);
    const auto *sog1_489 = buffer.data(sog1 + 489);
    const auto *sog1_490 = buffer.data(sog1 + 490);

    const auto *soh_597 = buffer.data(soh + 597);
    const auto *soh_602 = buffer.data(soh + 602);
    const auto *soh_603 = buffer.data(soh + 603);
    const auto *soh_604 = buffer.data(soh + 604);
    const auto *soh_605 = buffer.data(soh + 605);
    const auto *soh_606 = buffer.data(soh + 606);
    const auto *soh_607 = buffer.data(soh + 607);
    const auto *soh_608 = buffer.data(soh + 608);
    const auto *soh_609 = buffer.data(soh + 609);
    const auto *soh_611 = buffer.data(soh + 611);
    const auto *soh_612 = buffer.data(soh + 612);
    const auto *soh_614 = buffer.data(soh + 614);
    const auto *soh_615 = buffer.data(soh + 615);
    const auto *soh_618 = buffer.data(soh + 618);
    const auto *soh_623 = buffer.data(soh + 623);
    const auto *soh_624 = buffer.data(soh + 624);
    const auto *soh_625 = buffer.data(soh + 625);
    const auto *soh_626 = buffer.data(soh + 626);
    const auto *soh_627 = buffer.data(soh + 627);
    const auto *soh_628 = buffer.data(soh + 628);
    const auto *soh_629 = buffer.data(soh + 629);
    const auto *soh_630 = buffer.data(soh + 630);
    const auto *soh_632 = buffer.data(soh + 632);
    const auto *soh_633 = buffer.data(soh + 633);
    const auto *soh_635 = buffer.data(soh + 635);
    const auto *soh_636 = buffer.data(soh + 636);
    const auto *soh_639 = buffer.data(soh + 639);
    const auto *soh_640 = buffer.data(soh + 640);
    const auto *soh_642 = buffer.data(soh + 642);
    const auto *soh_644 = buffer.data(soh + 644);
    const auto *soh_645 = buffer.data(soh + 645);
    const auto *soh_646 = buffer.data(soh + 646);
    const auto *soh_647 = buffer.data(soh + 647);
    const auto *soh_648 = buffer.data(soh + 648);
    const auto *soh_649 = buffer.data(soh + 649);
    const auto *soh_650 = buffer.data(soh + 650);
    const auto *soh_651 = buffer.data(soh + 651);
    const auto *soh_653 = buffer.data(soh + 653);
    const auto *soh_654 = buffer.data(soh + 654);
    const auto *soh_656 = buffer.data(soh + 656);
    const auto *soh_657 = buffer.data(soh + 657);
    const auto *soh_660 = buffer.data(soh + 660);
    const auto *soh_661 = buffer.data(soh + 661);
    const auto *soh_663 = buffer.data(soh + 663);
    const auto *soh_665 = buffer.data(soh + 665);
    const auto *soh_666 = buffer.data(soh + 666);
    const auto *soh_667 = buffer.data(soh + 667);
    const auto *soh_668 = buffer.data(soh + 668);
    const auto *soh_669 = buffer.data(soh + 669);
    const auto *soh_670 = buffer.data(soh + 670);
    const auto *soh_671 = buffer.data(soh + 671);
    const auto *soh_672 = buffer.data(soh + 672);
    const auto *soh_674 = buffer.data(soh + 674);
    const auto *soh_675 = buffer.data(soh + 675);
    const auto *soh_677 = buffer.data(soh + 677);
    const auto *soh_678 = buffer.data(soh + 678);
    const auto *soh_681 = buffer.data(soh + 681);
    const auto *soh_682 = buffer.data(soh + 682);

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pc_x, pc_y, snh_450, snh_602, snh_603, \
                         snh_604, sog0_434, sog1_434, soh_597, soh_602, soh_603, \
                         soh_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_18 * snh_450[k]
                   + f_3 * pc_y[k] * soh_597[k];

        t_798[k] = f_14 * snh_602[k]
                   + f_8 * sog0_434[k]
                   - f_9 * sog1_434[k]
                   + f_3 * pc_x[k] * soh_602[k];

        t_799[k] = f_14 * snh_603[k]
                   + f_3 * pc_x[k] * soh_603[k];

        t_800[k] = f_14 * snh_604[k]
                   + f_3 * pc_x[k] * soh_604[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pc_x, snh_605, snh_606, snh_607, snh_608, \
                         soh_605, soh_606, soh_607, soh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_14 * snh_605[k]
                   + f_3 * pc_x[k] * soh_605[k];

        t_802[k] = f_14 * snh_606[k]
                   + f_3 * pc_x[k] * soh_606[k];

        t_803[k] = f_14 * snh_607[k]
                   + f_3 * pc_x[k] * soh_607[k];

        t_804[k] = f_14 * snh_608[k]
                   + f_3 * pc_x[k] * soh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pc_y, pc_z, snh_456, snh_458, sog0_430, \
                         sog0_432, sog1_430, sog1_432, soh_603, \
                         soh_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_18 * snh_456[k]
                   + f_1 * sog0_430[k]
                   - f_2 * sog1_430[k]
                   + f_3 * pc_y[k] * soh_603[k];

        t_806[k] = f_3 * pc_z[k] * soh_603[k];

        t_807[k] = f_18 * snh_458[k]
                   + f_4 * sog0_432[k]
                   - f_5 * sog1_432[k]
                   + f_3 * pc_y[k] * soh_605[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pc_y, pc_z, snh_459, snh_460, snh_461, \
                         sog0_433, sog0_434, sog1_433, sog1_434, soh_606, soh_607, \
                         soh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_18 * snh_459[k]
                   + f_6 * sog0_433[k]
                   - f_7 * sog1_433[k]
                   + f_3 * pc_y[k] * soh_606[k];

        t_809[k] = f_18 * snh_460[k]
                   + f_8 * sog0_434[k]
                   - f_9 * sog1_434[k]
                   + f_3 * pc_y[k] * soh_607[k];

        t_810[k] = f_18 * snh_461[k]
                   + f_3 * pc_y[k] * soh_608[k];

        t_811[k] = f_1 * sog0_434[k]
                   - f_2 * sog1_434[k]
                   + f_3 * pc_z[k] * soh_608[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pb_z, pc_y, pc_z, sni0_588, sni0_591, \
                         snh_441, snh_462, sni1_588, sni1_591, \
                         soh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pb_z[k] * sni0_588[k]
                   - f_10 * pc_z[k] * sni1_588[k];

        t_813[k] = f_19 * snh_462[k]
                   + f_3 * pc_y[k] * soh_609[k];

        t_814[k] = f_11 * snh_441[k]
                   + f_3 * pc_z[k] * soh_609[k];

        t_815[k] = pb_z[k] * sni0_591[k]
                   - f_10 * pc_z[k] * sni1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pb_z, pc_x, pc_y, pc_z, sni0_594, snh_464, \
                         snh_614, sni1_594, sog0_440, sog1_440, soh_611, \
                         soh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_19 * snh_464[k]
                   + f_3 * pc_y[k] * soh_611[k];

        t_817[k] = f_14 * snh_614[k]
                   + f_4 * sog0_440[k]
                   - f_5 * sog1_440[k]
                   + f_3 * pc_x[k] * soh_614[k];

        t_818[k] = pb_z[k] * sni0_594[k]
                   - f_10 * pc_z[k] * sni1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_y, pc_z, snh_444, snh_467, snh_618, \
                         sog0_444, sog1_444, soh_612, soh_614, \
                         soh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * snh_444[k]
                   + f_3 * pc_z[k] * soh_612[k];

        t_820[k] = f_19 * snh_467[k]
                   + f_3 * pc_y[k] * soh_614[k];

        t_821[k] = f_14 * snh_618[k]
                   + f_6 * sog0_444[k]
                   - f_7 * sog1_444[k]
                   + f_3 * pc_x[k] * soh_618[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pb_z, pc_y, pc_z, sni0_598, sni0_600, \
                         snh_447, snh_448, snh_471, sni1_598, sni1_600, soh_615, \
                         soh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pb_z[k] * sni0_598[k]
                   - f_10 * pc_z[k] * sni1_598[k];

        t_823[k] = f_11 * snh_447[k]
                   + f_3 * pc_z[k] * soh_615[k];

        t_824[k] = pb_z[k] * sni0_600[k]
                   + f_12 * snh_448[k]
                   - f_10 * pc_z[k] * sni1_600[k];

        t_825[k] = f_19 * snh_471[k]
                   + f_3 * pc_y[k] * soh_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pc_x, snh_623, snh_624, snh_625, snh_626, \
                         sog0_449, sog1_449, soh_623, soh_624, soh_625, \
                         soh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_14 * snh_623[k]
                   + f_8 * sog0_449[k]
                   - f_9 * sog1_449[k]
                   + f_3 * pc_x[k] * soh_623[k];

        t_827[k] = f_14 * snh_624[k]
                   + f_3 * pc_x[k] * soh_624[k];

        t_828[k] = f_14 * snh_625[k]
                   + f_3 * pc_x[k] * soh_625[k];

        t_829[k] = f_14 * snh_626[k]
                   + f_3 * pc_x[k] * soh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pb_z, pc_x, pc_z, sni0_609, snh_627, \
                         snh_628, snh_629, sni1_609, soh_627, soh_628, \
                         soh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_14 * snh_627[k]
                   + f_3 * pc_x[k] * soh_627[k];

        t_831[k] = f_14 * snh_628[k]
                   + f_3 * pc_x[k] * soh_628[k];

        t_832[k] = f_14 * snh_629[k]
                   + f_3 * pc_x[k] * soh_629[k];

        t_833[k] = pb_z[k] * sni0_609[k]
                   - f_10 * pc_z[k] * sni1_609[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pc_z, snh_456, snh_479, snh_480, sog0_447, \
                         sog0_448, sog1_447, sog1_448, soh_624, soh_626, \
                         soh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_11 * snh_456[k]
                   + f_3 * pc_z[k] * soh_624[k];

        t_835[k] = f_19 * snh_479[k]
                   + f_4 * sog0_447[k]
                   - f_5 * sog1_447[k]
                   + f_3 * pc_y[k] * soh_626[k];

        t_836[k] = f_19 * snh_480[k]
                   + f_6 * sog0_448[k]
                   - f_7 * sog1_448[k]
                   + f_3 * pc_y[k] * soh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, snh_461, snh_481, snh_482, sog0_449, \
                         sog1_449, soh_628, soh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_19 * snh_481[k]
                   + f_8 * sog0_449[k]
                   - f_9 * sog1_449[k]
                   + f_3 * pc_y[k] * soh_628[k];

        t_838[k] = f_19 * snh_482[k]
                   + f_3 * pc_y[k] * soh_629[k];

        t_839[k] = f_11 * snh_461[k]
                   + f_1 * sog0_449[k]
                   - f_2 * sog1_449[k]
                   + f_3 * pc_z[k] * soh_629[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, pc_x, pc_y, pc_z, snh_462, snh_483, snh_630, \
                         sog0_450, sog1_450, soh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_14 * snh_630[k]
                   + f_1 * sog0_450[k]
                   - f_2 * sog1_450[k]
                   + f_3 * pc_x[k] * soh_630[k];

        t_841[k] = f_20 * snh_483[k]
                   + f_3 * pc_y[k] * soh_630[k];

        t_842[k] = f_12 * snh_462[k]
                   + f_3 * pc_z[k] * soh_630[k];
    }

#pragma omp simd aligned(t_843, t_844, t_845, pc_x, pc_y, snh_485, snh_633, snh_635, sog0_453, \
                         sog0_455, sog1_453, sog1_455, soh_632, soh_633, \
                         soh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_843[k] = f_14 * snh_633[k]
                   + f_4 * sog0_453[k]
                   - f_5 * sog1_453[k]
                   + f_3 * pc_x[k] * soh_633[k];

        t_844[k] = f_20 * snh_485[k]
                   + f_3 * pc_y[k] * soh_632[k];

        t_845[k] = f_14 * snh_635[k]
                   + f_4 * sog0_455[k]
                   - f_5 * sog1_455[k]
                   + f_3 * pc_x[k] * soh_635[k];
    }

#pragma omp simd aligned(t_846, t_847, t_848, pc_x, pc_y, pc_z, snh_465, snh_488, snh_636, \
                         sog0_456, sog1_456, soh_633, soh_635, \
                         soh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_846[k] = f_14 * snh_636[k]
                   + f_6 * sog0_456[k]
                   - f_7 * sog1_456[k]
                   + f_3 * pc_x[k] * soh_636[k];

        t_847[k] = f_12 * snh_465[k]
                   + f_3 * pc_z[k] * soh_633[k];

        t_848[k] = f_20 * snh_488[k]
                   + f_3 * pc_y[k] * soh_635[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, pc_x, pc_z, snh_468, snh_639, snh_640, sog0_459, \
                         sog0_460, sog1_459, sog1_460, soh_636, soh_639, \
                         soh_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_14 * snh_639[k]
                   + f_6 * sog0_459[k]
                   - f_7 * sog1_459[k]
                   + f_3 * pc_x[k] * soh_639[k];

        t_850[k] = f_14 * snh_640[k]
                   + f_8 * sog0_460[k]
                   - f_9 * sog1_460[k]
                   + f_3 * pc_x[k] * soh_640[k];

        t_851[k] = f_12 * snh_468[k]
                   + f_3 * pc_z[k] * soh_636[k];
    }

#pragma omp simd aligned(t_852, t_853, t_854, pc_x, pc_y, snh_492, snh_642, snh_644, sog0_462, \
                         sog0_464, sog1_462, sog1_464, soh_639, soh_642, \
                         soh_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_852[k] = f_14 * snh_642[k]
                   + f_8 * sog0_462[k]
                   - f_9 * sog1_462[k]
                   + f_3 * pc_x[k] * soh_642[k];

        t_853[k] = f_20 * snh_492[k]
                   + f_3 * pc_y[k] * soh_639[k];

        t_854[k] = f_14 * snh_644[k]
                   + f_8 * sog0_464[k]
                   - f_9 * sog1_464[k]
                   + f_3 * pc_x[k] * soh_644[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, t_858, t_859, pc_x, snh_645, snh_646, snh_647, \
                         snh_648, snh_649, soh_645, soh_646, soh_647, soh_648, \
                         soh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_14 * snh_645[k]
                   + f_3 * pc_x[k] * soh_645[k];

        t_856[k] = f_14 * snh_646[k]
                   + f_3 * pc_x[k] * soh_646[k];

        t_857[k] = f_14 * snh_647[k]
                   + f_3 * pc_x[k] * soh_647[k];

        t_858[k] = f_14 * snh_648[k]
                   + f_3 * pc_x[k] * soh_648[k];

        t_859[k] = f_14 * snh_649[k]
                   + f_3 * pc_x[k] * soh_649[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_y, pc_z, snh_477, snh_498, snh_650, \
                         sog0_460, sog1_460, soh_645, soh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_14 * snh_650[k]
                   + f_3 * pc_x[k] * soh_650[k];

        t_861[k] = f_20 * snh_498[k]
                   + f_1 * sog0_460[k]
                   - f_2 * sog1_460[k]
                   + f_3 * pc_y[k] * soh_645[k];

        t_862[k] = f_12 * snh_477[k]
                   + f_3 * pc_z[k] * soh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, snh_500, snh_501, snh_502, sog0_462, \
                         sog0_463, sog0_464, sog1_462, sog1_463, sog1_464, soh_647, soh_648, \
                         soh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_20 * snh_500[k]
                   + f_4 * sog0_462[k]
                   - f_5 * sog1_462[k]
                   + f_3 * pc_y[k] * soh_647[k];

        t_864[k] = f_20 * snh_501[k]
                   + f_6 * sog0_463[k]
                   - f_7 * sog1_463[k]
                   + f_3 * pc_y[k] * soh_648[k];

        t_865[k] = f_20 * snh_502[k]
                   + f_8 * sog0_464[k]
                   - f_9 * sog1_464[k]
                   + f_3 * pc_y[k] * soh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, snh_482, snh_503, snh_651, \
                         sog0_464, sog0_465, sog1_464, sog1_465, soh_650, \
                         soh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_20 * snh_503[k]
                   + f_3 * pc_y[k] * soh_650[k];

        t_867[k] = f_12 * snh_482[k]
                   + f_1 * sog0_464[k]
                   - f_2 * sog1_464[k]
                   + f_3 * pc_z[k] * soh_650[k];

        t_868[k] = f_14 * snh_651[k]
                   + f_1 * sog0_465[k]
                   - f_2 * sog1_465[k]
                   + f_3 * pc_x[k] * soh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pc_x, pc_y, pc_z, snh_483, snh_504, \
                         snh_506, snh_654, sog0_468, sog1_468, soh_651, soh_653, \
                         soh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_14 * snh_504[k]
                   + f_3 * pc_y[k] * soh_651[k];

        t_870[k] = f_13 * snh_483[k]
                   + f_3 * pc_z[k] * soh_651[k];

        t_871[k] = f_14 * snh_654[k]
                   + f_4 * sog0_468[k]
                   - f_5 * sog1_468[k]
                   + f_3 * pc_x[k] * soh_654[k];

        t_872[k] = f_14 * snh_506[k]
                   + f_3 * pc_y[k] * soh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pc_x, pc_z, snh_486, snh_656, snh_657, sog0_470, \
                         sog0_471, sog1_470, sog1_471, soh_654, soh_656, \
                         soh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_14 * snh_656[k]
                   + f_4 * sog0_470[k]
                   - f_5 * sog1_470[k]
                   + f_3 * pc_x[k] * soh_656[k];

        t_874[k] = f_14 * snh_657[k]
                   + f_6 * sog0_471[k]
                   - f_7 * sog1_471[k]
                   + f_3 * pc_x[k] * soh_657[k];

        t_875[k] = f_13 * snh_486[k]
                   + f_3 * pc_z[k] * soh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pc_x, pc_y, snh_509, snh_660, snh_661, sog0_474, \
                         sog0_475, sog1_474, sog1_475, soh_656, soh_660, \
                         soh_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * snh_509[k]
                   + f_3 * pc_y[k] * soh_656[k];

        t_877[k] = f_14 * snh_660[k]
                   + f_6 * sog0_474[k]
                   - f_7 * sog1_474[k]
                   + f_3 * pc_x[k] * soh_660[k];

        t_878[k] = f_14 * snh_661[k]
                   + f_8 * sog0_475[k]
                   - f_9 * sog1_475[k]
                   + f_3 * pc_x[k] * soh_661[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pc_x, pc_y, pc_z, snh_489, snh_513, snh_663, \
                         sog0_477, sog1_477, soh_657, soh_660, \
                         soh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * snh_489[k]
                   + f_3 * pc_z[k] * soh_657[k];

        t_880[k] = f_14 * snh_663[k]
                   + f_8 * sog0_477[k]
                   - f_9 * sog1_477[k]
                   + f_3 * pc_x[k] * soh_663[k];

        t_881[k] = f_14 * snh_513[k]
                   + f_3 * pc_y[k] * soh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pc_x, snh_665, snh_666, snh_667, snh_668, \
                         sog0_479, sog1_479, soh_665, soh_666, soh_667, \
                         soh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = f_14 * snh_665[k]
                   + f_8 * sog0_479[k]
                   - f_9 * sog1_479[k]
                   + f_3 * pc_x[k] * soh_665[k];

        t_883[k] = f_14 * snh_666[k]
                   + f_3 * pc_x[k] * soh_666[k];

        t_884[k] = f_14 * snh_667[k]
                   + f_3 * pc_x[k] * soh_667[k];

        t_885[k] = f_14 * snh_668[k]
                   + f_3 * pc_x[k] * soh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pc_x, pc_y, snh_519, snh_669, snh_670, \
                         snh_671, sog0_475, sog1_475, soh_666, soh_669, soh_670, \
                         soh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_14 * snh_669[k]
                   + f_3 * pc_x[k] * soh_669[k];

        t_887[k] = f_14 * snh_670[k]
                   + f_3 * pc_x[k] * soh_670[k];

        t_888[k] = f_14 * snh_671[k]
                   + f_3 * pc_x[k] * soh_671[k];

        t_889[k] = f_14 * snh_519[k]
                   + f_1 * sog0_475[k]
                   - f_2 * sog1_475[k]
                   + f_3 * pc_y[k] * soh_666[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, pc_y, pc_z, snh_498, snh_521, snh_522, sog0_477, \
                         sog0_478, sog1_477, sog1_478, soh_666, soh_668, \
                         soh_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * snh_498[k]
                   + f_3 * pc_z[k] * soh_666[k];

        t_891[k] = f_14 * snh_521[k]
                   + f_4 * sog0_477[k]
                   - f_5 * sog1_477[k]
                   + f_3 * pc_y[k] * soh_668[k];

        t_892[k] = f_14 * snh_522[k]
                   + f_6 * sog0_478[k]
                   - f_7 * sog1_478[k]
                   + f_3 * pc_y[k] * soh_669[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, pc_z, snh_503, snh_523, snh_524, sog0_479, \
                         sog1_479, soh_670, soh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_14 * snh_523[k]
                   + f_8 * sog0_479[k]
                   - f_9 * sog1_479[k]
                   + f_3 * pc_y[k] * soh_670[k];

        t_894[k] = f_14 * snh_524[k]
                   + f_3 * pc_y[k] * soh_671[k];

        t_895[k] = f_13 * snh_503[k]
                   + f_1 * sog0_479[k]
                   - f_2 * sog1_479[k]
                   + f_3 * pc_z[k] * soh_671[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_x, pc_y, pc_z, snh_504, snh_525, snh_672, \
                         sog0_480, sog1_480, soh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_14 * snh_672[k]
                   + f_1 * sog0_480[k]
                   - f_2 * sog1_480[k]
                   + f_3 * pc_x[k] * soh_672[k];

        t_897[k] = f_13 * snh_525[k]
                   + f_3 * pc_y[k] * soh_672[k];

        t_898[k] = f_14 * snh_504[k]
                   + f_3 * pc_z[k] * soh_672[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_x, pc_y, snh_527, snh_675, snh_677, sog0_483, \
                         sog0_485, sog1_483, sog1_485, soh_674, soh_675, \
                         soh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_14 * snh_675[k]
                   + f_4 * sog0_483[k]
                   - f_5 * sog1_483[k]
                   + f_3 * pc_x[k] * soh_675[k];

        t_900[k] = f_13 * snh_527[k]
                   + f_3 * pc_y[k] * soh_674[k];

        t_901[k] = f_14 * snh_677[k]
                   + f_4 * sog0_485[k]
                   - f_5 * sog1_485[k]
                   + f_3 * pc_x[k] * soh_677[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pc_x, pc_y, pc_z, snh_507, snh_530, snh_678, \
                         sog0_486, sog1_486, soh_675, soh_677, \
                         soh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_14 * snh_678[k]
                   + f_6 * sog0_486[k]
                   - f_7 * sog1_486[k]
                   + f_3 * pc_x[k] * soh_678[k];

        t_903[k] = f_14 * snh_507[k]
                   + f_3 * pc_z[k] * soh_675[k];

        t_904[k] = f_13 * snh_530[k]
                   + f_3 * pc_y[k] * soh_677[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pc_x, pc_z, snh_510, snh_681, snh_682, sog0_489, \
                         sog0_490, sog1_489, sog1_490, soh_678, soh_681, \
                         soh_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = f_14 * snh_681[k]
                   + f_6 * sog0_489[k]
                   - f_7 * sog1_489[k]
                   + f_3 * pc_x[k] * soh_681[k];

        t_906[k] = f_14 * snh_682[k]
                   + f_8 * sog0_490[k]
                   - f_9 * sog1_490[k]
                   + f_3 * pc_x[k] * soh_682[k];

        t_907[k] = f_14 * snh_510[k]
                   + f_3 * pc_z[k] * soh_678[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_756 = buffer.data(sni0 + 756);
    const auto *sni0_759 = buffer.data(sni0 + 759);
    const auto *sni0_761 = buffer.data(sni0 + 761);
    const auto *sni0_762 = buffer.data(sni0 + 762);
    const auto *sni0_765 = buffer.data(sni0 + 765);
    const auto *sni0_766 = buffer.data(sni0 + 766);
    const auto *sni0_768 = buffer.data(sni0 + 768);
    const auto *sni0_770 = buffer.data(sni0 + 770);
    const auto *sni0_783 = buffer.data(sni0 + 783);

    const auto *snh_519 = buffer.data(snh + 519);
    const auto *snh_524 = buffer.data(snh + 524);
    const auto *snh_525 = buffer.data(snh + 525);
    const auto *snh_528 = buffer.data(snh + 528);
    const auto *snh_531 = buffer.data(snh + 531);
    const auto *snh_534 = buffer.data(snh + 534);
    const auto *snh_540 = buffer.data(snh + 540);
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
    const auto *snh_563 = buffer.data(snh + 563);
    const auto *snh_564 = buffer.data(snh + 564);
    const auto *snh_565 = buffer.data(snh + 565);
    const auto *snh_566 = buffer.data(snh + 566);
    const auto *snh_567 = buffer.data(snh + 567);
    const auto *snh_568 = buffer.data(snh + 568);
    const auto *snh_569 = buffer.data(snh + 569);
    const auto *snh_570 = buffer.data(snh + 570);
    const auto *snh_572 = buffer.data(snh + 572);
    const auto *snh_573 = buffer.data(snh + 573);
    const auto *snh_575 = buffer.data(snh + 575);
    const auto *snh_576 = buffer.data(snh + 576);
    const auto *snh_582 = buffer.data(snh + 582);
    const auto *snh_584 = buffer.data(snh + 584);
    const auto *snh_585 = buffer.data(snh + 585);
    const auto *snh_586 = buffer.data(snh + 586);
    const auto *snh_587 = buffer.data(snh + 587);
    const auto *snh_588 = buffer.data(snh + 588);
    const auto *snh_590 = buffer.data(snh + 590);
    const auto *snh_593 = buffer.data(snh + 593);
    const auto *snh_684 = buffer.data(snh + 684);
    const auto *snh_686 = buffer.data(snh + 686);
    const auto *snh_687 = buffer.data(snh + 687);
    const auto *snh_688 = buffer.data(snh + 688);
    const auto *snh_689 = buffer.data(snh + 689);
    const auto *snh_690 = buffer.data(snh + 690);
    const auto *snh_691 = buffer.data(snh + 691);
    const auto *snh_692 = buffer.data(snh + 692);
    const auto *snh_693 = buffer.data(snh + 693);
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
    const auto *snh_729 = buffer.data(snh + 729);
    const auto *snh_730 = buffer.data(snh + 730);
    const auto *snh_731 = buffer.data(snh + 731);
    const auto *snh_732 = buffer.data(snh + 732);
    const auto *snh_733 = buffer.data(snh + 733);
    const auto *snh_734 = buffer.data(snh + 734);
    const auto *snh_735 = buffer.data(snh + 735);
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
    const auto *snh_759 = buffer.data(snh + 759);
    const auto *snh_761 = buffer.data(snh + 761);
    const auto *snh_762 = buffer.data(snh + 762);
    const auto *snh_765 = buffer.data(snh + 765);
    const auto *snh_766 = buffer.data(snh + 766);
    const auto *snh_768 = buffer.data(snh + 768);

    const auto *sni1_756 = buffer.data(sni1 + 756);
    const auto *sni1_759 = buffer.data(sni1 + 759);
    const auto *sni1_761 = buffer.data(sni1 + 761);
    const auto *sni1_762 = buffer.data(sni1 + 762);
    const auto *sni1_765 = buffer.data(sni1 + 765);
    const auto *sni1_766 = buffer.data(sni1 + 766);
    const auto *sni1_768 = buffer.data(sni1 + 768);
    const auto *sni1_770 = buffer.data(sni1 + 770);
    const auto *sni1_783 = buffer.data(sni1 + 783);

    const auto *sog0_490 = buffer.data(sog0 + 490);
    const auto *sog0_492 = buffer.data(sog0 + 492);
    const auto *sog0_493 = buffer.data(sog0 + 493);
    const auto *sog0_494 = buffer.data(sog0 + 494);
    const auto *sog0_495 = buffer.data(sog0 + 495);
    const auto *sog0_498 = buffer.data(sog0 + 498);
    const auto *sog0_500 = buffer.data(sog0 + 500);
    const auto *sog0_501 = buffer.data(sog0 + 501);
    const auto *sog0_504 = buffer.data(sog0 + 504);
    const auto *sog0_505 = buffer.data(sog0 + 505);
    const auto *sog0_507 = buffer.data(sog0 + 507);
    const auto *sog0_508 = buffer.data(sog0 + 508);
    const auto *sog0_509 = buffer.data(sog0 + 509);
    const auto *sog0_520 = buffer.data(sog0 + 520);
    const auto *sog0_522 = buffer.data(sog0 + 522);
    const auto *sog0_523 = buffer.data(sog0 + 523);
    const auto *sog0_524 = buffer.data(sog0 + 524);
    const auto *sog0_525 = buffer.data(sog0 + 525);
    const auto *sog0_528 = buffer.data(sog0 + 528);
    const auto *sog0_530 = buffer.data(sog0 + 530);
    const auto *sog0_531 = buffer.data(sog0 + 531);
    const auto *sog0_534 = buffer.data(sog0 + 534);
    const auto *sog0_535 = buffer.data(sog0 + 535);
    const auto *sog0_537 = buffer.data(sog0 + 537);
    const auto *sog0_538 = buffer.data(sog0 + 538);
    const auto *sog0_539 = buffer.data(sog0 + 539);
    const auto *sog0_540 = buffer.data(sog0 + 540);
    const auto *sog0_543 = buffer.data(sog0 + 543);
    const auto *sog0_545 = buffer.data(sog0 + 545);
    const auto *sog0_546 = buffer.data(sog0 + 546);
    const auto *sog0_549 = buffer.data(sog0 + 549);
    const auto *sog0_550 = buffer.data(sog0 + 550);
    const auto *sog0_552 = buffer.data(sog0 + 552);

    const auto *sog1_490 = buffer.data(sog1 + 490);
    const auto *sog1_492 = buffer.data(sog1 + 492);
    const auto *sog1_493 = buffer.data(sog1 + 493);
    const auto *sog1_494 = buffer.data(sog1 + 494);
    const auto *sog1_495 = buffer.data(sog1 + 495);
    const auto *sog1_498 = buffer.data(sog1 + 498);
    const auto *sog1_500 = buffer.data(sog1 + 500);
    const auto *sog1_501 = buffer.data(sog1 + 501);
    const auto *sog1_504 = buffer.data(sog1 + 504);
    const auto *sog1_505 = buffer.data(sog1 + 505);
    const auto *sog1_507 = buffer.data(sog1 + 507);
    const auto *sog1_508 = buffer.data(sog1 + 508);
    const auto *sog1_509 = buffer.data(sog1 + 509);
    const auto *sog1_520 = buffer.data(sog1 + 520);
    const auto *sog1_522 = buffer.data(sog1 + 522);
    const auto *sog1_523 = buffer.data(sog1 + 523);
    const auto *sog1_524 = buffer.data(sog1 + 524);
    const auto *sog1_525 = buffer.data(sog1 + 525);
    const auto *sog1_528 = buffer.data(sog1 + 528);
    const auto *sog1_530 = buffer.data(sog1 + 530);
    const auto *sog1_531 = buffer.data(sog1 + 531);
    const auto *sog1_534 = buffer.data(sog1 + 534);
    const auto *sog1_535 = buffer.data(sog1 + 535);
    const auto *sog1_537 = buffer.data(sog1 + 537);
    const auto *sog1_538 = buffer.data(sog1 + 538);
    const auto *sog1_539 = buffer.data(sog1 + 539);
    const auto *sog1_540 = buffer.data(sog1 + 540);
    const auto *sog1_543 = buffer.data(sog1 + 543);
    const auto *sog1_545 = buffer.data(sog1 + 545);
    const auto *sog1_546 = buffer.data(sog1 + 546);
    const auto *sog1_549 = buffer.data(sog1 + 549);
    const auto *sog1_550 = buffer.data(sog1 + 550);
    const auto *sog1_552 = buffer.data(sog1 + 552);

    const auto *soh_681 = buffer.data(soh + 681);
    const auto *soh_684 = buffer.data(soh + 684);
    const auto *soh_686 = buffer.data(soh + 686);
    const auto *soh_687 = buffer.data(soh + 687);
    const auto *soh_688 = buffer.data(soh + 688);
    const auto *soh_689 = buffer.data(soh + 689);
    const auto *soh_690 = buffer.data(soh + 690);
    const auto *soh_691 = buffer.data(soh + 691);
    const auto *soh_692 = buffer.data(soh + 692);
    const auto *soh_693 = buffer.data(soh + 693);
    const auto *soh_695 = buffer.data(soh + 695);
    const auto *soh_696 = buffer.data(soh + 696);
    const auto *soh_698 = buffer.data(soh + 698);
    const auto *soh_699 = buffer.data(soh + 699);
    const auto *soh_702 = buffer.data(soh + 702);
    const auto *soh_703 = buffer.data(soh + 703);
    const auto *soh_705 = buffer.data(soh + 705);
    const auto *soh_707 = buffer.data(soh + 707);
    const auto *soh_708 = buffer.data(soh + 708);
    const auto *soh_709 = buffer.data(soh + 709);
    const auto *soh_710 = buffer.data(soh + 710);
    const auto *soh_711 = buffer.data(soh + 711);
    const auto *soh_712 = buffer.data(soh + 712);
    const auto *soh_713 = buffer.data(soh + 713);
    const auto *soh_714 = buffer.data(soh + 714);
    const auto *soh_716 = buffer.data(soh + 716);
    const auto *soh_717 = buffer.data(soh + 717);
    const auto *soh_719 = buffer.data(soh + 719);
    const auto *soh_720 = buffer.data(soh + 720);
    const auto *soh_723 = buffer.data(soh + 723);
    const auto *soh_729 = buffer.data(soh + 729);
    const auto *soh_730 = buffer.data(soh + 730);
    const auto *soh_731 = buffer.data(soh + 731);
    const auto *soh_732 = buffer.data(soh + 732);
    const auto *soh_733 = buffer.data(soh + 733);
    const auto *soh_734 = buffer.data(soh + 734);
    const auto *soh_735 = buffer.data(soh + 735);
    const auto *soh_737 = buffer.data(soh + 737);
    const auto *soh_738 = buffer.data(soh + 738);
    const auto *soh_740 = buffer.data(soh + 740);
    const auto *soh_741 = buffer.data(soh + 741);
    const auto *soh_744 = buffer.data(soh + 744);
    const auto *soh_745 = buffer.data(soh + 745);
    const auto *soh_747 = buffer.data(soh + 747);
    const auto *soh_749 = buffer.data(soh + 749);
    const auto *soh_750 = buffer.data(soh + 750);
    const auto *soh_751 = buffer.data(soh + 751);
    const auto *soh_752 = buffer.data(soh + 752);
    const auto *soh_753 = buffer.data(soh + 753);
    const auto *soh_754 = buffer.data(soh + 754);
    const auto *soh_755 = buffer.data(soh + 755);
    const auto *soh_756 = buffer.data(soh + 756);
    const auto *soh_758 = buffer.data(soh + 758);
    const auto *soh_759 = buffer.data(soh + 759);
    const auto *soh_761 = buffer.data(soh + 761);
    const auto *soh_762 = buffer.data(soh + 762);
    const auto *soh_765 = buffer.data(soh + 765);
    const auto *soh_766 = buffer.data(soh + 766);
    const auto *soh_768 = buffer.data(soh + 768);

#pragma omp simd aligned(t_908, t_909, t_910, pc_x, pc_y, snh_534, snh_684, snh_686, sog0_492, \
                         sog0_494, sog1_492, sog1_494, soh_681, soh_684, \
                         soh_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_14 * snh_684[k]
                   + f_8 * sog0_492[k]
                   - f_9 * sog1_492[k]
                   + f_3 * pc_x[k] * soh_684[k];

        t_909[k] = f_13 * snh_534[k]
                   + f_3 * pc_y[k] * soh_681[k];

        t_910[k] = f_14 * snh_686[k]
                   + f_8 * sog0_494[k]
                   - f_9 * sog1_494[k]
                   + f_3 * pc_x[k] * soh_686[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pc_x, snh_687, snh_688, snh_689, \
                         snh_690, snh_691, soh_687, soh_688, soh_689, soh_690, \
                         soh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_14 * snh_687[k]
                   + f_3 * pc_x[k] * soh_687[k];

        t_912[k] = f_14 * snh_688[k]
                   + f_3 * pc_x[k] * soh_688[k];

        t_913[k] = f_14 * snh_689[k]
                   + f_3 * pc_x[k] * soh_689[k];

        t_914[k] = f_14 * snh_690[k]
                   + f_3 * pc_x[k] * soh_690[k];

        t_915[k] = f_14 * snh_691[k]
                   + f_3 * pc_x[k] * soh_691[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, pc_x, pc_y, pc_z, snh_519, snh_540, snh_692, \
                         sog0_490, sog1_490, soh_687, soh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_14 * snh_692[k]
                   + f_3 * pc_x[k] * soh_692[k];

        t_917[k] = f_13 * snh_540[k]
                   + f_1 * sog0_490[k]
                   - f_2 * sog1_490[k]
                   + f_3 * pc_y[k] * soh_687[k];

        t_918[k] = f_14 * snh_519[k]
                   + f_3 * pc_z[k] * soh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, snh_542, snh_543, snh_544, sog0_492, \
                         sog0_493, sog0_494, sog1_492, sog1_493, sog1_494, soh_689, soh_690, \
                         soh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * snh_542[k]
                   + f_4 * sog0_492[k]
                   - f_5 * sog1_492[k]
                   + f_3 * pc_y[k] * soh_689[k];

        t_920[k] = f_13 * snh_543[k]
                   + f_6 * sog0_493[k]
                   - f_7 * sog1_493[k]
                   + f_3 * pc_y[k] * soh_690[k];

        t_921[k] = f_13 * snh_544[k]
                   + f_8 * sog0_494[k]
                   - f_9 * sog1_494[k]
                   + f_3 * pc_y[k] * soh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, snh_524, snh_545, snh_693, \
                         sog0_494, sog0_495, sog1_494, sog1_495, soh_692, \
                         soh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * snh_545[k]
                   + f_3 * pc_y[k] * soh_692[k];

        t_923[k] = f_14 * snh_524[k]
                   + f_1 * sog0_494[k]
                   - f_2 * sog1_494[k]
                   + f_3 * pc_z[k] * soh_692[k];

        t_924[k] = f_14 * snh_693[k]
                   + f_1 * sog0_495[k]
                   - f_2 * sog1_495[k]
                   + f_3 * pc_x[k] * soh_693[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, t_928, pc_x, pc_y, pc_z, snh_525, snh_546, \
                         snh_548, snh_696, sog0_498, sog1_498, soh_693, soh_695, \
                         soh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_12 * snh_546[k]
                   + f_3 * pc_y[k] * soh_693[k];

        t_926[k] = f_20 * snh_525[k]
                   + f_3 * pc_z[k] * soh_693[k];

        t_927[k] = f_14 * snh_696[k]
                   + f_4 * sog0_498[k]
                   - f_5 * sog1_498[k]
                   + f_3 * pc_x[k] * soh_696[k];

        t_928[k] = f_12 * snh_548[k]
                   + f_3 * pc_y[k] * soh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pc_x, pc_z, snh_528, snh_698, snh_699, sog0_500, \
                         sog0_501, sog1_500, sog1_501, soh_696, soh_698, \
                         soh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_14 * snh_698[k]
                   + f_4 * sog0_500[k]
                   - f_5 * sog1_500[k]
                   + f_3 * pc_x[k] * soh_698[k];

        t_930[k] = f_14 * snh_699[k]
                   + f_6 * sog0_501[k]
                   - f_7 * sog1_501[k]
                   + f_3 * pc_x[k] * soh_699[k];

        t_931[k] = f_20 * snh_528[k]
                   + f_3 * pc_z[k] * soh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pc_x, pc_y, snh_551, snh_702, snh_703, sog0_504, \
                         sog0_505, sog1_504, sog1_505, soh_698, soh_702, \
                         soh_703 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * snh_551[k]
                   + f_3 * pc_y[k] * soh_698[k];

        t_933[k] = f_14 * snh_702[k]
                   + f_6 * sog0_504[k]
                   - f_7 * sog1_504[k]
                   + f_3 * pc_x[k] * soh_702[k];

        t_934[k] = f_14 * snh_703[k]
                   + f_8 * sog0_505[k]
                   - f_9 * sog1_505[k]
                   + f_3 * pc_x[k] * soh_703[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pc_x, pc_y, pc_z, snh_531, snh_555, snh_705, \
                         sog0_507, sog1_507, soh_699, soh_702, \
                         soh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_20 * snh_531[k]
                   + f_3 * pc_z[k] * soh_699[k];

        t_936[k] = f_14 * snh_705[k]
                   + f_8 * sog0_507[k]
                   - f_9 * sog1_507[k]
                   + f_3 * pc_x[k] * soh_705[k];

        t_937[k] = f_12 * snh_555[k]
                   + f_3 * pc_y[k] * soh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pc_x, snh_707, snh_708, snh_709, snh_710, \
                         sog0_509, sog1_509, soh_707, soh_708, soh_709, \
                         soh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = f_14 * snh_707[k]
                   + f_8 * sog0_509[k]
                   - f_9 * sog1_509[k]
                   + f_3 * pc_x[k] * soh_707[k];

        t_939[k] = f_14 * snh_708[k]
                   + f_3 * pc_x[k] * soh_708[k];

        t_940[k] = f_14 * snh_709[k]
                   + f_3 * pc_x[k] * soh_709[k];

        t_941[k] = f_14 * snh_710[k]
                   + f_3 * pc_x[k] * soh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pc_x, pc_y, snh_561, snh_711, snh_712, \
                         snh_713, sog0_505, sog1_505, soh_708, soh_711, soh_712, \
                         soh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_14 * snh_711[k]
                   + f_3 * pc_x[k] * soh_711[k];

        t_943[k] = f_14 * snh_712[k]
                   + f_3 * pc_x[k] * soh_712[k];

        t_944[k] = f_14 * snh_713[k]
                   + f_3 * pc_x[k] * soh_713[k];

        t_945[k] = f_12 * snh_561[k]
                   + f_1 * sog0_505[k]
                   - f_2 * sog1_505[k]
                   + f_3 * pc_y[k] * soh_708[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, pc_y, pc_z, snh_540, snh_563, snh_564, sog0_507, \
                         sog0_508, sog1_507, sog1_508, soh_708, soh_710, \
                         soh_711 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_20 * snh_540[k]
                   + f_3 * pc_z[k] * soh_708[k];

        t_947[k] = f_12 * snh_563[k]
                   + f_4 * sog0_507[k]
                   - f_5 * sog1_507[k]
                   + f_3 * pc_y[k] * soh_710[k];

        t_948[k] = f_12 * snh_564[k]
                   + f_6 * sog0_508[k]
                   - f_7 * sog1_508[k]
                   + f_3 * pc_y[k] * soh_711[k];
    }

#pragma omp simd aligned(t_949, t_950, t_951, t_952, pb_y, pc_y, pc_z, sni0_756, snh_545, \
                         snh_565, snh_566, sni1_756, sog0_509, sog1_509, soh_712, \
                         soh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_949[k] = f_12 * snh_565[k]
                   + f_8 * sog0_509[k]
                   - f_9 * sog1_509[k]
                   + f_3 * pc_y[k] * soh_712[k];

        t_950[k] = f_12 * snh_566[k]
                   + f_3 * pc_y[k] * soh_713[k];

        t_951[k] = f_20 * snh_545[k]
                   + f_1 * sog0_509[k]
                   - f_2 * sog1_509[k]
                   + f_3 * pc_z[k] * soh_713[k];

        t_952[k] = pb_y[k] * sni0_756[k]
                   - f_10 * pc_y[k] * sni1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, t_956, pb_y, pc_y, pc_z, sni0_759, snh_546, \
                         snh_567, snh_568, snh_569, sni1_759, soh_714, \
                         soh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_11 * snh_567[k]
                   + f_3 * pc_y[k] * soh_714[k];

        t_954[k] = f_19 * snh_546[k]
                   + f_3 * pc_z[k] * soh_714[k];

        t_955[k] = pb_y[k] * sni0_759[k]
                   + f_12 * snh_568[k]
                   - f_10 * pc_y[k] * sni1_759[k];

        t_956[k] = f_11 * snh_569[k]
                   + f_3 * pc_y[k] * soh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, t_960, pb_y, pc_y, pc_z, sni0_761, sni0_762, \
                         snh_549, snh_570, snh_572, sni1_761, sni1_762, soh_717, \
                         soh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pb_y[k] * sni0_761[k]
                   - f_10 * pc_y[k] * sni1_761[k];

        t_958[k] = pb_y[k] * sni0_762[k]
                   + f_13 * snh_570[k]
                   - f_10 * pc_y[k] * sni1_762[k];

        t_959[k] = f_19 * snh_549[k]
                   + f_3 * pc_z[k] * soh_717[k];

        t_960[k] = f_11 * snh_572[k]
                   + f_3 * pc_y[k] * soh_719[k];
    }

#pragma omp simd aligned(t_961, t_962, t_963, pb_y, pc_y, pc_z, sni0_765, sni0_766, snh_552, \
                         snh_573, sni1_765, sni1_766, soh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_961[k] = pb_y[k] * sni0_765[k]
                   - f_10 * pc_y[k] * sni1_765[k];

        t_962[k] = pb_y[k] * sni0_766[k]
                   + f_14 * snh_573[k]
                   - f_10 * pc_y[k] * sni1_766[k];

        t_963[k] = f_19 * snh_552[k]
                   + f_3 * pc_z[k] * soh_720[k];
    }

#pragma omp simd aligned(t_964, t_965, t_966, t_967, pb_y, pc_x, pc_y, sni0_768, sni0_770, \
                         snh_575, snh_576, snh_729, sni1_768, sni1_770, soh_723, \
                         soh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_964[k] = pb_y[k] * sni0_768[k]
                   + f_12 * snh_575[k]
                   - f_10 * pc_y[k] * sni1_768[k];

        t_965[k] = f_11 * snh_576[k]
                   + f_3 * pc_y[k] * soh_723[k];

        t_966[k] = pb_y[k] * sni0_770[k]
                   - f_10 * pc_y[k] * sni1_770[k];

        t_967[k] = f_14 * snh_729[k]
                   + f_3 * pc_x[k] * soh_729[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, t_971, t_972, pc_x, snh_730, snh_731, snh_732, \
                         snh_733, snh_734, soh_730, soh_731, soh_732, soh_733, \
                         soh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = f_14 * snh_730[k]
                   + f_3 * pc_x[k] * soh_730[k];

        t_969[k] = f_14 * snh_731[k]
                   + f_3 * pc_x[k] * soh_731[k];

        t_970[k] = f_14 * snh_732[k]
                   + f_3 * pc_x[k] * soh_732[k];

        t_971[k] = f_14 * snh_733[k]
                   + f_3 * pc_x[k] * soh_733[k];

        t_972[k] = f_14 * snh_734[k]
                   + f_3 * pc_x[k] * soh_734[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, pc_y, pc_z, snh_561, snh_582, snh_584, sog0_520, \
                         sog0_522, sog1_520, sog1_522, soh_729, \
                         soh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_11 * snh_582[k]
                   + f_1 * sog0_520[k]
                   - f_2 * sog1_520[k]
                   + f_3 * pc_y[k] * soh_729[k];

        t_974[k] = f_19 * snh_561[k]
                   + f_3 * pc_z[k] * soh_729[k];

        t_975[k] = f_11 * snh_584[k]
                   + f_4 * sog0_522[k]
                   - f_5 * sog1_522[k]
                   + f_3 * pc_y[k] * soh_731[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, pc_y, snh_585, snh_586, snh_587, sog0_523, \
                         sog0_524, sog1_523, sog1_524, soh_732, soh_733, \
                         soh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_11 * snh_585[k]
                   + f_6 * sog0_523[k]
                   - f_7 * sog1_523[k]
                   + f_3 * pc_y[k] * soh_732[k];

        t_977[k] = f_11 * snh_586[k]
                   + f_8 * sog0_524[k]
                   - f_9 * sog1_524[k]
                   + f_3 * pc_y[k] * soh_733[k];

        t_978[k] = f_11 * snh_587[k]
                   + f_3 * pc_y[k] * soh_734[k];
    }

#pragma omp simd aligned(t_979, t_980, t_981, t_982, pb_y, pc_x, pc_y, pc_z, sni0_783, \
                         snh_567, snh_735, sni1_783, sog0_525, sog1_525, \
                         soh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_979[k] = pb_y[k] * sni0_783[k]
                   - f_10 * pc_y[k] * sni1_783[k];

        t_980[k] = f_14 * snh_735[k]
                   + f_1 * sog0_525[k]
                   - f_2 * sog1_525[k]
                   + f_3 * pc_x[k] * soh_735[k];

        t_981[k] = f_3 * pc_y[k] * soh_735[k];

        t_982[k] = f_18 * snh_567[k]
                   + f_3 * pc_z[k] * soh_735[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pc_x, pc_y, snh_738, snh_740, sog0_528, \
                         sog0_530, sog1_528, sog1_530, soh_737, soh_738, \
                         soh_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_14 * snh_738[k]
                   + f_4 * sog0_528[k]
                   - f_5 * sog1_528[k]
                   + f_3 * pc_x[k] * soh_738[k];

        t_984[k] = f_3 * pc_y[k] * soh_737[k];

        t_985[k] = f_14 * snh_740[k]
                   + f_4 * sog0_530[k]
                   - f_5 * sog1_530[k]
                   + f_3 * pc_x[k] * soh_740[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, pc_x, pc_y, pc_z, snh_570, snh_741, sog0_531, \
                         sog1_531, soh_738, soh_740, soh_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_14 * snh_741[k]
                   + f_6 * sog0_531[k]
                   - f_7 * sog1_531[k]
                   + f_3 * pc_x[k] * soh_741[k];

        t_987[k] = f_18 * snh_570[k]
                   + f_3 * pc_z[k] * soh_738[k];

        t_988[k] = f_3 * pc_y[k] * soh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, pc_z, snh_573, snh_744, snh_745, sog0_534, \
                         sog0_535, sog1_534, sog1_535, soh_741, soh_744, \
                         soh_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_14 * snh_744[k]
                   + f_6 * sog0_534[k]
                   - f_7 * sog1_534[k]
                   + f_3 * pc_x[k] * soh_744[k];

        t_990[k] = f_14 * snh_745[k]
                   + f_8 * sog0_535[k]
                   - f_9 * sog1_535[k]
                   + f_3 * pc_x[k] * soh_745[k];

        t_991[k] = f_18 * snh_573[k]
                   + f_3 * pc_z[k] * soh_741[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, pc_x, pc_y, snh_747, snh_749, sog0_537, \
                         sog0_539, sog1_537, sog1_539, soh_744, soh_747, \
                         soh_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_14 * snh_747[k]
                   + f_8 * sog0_537[k]
                   - f_9 * sog1_537[k]
                   + f_3 * pc_x[k] * soh_747[k];

        t_993[k] = f_3 * pc_y[k] * soh_744[k];

        t_994[k] = f_14 * snh_749[k]
                   + f_8 * sog0_539[k]
                   - f_9 * sog1_539[k]
                   + f_3 * pc_x[k] * soh_749[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, t_999, pc_x, snh_750, snh_751, snh_752, \
                         snh_753, snh_754, soh_750, soh_751, soh_752, soh_753, \
                         soh_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_14 * snh_750[k]
                   + f_3 * pc_x[k] * soh_750[k];

        t_996[k] = f_14 * snh_751[k]
                   + f_3 * pc_x[k] * soh_751[k];

        t_997[k] = f_14 * snh_752[k]
                   + f_3 * pc_x[k] * soh_752[k];

        t_998[k] = f_14 * snh_753[k]
                   + f_3 * pc_x[k] * soh_753[k];

        t_999[k] = f_14 * snh_754[k]
                   + f_3 * pc_x[k] * soh_754[k];
    }

#pragma omp simd aligned(t_1000, t_1001, t_1002, t_1003, pc_x, pc_y, pc_z, snh_582, snh_755, \
                         sog0_535, sog0_537, sog1_535, sog1_537, soh_750, soh_752, \
                         soh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1000[k] = f_14 * snh_755[k]
                    + f_3 * pc_x[k] * soh_755[k];

        t_1001[k] = f_1 * sog0_535[k]
                    - f_2 * sog1_535[k]
                    + f_3 * pc_y[k] * soh_750[k];

        t_1002[k] = f_18 * snh_582[k]
                    + f_3 * pc_z[k] * soh_750[k];

        t_1003[k] = f_4 * sog0_537[k]
                    - f_5 * sog1_537[k]
                    + f_3 * pc_y[k] * soh_752[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, t_1007, pc_y, pc_z, snh_587, sog0_538, \
                         sog0_539, sog1_538, sog1_539, soh_753, soh_754, \
                         soh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_6 * sog0_538[k]
                    - f_7 * sog1_538[k]
                    + f_3 * pc_y[k] * soh_753[k];

        t_1005[k] = f_8 * sog0_539[k]
                    - f_9 * sog1_539[k]
                    + f_3 * pc_y[k] * soh_754[k];

        t_1006[k] = f_3 * pc_y[k] * soh_755[k];

        t_1007[k] = f_18 * snh_587[k]
                    + f_1 * sog0_539[k]
                    - f_2 * sog1_539[k]
                    + f_3 * pc_z[k] * soh_755[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, t_1011, pc_x, pc_y, pc_z, snh_588, snh_756, \
                         snh_759, sog0_540, sog0_543, sog1_540, sog1_543, soh_756, \
                         soh_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_13 * snh_756[k]
                    + f_1 * sog0_540[k]
                    - f_2 * sog1_540[k]
                    + f_3 * pc_x[k] * soh_756[k];

        t_1009[k] = f_17 * snh_588[k]
                    + f_3 * pc_y[k] * soh_756[k];

        t_1010[k] = f_3 * pc_z[k] * soh_756[k];

        t_1011[k] = f_13 * snh_759[k]
                    + f_4 * sog0_543[k]
                    - f_5 * sog1_543[k]
                    + f_3 * pc_x[k] * soh_759[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, pc_x, pc_y, snh_590, snh_761, snh_762, \
                         sog0_545, sog0_546, sog1_545, sog1_546, soh_758, soh_761, \
                         soh_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_17 * snh_590[k]
                    + f_3 * pc_y[k] * soh_758[k];

        t_1013[k] = f_13 * snh_761[k]
                    + f_4 * sog0_545[k]
                    - f_5 * sog1_545[k]
                    + f_3 * pc_x[k] * soh_761[k];

        t_1014[k] = f_13 * snh_762[k]
                    + f_6 * sog0_546[k]
                    - f_7 * sog1_546[k]
                    + f_3 * pc_x[k] * soh_762[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, pc_x, pc_y, pc_z, snh_593, snh_765, sog0_549, \
                         sog1_549, soh_759, soh_761, soh_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = f_3 * pc_z[k] * soh_759[k];

        t_1016[k] = f_17 * snh_593[k]
                    + f_3 * pc_y[k] * soh_761[k];

        t_1017[k] = f_13 * snh_765[k]
                    + f_6 * sog0_549[k]
                    - f_7 * sog1_549[k]
                    + f_3 * pc_x[k] * soh_765[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, pc_x, pc_z, snh_766, snh_768, sog0_550, \
                         sog0_552, sog1_550, sog1_552, soh_762, soh_766, \
                         soh_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_13 * snh_766[k]
                    + f_8 * sog0_550[k]
                    - f_9 * sog1_550[k]
                    + f_3 * pc_x[k] * soh_766[k];

        t_1019[k] = f_3 * pc_z[k] * soh_762[k];

        t_1020[k] = f_13 * snh_768[k]
                    + f_8 * sog0_552[k]
                    - f_9 * sog1_552[k]
                    + f_3 * pc_x[k] * soh_768[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sni0,
                                                          const size_t snh, const size_t sni1,
                                                          const size_t sog0, const size_t sog1,
                                                          const size_t soh, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_784 = buffer.data(sni0 + 784);
    const auto *sni0_787 = buffer.data(sni0 + 787);
    const auto *sni0_790 = buffer.data(sni0 + 790);
    const auto *sni0_794 = buffer.data(sni0 + 794);
    const auto *sni0_796 = buffer.data(sni0 + 796);
    const auto *sni0_805 = buffer.data(sni0 + 805);

    const auto *snh_588 = buffer.data(snh + 588);
    const auto *snh_591 = buffer.data(snh + 591);
    const auto *snh_594 = buffer.data(snh + 594);
    const auto *snh_595 = buffer.data(snh + 595);
    const auto *snh_597 = buffer.data(snh + 597);
    const auto *snh_603 = buffer.data(snh + 603);
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
    const auto *snh_624 = buffer.data(snh + 624);
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
    const auto *snh_645 = buffer.data(snh + 645);
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
    const auto *snh_666 = buffer.data(snh + 666);
    const auto *snh_668 = buffer.data(snh + 668);
    const auto *snh_669 = buffer.data(snh + 669);
    const auto *snh_670 = buffer.data(snh + 670);
    const auto *snh_671 = buffer.data(snh + 671);
    const auto *snh_672 = buffer.data(snh + 672);
    const auto *snh_674 = buffer.data(snh + 674);
    const auto *snh_677 = buffer.data(snh + 677);
    const auto *snh_770 = buffer.data(snh + 770);
    const auto *snh_771 = buffer.data(snh + 771);
    const auto *snh_772 = buffer.data(snh + 772);
    const auto *snh_773 = buffer.data(snh + 773);
    const auto *snh_774 = buffer.data(snh + 774);
    const auto *snh_775 = buffer.data(snh + 775);
    const auto *snh_776 = buffer.data(snh + 776);
    const auto *snh_782 = buffer.data(snh + 782);
    const auto *snh_786 = buffer.data(snh + 786);
    const auto *snh_791 = buffer.data(snh + 791);
    const auto *snh_792 = buffer.data(snh + 792);
    const auto *snh_793 = buffer.data(snh + 793);
    const auto *snh_794 = buffer.data(snh + 794);
    const auto *snh_795 = buffer.data(snh + 795);
    const auto *snh_796 = buffer.data(snh + 796);
    const auto *snh_797 = buffer.data(snh + 797);
    const auto *snh_798 = buffer.data(snh + 798);
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
    const auto *snh_843 = buffer.data(snh + 843);
    const auto *snh_845 = buffer.data(snh + 845);
    const auto *snh_846 = buffer.data(snh + 846);
    const auto *snh_849 = buffer.data(snh + 849);
    const auto *snh_850 = buffer.data(snh + 850);

    const auto *sni1_784 = buffer.data(sni1 + 784);
    const auto *sni1_787 = buffer.data(sni1 + 787);
    const auto *sni1_790 = buffer.data(sni1 + 790);
    const auto *sni1_794 = buffer.data(sni1 + 794);
    const auto *sni1_796 = buffer.data(sni1 + 796);
    const auto *sni1_805 = buffer.data(sni1 + 805);

    const auto *sog0_550 = buffer.data(sog0 + 550);
    const auto *sog0_552 = buffer.data(sog0 + 552);
    const auto *sog0_553 = buffer.data(sog0 + 553);
    const auto *sog0_554 = buffer.data(sog0 + 554);
    const auto *sog0_560 = buffer.data(sog0 + 560);
    const auto *sog0_564 = buffer.data(sog0 + 564);
    const auto *sog0_567 = buffer.data(sog0 + 567);
    const auto *sog0_568 = buffer.data(sog0 + 568);
    const auto *sog0_569 = buffer.data(sog0 + 569);
    const auto *sog0_570 = buffer.data(sog0 + 570);
    const auto *sog0_573 = buffer.data(sog0 + 573);
    const auto *sog0_575 = buffer.data(sog0 + 575);
    const auto *sog0_576 = buffer.data(sog0 + 576);
    const auto *sog0_579 = buffer.data(sog0 + 579);
    const auto *sog0_580 = buffer.data(sog0 + 580);
    const auto *sog0_582 = buffer.data(sog0 + 582);
    const auto *sog0_583 = buffer.data(sog0 + 583);
    const auto *sog0_584 = buffer.data(sog0 + 584);
    const auto *sog0_585 = buffer.data(sog0 + 585);
    const auto *sog0_588 = buffer.data(sog0 + 588);
    const auto *sog0_590 = buffer.data(sog0 + 590);
    const auto *sog0_591 = buffer.data(sog0 + 591);
    const auto *sog0_594 = buffer.data(sog0 + 594);
    const auto *sog0_595 = buffer.data(sog0 + 595);
    const auto *sog0_597 = buffer.data(sog0 + 597);
    const auto *sog0_598 = buffer.data(sog0 + 598);
    const auto *sog0_599 = buffer.data(sog0 + 599);
    const auto *sog0_600 = buffer.data(sog0 + 600);
    const auto *sog0_603 = buffer.data(sog0 + 603);
    const auto *sog0_605 = buffer.data(sog0 + 605);
    const auto *sog0_606 = buffer.data(sog0 + 606);
    const auto *sog0_609 = buffer.data(sog0 + 609);
    const auto *sog0_610 = buffer.data(sog0 + 610);

    const auto *sog1_550 = buffer.data(sog1 + 550);
    const auto *sog1_552 = buffer.data(sog1 + 552);
    const auto *sog1_553 = buffer.data(sog1 + 553);
    const auto *sog1_554 = buffer.data(sog1 + 554);
    const auto *sog1_560 = buffer.data(sog1 + 560);
    const auto *sog1_564 = buffer.data(sog1 + 564);
    const auto *sog1_567 = buffer.data(sog1 + 567);
    const auto *sog1_568 = buffer.data(sog1 + 568);
    const auto *sog1_569 = buffer.data(sog1 + 569);
    const auto *sog1_570 = buffer.data(sog1 + 570);
    const auto *sog1_573 = buffer.data(sog1 + 573);
    const auto *sog1_575 = buffer.data(sog1 + 575);
    const auto *sog1_576 = buffer.data(sog1 + 576);
    const auto *sog1_579 = buffer.data(sog1 + 579);
    const auto *sog1_580 = buffer.data(sog1 + 580);
    const auto *sog1_582 = buffer.data(sog1 + 582);
    const auto *sog1_583 = buffer.data(sog1 + 583);
    const auto *sog1_584 = buffer.data(sog1 + 584);
    const auto *sog1_585 = buffer.data(sog1 + 585);
    const auto *sog1_588 = buffer.data(sog1 + 588);
    const auto *sog1_590 = buffer.data(sog1 + 590);
    const auto *sog1_591 = buffer.data(sog1 + 591);
    const auto *sog1_594 = buffer.data(sog1 + 594);
    const auto *sog1_595 = buffer.data(sog1 + 595);
    const auto *sog1_597 = buffer.data(sog1 + 597);
    const auto *sog1_598 = buffer.data(sog1 + 598);
    const auto *sog1_599 = buffer.data(sog1 + 599);
    const auto *sog1_600 = buffer.data(sog1 + 600);
    const auto *sog1_603 = buffer.data(sog1 + 603);
    const auto *sog1_605 = buffer.data(sog1 + 605);
    const auto *sog1_606 = buffer.data(sog1 + 606);
    const auto *sog1_609 = buffer.data(sog1 + 609);
    const auto *sog1_610 = buffer.data(sog1 + 610);

    const auto *soh_765 = buffer.data(soh + 765);
    const auto *soh_770 = buffer.data(soh + 770);
    const auto *soh_771 = buffer.data(soh + 771);
    const auto *soh_772 = buffer.data(soh + 772);
    const auto *soh_773 = buffer.data(soh + 773);
    const auto *soh_774 = buffer.data(soh + 774);
    const auto *soh_775 = buffer.data(soh + 775);
    const auto *soh_776 = buffer.data(soh + 776);
    const auto *soh_777 = buffer.data(soh + 777);
    const auto *soh_779 = buffer.data(soh + 779);
    const auto *soh_780 = buffer.data(soh + 780);
    const auto *soh_782 = buffer.data(soh + 782);
    const auto *soh_783 = buffer.data(soh + 783);
    const auto *soh_786 = buffer.data(soh + 786);
    const auto *soh_791 = buffer.data(soh + 791);
    const auto *soh_792 = buffer.data(soh + 792);
    const auto *soh_793 = buffer.data(soh + 793);
    const auto *soh_794 = buffer.data(soh + 794);
    const auto *soh_795 = buffer.data(soh + 795);
    const auto *soh_796 = buffer.data(soh + 796);
    const auto *soh_797 = buffer.data(soh + 797);
    const auto *soh_798 = buffer.data(soh + 798);
    const auto *soh_800 = buffer.data(soh + 800);
    const auto *soh_801 = buffer.data(soh + 801);
    const auto *soh_803 = buffer.data(soh + 803);
    const auto *soh_804 = buffer.data(soh + 804);
    const auto *soh_807 = buffer.data(soh + 807);
    const auto *soh_808 = buffer.data(soh + 808);
    const auto *soh_810 = buffer.data(soh + 810);
    const auto *soh_812 = buffer.data(soh + 812);
    const auto *soh_813 = buffer.data(soh + 813);
    const auto *soh_814 = buffer.data(soh + 814);
    const auto *soh_815 = buffer.data(soh + 815);
    const auto *soh_816 = buffer.data(soh + 816);
    const auto *soh_817 = buffer.data(soh + 817);
    const auto *soh_818 = buffer.data(soh + 818);
    const auto *soh_819 = buffer.data(soh + 819);
    const auto *soh_821 = buffer.data(soh + 821);
    const auto *soh_822 = buffer.data(soh + 822);
    const auto *soh_824 = buffer.data(soh + 824);
    const auto *soh_825 = buffer.data(soh + 825);
    const auto *soh_828 = buffer.data(soh + 828);
    const auto *soh_829 = buffer.data(soh + 829);
    const auto *soh_831 = buffer.data(soh + 831);
    const auto *soh_833 = buffer.data(soh + 833);
    const auto *soh_834 = buffer.data(soh + 834);
    const auto *soh_835 = buffer.data(soh + 835);
    const auto *soh_836 = buffer.data(soh + 836);
    const auto *soh_837 = buffer.data(soh + 837);
    const auto *soh_838 = buffer.data(soh + 838);
    const auto *soh_839 = buffer.data(soh + 839);
    const auto *soh_840 = buffer.data(soh + 840);
    const auto *soh_842 = buffer.data(soh + 842);
    const auto *soh_843 = buffer.data(soh + 843);
    const auto *soh_845 = buffer.data(soh + 845);
    const auto *soh_846 = buffer.data(soh + 846);
    const auto *soh_849 = buffer.data(soh + 849);
    const auto *soh_850 = buffer.data(soh + 850);

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, pc_x, pc_y, snh_597, snh_770, \
                         snh_771, snh_772, sog0_554, sog1_554, soh_765, soh_770, soh_771, \
                         soh_772 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = f_17 * snh_597[k]
                    + f_3 * pc_y[k] * soh_765[k];

        t_1022[k] = f_13 * snh_770[k]
                    + f_8 * sog0_554[k]
                    - f_9 * sog1_554[k]
                    + f_3 * pc_x[k] * soh_770[k];

        t_1023[k] = f_13 * snh_771[k]
                    + f_3 * pc_x[k] * soh_771[k];

        t_1024[k] = f_13 * snh_772[k]
                    + f_3 * pc_x[k] * soh_772[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pc_x, snh_773, snh_774, snh_775, \
                         snh_776, soh_773, soh_774, soh_775, soh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_13 * snh_773[k]
                    + f_3 * pc_x[k] * soh_773[k];

        t_1026[k] = f_13 * snh_774[k]
                    + f_3 * pc_x[k] * soh_774[k];

        t_1027[k] = f_13 * snh_775[k]
                    + f_3 * pc_x[k] * soh_775[k];

        t_1028[k] = f_13 * snh_776[k]
                    + f_3 * pc_x[k] * soh_776[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, pc_y, pc_z, snh_603, snh_605, sog0_550, \
                         sog0_552, sog1_550, sog1_552, soh_771, \
                         soh_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_17 * snh_603[k]
                    + f_1 * sog0_550[k]
                    - f_2 * sog1_550[k]
                    + f_3 * pc_y[k] * soh_771[k];

        t_1030[k] = f_3 * pc_z[k] * soh_771[k];

        t_1031[k] = f_17 * snh_605[k]
                    + f_4 * sog0_552[k]
                    - f_5 * sog1_552[k]
                    + f_3 * pc_y[k] * soh_773[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, pc_y, pc_z, snh_606, snh_607, \
                         snh_608, sog0_553, sog0_554, sog1_553, sog1_554, soh_774, soh_775, \
                         soh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_17 * snh_606[k]
                    + f_6 * sog0_553[k]
                    - f_7 * sog1_553[k]
                    + f_3 * pc_y[k] * soh_774[k];

        t_1033[k] = f_17 * snh_607[k]
                    + f_8 * sog0_554[k]
                    - f_9 * sog1_554[k]
                    + f_3 * pc_y[k] * soh_775[k];

        t_1034[k] = f_17 * snh_608[k]
                    + f_3 * pc_y[k] * soh_776[k];

        t_1035[k] = f_1 * sog0_554[k]
                    - f_2 * sog1_554[k]
                    + f_3 * pc_z[k] * soh_776[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pb_z, pc_y, pc_z, sni0_784, sni0_787, \
                         snh_588, snh_609, sni1_784, sni1_787, \
                         soh_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pb_z[k] * sni0_784[k]
                    - f_10 * pc_z[k] * sni1_784[k];

        t_1037[k] = f_18 * snh_609[k]
                    + f_3 * pc_y[k] * soh_777[k];

        t_1038[k] = f_11 * snh_588[k]
                    + f_3 * pc_z[k] * soh_777[k];

        t_1039[k] = pb_z[k] * sni0_787[k]
                    - f_10 * pc_z[k] * sni1_787[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pb_z, pc_x, pc_y, pc_z, sni0_790, snh_611, \
                         snh_782, sni1_790, sog0_560, sog1_560, soh_779, \
                         soh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_18 * snh_611[k]
                    + f_3 * pc_y[k] * soh_779[k];

        t_1041[k] = f_13 * snh_782[k]
                    + f_4 * sog0_560[k]
                    - f_5 * sog1_560[k]
                    + f_3 * pc_x[k] * soh_782[k];

        t_1042[k] = pb_z[k] * sni0_790[k]
                    - f_10 * pc_z[k] * sni1_790[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, pc_z, snh_591, snh_614, snh_786, \
                         sog0_564, sog1_564, soh_780, soh_782, \
                         soh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_11 * snh_591[k]
                    + f_3 * pc_z[k] * soh_780[k];

        t_1044[k] = f_18 * snh_614[k]
                    + f_3 * pc_y[k] * soh_782[k];

        t_1045[k] = f_13 * snh_786[k]
                    + f_6 * sog0_564[k]
                    - f_7 * sog1_564[k]
                    + f_3 * pc_x[k] * soh_786[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, t_1049, pb_z, pc_y, pc_z, sni0_794, sni0_796, \
                         snh_594, snh_595, snh_618, sni1_794, sni1_796, soh_783, \
                         soh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pb_z[k] * sni0_794[k]
                    - f_10 * pc_z[k] * sni1_794[k];

        t_1047[k] = f_11 * snh_594[k]
                    + f_3 * pc_z[k] * soh_783[k];

        t_1048[k] = pb_z[k] * sni0_796[k]
                    + f_12 * snh_595[k]
                    - f_10 * pc_z[k] * sni1_796[k];

        t_1049[k] = f_18 * snh_618[k]
                    + f_3 * pc_y[k] * soh_786[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, t_1053, pc_x, snh_791, snh_792, snh_793, \
                         snh_794, sog0_569, sog1_569, soh_791, soh_792, soh_793, \
                         soh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_13 * snh_791[k]
                    + f_8 * sog0_569[k]
                    - f_9 * sog1_569[k]
                    + f_3 * pc_x[k] * soh_791[k];

        t_1051[k] = f_13 * snh_792[k]
                    + f_3 * pc_x[k] * soh_792[k];

        t_1052[k] = f_13 * snh_793[k]
                    + f_3 * pc_x[k] * soh_793[k];

        t_1053[k] = f_13 * snh_794[k]
                    + f_3 * pc_x[k] * soh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, pb_z, pc_x, pc_z, sni0_805, snh_795, \
                         snh_796, snh_797, sni1_805, soh_795, soh_796, \
                         soh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_13 * snh_795[k]
                    + f_3 * pc_x[k] * soh_795[k];

        t_1055[k] = f_13 * snh_796[k]
                    + f_3 * pc_x[k] * soh_796[k];

        t_1056[k] = f_13 * snh_797[k]
                    + f_3 * pc_x[k] * soh_797[k];

        t_1057[k] = pb_z[k] * sni0_805[k]
                    - f_10 * pc_z[k] * sni1_805[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_y, pc_z, snh_603, snh_626, snh_627, \
                         sog0_567, sog0_568, sog1_567, sog1_568, soh_792, soh_794, \
                         soh_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_11 * snh_603[k]
                    + f_3 * pc_z[k] * soh_792[k];

        t_1059[k] = f_18 * snh_626[k]
                    + f_4 * sog0_567[k]
                    - f_5 * sog1_567[k]
                    + f_3 * pc_y[k] * soh_794[k];

        t_1060[k] = f_18 * snh_627[k]
                    + f_6 * sog0_568[k]
                    - f_7 * sog1_568[k]
                    + f_3 * pc_y[k] * soh_795[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, pc_y, pc_z, snh_608, snh_628, snh_629, \
                         sog0_569, sog1_569, soh_796, soh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_18 * snh_628[k]
                    + f_8 * sog0_569[k]
                    - f_9 * sog1_569[k]
                    + f_3 * pc_y[k] * soh_796[k];

        t_1062[k] = f_18 * snh_629[k]
                    + f_3 * pc_y[k] * soh_797[k];

        t_1063[k] = f_11 * snh_608[k]
                    + f_1 * sog0_569[k]
                    - f_2 * sog1_569[k]
                    + f_3 * pc_z[k] * soh_797[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, pc_x, pc_y, pc_z, snh_609, snh_630, snh_798, \
                         sog0_570, sog1_570, soh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = f_13 * snh_798[k]
                    + f_1 * sog0_570[k]
                    - f_2 * sog1_570[k]
                    + f_3 * pc_x[k] * soh_798[k];

        t_1065[k] = f_19 * snh_630[k]
                    + f_3 * pc_y[k] * soh_798[k];

        t_1066[k] = f_12 * snh_609[k]
                    + f_3 * pc_z[k] * soh_798[k];
    }

#pragma omp simd aligned(t_1067, t_1068, t_1069, pc_x, pc_y, snh_632, snh_801, snh_803, \
                         sog0_573, sog0_575, sog1_573, sog1_575, soh_800, soh_801, \
                         soh_803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1067[k] = f_13 * snh_801[k]
                    + f_4 * sog0_573[k]
                    - f_5 * sog1_573[k]
                    + f_3 * pc_x[k] * soh_801[k];

        t_1068[k] = f_19 * snh_632[k]
                    + f_3 * pc_y[k] * soh_800[k];

        t_1069[k] = f_13 * snh_803[k]
                    + f_4 * sog0_575[k]
                    - f_5 * sog1_575[k]
                    + f_3 * pc_x[k] * soh_803[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, snh_612, snh_635, snh_804, \
                         sog0_576, sog1_576, soh_801, soh_803, \
                         soh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_13 * snh_804[k]
                    + f_6 * sog0_576[k]
                    - f_7 * sog1_576[k]
                    + f_3 * pc_x[k] * soh_804[k];

        t_1071[k] = f_12 * snh_612[k]
                    + f_3 * pc_z[k] * soh_801[k];

        t_1072[k] = f_19 * snh_635[k]
                    + f_3 * pc_y[k] * soh_803[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_x, pc_z, snh_615, snh_807, snh_808, \
                         sog0_579, sog0_580, sog1_579, sog1_580, soh_804, soh_807, \
                         soh_808 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_13 * snh_807[k]
                    + f_6 * sog0_579[k]
                    - f_7 * sog1_579[k]
                    + f_3 * pc_x[k] * soh_807[k];

        t_1074[k] = f_13 * snh_808[k]
                    + f_8 * sog0_580[k]
                    - f_9 * sog1_580[k]
                    + f_3 * pc_x[k] * soh_808[k];

        t_1075[k] = f_12 * snh_615[k]
                    + f_3 * pc_z[k] * soh_804[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_x, pc_y, snh_639, snh_810, snh_812, \
                         sog0_582, sog0_584, sog1_582, sog1_584, soh_807, soh_810, \
                         soh_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_13 * snh_810[k]
                    + f_8 * sog0_582[k]
                    - f_9 * sog1_582[k]
                    + f_3 * pc_x[k] * soh_810[k];

        t_1077[k] = f_19 * snh_639[k]
                    + f_3 * pc_y[k] * soh_807[k];

        t_1078[k] = f_13 * snh_812[k]
                    + f_8 * sog0_584[k]
                    - f_9 * sog1_584[k]
                    + f_3 * pc_x[k] * soh_812[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, t_1083, pc_x, snh_813, snh_814, \
                         snh_815, snh_816, snh_817, soh_813, soh_814, soh_815, soh_816, \
                         soh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_13 * snh_813[k]
                    + f_3 * pc_x[k] * soh_813[k];

        t_1080[k] = f_13 * snh_814[k]
                    + f_3 * pc_x[k] * soh_814[k];

        t_1081[k] = f_13 * snh_815[k]
                    + f_3 * pc_x[k] * soh_815[k];

        t_1082[k] = f_13 * snh_816[k]
                    + f_3 * pc_x[k] * soh_816[k];

        t_1083[k] = f_13 * snh_817[k]
                    + f_3 * pc_x[k] * soh_817[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, pc_x, pc_y, pc_z, snh_624, snh_645, snh_818, \
                         sog0_580, sog1_580, soh_813, soh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_13 * snh_818[k]
                    + f_3 * pc_x[k] * soh_818[k];

        t_1085[k] = f_19 * snh_645[k]
                    + f_1 * sog0_580[k]
                    - f_2 * sog1_580[k]
                    + f_3 * pc_y[k] * soh_813[k];

        t_1086[k] = f_12 * snh_624[k]
                    + f_3 * pc_z[k] * soh_813[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, snh_647, snh_648, snh_649, sog0_582, \
                         sog0_583, sog0_584, sog1_582, sog1_583, sog1_584, soh_815, soh_816, \
                         soh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_19 * snh_647[k]
                    + f_4 * sog0_582[k]
                    - f_5 * sog1_582[k]
                    + f_3 * pc_y[k] * soh_815[k];

        t_1088[k] = f_19 * snh_648[k]
                    + f_6 * sog0_583[k]
                    - f_7 * sog1_583[k]
                    + f_3 * pc_y[k] * soh_816[k];

        t_1089[k] = f_19 * snh_649[k]
                    + f_8 * sog0_584[k]
                    - f_9 * sog1_584[k]
                    + f_3 * pc_y[k] * soh_817[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, snh_629, snh_650, snh_819, \
                         sog0_584, sog0_585, sog1_584, sog1_585, soh_818, \
                         soh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_19 * snh_650[k]
                    + f_3 * pc_y[k] * soh_818[k];

        t_1091[k] = f_12 * snh_629[k]
                    + f_1 * sog0_584[k]
                    - f_2 * sog1_584[k]
                    + f_3 * pc_z[k] * soh_818[k];

        t_1092[k] = f_13 * snh_819[k]
                    + f_1 * sog0_585[k]
                    - f_2 * sog1_585[k]
                    + f_3 * pc_x[k] * soh_819[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, t_1096, pc_x, pc_y, pc_z, snh_630, snh_651, \
                         snh_653, snh_822, sog0_588, sog1_588, soh_819, soh_821, \
                         soh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_20 * snh_651[k]
                    + f_3 * pc_y[k] * soh_819[k];

        t_1094[k] = f_13 * snh_630[k]
                    + f_3 * pc_z[k] * soh_819[k];

        t_1095[k] = f_13 * snh_822[k]
                    + f_4 * sog0_588[k]
                    - f_5 * sog1_588[k]
                    + f_3 * pc_x[k] * soh_822[k];

        t_1096[k] = f_20 * snh_653[k]
                    + f_3 * pc_y[k] * soh_821[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_z, snh_633, snh_824, snh_825, \
                         sog0_590, sog0_591, sog1_590, sog1_591, soh_822, soh_824, \
                         soh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_13 * snh_824[k]
                    + f_4 * sog0_590[k]
                    - f_5 * sog1_590[k]
                    + f_3 * pc_x[k] * soh_824[k];

        t_1098[k] = f_13 * snh_825[k]
                    + f_6 * sog0_591[k]
                    - f_7 * sog1_591[k]
                    + f_3 * pc_x[k] * soh_825[k];

        t_1099[k] = f_13 * snh_633[k]
                    + f_3 * pc_z[k] * soh_822[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_y, snh_656, snh_828, snh_829, \
                         sog0_594, sog0_595, sog1_594, sog1_595, soh_824, soh_828, \
                         soh_829 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_20 * snh_656[k]
                    + f_3 * pc_y[k] * soh_824[k];

        t_1101[k] = f_13 * snh_828[k]
                    + f_6 * sog0_594[k]
                    - f_7 * sog1_594[k]
                    + f_3 * pc_x[k] * soh_828[k];

        t_1102[k] = f_13 * snh_829[k]
                    + f_8 * sog0_595[k]
                    - f_9 * sog1_595[k]
                    + f_3 * pc_x[k] * soh_829[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, pc_y, pc_z, snh_636, snh_660, snh_831, \
                         sog0_597, sog1_597, soh_825, soh_828, \
                         soh_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_13 * snh_636[k]
                    + f_3 * pc_z[k] * soh_825[k];

        t_1104[k] = f_13 * snh_831[k]
                    + f_8 * sog0_597[k]
                    - f_9 * sog1_597[k]
                    + f_3 * pc_x[k] * soh_831[k];

        t_1105[k] = f_20 * snh_660[k]
                    + f_3 * pc_y[k] * soh_828[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, snh_833, snh_834, snh_835, \
                         snh_836, sog0_599, sog1_599, soh_833, soh_834, soh_835, \
                         soh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_13 * snh_833[k]
                    + f_8 * sog0_599[k]
                    - f_9 * sog1_599[k]
                    + f_3 * pc_x[k] * soh_833[k];

        t_1107[k] = f_13 * snh_834[k]
                    + f_3 * pc_x[k] * soh_834[k];

        t_1108[k] = f_13 * snh_835[k]
                    + f_3 * pc_x[k] * soh_835[k];

        t_1109[k] = f_13 * snh_836[k]
                    + f_3 * pc_x[k] * soh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, pc_x, pc_y, snh_666, snh_837, \
                         snh_838, snh_839, sog0_595, sog1_595, soh_834, soh_837, soh_838, \
                         soh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_13 * snh_837[k]
                    + f_3 * pc_x[k] * soh_837[k];

        t_1111[k] = f_13 * snh_838[k]
                    + f_3 * pc_x[k] * soh_838[k];

        t_1112[k] = f_13 * snh_839[k]
                    + f_3 * pc_x[k] * soh_839[k];

        t_1113[k] = f_20 * snh_666[k]
                    + f_1 * sog0_595[k]
                    - f_2 * sog1_595[k]
                    + f_3 * pc_y[k] * soh_834[k];
    }

#pragma omp simd aligned(t_1114, t_1115, t_1116, pc_y, pc_z, snh_645, snh_668, snh_669, \
                         sog0_597, sog0_598, sog1_597, sog1_598, soh_834, soh_836, \
                         soh_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1114[k] = f_13 * snh_645[k]
                    + f_3 * pc_z[k] * soh_834[k];

        t_1115[k] = f_20 * snh_668[k]
                    + f_4 * sog0_597[k]
                    - f_5 * sog1_597[k]
                    + f_3 * pc_y[k] * soh_836[k];

        t_1116[k] = f_20 * snh_669[k]
                    + f_6 * sog0_598[k]
                    - f_7 * sog1_598[k]
                    + f_3 * pc_y[k] * soh_837[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, pc_y, pc_z, snh_650, snh_670, snh_671, \
                         sog0_599, sog1_599, soh_838, soh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_20 * snh_670[k]
                    + f_8 * sog0_599[k]
                    - f_9 * sog1_599[k]
                    + f_3 * pc_y[k] * soh_838[k];

        t_1118[k] = f_20 * snh_671[k]
                    + f_3 * pc_y[k] * soh_839[k];

        t_1119[k] = f_13 * snh_650[k]
                    + f_1 * sog0_599[k]
                    - f_2 * sog1_599[k]
                    + f_3 * pc_z[k] * soh_839[k];
    }

#pragma omp simd aligned(t_1120, t_1121, t_1122, pc_x, pc_y, pc_z, snh_651, snh_672, snh_840, \
                         sog0_600, sog1_600, soh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1120[k] = f_13 * snh_840[k]
                    + f_1 * sog0_600[k]
                    - f_2 * sog1_600[k]
                    + f_3 * pc_x[k] * soh_840[k];

        t_1121[k] = f_14 * snh_672[k]
                    + f_3 * pc_y[k] * soh_840[k];

        t_1122[k] = f_14 * snh_651[k]
                    + f_3 * pc_z[k] * soh_840[k];
    }

#pragma omp simd aligned(t_1123, t_1124, t_1125, pc_x, pc_y, snh_674, snh_843, snh_845, \
                         sog0_603, sog0_605, sog1_603, sog1_605, soh_842, soh_843, \
                         soh_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_13 * snh_843[k]
                    + f_4 * sog0_603[k]
                    - f_5 * sog1_603[k]
                    + f_3 * pc_x[k] * soh_843[k];

        t_1124[k] = f_14 * snh_674[k]
                    + f_3 * pc_y[k] * soh_842[k];

        t_1125[k] = f_13 * snh_845[k]
                    + f_4 * sog0_605[k]
                    - f_5 * sog1_605[k]
                    + f_3 * pc_x[k] * soh_845[k];
    }

#pragma omp simd aligned(t_1126, t_1127, t_1128, pc_x, pc_y, pc_z, snh_654, snh_677, snh_846, \
                         sog0_606, sog1_606, soh_843, soh_845, \
                         soh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1126[k] = f_13 * snh_846[k]
                    + f_6 * sog0_606[k]
                    - f_7 * sog1_606[k]
                    + f_3 * pc_x[k] * soh_846[k];

        t_1127[k] = f_14 * snh_654[k]
                    + f_3 * pc_z[k] * soh_843[k];

        t_1128[k] = f_14 * snh_677[k]
                    + f_3 * pc_y[k] * soh_845[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pc_x, pc_z, snh_657, snh_849, snh_850, \
                         sog0_609, sog0_610, sog1_609, sog1_610, soh_846, soh_849, \
                         soh_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_13 * snh_849[k]
                    + f_6 * sog0_609[k]
                    - f_7 * sog1_609[k]
                    + f_3 * pc_x[k] * soh_849[k];

        t_1130[k] = f_13 * snh_850[k]
                    + f_8 * sog0_610[k]
                    - f_9 * sog1_610[k]
                    + f_3 * pc_x[k] * soh_850[k];

        t_1131[k] = f_14 * snh_657[k]
                    + f_3 * pc_z[k] * soh_846[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t sog0, const size_t sog1,
                                                           const size_t soh, const size_t ncols,
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
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *sni0_980 = buffer.data(sni0 + 980);
    const auto *sni0_983 = buffer.data(sni0 + 983);
    const auto *sni0_985 = buffer.data(sni0 + 985);
    const auto *sni0_986 = buffer.data(sni0 + 986);
    const auto *sni0_989 = buffer.data(sni0 + 989);
    const auto *sni0_990 = buffer.data(sni0 + 990);
    const auto *sni0_992 = buffer.data(sni0 + 992);
    const auto *sni0_994 = buffer.data(sni0 + 994);
    const auto *sni0_1007 = buffer.data(sni0 + 1007);

    const auto *snh_666 = buffer.data(snh + 666);
    const auto *snh_671 = buffer.data(snh + 671);
    const auto *snh_672 = buffer.data(snh + 672);
    const auto *snh_675 = buffer.data(snh + 675);
    const auto *snh_678 = buffer.data(snh + 678);
    const auto *snh_681 = buffer.data(snh + 681);
    const auto *snh_687 = buffer.data(snh + 687);
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
    const auto *snh_708 = buffer.data(snh + 708);
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
    const auto *snh_731 = buffer.data(snh + 731);
    const auto *snh_732 = buffer.data(snh + 732);
    const auto *snh_733 = buffer.data(snh + 733);
    const auto *snh_734 = buffer.data(snh + 734);
    const auto *snh_735 = buffer.data(snh + 735);
    const auto *snh_736 = buffer.data(snh + 736);
    const auto *snh_737 = buffer.data(snh + 737);
    const auto *snh_738 = buffer.data(snh + 738);
    const auto *snh_740 = buffer.data(snh + 740);
    const auto *snh_741 = buffer.data(snh + 741);
    const auto *snh_743 = buffer.data(snh + 743);
    const auto *snh_744 = buffer.data(snh + 744);
    const auto *snh_750 = buffer.data(snh + 750);
    const auto *snh_752 = buffer.data(snh + 752);
    const auto *snh_753 = buffer.data(snh + 753);
    const auto *snh_754 = buffer.data(snh + 754);
    const auto *snh_755 = buffer.data(snh + 755);
    const auto *snh_852 = buffer.data(snh + 852);
    const auto *snh_854 = buffer.data(snh + 854);
    const auto *snh_855 = buffer.data(snh + 855);
    const auto *snh_856 = buffer.data(snh + 856);
    const auto *snh_857 = buffer.data(snh + 857);
    const auto *snh_858 = buffer.data(snh + 858);
    const auto *snh_859 = buffer.data(snh + 859);
    const auto *snh_860 = buffer.data(snh + 860);
    const auto *snh_861 = buffer.data(snh + 861);
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
    const auto *snh_918 = buffer.data(snh + 918);
    const auto *snh_919 = buffer.data(snh + 919);
    const auto *snh_920 = buffer.data(snh + 920);
    const auto *snh_921 = buffer.data(snh + 921);
    const auto *snh_922 = buffer.data(snh + 922);
    const auto *snh_923 = buffer.data(snh + 923);
    const auto *snh_924 = buffer.data(snh + 924);
    const auto *snh_927 = buffer.data(snh + 927);
    const auto *snh_929 = buffer.data(snh + 929);
    const auto *snh_930 = buffer.data(snh + 930);

    const auto *sni1_980 = buffer.data(sni1 + 980);
    const auto *sni1_983 = buffer.data(sni1 + 983);
    const auto *sni1_985 = buffer.data(sni1 + 985);
    const auto *sni1_986 = buffer.data(sni1 + 986);
    const auto *sni1_989 = buffer.data(sni1 + 989);
    const auto *sni1_990 = buffer.data(sni1 + 990);
    const auto *sni1_992 = buffer.data(sni1 + 992);
    const auto *sni1_994 = buffer.data(sni1 + 994);
    const auto *sni1_1007 = buffer.data(sni1 + 1007);

    const auto *sog0_610 = buffer.data(sog0 + 610);
    const auto *sog0_612 = buffer.data(sog0 + 612);
    const auto *sog0_613 = buffer.data(sog0 + 613);
    const auto *sog0_614 = buffer.data(sog0 + 614);
    const auto *sog0_615 = buffer.data(sog0 + 615);
    const auto *sog0_618 = buffer.data(sog0 + 618);
    const auto *sog0_620 = buffer.data(sog0 + 620);
    const auto *sog0_621 = buffer.data(sog0 + 621);
    const auto *sog0_624 = buffer.data(sog0 + 624);
    const auto *sog0_625 = buffer.data(sog0 + 625);
    const auto *sog0_627 = buffer.data(sog0 + 627);
    const auto *sog0_628 = buffer.data(sog0 + 628);
    const auto *sog0_629 = buffer.data(sog0 + 629);
    const auto *sog0_630 = buffer.data(sog0 + 630);
    const auto *sog0_633 = buffer.data(sog0 + 633);
    const auto *sog0_635 = buffer.data(sog0 + 635);
    const auto *sog0_636 = buffer.data(sog0 + 636);
    const auto *sog0_639 = buffer.data(sog0 + 639);
    const auto *sog0_640 = buffer.data(sog0 + 640);
    const auto *sog0_642 = buffer.data(sog0 + 642);
    const auto *sog0_643 = buffer.data(sog0 + 643);
    const auto *sog0_644 = buffer.data(sog0 + 644);
    const auto *sog0_655 = buffer.data(sog0 + 655);
    const auto *sog0_657 = buffer.data(sog0 + 657);
    const auto *sog0_658 = buffer.data(sog0 + 658);
    const auto *sog0_659 = buffer.data(sog0 + 659);
    const auto *sog0_660 = buffer.data(sog0 + 660);
    const auto *sog0_663 = buffer.data(sog0 + 663);
    const auto *sog0_665 = buffer.data(sog0 + 665);
    const auto *sog0_666 = buffer.data(sog0 + 666);

    const auto *sog1_610 = buffer.data(sog1 + 610);
    const auto *sog1_612 = buffer.data(sog1 + 612);
    const auto *sog1_613 = buffer.data(sog1 + 613);
    const auto *sog1_614 = buffer.data(sog1 + 614);
    const auto *sog1_615 = buffer.data(sog1 + 615);
    const auto *sog1_618 = buffer.data(sog1 + 618);
    const auto *sog1_620 = buffer.data(sog1 + 620);
    const auto *sog1_621 = buffer.data(sog1 + 621);
    const auto *sog1_624 = buffer.data(sog1 + 624);
    const auto *sog1_625 = buffer.data(sog1 + 625);
    const auto *sog1_627 = buffer.data(sog1 + 627);
    const auto *sog1_628 = buffer.data(sog1 + 628);
    const auto *sog1_629 = buffer.data(sog1 + 629);
    const auto *sog1_630 = buffer.data(sog1 + 630);
    const auto *sog1_633 = buffer.data(sog1 + 633);
    const auto *sog1_635 = buffer.data(sog1 + 635);
    const auto *sog1_636 = buffer.data(sog1 + 636);
    const auto *sog1_639 = buffer.data(sog1 + 639);
    const auto *sog1_640 = buffer.data(sog1 + 640);
    const auto *sog1_642 = buffer.data(sog1 + 642);
    const auto *sog1_643 = buffer.data(sog1 + 643);
    const auto *sog1_644 = buffer.data(sog1 + 644);
    const auto *sog1_655 = buffer.data(sog1 + 655);
    const auto *sog1_657 = buffer.data(sog1 + 657);
    const auto *sog1_658 = buffer.data(sog1 + 658);
    const auto *sog1_659 = buffer.data(sog1 + 659);
    const auto *sog1_660 = buffer.data(sog1 + 660);
    const auto *sog1_663 = buffer.data(sog1 + 663);
    const auto *sog1_665 = buffer.data(sog1 + 665);
    const auto *sog1_666 = buffer.data(sog1 + 666);

    const auto *soh_849 = buffer.data(soh + 849);
    const auto *soh_852 = buffer.data(soh + 852);
    const auto *soh_854 = buffer.data(soh + 854);
    const auto *soh_855 = buffer.data(soh + 855);
    const auto *soh_856 = buffer.data(soh + 856);
    const auto *soh_857 = buffer.data(soh + 857);
    const auto *soh_858 = buffer.data(soh + 858);
    const auto *soh_859 = buffer.data(soh + 859);
    const auto *soh_860 = buffer.data(soh + 860);
    const auto *soh_861 = buffer.data(soh + 861);
    const auto *soh_863 = buffer.data(soh + 863);
    const auto *soh_864 = buffer.data(soh + 864);
    const auto *soh_866 = buffer.data(soh + 866);
    const auto *soh_867 = buffer.data(soh + 867);
    const auto *soh_870 = buffer.data(soh + 870);
    const auto *soh_871 = buffer.data(soh + 871);
    const auto *soh_873 = buffer.data(soh + 873);
    const auto *soh_875 = buffer.data(soh + 875);
    const auto *soh_876 = buffer.data(soh + 876);
    const auto *soh_877 = buffer.data(soh + 877);
    const auto *soh_878 = buffer.data(soh + 878);
    const auto *soh_879 = buffer.data(soh + 879);
    const auto *soh_880 = buffer.data(soh + 880);
    const auto *soh_881 = buffer.data(soh + 881);
    const auto *soh_882 = buffer.data(soh + 882);
    const auto *soh_884 = buffer.data(soh + 884);
    const auto *soh_885 = buffer.data(soh + 885);
    const auto *soh_887 = buffer.data(soh + 887);
    const auto *soh_888 = buffer.data(soh + 888);
    const auto *soh_891 = buffer.data(soh + 891);
    const auto *soh_892 = buffer.data(soh + 892);
    const auto *soh_894 = buffer.data(soh + 894);
    const auto *soh_896 = buffer.data(soh + 896);
    const auto *soh_897 = buffer.data(soh + 897);
    const auto *soh_898 = buffer.data(soh + 898);
    const auto *soh_899 = buffer.data(soh + 899);
    const auto *soh_900 = buffer.data(soh + 900);
    const auto *soh_901 = buffer.data(soh + 901);
    const auto *soh_902 = buffer.data(soh + 902);
    const auto *soh_903 = buffer.data(soh + 903);
    const auto *soh_905 = buffer.data(soh + 905);
    const auto *soh_906 = buffer.data(soh + 906);
    const auto *soh_908 = buffer.data(soh + 908);
    const auto *soh_909 = buffer.data(soh + 909);
    const auto *soh_912 = buffer.data(soh + 912);
    const auto *soh_918 = buffer.data(soh + 918);
    const auto *soh_919 = buffer.data(soh + 919);
    const auto *soh_920 = buffer.data(soh + 920);
    const auto *soh_921 = buffer.data(soh + 921);
    const auto *soh_922 = buffer.data(soh + 922);
    const auto *soh_923 = buffer.data(soh + 923);
    const auto *soh_924 = buffer.data(soh + 924);
    const auto *soh_926 = buffer.data(soh + 926);
    const auto *soh_927 = buffer.data(soh + 927);
    const auto *soh_929 = buffer.data(soh + 929);
    const auto *soh_930 = buffer.data(soh + 930);

#pragma omp simd aligned(t_1132, t_1133, t_1134, pc_x, pc_y, snh_681, snh_852, snh_854, \
                         sog0_612, sog0_614, sog1_612, sog1_614, soh_849, soh_852, \
                         soh_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_13 * snh_852[k]
                    + f_8 * sog0_612[k]
                    - f_9 * sog1_612[k]
                    + f_3 * pc_x[k] * soh_852[k];

        t_1133[k] = f_14 * snh_681[k]
                    + f_3 * pc_y[k] * soh_849[k];

        t_1134[k] = f_13 * snh_854[k]
                    + f_8 * sog0_614[k]
                    - f_9 * sog1_614[k]
                    + f_3 * pc_x[k] * soh_854[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, t_1138, t_1139, pc_x, snh_855, snh_856, \
                         snh_857, snh_858, snh_859, soh_855, soh_856, soh_857, soh_858, \
                         soh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_13 * snh_855[k]
                    + f_3 * pc_x[k] * soh_855[k];

        t_1136[k] = f_13 * snh_856[k]
                    + f_3 * pc_x[k] * soh_856[k];

        t_1137[k] = f_13 * snh_857[k]
                    + f_3 * pc_x[k] * soh_857[k];

        t_1138[k] = f_13 * snh_858[k]
                    + f_3 * pc_x[k] * soh_858[k];

        t_1139[k] = f_13 * snh_859[k]
                    + f_3 * pc_x[k] * soh_859[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, snh_666, snh_687, snh_860, \
                         sog0_610, sog1_610, soh_855, soh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_13 * snh_860[k]
                    + f_3 * pc_x[k] * soh_860[k];

        t_1141[k] = f_14 * snh_687[k]
                    + f_1 * sog0_610[k]
                    - f_2 * sog1_610[k]
                    + f_3 * pc_y[k] * soh_855[k];

        t_1142[k] = f_14 * snh_666[k]
                    + f_3 * pc_z[k] * soh_855[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_y, snh_689, snh_690, snh_691, sog0_612, \
                         sog0_613, sog0_614, sog1_612, sog1_613, sog1_614, soh_857, soh_858, \
                         soh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_14 * snh_689[k]
                    + f_4 * sog0_612[k]
                    - f_5 * sog1_612[k]
                    + f_3 * pc_y[k] * soh_857[k];

        t_1144[k] = f_14 * snh_690[k]
                    + f_6 * sog0_613[k]
                    - f_7 * sog1_613[k]
                    + f_3 * pc_y[k] * soh_858[k];

        t_1145[k] = f_14 * snh_691[k]
                    + f_8 * sog0_614[k]
                    - f_9 * sog1_614[k]
                    + f_3 * pc_y[k] * soh_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_y, pc_z, snh_671, snh_692, snh_861, \
                         sog0_614, sog0_615, sog1_614, sog1_615, soh_860, \
                         soh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * snh_692[k]
                    + f_3 * pc_y[k] * soh_860[k];

        t_1147[k] = f_14 * snh_671[k]
                    + f_1 * sog0_614[k]
                    - f_2 * sog1_614[k]
                    + f_3 * pc_z[k] * soh_860[k];

        t_1148[k] = f_13 * snh_861[k]
                    + f_1 * sog0_615[k]
                    - f_2 * sog1_615[k]
                    + f_3 * pc_x[k] * soh_861[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, t_1152, pc_x, pc_y, pc_z, snh_672, snh_693, \
                         snh_695, snh_864, sog0_618, sog1_618, soh_861, soh_863, \
                         soh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_13 * snh_693[k]
                    + f_3 * pc_y[k] * soh_861[k];

        t_1150[k] = f_20 * snh_672[k]
                    + f_3 * pc_z[k] * soh_861[k];

        t_1151[k] = f_13 * snh_864[k]
                    + f_4 * sog0_618[k]
                    - f_5 * sog1_618[k]
                    + f_3 * pc_x[k] * soh_864[k];

        t_1152[k] = f_13 * snh_695[k]
                    + f_3 * pc_y[k] * soh_863[k];
    }

#pragma omp simd aligned(t_1153, t_1154, t_1155, pc_x, pc_z, snh_675, snh_866, snh_867, \
                         sog0_620, sog0_621, sog1_620, sog1_621, soh_864, soh_866, \
                         soh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1153[k] = f_13 * snh_866[k]
                    + f_4 * sog0_620[k]
                    - f_5 * sog1_620[k]
                    + f_3 * pc_x[k] * soh_866[k];

        t_1154[k] = f_13 * snh_867[k]
                    + f_6 * sog0_621[k]
                    - f_7 * sog1_621[k]
                    + f_3 * pc_x[k] * soh_867[k];

        t_1155[k] = f_20 * snh_675[k]
                    + f_3 * pc_z[k] * soh_864[k];
    }

#pragma omp simd aligned(t_1156, t_1157, t_1158, pc_x, pc_y, snh_698, snh_870, snh_871, \
                         sog0_624, sog0_625, sog1_624, sog1_625, soh_866, soh_870, \
                         soh_871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1156[k] = f_13 * snh_698[k]
                    + f_3 * pc_y[k] * soh_866[k];

        t_1157[k] = f_13 * snh_870[k]
                    + f_6 * sog0_624[k]
                    - f_7 * sog1_624[k]
                    + f_3 * pc_x[k] * soh_870[k];

        t_1158[k] = f_13 * snh_871[k]
                    + f_8 * sog0_625[k]
                    - f_9 * sog1_625[k]
                    + f_3 * pc_x[k] * soh_871[k];
    }

#pragma omp simd aligned(t_1159, t_1160, t_1161, pc_x, pc_y, pc_z, snh_678, snh_702, snh_873, \
                         sog0_627, sog1_627, soh_867, soh_870, \
                         soh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1159[k] = f_20 * snh_678[k]
                    + f_3 * pc_z[k] * soh_867[k];

        t_1160[k] = f_13 * snh_873[k]
                    + f_8 * sog0_627[k]
                    - f_9 * sog1_627[k]
                    + f_3 * pc_x[k] * soh_873[k];

        t_1161[k] = f_13 * snh_702[k]
                    + f_3 * pc_y[k] * soh_870[k];
    }

#pragma omp simd aligned(t_1162, t_1163, t_1164, t_1165, pc_x, snh_875, snh_876, snh_877, \
                         snh_878, sog0_629, sog1_629, soh_875, soh_876, soh_877, \
                         soh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1162[k] = f_13 * snh_875[k]
                    + f_8 * sog0_629[k]
                    - f_9 * sog1_629[k]
                    + f_3 * pc_x[k] * soh_875[k];

        t_1163[k] = f_13 * snh_876[k]
                    + f_3 * pc_x[k] * soh_876[k];

        t_1164[k] = f_13 * snh_877[k]
                    + f_3 * pc_x[k] * soh_877[k];

        t_1165[k] = f_13 * snh_878[k]
                    + f_3 * pc_x[k] * soh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, pc_x, pc_y, snh_708, snh_879, \
                         snh_880, snh_881, sog0_625, sog1_625, soh_876, soh_879, soh_880, \
                         soh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_13 * snh_879[k]
                    + f_3 * pc_x[k] * soh_879[k];

        t_1167[k] = f_13 * snh_880[k]
                    + f_3 * pc_x[k] * soh_880[k];

        t_1168[k] = f_13 * snh_881[k]
                    + f_3 * pc_x[k] * soh_881[k];

        t_1169[k] = f_13 * snh_708[k]
                    + f_1 * sog0_625[k]
                    - f_2 * sog1_625[k]
                    + f_3 * pc_y[k] * soh_876[k];
    }

#pragma omp simd aligned(t_1170, t_1171, t_1172, pc_y, pc_z, snh_687, snh_710, snh_711, \
                         sog0_627, sog0_628, sog1_627, sog1_628, soh_876, soh_878, \
                         soh_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1170[k] = f_20 * snh_687[k]
                    + f_3 * pc_z[k] * soh_876[k];

        t_1171[k] = f_13 * snh_710[k]
                    + f_4 * sog0_627[k]
                    - f_5 * sog1_627[k]
                    + f_3 * pc_y[k] * soh_878[k];

        t_1172[k] = f_13 * snh_711[k]
                    + f_6 * sog0_628[k]
                    - f_7 * sog1_628[k]
                    + f_3 * pc_y[k] * soh_879[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pc_y, pc_z, snh_692, snh_712, snh_713, \
                         sog0_629, sog1_629, soh_880, soh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_13 * snh_712[k]
                    + f_8 * sog0_629[k]
                    - f_9 * sog1_629[k]
                    + f_3 * pc_y[k] * soh_880[k];

        t_1174[k] = f_13 * snh_713[k]
                    + f_3 * pc_y[k] * soh_881[k];

        t_1175[k] = f_20 * snh_692[k]
                    + f_1 * sog0_629[k]
                    - f_2 * sog1_629[k]
                    + f_3 * pc_z[k] * soh_881[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, snh_693, snh_714, snh_882, \
                         sog0_630, sog1_630, soh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_13 * snh_882[k]
                    + f_1 * sog0_630[k]
                    - f_2 * sog1_630[k]
                    + f_3 * pc_x[k] * soh_882[k];

        t_1177[k] = f_12 * snh_714[k]
                    + f_3 * pc_y[k] * soh_882[k];

        t_1178[k] = f_19 * snh_693[k]
                    + f_3 * pc_z[k] * soh_882[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pc_x, pc_y, snh_716, snh_885, snh_887, \
                         sog0_633, sog0_635, sog1_633, sog1_635, soh_884, soh_885, \
                         soh_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = f_13 * snh_885[k]
                    + f_4 * sog0_633[k]
                    - f_5 * sog1_633[k]
                    + f_3 * pc_x[k] * soh_885[k];

        t_1180[k] = f_12 * snh_716[k]
                    + f_3 * pc_y[k] * soh_884[k];

        t_1181[k] = f_13 * snh_887[k]
                    + f_4 * sog0_635[k]
                    - f_5 * sog1_635[k]
                    + f_3 * pc_x[k] * soh_887[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pc_x, pc_y, pc_z, snh_696, snh_719, snh_888, \
                         sog0_636, sog1_636, soh_885, soh_887, \
                         soh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_13 * snh_888[k]
                    + f_6 * sog0_636[k]
                    - f_7 * sog1_636[k]
                    + f_3 * pc_x[k] * soh_888[k];

        t_1183[k] = f_19 * snh_696[k]
                    + f_3 * pc_z[k] * soh_885[k];

        t_1184[k] = f_12 * snh_719[k]
                    + f_3 * pc_y[k] * soh_887[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, snh_699, snh_891, snh_892, \
                         sog0_639, sog0_640, sog1_639, sog1_640, soh_888, soh_891, \
                         soh_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_13 * snh_891[k]
                    + f_6 * sog0_639[k]
                    - f_7 * sog1_639[k]
                    + f_3 * pc_x[k] * soh_891[k];

        t_1186[k] = f_13 * snh_892[k]
                    + f_8 * sog0_640[k]
                    - f_9 * sog1_640[k]
                    + f_3 * pc_x[k] * soh_892[k];

        t_1187[k] = f_19 * snh_699[k]
                    + f_3 * pc_z[k] * soh_888[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pc_x, pc_y, snh_723, snh_894, snh_896, \
                         sog0_642, sog0_644, sog1_642, sog1_644, soh_891, soh_894, \
                         soh_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_13 * snh_894[k]
                    + f_8 * sog0_642[k]
                    - f_9 * sog1_642[k]
                    + f_3 * pc_x[k] * soh_894[k];

        t_1189[k] = f_12 * snh_723[k]
                    + f_3 * pc_y[k] * soh_891[k];

        t_1190[k] = f_13 * snh_896[k]
                    + f_8 * sog0_644[k]
                    - f_9 * sog1_644[k]
                    + f_3 * pc_x[k] * soh_896[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, t_1194, t_1195, pc_x, snh_897, snh_898, \
                         snh_899, snh_900, snh_901, soh_897, soh_898, soh_899, soh_900, \
                         soh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_13 * snh_897[k]
                    + f_3 * pc_x[k] * soh_897[k];

        t_1192[k] = f_13 * snh_898[k]
                    + f_3 * pc_x[k] * soh_898[k];

        t_1193[k] = f_13 * snh_899[k]
                    + f_3 * pc_x[k] * soh_899[k];

        t_1194[k] = f_13 * snh_900[k]
                    + f_3 * pc_x[k] * soh_900[k];

        t_1195[k] = f_13 * snh_901[k]
                    + f_3 * pc_x[k] * soh_901[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, snh_708, snh_729, snh_902, \
                         sog0_640, sog1_640, soh_897, soh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_13 * snh_902[k]
                    + f_3 * pc_x[k] * soh_902[k];

        t_1197[k] = f_12 * snh_729[k]
                    + f_1 * sog0_640[k]
                    - f_2 * sog1_640[k]
                    + f_3 * pc_y[k] * soh_897[k];

        t_1198[k] = f_19 * snh_708[k]
                    + f_3 * pc_z[k] * soh_897[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_y, snh_731, snh_732, snh_733, sog0_642, \
                         sog0_643, sog0_644, sog1_642, sog1_643, sog1_644, soh_899, soh_900, \
                         soh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_12 * snh_731[k]
                    + f_4 * sog0_642[k]
                    - f_5 * sog1_642[k]
                    + f_3 * pc_y[k] * soh_899[k];

        t_1200[k] = f_12 * snh_732[k]
                    + f_6 * sog0_643[k]
                    - f_7 * sog1_643[k]
                    + f_3 * pc_y[k] * soh_900[k];

        t_1201[k] = f_12 * snh_733[k]
                    + f_8 * sog0_644[k]
                    - f_9 * sog1_644[k]
                    + f_3 * pc_y[k] * soh_901[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, t_1205, pb_y, pc_y, pc_z, sni0_980, snh_713, \
                         snh_734, snh_735, sni1_980, sog0_644, sog1_644, soh_902, \
                         soh_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * snh_734[k]
                    + f_3 * pc_y[k] * soh_902[k];

        t_1203[k] = f_19 * snh_713[k]
                    + f_1 * sog0_644[k]
                    - f_2 * sog1_644[k]
                    + f_3 * pc_z[k] * soh_902[k];

        t_1204[k] = pb_y[k] * sni0_980[k]
                    - f_10 * pc_y[k] * sni1_980[k];

        t_1205[k] = f_11 * snh_735[k]
                    + f_3 * pc_y[k] * soh_903[k];
    }

#pragma omp simd aligned(t_1206, t_1207, t_1208, t_1209, pb_y, pc_y, pc_z, sni0_983, sni0_985, \
                         snh_714, snh_736, snh_737, sni1_983, sni1_985, soh_903, \
                         soh_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1206[k] = f_18 * snh_714[k]
                    + f_3 * pc_z[k] * soh_903[k];

        t_1207[k] = pb_y[k] * sni0_983[k]
                    + f_12 * snh_736[k]
                    - f_10 * pc_y[k] * sni1_983[k];

        t_1208[k] = f_11 * snh_737[k]
                    + f_3 * pc_y[k] * soh_905[k];

        t_1209[k] = pb_y[k] * sni0_985[k]
                    - f_10 * pc_y[k] * sni1_985[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, t_1213, pb_y, pc_y, pc_z, sni0_986, sni0_989, \
                         snh_717, snh_738, snh_740, sni1_986, sni1_989, soh_906, \
                         soh_908 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pb_y[k] * sni0_986[k]
                    + f_13 * snh_738[k]
                    - f_10 * pc_y[k] * sni1_986[k];

        t_1211[k] = f_18 * snh_717[k]
                    + f_3 * pc_z[k] * soh_906[k];

        t_1212[k] = f_11 * snh_740[k]
                    + f_3 * pc_y[k] * soh_908[k];

        t_1213[k] = pb_y[k] * sni0_989[k]
                    - f_10 * pc_y[k] * sni1_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pb_y, pc_y, pc_z, sni0_990, sni0_992, \
                         snh_720, snh_741, snh_743, sni1_990, sni1_992, \
                         soh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pb_y[k] * sni0_990[k]
                    + f_14 * snh_741[k]
                    - f_10 * pc_y[k] * sni1_990[k];

        t_1215[k] = f_18 * snh_720[k]
                    + f_3 * pc_z[k] * soh_909[k];

        t_1216[k] = pb_y[k] * sni0_992[k]
                    + f_12 * snh_743[k]
                    - f_10 * pc_y[k] * sni1_992[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, pb_y, pc_x, pc_y, sni0_994, snh_744, \
                         snh_918, snh_919, sni1_994, soh_912, soh_918, \
                         soh_919 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_11 * snh_744[k]
                    + f_3 * pc_y[k] * soh_912[k];

        t_1218[k] = pb_y[k] * sni0_994[k]
                    - f_10 * pc_y[k] * sni1_994[k];

        t_1219[k] = f_13 * snh_918[k]
                    + f_3 * pc_x[k] * soh_918[k];

        t_1220[k] = f_13 * snh_919[k]
                    + f_3 * pc_x[k] * soh_919[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, pc_x, snh_920, snh_921, snh_922, \
                         snh_923, soh_920, soh_921, soh_922, soh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_13 * snh_920[k]
                    + f_3 * pc_x[k] * soh_920[k];

        t_1222[k] = f_13 * snh_921[k]
                    + f_3 * pc_x[k] * soh_921[k];

        t_1223[k] = f_13 * snh_922[k]
                    + f_3 * pc_x[k] * soh_922[k];

        t_1224[k] = f_13 * snh_923[k]
                    + f_3 * pc_x[k] * soh_923[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, pc_y, pc_z, snh_729, snh_750, snh_752, \
                         sog0_655, sog0_657, sog1_655, sog1_657, soh_918, \
                         soh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_11 * snh_750[k]
                    + f_1 * sog0_655[k]
                    - f_2 * sog1_655[k]
                    + f_3 * pc_y[k] * soh_918[k];

        t_1226[k] = f_18 * snh_729[k]
                    + f_3 * pc_z[k] * soh_918[k];

        t_1227[k] = f_11 * snh_752[k]
                    + f_4 * sog0_657[k]
                    - f_5 * sog1_657[k]
                    + f_3 * pc_y[k] * soh_920[k];
    }

#pragma omp simd aligned(t_1228, t_1229, t_1230, pc_y, snh_753, snh_754, snh_755, sog0_658, \
                         sog0_659, sog1_658, sog1_659, soh_921, soh_922, \
                         soh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1228[k] = f_11 * snh_753[k]
                    + f_6 * sog0_658[k]
                    - f_7 * sog1_658[k]
                    + f_3 * pc_y[k] * soh_921[k];

        t_1229[k] = f_11 * snh_754[k]
                    + f_8 * sog0_659[k]
                    - f_9 * sog1_659[k]
                    + f_3 * pc_y[k] * soh_922[k];

        t_1230[k] = f_11 * snh_755[k]
                    + f_3 * pc_y[k] * soh_923[k];
    }

#pragma omp simd aligned(t_1231, t_1232, t_1233, t_1234, pb_y, pc_x, pc_y, pc_z, sni0_1007, \
                         snh_735, snh_924, sni1_1007, sog0_660, sog1_660, \
                         soh_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1231[k] = pb_y[k] * sni0_1007[k]
                    - f_10 * pc_y[k] * sni1_1007[k];

        t_1232[k] = f_13 * snh_924[k]
                    + f_1 * sog0_660[k]
                    - f_2 * sog1_660[k]
                    + f_3 * pc_x[k] * soh_924[k];

        t_1233[k] = f_3 * pc_y[k] * soh_924[k];

        t_1234[k] = f_17 * snh_735[k]
                    + f_3 * pc_z[k] * soh_924[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pc_x, pc_y, snh_927, snh_929, sog0_663, \
                         sog0_665, sog1_663, sog1_665, soh_926, soh_927, \
                         soh_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = f_13 * snh_927[k]
                    + f_4 * sog0_663[k]
                    - f_5 * sog1_663[k]
                    + f_3 * pc_x[k] * soh_927[k];

        t_1236[k] = f_3 * pc_y[k] * soh_926[k];

        t_1237[k] = f_13 * snh_929[k]
                    + f_4 * sog0_665[k]
                    - f_5 * sog1_665[k]
                    + f_3 * pc_x[k] * soh_929[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, pc_x, pc_y, pc_z, snh_738, snh_930, sog0_666, \
                         sog1_666, soh_927, soh_929, soh_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_13 * snh_930[k]
                    + f_6 * sog0_666[k]
                    - f_7 * sog1_666[k]
                    + f_3 * pc_x[k] * soh_930[k];

        t_1239[k] = f_17 * snh_738[k]
                    + f_3 * pc_z[k] * soh_927[k];

        t_1240[k] = f_3 * pc_y[k] * soh_929[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t sog0, const size_t sog1,
                                                           const size_t soh, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sni0_1008 = buffer.data(sni0 + 1008);
    const auto *sni0_1011 = buffer.data(sni0 + 1011);
    const auto *sni0_1014 = buffer.data(sni0 + 1014);
    const auto *sni0_1018 = buffer.data(sni0 + 1018);
    const auto *sni0_1020 = buffer.data(sni0 + 1020);
    const auto *sni0_1029 = buffer.data(sni0 + 1029);

    const auto *snh_741 = buffer.data(snh + 741);
    const auto *snh_750 = buffer.data(snh + 750);
    const auto *snh_755 = buffer.data(snh + 755);
    const auto *snh_756 = buffer.data(snh + 756);
    const auto *snh_758 = buffer.data(snh + 758);
    const auto *snh_759 = buffer.data(snh + 759);
    const auto *snh_761 = buffer.data(snh + 761);
    const auto *snh_762 = buffer.data(snh + 762);
    const auto *snh_763 = buffer.data(snh + 763);
    const auto *snh_765 = buffer.data(snh + 765);
    const auto *snh_771 = buffer.data(snh + 771);
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
    const auto *snh_792 = buffer.data(snh + 792);
    const auto *snh_794 = buffer.data(snh + 794);
    const auto *snh_795 = buffer.data(snh + 795);
    const auto *snh_796 = buffer.data(snh + 796);
    const auto *snh_797 = buffer.data(snh + 797);
    const auto *snh_798 = buffer.data(snh + 798);
    const auto *snh_800 = buffer.data(snh + 800);
    const auto *snh_801 = buffer.data(snh + 801);
    const auto *snh_803 = buffer.data(snh + 803);
    const auto *snh_807 = buffer.data(snh + 807);
    const auto *snh_813 = buffer.data(snh + 813);
    const auto *snh_815 = buffer.data(snh + 815);
    const auto *snh_816 = buffer.data(snh + 816);
    const auto *snh_817 = buffer.data(snh + 817);
    const auto *snh_818 = buffer.data(snh + 818);
    const auto *snh_819 = buffer.data(snh + 819);
    const auto *snh_821 = buffer.data(snh + 821);
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
    const auto *snh_948 = buffer.data(snh + 948);
    const auto *snh_950 = buffer.data(snh + 950);
    const auto *snh_951 = buffer.data(snh + 951);
    const auto *snh_954 = buffer.data(snh + 954);
    const auto *snh_955 = buffer.data(snh + 955);
    const auto *snh_957 = buffer.data(snh + 957);
    const auto *snh_959 = buffer.data(snh + 959);
    const auto *snh_960 = buffer.data(snh + 960);
    const auto *snh_961 = buffer.data(snh + 961);
    const auto *snh_962 = buffer.data(snh + 962);
    const auto *snh_963 = buffer.data(snh + 963);
    const auto *snh_964 = buffer.data(snh + 964);
    const auto *snh_965 = buffer.data(snh + 965);
    const auto *snh_971 = buffer.data(snh + 971);
    const auto *snh_975 = buffer.data(snh + 975);
    const auto *snh_980 = buffer.data(snh + 980);
    const auto *snh_981 = buffer.data(snh + 981);
    const auto *snh_982 = buffer.data(snh + 982);
    const auto *snh_983 = buffer.data(snh + 983);
    const auto *snh_984 = buffer.data(snh + 984);
    const auto *snh_985 = buffer.data(snh + 985);
    const auto *snh_986 = buffer.data(snh + 986);
    const auto *snh_987 = buffer.data(snh + 987);
    const auto *snh_990 = buffer.data(snh + 990);
    const auto *snh_992 = buffer.data(snh + 992);
    const auto *snh_993 = buffer.data(snh + 993);
    const auto *snh_996 = buffer.data(snh + 996);
    const auto *snh_997 = buffer.data(snh + 997);
    const auto *snh_999 = buffer.data(snh + 999);
    const auto *snh_1001 = buffer.data(snh + 1001);
    const auto *snh_1002 = buffer.data(snh + 1002);
    const auto *snh_1003 = buffer.data(snh + 1003);
    const auto *snh_1004 = buffer.data(snh + 1004);
    const auto *snh_1005 = buffer.data(snh + 1005);
    const auto *snh_1006 = buffer.data(snh + 1006);
    const auto *snh_1007 = buffer.data(snh + 1007);
    const auto *snh_1008 = buffer.data(snh + 1008);
    const auto *snh_1011 = buffer.data(snh + 1011);
    const auto *snh_1013 = buffer.data(snh + 1013);
    const auto *snh_1014 = buffer.data(snh + 1014);

    const auto *sni1_1008 = buffer.data(sni1 + 1008);
    const auto *sni1_1011 = buffer.data(sni1 + 1011);
    const auto *sni1_1014 = buffer.data(sni1 + 1014);
    const auto *sni1_1018 = buffer.data(sni1 + 1018);
    const auto *sni1_1020 = buffer.data(sni1 + 1020);
    const auto *sni1_1029 = buffer.data(sni1 + 1029);

    const auto *sog0_669 = buffer.data(sog0 + 669);
    const auto *sog0_670 = buffer.data(sog0 + 670);
    const auto *sog0_672 = buffer.data(sog0 + 672);
    const auto *sog0_673 = buffer.data(sog0 + 673);
    const auto *sog0_674 = buffer.data(sog0 + 674);
    const auto *sog0_675 = buffer.data(sog0 + 675);
    const auto *sog0_678 = buffer.data(sog0 + 678);
    const auto *sog0_680 = buffer.data(sog0 + 680);
    const auto *sog0_681 = buffer.data(sog0 + 681);
    const auto *sog0_684 = buffer.data(sog0 + 684);
    const auto *sog0_685 = buffer.data(sog0 + 685);
    const auto *sog0_687 = buffer.data(sog0 + 687);
    const auto *sog0_688 = buffer.data(sog0 + 688);
    const auto *sog0_689 = buffer.data(sog0 + 689);
    const auto *sog0_695 = buffer.data(sog0 + 695);
    const auto *sog0_699 = buffer.data(sog0 + 699);
    const auto *sog0_702 = buffer.data(sog0 + 702);
    const auto *sog0_703 = buffer.data(sog0 + 703);
    const auto *sog0_704 = buffer.data(sog0 + 704);
    const auto *sog0_705 = buffer.data(sog0 + 705);
    const auto *sog0_708 = buffer.data(sog0 + 708);
    const auto *sog0_710 = buffer.data(sog0 + 710);
    const auto *sog0_711 = buffer.data(sog0 + 711);
    const auto *sog0_714 = buffer.data(sog0 + 714);
    const auto *sog0_715 = buffer.data(sog0 + 715);
    const auto *sog0_717 = buffer.data(sog0 + 717);
    const auto *sog0_718 = buffer.data(sog0 + 718);
    const auto *sog0_719 = buffer.data(sog0 + 719);
    const auto *sog0_720 = buffer.data(sog0 + 720);
    const auto *sog0_723 = buffer.data(sog0 + 723);
    const auto *sog0_725 = buffer.data(sog0 + 725);
    const auto *sog0_726 = buffer.data(sog0 + 726);

    const auto *sog1_669 = buffer.data(sog1 + 669);
    const auto *sog1_670 = buffer.data(sog1 + 670);
    const auto *sog1_672 = buffer.data(sog1 + 672);
    const auto *sog1_673 = buffer.data(sog1 + 673);
    const auto *sog1_674 = buffer.data(sog1 + 674);
    const auto *sog1_675 = buffer.data(sog1 + 675);
    const auto *sog1_678 = buffer.data(sog1 + 678);
    const auto *sog1_680 = buffer.data(sog1 + 680);
    const auto *sog1_681 = buffer.data(sog1 + 681);
    const auto *sog1_684 = buffer.data(sog1 + 684);
    const auto *sog1_685 = buffer.data(sog1 + 685);
    const auto *sog1_687 = buffer.data(sog1 + 687);
    const auto *sog1_688 = buffer.data(sog1 + 688);
    const auto *sog1_689 = buffer.data(sog1 + 689);
    const auto *sog1_695 = buffer.data(sog1 + 695);
    const auto *sog1_699 = buffer.data(sog1 + 699);
    const auto *sog1_702 = buffer.data(sog1 + 702);
    const auto *sog1_703 = buffer.data(sog1 + 703);
    const auto *sog1_704 = buffer.data(sog1 + 704);
    const auto *sog1_705 = buffer.data(sog1 + 705);
    const auto *sog1_708 = buffer.data(sog1 + 708);
    const auto *sog1_710 = buffer.data(sog1 + 710);
    const auto *sog1_711 = buffer.data(sog1 + 711);
    const auto *sog1_714 = buffer.data(sog1 + 714);
    const auto *sog1_715 = buffer.data(sog1 + 715);
    const auto *sog1_717 = buffer.data(sog1 + 717);
    const auto *sog1_718 = buffer.data(sog1 + 718);
    const auto *sog1_719 = buffer.data(sog1 + 719);
    const auto *sog1_720 = buffer.data(sog1 + 720);
    const auto *sog1_723 = buffer.data(sog1 + 723);
    const auto *sog1_725 = buffer.data(sog1 + 725);
    const auto *sog1_726 = buffer.data(sog1 + 726);

    const auto *soh_930 = buffer.data(soh + 930);
    const auto *soh_933 = buffer.data(soh + 933);
    const auto *soh_934 = buffer.data(soh + 934);
    const auto *soh_936 = buffer.data(soh + 936);
    const auto *soh_938 = buffer.data(soh + 938);
    const auto *soh_939 = buffer.data(soh + 939);
    const auto *soh_940 = buffer.data(soh + 940);
    const auto *soh_941 = buffer.data(soh + 941);
    const auto *soh_942 = buffer.data(soh + 942);
    const auto *soh_943 = buffer.data(soh + 943);
    const auto *soh_944 = buffer.data(soh + 944);
    const auto *soh_945 = buffer.data(soh + 945);
    const auto *soh_947 = buffer.data(soh + 947);
    const auto *soh_948 = buffer.data(soh + 948);
    const auto *soh_950 = buffer.data(soh + 950);
    const auto *soh_951 = buffer.data(soh + 951);
    const auto *soh_954 = buffer.data(soh + 954);
    const auto *soh_955 = buffer.data(soh + 955);
    const auto *soh_957 = buffer.data(soh + 957);
    const auto *soh_959 = buffer.data(soh + 959);
    const auto *soh_960 = buffer.data(soh + 960);
    const auto *soh_961 = buffer.data(soh + 961);
    const auto *soh_962 = buffer.data(soh + 962);
    const auto *soh_963 = buffer.data(soh + 963);
    const auto *soh_964 = buffer.data(soh + 964);
    const auto *soh_965 = buffer.data(soh + 965);
    const auto *soh_966 = buffer.data(soh + 966);
    const auto *soh_968 = buffer.data(soh + 968);
    const auto *soh_969 = buffer.data(soh + 969);
    const auto *soh_971 = buffer.data(soh + 971);
    const auto *soh_972 = buffer.data(soh + 972);
    const auto *soh_975 = buffer.data(soh + 975);
    const auto *soh_980 = buffer.data(soh + 980);
    const auto *soh_981 = buffer.data(soh + 981);
    const auto *soh_982 = buffer.data(soh + 982);
    const auto *soh_983 = buffer.data(soh + 983);
    const auto *soh_984 = buffer.data(soh + 984);
    const auto *soh_985 = buffer.data(soh + 985);
    const auto *soh_986 = buffer.data(soh + 986);
    const auto *soh_987 = buffer.data(soh + 987);
    const auto *soh_989 = buffer.data(soh + 989);
    const auto *soh_990 = buffer.data(soh + 990);
    const auto *soh_992 = buffer.data(soh + 992);
    const auto *soh_993 = buffer.data(soh + 993);
    const auto *soh_996 = buffer.data(soh + 996);
    const auto *soh_997 = buffer.data(soh + 997);
    const auto *soh_999 = buffer.data(soh + 999);
    const auto *soh_1001 = buffer.data(soh + 1001);
    const auto *soh_1002 = buffer.data(soh + 1002);
    const auto *soh_1003 = buffer.data(soh + 1003);
    const auto *soh_1004 = buffer.data(soh + 1004);
    const auto *soh_1005 = buffer.data(soh + 1005);
    const auto *soh_1006 = buffer.data(soh + 1006);
    const auto *soh_1007 = buffer.data(soh + 1007);
    const auto *soh_1008 = buffer.data(soh + 1008);
    const auto *soh_1010 = buffer.data(soh + 1010);
    const auto *soh_1011 = buffer.data(soh + 1011);
    const auto *soh_1013 = buffer.data(soh + 1013);
    const auto *soh_1014 = buffer.data(soh + 1014);

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, pc_z, snh_741, snh_933, snh_934, \
                         sog0_669, sog0_670, sog1_669, sog1_670, soh_930, soh_933, \
                         soh_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_13 * snh_933[k]
                    + f_6 * sog0_669[k]
                    - f_7 * sog1_669[k]
                    + f_3 * pc_x[k] * soh_933[k];

        t_1242[k] = f_13 * snh_934[k]
                    + f_8 * sog0_670[k]
                    - f_9 * sog1_670[k]
                    + f_3 * pc_x[k] * soh_934[k];

        t_1243[k] = f_17 * snh_741[k]
                    + f_3 * pc_z[k] * soh_930[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, pc_x, pc_y, snh_936, snh_938, sog0_672, \
                         sog0_674, sog1_672, sog1_674, soh_933, soh_936, \
                         soh_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_13 * snh_936[k]
                    + f_8 * sog0_672[k]
                    - f_9 * sog1_672[k]
                    + f_3 * pc_x[k] * soh_936[k];

        t_1245[k] = f_3 * pc_y[k] * soh_933[k];

        t_1246[k] = f_13 * snh_938[k]
                    + f_8 * sog0_674[k]
                    - f_9 * sog1_674[k]
                    + f_3 * pc_x[k] * soh_938[k];
    }

#pragma omp simd aligned(t_1247, t_1248, t_1249, t_1250, t_1251, pc_x, snh_939, snh_940, \
                         snh_941, snh_942, snh_943, soh_939, soh_940, soh_941, soh_942, \
                         soh_943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1247[k] = f_13 * snh_939[k]
                    + f_3 * pc_x[k] * soh_939[k];

        t_1248[k] = f_13 * snh_940[k]
                    + f_3 * pc_x[k] * soh_940[k];

        t_1249[k] = f_13 * snh_941[k]
                    + f_3 * pc_x[k] * soh_941[k];

        t_1250[k] = f_13 * snh_942[k]
                    + f_3 * pc_x[k] * soh_942[k];

        t_1251[k] = f_13 * snh_943[k]
                    + f_3 * pc_x[k] * soh_943[k];
    }

#pragma omp simd aligned(t_1252, t_1253, t_1254, t_1255, pc_x, pc_y, pc_z, snh_750, snh_944, \
                         sog0_670, sog0_672, sog1_670, sog1_672, soh_939, soh_941, \
                         soh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1252[k] = f_13 * snh_944[k]
                    + f_3 * pc_x[k] * soh_944[k];

        t_1253[k] = f_1 * sog0_670[k]
                    - f_2 * sog1_670[k]
                    + f_3 * pc_y[k] * soh_939[k];

        t_1254[k] = f_17 * snh_750[k]
                    + f_3 * pc_z[k] * soh_939[k];

        t_1255[k] = f_4 * sog0_672[k]
                    - f_5 * sog1_672[k]
                    + f_3 * pc_y[k] * soh_941[k];
    }

#pragma omp simd aligned(t_1256, t_1257, t_1258, t_1259, pc_y, pc_z, snh_755, sog0_673, \
                         sog0_674, sog1_673, sog1_674, soh_942, soh_943, \
                         soh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1256[k] = f_6 * sog0_673[k]
                    - f_7 * sog1_673[k]
                    + f_3 * pc_y[k] * soh_942[k];

        t_1257[k] = f_8 * sog0_674[k]
                    - f_9 * sog1_674[k]
                    + f_3 * pc_y[k] * soh_943[k];

        t_1258[k] = f_3 * pc_y[k] * soh_944[k];

        t_1259[k] = f_17 * snh_755[k]
                    + f_1 * sog0_674[k]
                    - f_2 * sog1_674[k]
                    + f_3 * pc_z[k] * soh_944[k];
    }

#pragma omp simd aligned(t_1260, t_1261, t_1262, t_1263, pc_x, pc_y, pc_z, snh_756, snh_945, \
                         snh_948, sog0_675, sog0_678, sog1_675, sog1_678, soh_945, \
                         soh_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1260[k] = f_12 * snh_945[k]
                    + f_1 * sog0_675[k]
                    - f_2 * sog1_675[k]
                    + f_3 * pc_x[k] * soh_945[k];

        t_1261[k] = f_16 * snh_756[k]
                    + f_3 * pc_y[k] * soh_945[k];

        t_1262[k] = f_3 * pc_z[k] * soh_945[k];

        t_1263[k] = f_12 * snh_948[k]
                    + f_4 * sog0_678[k]
                    - f_5 * sog1_678[k]
                    + f_3 * pc_x[k] * soh_948[k];
    }

#pragma omp simd aligned(t_1264, t_1265, t_1266, pc_x, pc_y, snh_758, snh_950, snh_951, \
                         sog0_680, sog0_681, sog1_680, sog1_681, soh_947, soh_950, \
                         soh_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1264[k] = f_16 * snh_758[k]
                    + f_3 * pc_y[k] * soh_947[k];

        t_1265[k] = f_12 * snh_950[k]
                    + f_4 * sog0_680[k]
                    - f_5 * sog1_680[k]
                    + f_3 * pc_x[k] * soh_950[k];

        t_1266[k] = f_12 * snh_951[k]
                    + f_6 * sog0_681[k]
                    - f_7 * sog1_681[k]
                    + f_3 * pc_x[k] * soh_951[k];
    }

#pragma omp simd aligned(t_1267, t_1268, t_1269, pc_x, pc_y, pc_z, snh_761, snh_954, sog0_684, \
                         sog1_684, soh_948, soh_950, soh_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1267[k] = f_3 * pc_z[k] * soh_948[k];

        t_1268[k] = f_16 * snh_761[k]
                    + f_3 * pc_y[k] * soh_950[k];

        t_1269[k] = f_12 * snh_954[k]
                    + f_6 * sog0_684[k]
                    - f_7 * sog1_684[k]
                    + f_3 * pc_x[k] * soh_954[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, pc_x, pc_z, snh_955, snh_957, sog0_685, \
                         sog0_687, sog1_685, sog1_687, soh_951, soh_955, \
                         soh_957 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_12 * snh_955[k]
                    + f_8 * sog0_685[k]
                    - f_9 * sog1_685[k]
                    + f_3 * pc_x[k] * soh_955[k];

        t_1271[k] = f_3 * pc_z[k] * soh_951[k];

        t_1272[k] = f_12 * snh_957[k]
                    + f_8 * sog0_687[k]
                    - f_9 * sog1_687[k]
                    + f_3 * pc_x[k] * soh_957[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, pc_x, pc_y, snh_765, snh_959, \
                         snh_960, snh_961, sog0_689, sog1_689, soh_954, soh_959, soh_960, \
                         soh_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_16 * snh_765[k]
                    + f_3 * pc_y[k] * soh_954[k];

        t_1274[k] = f_12 * snh_959[k]
                    + f_8 * sog0_689[k]
                    - f_9 * sog1_689[k]
                    + f_3 * pc_x[k] * soh_959[k];

        t_1275[k] = f_12 * snh_960[k]
                    + f_3 * pc_x[k] * soh_960[k];

        t_1276[k] = f_12 * snh_961[k]
                    + f_3 * pc_x[k] * soh_961[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, t_1280, pc_x, snh_962, snh_963, snh_964, \
                         snh_965, soh_962, soh_963, soh_964, soh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = f_12 * snh_962[k]
                    + f_3 * pc_x[k] * soh_962[k];

        t_1278[k] = f_12 * snh_963[k]
                    + f_3 * pc_x[k] * soh_963[k];

        t_1279[k] = f_12 * snh_964[k]
                    + f_3 * pc_x[k] * soh_964[k];

        t_1280[k] = f_12 * snh_965[k]
                    + f_3 * pc_x[k] * soh_965[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pc_y, pc_z, snh_771, snh_773, sog0_685, \
                         sog0_687, sog1_685, sog1_687, soh_960, \
                         soh_962 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_16 * snh_771[k]
                    + f_1 * sog0_685[k]
                    - f_2 * sog1_685[k]
                    + f_3 * pc_y[k] * soh_960[k];

        t_1282[k] = f_3 * pc_z[k] * soh_960[k];

        t_1283[k] = f_16 * snh_773[k]
                    + f_4 * sog0_687[k]
                    - f_5 * sog1_687[k]
                    + f_3 * pc_y[k] * soh_962[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pc_y, pc_z, snh_774, snh_775, \
                         snh_776, sog0_688, sog0_689, sog1_688, sog1_689, soh_963, soh_964, \
                         soh_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_16 * snh_774[k]
                    + f_6 * sog0_688[k]
                    - f_7 * sog1_688[k]
                    + f_3 * pc_y[k] * soh_963[k];

        t_1285[k] = f_16 * snh_775[k]
                    + f_8 * sog0_689[k]
                    - f_9 * sog1_689[k]
                    + f_3 * pc_y[k] * soh_964[k];

        t_1286[k] = f_16 * snh_776[k]
                    + f_3 * pc_y[k] * soh_965[k];

        t_1287[k] = f_1 * sog0_689[k]
                    - f_2 * sog1_689[k]
                    + f_3 * pc_z[k] * soh_965[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, pb_z, pc_y, pc_z, sni0_1008, \
                         sni0_1011, snh_756, snh_777, sni1_1008, sni1_1011, \
                         soh_966 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = pb_z[k] * sni0_1008[k]
                    - f_10 * pc_z[k] * sni1_1008[k];

        t_1289[k] = f_17 * snh_777[k]
                    + f_3 * pc_y[k] * soh_966[k];

        t_1290[k] = f_11 * snh_756[k]
                    + f_3 * pc_z[k] * soh_966[k];

        t_1291[k] = pb_z[k] * sni0_1011[k]
                    - f_10 * pc_z[k] * sni1_1011[k];
    }

#pragma omp simd aligned(t_1292, t_1293, t_1294, pb_z, pc_x, pc_y, pc_z, sni0_1014, snh_779, \
                         snh_971, sni1_1014, sog0_695, sog1_695, soh_968, \
                         soh_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1292[k] = f_17 * snh_779[k]
                    + f_3 * pc_y[k] * soh_968[k];

        t_1293[k] = f_12 * snh_971[k]
                    + f_4 * sog0_695[k]
                    - f_5 * sog1_695[k]
                    + f_3 * pc_x[k] * soh_971[k];

        t_1294[k] = pb_z[k] * sni0_1014[k]
                    - f_10 * pc_z[k] * sni1_1014[k];
    }

#pragma omp simd aligned(t_1295, t_1296, t_1297, pc_x, pc_y, pc_z, snh_759, snh_782, snh_975, \
                         sog0_699, sog1_699, soh_969, soh_971, \
                         soh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1295[k] = f_11 * snh_759[k]
                    + f_3 * pc_z[k] * soh_969[k];

        t_1296[k] = f_17 * snh_782[k]
                    + f_3 * pc_y[k] * soh_971[k];

        t_1297[k] = f_12 * snh_975[k]
                    + f_6 * sog0_699[k]
                    - f_7 * sog1_699[k]
                    + f_3 * pc_x[k] * soh_975[k];
    }

#pragma omp simd aligned(t_1298, t_1299, t_1300, t_1301, pb_z, pc_y, pc_z, sni0_1018, \
                         sni0_1020, snh_762, snh_763, snh_786, sni1_1018, sni1_1020, soh_972, \
                         soh_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1298[k] = pb_z[k] * sni0_1018[k]
                    - f_10 * pc_z[k] * sni1_1018[k];

        t_1299[k] = f_11 * snh_762[k]
                    + f_3 * pc_z[k] * soh_972[k];

        t_1300[k] = pb_z[k] * sni0_1020[k]
                    + f_12 * snh_763[k]
                    - f_10 * pc_z[k] * sni1_1020[k];

        t_1301[k] = f_17 * snh_786[k]
                    + f_3 * pc_y[k] * soh_975[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, t_1305, pc_x, snh_980, snh_981, snh_982, \
                         snh_983, sog0_704, sog1_704, soh_980, soh_981, soh_982, \
                         soh_983 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_12 * snh_980[k]
                    + f_8 * sog0_704[k]
                    - f_9 * sog1_704[k]
                    + f_3 * pc_x[k] * soh_980[k];

        t_1303[k] = f_12 * snh_981[k]
                    + f_3 * pc_x[k] * soh_981[k];

        t_1304[k] = f_12 * snh_982[k]
                    + f_3 * pc_x[k] * soh_982[k];

        t_1305[k] = f_12 * snh_983[k]
                    + f_3 * pc_x[k] * soh_983[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, pb_z, pc_x, pc_z, sni0_1029, snh_984, \
                         snh_985, snh_986, sni1_1029, soh_984, soh_985, \
                         soh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_12 * snh_984[k]
                    + f_3 * pc_x[k] * soh_984[k];

        t_1307[k] = f_12 * snh_985[k]
                    + f_3 * pc_x[k] * soh_985[k];

        t_1308[k] = f_12 * snh_986[k]
                    + f_3 * pc_x[k] * soh_986[k];

        t_1309[k] = pb_z[k] * sni0_1029[k]
                    - f_10 * pc_z[k] * sni1_1029[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, pc_y, pc_z, snh_771, snh_794, snh_795, \
                         sog0_702, sog0_703, sog1_702, sog1_703, soh_981, soh_983, \
                         soh_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_11 * snh_771[k]
                    + f_3 * pc_z[k] * soh_981[k];

        t_1311[k] = f_17 * snh_794[k]
                    + f_4 * sog0_702[k]
                    - f_5 * sog1_702[k]
                    + f_3 * pc_y[k] * soh_983[k];

        t_1312[k] = f_17 * snh_795[k]
                    + f_6 * sog0_703[k]
                    - f_7 * sog1_703[k]
                    + f_3 * pc_y[k] * soh_984[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pc_y, pc_z, snh_776, snh_796, snh_797, \
                         sog0_704, sog1_704, soh_985, soh_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_17 * snh_796[k]
                    + f_8 * sog0_704[k]
                    - f_9 * sog1_704[k]
                    + f_3 * pc_y[k] * soh_985[k];

        t_1314[k] = f_17 * snh_797[k]
                    + f_3 * pc_y[k] * soh_986[k];

        t_1315[k] = f_11 * snh_776[k]
                    + f_1 * sog0_704[k]
                    - f_2 * sog1_704[k]
                    + f_3 * pc_z[k] * soh_986[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, pc_x, pc_y, pc_z, snh_777, snh_798, snh_987, \
                         sog0_705, sog1_705, soh_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_12 * snh_987[k]
                    + f_1 * sog0_705[k]
                    - f_2 * sog1_705[k]
                    + f_3 * pc_x[k] * soh_987[k];

        t_1317[k] = f_18 * snh_798[k]
                    + f_3 * pc_y[k] * soh_987[k];

        t_1318[k] = f_12 * snh_777[k]
                    + f_3 * pc_z[k] * soh_987[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pc_x, pc_y, snh_800, snh_990, snh_992, \
                         sog0_708, sog0_710, sog1_708, sog1_710, soh_989, soh_990, \
                         soh_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = f_12 * snh_990[k]
                    + f_4 * sog0_708[k]
                    - f_5 * sog1_708[k]
                    + f_3 * pc_x[k] * soh_990[k];

        t_1320[k] = f_18 * snh_800[k]
                    + f_3 * pc_y[k] * soh_989[k];

        t_1321[k] = f_12 * snh_992[k]
                    + f_4 * sog0_710[k]
                    - f_5 * sog1_710[k]
                    + f_3 * pc_x[k] * soh_992[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pc_x, pc_y, pc_z, snh_780, snh_803, snh_993, \
                         sog0_711, sog1_711, soh_990, soh_992, \
                         soh_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_12 * snh_993[k]
                    + f_6 * sog0_711[k]
                    - f_7 * sog1_711[k]
                    + f_3 * pc_x[k] * soh_993[k];

        t_1323[k] = f_12 * snh_780[k]
                    + f_3 * pc_z[k] * soh_990[k];

        t_1324[k] = f_18 * snh_803[k]
                    + f_3 * pc_y[k] * soh_992[k];
    }

#pragma omp simd aligned(t_1325, t_1326, t_1327, pc_x, pc_z, snh_783, snh_996, snh_997, \
                         sog0_714, sog0_715, sog1_714, sog1_715, soh_993, soh_996, \
                         soh_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = f_12 * snh_996[k]
                    + f_6 * sog0_714[k]
                    - f_7 * sog1_714[k]
                    + f_3 * pc_x[k] * soh_996[k];

        t_1326[k] = f_12 * snh_997[k]
                    + f_8 * sog0_715[k]
                    - f_9 * sog1_715[k]
                    + f_3 * pc_x[k] * soh_997[k];

        t_1327[k] = f_12 * snh_783[k]
                    + f_3 * pc_z[k] * soh_993[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pc_x, pc_y, snh_807, snh_999, snh_1001, \
                         sog0_717, sog0_719, sog1_717, sog1_719, soh_996, soh_999, \
                         soh_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = f_12 * snh_999[k]
                    + f_8 * sog0_717[k]
                    - f_9 * sog1_717[k]
                    + f_3 * pc_x[k] * soh_999[k];

        t_1329[k] = f_18 * snh_807[k]
                    + f_3 * pc_y[k] * soh_996[k];

        t_1330[k] = f_12 * snh_1001[k]
                    + f_8 * sog0_719[k]
                    - f_9 * sog1_719[k]
                    + f_3 * pc_x[k] * soh_1001[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, t_1335, pc_x, snh_1002, snh_1003, \
                         snh_1004, snh_1005, snh_1006, soh_1002, soh_1003, soh_1004, soh_1005, \
                         soh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_12 * snh_1002[k]
                    + f_3 * pc_x[k] * soh_1002[k];

        t_1332[k] = f_12 * snh_1003[k]
                    + f_3 * pc_x[k] * soh_1003[k];

        t_1333[k] = f_12 * snh_1004[k]
                    + f_3 * pc_x[k] * soh_1004[k];

        t_1334[k] = f_12 * snh_1005[k]
                    + f_3 * pc_x[k] * soh_1005[k];

        t_1335[k] = f_12 * snh_1006[k]
                    + f_3 * pc_x[k] * soh_1006[k];
    }

#pragma omp simd aligned(t_1336, t_1337, t_1338, pc_x, pc_y, pc_z, snh_792, snh_813, snh_1007, \
                         sog0_715, sog1_715, soh_1002, soh_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1336[k] = f_12 * snh_1007[k]
                    + f_3 * pc_x[k] * soh_1007[k];

        t_1337[k] = f_18 * snh_813[k]
                    + f_1 * sog0_715[k]
                    - f_2 * sog1_715[k]
                    + f_3 * pc_y[k] * soh_1002[k];

        t_1338[k] = f_12 * snh_792[k]
                    + f_3 * pc_z[k] * soh_1002[k];
    }

#pragma omp simd aligned(t_1339, t_1340, t_1341, pc_y, snh_815, snh_816, snh_817, sog0_717, \
                         sog0_718, sog0_719, sog1_717, sog1_718, sog1_719, soh_1004, soh_1005, \
                         soh_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1339[k] = f_18 * snh_815[k]
                    + f_4 * sog0_717[k]
                    - f_5 * sog1_717[k]
                    + f_3 * pc_y[k] * soh_1004[k];

        t_1340[k] = f_18 * snh_816[k]
                    + f_6 * sog0_718[k]
                    - f_7 * sog1_718[k]
                    + f_3 * pc_y[k] * soh_1005[k];

        t_1341[k] = f_18 * snh_817[k]
                    + f_8 * sog0_719[k]
                    - f_9 * sog1_719[k]
                    + f_3 * pc_y[k] * soh_1006[k];
    }

#pragma omp simd aligned(t_1342, t_1343, t_1344, pc_x, pc_y, pc_z, snh_797, snh_818, snh_1008, \
                         sog0_719, sog0_720, sog1_719, sog1_720, soh_1007, \
                         soh_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1342[k] = f_18 * snh_818[k]
                    + f_3 * pc_y[k] * soh_1007[k];

        t_1343[k] = f_12 * snh_797[k]
                    + f_1 * sog0_719[k]
                    - f_2 * sog1_719[k]
                    + f_3 * pc_z[k] * soh_1007[k];

        t_1344[k] = f_12 * snh_1008[k]
                    + f_1 * sog0_720[k]
                    - f_2 * sog1_720[k]
                    + f_3 * pc_x[k] * soh_1008[k];
    }

#pragma omp simd aligned(t_1345, t_1346, t_1347, t_1348, pc_x, pc_y, pc_z, snh_798, snh_819, \
                         snh_821, snh_1011, sog0_723, sog1_723, soh_1008, soh_1010, \
                         soh_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1345[k] = f_19 * snh_819[k]
                    + f_3 * pc_y[k] * soh_1008[k];

        t_1346[k] = f_13 * snh_798[k]
                    + f_3 * pc_z[k] * soh_1008[k];

        t_1347[k] = f_12 * snh_1011[k]
                    + f_4 * sog0_723[k]
                    - f_5 * sog1_723[k]
                    + f_3 * pc_x[k] * soh_1011[k];

        t_1348[k] = f_19 * snh_821[k]
                    + f_3 * pc_y[k] * soh_1010[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, pc_x, pc_z, snh_801, snh_1013, snh_1014, \
                         sog0_725, sog0_726, sog1_725, sog1_726, soh_1011, soh_1013, \
                         soh_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_12 * snh_1013[k]
                    + f_4 * sog0_725[k]
                    - f_5 * sog1_725[k]
                    + f_3 * pc_x[k] * soh_1013[k];

        t_1350[k] = f_12 * snh_1014[k]
                    + f_6 * sog0_726[k]
                    - f_7 * sog1_726[k]
                    + f_3 * pc_x[k] * soh_1014[k];

        t_1351[k] = f_13 * snh_801[k]
                    + f_3 * pc_z[k] * soh_1011[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t snh, const size_t sog0,
                                                           const size_t sog1, const size_t soh,
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
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh_804 = buffer.data(snh + 804);
    const auto *snh_813 = buffer.data(snh + 813);
    const auto *snh_818 = buffer.data(snh + 818);
    const auto *snh_819 = buffer.data(snh + 819);
    const auto *snh_822 = buffer.data(snh + 822);
    const auto *snh_824 = buffer.data(snh + 824);
    const auto *snh_825 = buffer.data(snh + 825);
    const auto *snh_828 = buffer.data(snh + 828);
    const auto *snh_834 = buffer.data(snh + 834);
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
    const auto *snh_855 = buffer.data(snh + 855);
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
    const auto *snh_876 = buffer.data(snh + 876);
    const auto *snh_878 = buffer.data(snh + 878);
    const auto *snh_879 = buffer.data(snh + 879);
    const auto *snh_880 = buffer.data(snh + 880);
    const auto *snh_881 = buffer.data(snh + 881);
    const auto *snh_882 = buffer.data(snh + 882);
    const auto *snh_884 = buffer.data(snh + 884);
    const auto *snh_887 = buffer.data(snh + 887);
    const auto *snh_891 = buffer.data(snh + 891);
    const auto *snh_897 = buffer.data(snh + 897);
    const auto *snh_899 = buffer.data(snh + 899);
    const auto *snh_900 = buffer.data(snh + 900);
    const auto *snh_901 = buffer.data(snh + 901);
    const auto *snh_902 = buffer.data(snh + 902);
    const auto *snh_1017 = buffer.data(snh + 1017);
    const auto *snh_1018 = buffer.data(snh + 1018);
    const auto *snh_1020 = buffer.data(snh + 1020);
    const auto *snh_1022 = buffer.data(snh + 1022);
    const auto *snh_1023 = buffer.data(snh + 1023);
    const auto *snh_1024 = buffer.data(snh + 1024);
    const auto *snh_1025 = buffer.data(snh + 1025);
    const auto *snh_1026 = buffer.data(snh + 1026);
    const auto *snh_1027 = buffer.data(snh + 1027);
    const auto *snh_1028 = buffer.data(snh + 1028);
    const auto *snh_1029 = buffer.data(snh + 1029);
    const auto *snh_1032 = buffer.data(snh + 1032);
    const auto *snh_1034 = buffer.data(snh + 1034);
    const auto *snh_1035 = buffer.data(snh + 1035);
    const auto *snh_1038 = buffer.data(snh + 1038);
    const auto *snh_1039 = buffer.data(snh + 1039);
    const auto *snh_1041 = buffer.data(snh + 1041);
    const auto *snh_1043 = buffer.data(snh + 1043);
    const auto *snh_1044 = buffer.data(snh + 1044);
    const auto *snh_1045 = buffer.data(snh + 1045);
    const auto *snh_1046 = buffer.data(snh + 1046);
    const auto *snh_1047 = buffer.data(snh + 1047);
    const auto *snh_1048 = buffer.data(snh + 1048);
    const auto *snh_1049 = buffer.data(snh + 1049);
    const auto *snh_1050 = buffer.data(snh + 1050);
    const auto *snh_1053 = buffer.data(snh + 1053);
    const auto *snh_1055 = buffer.data(snh + 1055);
    const auto *snh_1056 = buffer.data(snh + 1056);
    const auto *snh_1059 = buffer.data(snh + 1059);
    const auto *snh_1060 = buffer.data(snh + 1060);
    const auto *snh_1062 = buffer.data(snh + 1062);
    const auto *snh_1064 = buffer.data(snh + 1064);
    const auto *snh_1065 = buffer.data(snh + 1065);
    const auto *snh_1066 = buffer.data(snh + 1066);
    const auto *snh_1067 = buffer.data(snh + 1067);
    const auto *snh_1068 = buffer.data(snh + 1068);
    const auto *snh_1069 = buffer.data(snh + 1069);
    const auto *snh_1070 = buffer.data(snh + 1070);
    const auto *snh_1071 = buffer.data(snh + 1071);
    const auto *snh_1074 = buffer.data(snh + 1074);
    const auto *snh_1076 = buffer.data(snh + 1076);
    const auto *snh_1077 = buffer.data(snh + 1077);
    const auto *snh_1080 = buffer.data(snh + 1080);
    const auto *snh_1081 = buffer.data(snh + 1081);
    const auto *snh_1083 = buffer.data(snh + 1083);
    const auto *snh_1085 = buffer.data(snh + 1085);
    const auto *snh_1086 = buffer.data(snh + 1086);
    const auto *snh_1087 = buffer.data(snh + 1087);
    const auto *snh_1088 = buffer.data(snh + 1088);
    const auto *snh_1089 = buffer.data(snh + 1089);
    const auto *snh_1090 = buffer.data(snh + 1090);
    const auto *snh_1091 = buffer.data(snh + 1091);
    const auto *snh_1092 = buffer.data(snh + 1092);

    const auto *sog0_729 = buffer.data(sog0 + 729);
    const auto *sog0_730 = buffer.data(sog0 + 730);
    const auto *sog0_732 = buffer.data(sog0 + 732);
    const auto *sog0_733 = buffer.data(sog0 + 733);
    const auto *sog0_734 = buffer.data(sog0 + 734);
    const auto *sog0_735 = buffer.data(sog0 + 735);
    const auto *sog0_738 = buffer.data(sog0 + 738);
    const auto *sog0_740 = buffer.data(sog0 + 740);
    const auto *sog0_741 = buffer.data(sog0 + 741);
    const auto *sog0_744 = buffer.data(sog0 + 744);
    const auto *sog0_745 = buffer.data(sog0 + 745);
    const auto *sog0_747 = buffer.data(sog0 + 747);
    const auto *sog0_748 = buffer.data(sog0 + 748);
    const auto *sog0_749 = buffer.data(sog0 + 749);
    const auto *sog0_750 = buffer.data(sog0 + 750);
    const auto *sog0_753 = buffer.data(sog0 + 753);
    const auto *sog0_755 = buffer.data(sog0 + 755);
    const auto *sog0_756 = buffer.data(sog0 + 756);
    const auto *sog0_759 = buffer.data(sog0 + 759);
    const auto *sog0_760 = buffer.data(sog0 + 760);
    const auto *sog0_762 = buffer.data(sog0 + 762);
    const auto *sog0_763 = buffer.data(sog0 + 763);
    const auto *sog0_764 = buffer.data(sog0 + 764);
    const auto *sog0_765 = buffer.data(sog0 + 765);
    const auto *sog0_768 = buffer.data(sog0 + 768);
    const auto *sog0_770 = buffer.data(sog0 + 770);
    const auto *sog0_771 = buffer.data(sog0 + 771);
    const auto *sog0_774 = buffer.data(sog0 + 774);
    const auto *sog0_775 = buffer.data(sog0 + 775);
    const auto *sog0_777 = buffer.data(sog0 + 777);
    const auto *sog0_778 = buffer.data(sog0 + 778);
    const auto *sog0_779 = buffer.data(sog0 + 779);
    const auto *sog0_780 = buffer.data(sog0 + 780);

    const auto *sog1_729 = buffer.data(sog1 + 729);
    const auto *sog1_730 = buffer.data(sog1 + 730);
    const auto *sog1_732 = buffer.data(sog1 + 732);
    const auto *sog1_733 = buffer.data(sog1 + 733);
    const auto *sog1_734 = buffer.data(sog1 + 734);
    const auto *sog1_735 = buffer.data(sog1 + 735);
    const auto *sog1_738 = buffer.data(sog1 + 738);
    const auto *sog1_740 = buffer.data(sog1 + 740);
    const auto *sog1_741 = buffer.data(sog1 + 741);
    const auto *sog1_744 = buffer.data(sog1 + 744);
    const auto *sog1_745 = buffer.data(sog1 + 745);
    const auto *sog1_747 = buffer.data(sog1 + 747);
    const auto *sog1_748 = buffer.data(sog1 + 748);
    const auto *sog1_749 = buffer.data(sog1 + 749);
    const auto *sog1_750 = buffer.data(sog1 + 750);
    const auto *sog1_753 = buffer.data(sog1 + 753);
    const auto *sog1_755 = buffer.data(sog1 + 755);
    const auto *sog1_756 = buffer.data(sog1 + 756);
    const auto *sog1_759 = buffer.data(sog1 + 759);
    const auto *sog1_760 = buffer.data(sog1 + 760);
    const auto *sog1_762 = buffer.data(sog1 + 762);
    const auto *sog1_763 = buffer.data(sog1 + 763);
    const auto *sog1_764 = buffer.data(sog1 + 764);
    const auto *sog1_765 = buffer.data(sog1 + 765);
    const auto *sog1_768 = buffer.data(sog1 + 768);
    const auto *sog1_770 = buffer.data(sog1 + 770);
    const auto *sog1_771 = buffer.data(sog1 + 771);
    const auto *sog1_774 = buffer.data(sog1 + 774);
    const auto *sog1_775 = buffer.data(sog1 + 775);
    const auto *sog1_777 = buffer.data(sog1 + 777);
    const auto *sog1_778 = buffer.data(sog1 + 778);
    const auto *sog1_779 = buffer.data(sog1 + 779);
    const auto *sog1_780 = buffer.data(sog1 + 780);

    const auto *soh_1013 = buffer.data(soh + 1013);
    const auto *soh_1014 = buffer.data(soh + 1014);
    const auto *soh_1017 = buffer.data(soh + 1017);
    const auto *soh_1018 = buffer.data(soh + 1018);
    const auto *soh_1020 = buffer.data(soh + 1020);
    const auto *soh_1022 = buffer.data(soh + 1022);
    const auto *soh_1023 = buffer.data(soh + 1023);
    const auto *soh_1024 = buffer.data(soh + 1024);
    const auto *soh_1025 = buffer.data(soh + 1025);
    const auto *soh_1026 = buffer.data(soh + 1026);
    const auto *soh_1027 = buffer.data(soh + 1027);
    const auto *soh_1028 = buffer.data(soh + 1028);
    const auto *soh_1029 = buffer.data(soh + 1029);
    const auto *soh_1031 = buffer.data(soh + 1031);
    const auto *soh_1032 = buffer.data(soh + 1032);
    const auto *soh_1034 = buffer.data(soh + 1034);
    const auto *soh_1035 = buffer.data(soh + 1035);
    const auto *soh_1038 = buffer.data(soh + 1038);
    const auto *soh_1039 = buffer.data(soh + 1039);
    const auto *soh_1041 = buffer.data(soh + 1041);
    const auto *soh_1043 = buffer.data(soh + 1043);
    const auto *soh_1044 = buffer.data(soh + 1044);
    const auto *soh_1045 = buffer.data(soh + 1045);
    const auto *soh_1046 = buffer.data(soh + 1046);
    const auto *soh_1047 = buffer.data(soh + 1047);
    const auto *soh_1048 = buffer.data(soh + 1048);
    const auto *soh_1049 = buffer.data(soh + 1049);
    const auto *soh_1050 = buffer.data(soh + 1050);
    const auto *soh_1052 = buffer.data(soh + 1052);
    const auto *soh_1053 = buffer.data(soh + 1053);
    const auto *soh_1055 = buffer.data(soh + 1055);
    const auto *soh_1056 = buffer.data(soh + 1056);
    const auto *soh_1059 = buffer.data(soh + 1059);
    const auto *soh_1060 = buffer.data(soh + 1060);
    const auto *soh_1062 = buffer.data(soh + 1062);
    const auto *soh_1064 = buffer.data(soh + 1064);
    const auto *soh_1065 = buffer.data(soh + 1065);
    const auto *soh_1066 = buffer.data(soh + 1066);
    const auto *soh_1067 = buffer.data(soh + 1067);
    const auto *soh_1068 = buffer.data(soh + 1068);
    const auto *soh_1069 = buffer.data(soh + 1069);
    const auto *soh_1070 = buffer.data(soh + 1070);
    const auto *soh_1071 = buffer.data(soh + 1071);
    const auto *soh_1073 = buffer.data(soh + 1073);
    const auto *soh_1074 = buffer.data(soh + 1074);
    const auto *soh_1076 = buffer.data(soh + 1076);
    const auto *soh_1077 = buffer.data(soh + 1077);
    const auto *soh_1080 = buffer.data(soh + 1080);
    const auto *soh_1081 = buffer.data(soh + 1081);
    const auto *soh_1083 = buffer.data(soh + 1083);
    const auto *soh_1085 = buffer.data(soh + 1085);
    const auto *soh_1086 = buffer.data(soh + 1086);
    const auto *soh_1087 = buffer.data(soh + 1087);
    const auto *soh_1088 = buffer.data(soh + 1088);
    const auto *soh_1089 = buffer.data(soh + 1089);
    const auto *soh_1090 = buffer.data(soh + 1090);
    const auto *soh_1091 = buffer.data(soh + 1091);
    const auto *soh_1092 = buffer.data(soh + 1092);

#pragma omp simd aligned(t_1352, t_1353, t_1354, pc_x, pc_y, snh_824, snh_1017, snh_1018, \
                         sog0_729, sog0_730, sog1_729, sog1_730, soh_1013, soh_1017, \
                         soh_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_19 * snh_824[k]
                    + f_3 * pc_y[k] * soh_1013[k];

        t_1353[k] = f_12 * snh_1017[k]
                    + f_6 * sog0_729[k]
                    - f_7 * sog1_729[k]
                    + f_3 * pc_x[k] * soh_1017[k];

        t_1354[k] = f_12 * snh_1018[k]
                    + f_8 * sog0_730[k]
                    - f_9 * sog1_730[k]
                    + f_3 * pc_x[k] * soh_1018[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pc_x, pc_y, pc_z, snh_804, snh_828, snh_1020, \
                         sog0_732, sog1_732, soh_1014, soh_1017, \
                         soh_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_13 * snh_804[k]
                    + f_3 * pc_z[k] * soh_1014[k];

        t_1356[k] = f_12 * snh_1020[k]
                    + f_8 * sog0_732[k]
                    - f_9 * sog1_732[k]
                    + f_3 * pc_x[k] * soh_1020[k];

        t_1357[k] = f_19 * snh_828[k]
                    + f_3 * pc_y[k] * soh_1017[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, t_1361, pc_x, snh_1022, snh_1023, snh_1024, \
                         snh_1025, sog0_734, sog1_734, soh_1022, soh_1023, soh_1024, \
                         soh_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_12 * snh_1022[k]
                    + f_8 * sog0_734[k]
                    - f_9 * sog1_734[k]
                    + f_3 * pc_x[k] * soh_1022[k];

        t_1359[k] = f_12 * snh_1023[k]
                    + f_3 * pc_x[k] * soh_1023[k];

        t_1360[k] = f_12 * snh_1024[k]
                    + f_3 * pc_x[k] * soh_1024[k];

        t_1361[k] = f_12 * snh_1025[k]
                    + f_3 * pc_x[k] * soh_1025[k];
    }

#pragma omp simd aligned(t_1362, t_1363, t_1364, t_1365, pc_x, pc_y, snh_834, snh_1026, \
                         snh_1027, snh_1028, sog0_730, sog1_730, soh_1023, soh_1026, soh_1027, \
                         soh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_12 * snh_1026[k]
                    + f_3 * pc_x[k] * soh_1026[k];

        t_1363[k] = f_12 * snh_1027[k]
                    + f_3 * pc_x[k] * soh_1027[k];

        t_1364[k] = f_12 * snh_1028[k]
                    + f_3 * pc_x[k] * soh_1028[k];

        t_1365[k] = f_19 * snh_834[k]
                    + f_1 * sog0_730[k]
                    - f_2 * sog1_730[k]
                    + f_3 * pc_y[k] * soh_1023[k];
    }

#pragma omp simd aligned(t_1366, t_1367, t_1368, pc_y, pc_z, snh_813, snh_836, snh_837, \
                         sog0_732, sog0_733, sog1_732, sog1_733, soh_1023, soh_1025, \
                         soh_1026 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1366[k] = f_13 * snh_813[k]
                    + f_3 * pc_z[k] * soh_1023[k];

        t_1367[k] = f_19 * snh_836[k]
                    + f_4 * sog0_732[k]
                    - f_5 * sog1_732[k]
                    + f_3 * pc_y[k] * soh_1025[k];

        t_1368[k] = f_19 * snh_837[k]
                    + f_6 * sog0_733[k]
                    - f_7 * sog1_733[k]
                    + f_3 * pc_y[k] * soh_1026[k];
    }

#pragma omp simd aligned(t_1369, t_1370, t_1371, pc_y, pc_z, snh_818, snh_838, snh_839, \
                         sog0_734, sog1_734, soh_1027, soh_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1369[k] = f_19 * snh_838[k]
                    + f_8 * sog0_734[k]
                    - f_9 * sog1_734[k]
                    + f_3 * pc_y[k] * soh_1027[k];

        t_1370[k] = f_19 * snh_839[k]
                    + f_3 * pc_y[k] * soh_1028[k];

        t_1371[k] = f_13 * snh_818[k]
                    + f_1 * sog0_734[k]
                    - f_2 * sog1_734[k]
                    + f_3 * pc_z[k] * soh_1028[k];
    }

#pragma omp simd aligned(t_1372, t_1373, t_1374, pc_x, pc_y, pc_z, snh_819, snh_840, snh_1029, \
                         sog0_735, sog1_735, soh_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1372[k] = f_12 * snh_1029[k]
                    + f_1 * sog0_735[k]
                    - f_2 * sog1_735[k]
                    + f_3 * pc_x[k] * soh_1029[k];

        t_1373[k] = f_20 * snh_840[k]
                    + f_3 * pc_y[k] * soh_1029[k];

        t_1374[k] = f_14 * snh_819[k]
                    + f_3 * pc_z[k] * soh_1029[k];
    }

#pragma omp simd aligned(t_1375, t_1376, t_1377, pc_x, pc_y, snh_842, snh_1032, snh_1034, \
                         sog0_738, sog0_740, sog1_738, sog1_740, soh_1031, soh_1032, \
                         soh_1034 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1375[k] = f_12 * snh_1032[k]
                    + f_4 * sog0_738[k]
                    - f_5 * sog1_738[k]
                    + f_3 * pc_x[k] * soh_1032[k];

        t_1376[k] = f_20 * snh_842[k]
                    + f_3 * pc_y[k] * soh_1031[k];

        t_1377[k] = f_12 * snh_1034[k]
                    + f_4 * sog0_740[k]
                    - f_5 * sog1_740[k]
                    + f_3 * pc_x[k] * soh_1034[k];
    }

#pragma omp simd aligned(t_1378, t_1379, t_1380, pc_x, pc_y, pc_z, snh_822, snh_845, snh_1035, \
                         sog0_741, sog1_741, soh_1032, soh_1034, \
                         soh_1035 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1378[k] = f_12 * snh_1035[k]
                    + f_6 * sog0_741[k]
                    - f_7 * sog1_741[k]
                    + f_3 * pc_x[k] * soh_1035[k];

        t_1379[k] = f_14 * snh_822[k]
                    + f_3 * pc_z[k] * soh_1032[k];

        t_1380[k] = f_20 * snh_845[k]
                    + f_3 * pc_y[k] * soh_1034[k];
    }

#pragma omp simd aligned(t_1381, t_1382, t_1383, pc_x, pc_z, snh_825, snh_1038, snh_1039, \
                         sog0_744, sog0_745, sog1_744, sog1_745, soh_1035, soh_1038, \
                         soh_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1381[k] = f_12 * snh_1038[k]
                    + f_6 * sog0_744[k]
                    - f_7 * sog1_744[k]
                    + f_3 * pc_x[k] * soh_1038[k];

        t_1382[k] = f_12 * snh_1039[k]
                    + f_8 * sog0_745[k]
                    - f_9 * sog1_745[k]
                    + f_3 * pc_x[k] * soh_1039[k];

        t_1383[k] = f_14 * snh_825[k]
                    + f_3 * pc_z[k] * soh_1035[k];
    }

#pragma omp simd aligned(t_1384, t_1385, t_1386, pc_x, pc_y, snh_849, snh_1041, snh_1043, \
                         sog0_747, sog0_749, sog1_747, sog1_749, soh_1038, soh_1041, \
                         soh_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1384[k] = f_12 * snh_1041[k]
                    + f_8 * sog0_747[k]
                    - f_9 * sog1_747[k]
                    + f_3 * pc_x[k] * soh_1041[k];

        t_1385[k] = f_20 * snh_849[k]
                    + f_3 * pc_y[k] * soh_1038[k];

        t_1386[k] = f_12 * snh_1043[k]
                    + f_8 * sog0_749[k]
                    - f_9 * sog1_749[k]
                    + f_3 * pc_x[k] * soh_1043[k];
    }

#pragma omp simd aligned(t_1387, t_1388, t_1389, t_1390, t_1391, pc_x, snh_1044, snh_1045, \
                         snh_1046, snh_1047, snh_1048, soh_1044, soh_1045, soh_1046, soh_1047, \
                         soh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1387[k] = f_12 * snh_1044[k]
                    + f_3 * pc_x[k] * soh_1044[k];

        t_1388[k] = f_12 * snh_1045[k]
                    + f_3 * pc_x[k] * soh_1045[k];

        t_1389[k] = f_12 * snh_1046[k]
                    + f_3 * pc_x[k] * soh_1046[k];

        t_1390[k] = f_12 * snh_1047[k]
                    + f_3 * pc_x[k] * soh_1047[k];

        t_1391[k] = f_12 * snh_1048[k]
                    + f_3 * pc_x[k] * soh_1048[k];
    }

#pragma omp simd aligned(t_1392, t_1393, t_1394, pc_x, pc_y, pc_z, snh_834, snh_855, snh_1049, \
                         sog0_745, sog1_745, soh_1044, soh_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1392[k] = f_12 * snh_1049[k]
                    + f_3 * pc_x[k] * soh_1049[k];

        t_1393[k] = f_20 * snh_855[k]
                    + f_1 * sog0_745[k]
                    - f_2 * sog1_745[k]
                    + f_3 * pc_y[k] * soh_1044[k];

        t_1394[k] = f_14 * snh_834[k]
                    + f_3 * pc_z[k] * soh_1044[k];
    }

#pragma omp simd aligned(t_1395, t_1396, t_1397, pc_y, snh_857, snh_858, snh_859, sog0_747, \
                         sog0_748, sog0_749, sog1_747, sog1_748, sog1_749, soh_1046, soh_1047, \
                         soh_1048 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1395[k] = f_20 * snh_857[k]
                    + f_4 * sog0_747[k]
                    - f_5 * sog1_747[k]
                    + f_3 * pc_y[k] * soh_1046[k];

        t_1396[k] = f_20 * snh_858[k]
                    + f_6 * sog0_748[k]
                    - f_7 * sog1_748[k]
                    + f_3 * pc_y[k] * soh_1047[k];

        t_1397[k] = f_20 * snh_859[k]
                    + f_8 * sog0_749[k]
                    - f_9 * sog1_749[k]
                    + f_3 * pc_y[k] * soh_1048[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, pc_x, pc_y, pc_z, snh_839, snh_860, snh_1050, \
                         sog0_749, sog0_750, sog1_749, sog1_750, soh_1049, \
                         soh_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_20 * snh_860[k]
                    + f_3 * pc_y[k] * soh_1049[k];

        t_1399[k] = f_14 * snh_839[k]
                    + f_1 * sog0_749[k]
                    - f_2 * sog1_749[k]
                    + f_3 * pc_z[k] * soh_1049[k];

        t_1400[k] = f_12 * snh_1050[k]
                    + f_1 * sog0_750[k]
                    - f_2 * sog1_750[k]
                    + f_3 * pc_x[k] * soh_1050[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, t_1404, pc_x, pc_y, pc_z, snh_840, snh_861, \
                         snh_863, snh_1053, sog0_753, sog1_753, soh_1050, soh_1052, \
                         soh_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = f_14 * snh_861[k]
                    + f_3 * pc_y[k] * soh_1050[k];

        t_1402[k] = f_20 * snh_840[k]
                    + f_3 * pc_z[k] * soh_1050[k];

        t_1403[k] = f_12 * snh_1053[k]
                    + f_4 * sog0_753[k]
                    - f_5 * sog1_753[k]
                    + f_3 * pc_x[k] * soh_1053[k];

        t_1404[k] = f_14 * snh_863[k]
                    + f_3 * pc_y[k] * soh_1052[k];
    }

#pragma omp simd aligned(t_1405, t_1406, t_1407, pc_x, pc_z, snh_843, snh_1055, snh_1056, \
                         sog0_755, sog0_756, sog1_755, sog1_756, soh_1053, soh_1055, \
                         soh_1056 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1405[k] = f_12 * snh_1055[k]
                    + f_4 * sog0_755[k]
                    - f_5 * sog1_755[k]
                    + f_3 * pc_x[k] * soh_1055[k];

        t_1406[k] = f_12 * snh_1056[k]
                    + f_6 * sog0_756[k]
                    - f_7 * sog1_756[k]
                    + f_3 * pc_x[k] * soh_1056[k];

        t_1407[k] = f_20 * snh_843[k]
                    + f_3 * pc_z[k] * soh_1053[k];
    }

#pragma omp simd aligned(t_1408, t_1409, t_1410, pc_x, pc_y, snh_866, snh_1059, snh_1060, \
                         sog0_759, sog0_760, sog1_759, sog1_760, soh_1055, soh_1059, \
                         soh_1060 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1408[k] = f_14 * snh_866[k]
                    + f_3 * pc_y[k] * soh_1055[k];

        t_1409[k] = f_12 * snh_1059[k]
                    + f_6 * sog0_759[k]
                    - f_7 * sog1_759[k]
                    + f_3 * pc_x[k] * soh_1059[k];

        t_1410[k] = f_12 * snh_1060[k]
                    + f_8 * sog0_760[k]
                    - f_9 * sog1_760[k]
                    + f_3 * pc_x[k] * soh_1060[k];
    }

#pragma omp simd aligned(t_1411, t_1412, t_1413, pc_x, pc_y, pc_z, snh_846, snh_870, snh_1062, \
                         sog0_762, sog1_762, soh_1056, soh_1059, \
                         soh_1062 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1411[k] = f_20 * snh_846[k]
                    + f_3 * pc_z[k] * soh_1056[k];

        t_1412[k] = f_12 * snh_1062[k]
                    + f_8 * sog0_762[k]
                    - f_9 * sog1_762[k]
                    + f_3 * pc_x[k] * soh_1062[k];

        t_1413[k] = f_14 * snh_870[k]
                    + f_3 * pc_y[k] * soh_1059[k];
    }

#pragma omp simd aligned(t_1414, t_1415, t_1416, t_1417, pc_x, snh_1064, snh_1065, snh_1066, \
                         snh_1067, sog0_764, sog1_764, soh_1064, soh_1065, soh_1066, \
                         soh_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1414[k] = f_12 * snh_1064[k]
                    + f_8 * sog0_764[k]
                    - f_9 * sog1_764[k]
                    + f_3 * pc_x[k] * soh_1064[k];

        t_1415[k] = f_12 * snh_1065[k]
                    + f_3 * pc_x[k] * soh_1065[k];

        t_1416[k] = f_12 * snh_1066[k]
                    + f_3 * pc_x[k] * soh_1066[k];

        t_1417[k] = f_12 * snh_1067[k]
                    + f_3 * pc_x[k] * soh_1067[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, t_1421, pc_x, pc_y, snh_876, snh_1068, \
                         snh_1069, snh_1070, sog0_760, sog1_760, soh_1065, soh_1068, soh_1069, \
                         soh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_12 * snh_1068[k]
                    + f_3 * pc_x[k] * soh_1068[k];

        t_1419[k] = f_12 * snh_1069[k]
                    + f_3 * pc_x[k] * soh_1069[k];

        t_1420[k] = f_12 * snh_1070[k]
                    + f_3 * pc_x[k] * soh_1070[k];

        t_1421[k] = f_14 * snh_876[k]
                    + f_1 * sog0_760[k]
                    - f_2 * sog1_760[k]
                    + f_3 * pc_y[k] * soh_1065[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, pc_y, pc_z, snh_855, snh_878, snh_879, \
                         sog0_762, sog0_763, sog1_762, sog1_763, soh_1065, soh_1067, \
                         soh_1068 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_20 * snh_855[k]
                    + f_3 * pc_z[k] * soh_1065[k];

        t_1423[k] = f_14 * snh_878[k]
                    + f_4 * sog0_762[k]
                    - f_5 * sog1_762[k]
                    + f_3 * pc_y[k] * soh_1067[k];

        t_1424[k] = f_14 * snh_879[k]
                    + f_6 * sog0_763[k]
                    - f_7 * sog1_763[k]
                    + f_3 * pc_y[k] * soh_1068[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, pc_y, pc_z, snh_860, snh_880, snh_881, \
                         sog0_764, sog1_764, soh_1069, soh_1070 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_14 * snh_880[k]
                    + f_8 * sog0_764[k]
                    - f_9 * sog1_764[k]
                    + f_3 * pc_y[k] * soh_1069[k];

        t_1426[k] = f_14 * snh_881[k]
                    + f_3 * pc_y[k] * soh_1070[k];

        t_1427[k] = f_20 * snh_860[k]
                    + f_1 * sog0_764[k]
                    - f_2 * sog1_764[k]
                    + f_3 * pc_z[k] * soh_1070[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, pc_x, pc_y, pc_z, snh_861, snh_882, snh_1071, \
                         sog0_765, sog1_765, soh_1071 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_12 * snh_1071[k]
                    + f_1 * sog0_765[k]
                    - f_2 * sog1_765[k]
                    + f_3 * pc_x[k] * soh_1071[k];

        t_1429[k] = f_13 * snh_882[k]
                    + f_3 * pc_y[k] * soh_1071[k];

        t_1430[k] = f_19 * snh_861[k]
                    + f_3 * pc_z[k] * soh_1071[k];
    }

#pragma omp simd aligned(t_1431, t_1432, t_1433, pc_x, pc_y, snh_884, snh_1074, snh_1076, \
                         sog0_768, sog0_770, sog1_768, sog1_770, soh_1073, soh_1074, \
                         soh_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1431[k] = f_12 * snh_1074[k]
                    + f_4 * sog0_768[k]
                    - f_5 * sog1_768[k]
                    + f_3 * pc_x[k] * soh_1074[k];

        t_1432[k] = f_13 * snh_884[k]
                    + f_3 * pc_y[k] * soh_1073[k];

        t_1433[k] = f_12 * snh_1076[k]
                    + f_4 * sog0_770[k]
                    - f_5 * sog1_770[k]
                    + f_3 * pc_x[k] * soh_1076[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, pc_x, pc_y, pc_z, snh_864, snh_887, snh_1077, \
                         sog0_771, sog1_771, soh_1074, soh_1076, \
                         soh_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = f_12 * snh_1077[k]
                    + f_6 * sog0_771[k]
                    - f_7 * sog1_771[k]
                    + f_3 * pc_x[k] * soh_1077[k];

        t_1435[k] = f_19 * snh_864[k]
                    + f_3 * pc_z[k] * soh_1074[k];

        t_1436[k] = f_13 * snh_887[k]
                    + f_3 * pc_y[k] * soh_1076[k];
    }

#pragma omp simd aligned(t_1437, t_1438, t_1439, pc_x, pc_z, snh_867, snh_1080, snh_1081, \
                         sog0_774, sog0_775, sog1_774, sog1_775, soh_1077, soh_1080, \
                         soh_1081 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1437[k] = f_12 * snh_1080[k]
                    + f_6 * sog0_774[k]
                    - f_7 * sog1_774[k]
                    + f_3 * pc_x[k] * soh_1080[k];

        t_1438[k] = f_12 * snh_1081[k]
                    + f_8 * sog0_775[k]
                    - f_9 * sog1_775[k]
                    + f_3 * pc_x[k] * soh_1081[k];

        t_1439[k] = f_19 * snh_867[k]
                    + f_3 * pc_z[k] * soh_1077[k];
    }

#pragma omp simd aligned(t_1440, t_1441, t_1442, pc_x, pc_y, snh_891, snh_1083, snh_1085, \
                         sog0_777, sog0_779, sog1_777, sog1_779, soh_1080, soh_1083, \
                         soh_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1440[k] = f_12 * snh_1083[k]
                    + f_8 * sog0_777[k]
                    - f_9 * sog1_777[k]
                    + f_3 * pc_x[k] * soh_1083[k];

        t_1441[k] = f_13 * snh_891[k]
                    + f_3 * pc_y[k] * soh_1080[k];

        t_1442[k] = f_12 * snh_1085[k]
                    + f_8 * sog0_779[k]
                    - f_9 * sog1_779[k]
                    + f_3 * pc_x[k] * soh_1085[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, t_1446, t_1447, pc_x, snh_1086, snh_1087, \
                         snh_1088, snh_1089, snh_1090, soh_1086, soh_1087, soh_1088, soh_1089, \
                         soh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = f_12 * snh_1086[k]
                    + f_3 * pc_x[k] * soh_1086[k];

        t_1444[k] = f_12 * snh_1087[k]
                    + f_3 * pc_x[k] * soh_1087[k];

        t_1445[k] = f_12 * snh_1088[k]
                    + f_3 * pc_x[k] * soh_1088[k];

        t_1446[k] = f_12 * snh_1089[k]
                    + f_3 * pc_x[k] * soh_1089[k];

        t_1447[k] = f_12 * snh_1090[k]
                    + f_3 * pc_x[k] * soh_1090[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pc_x, pc_y, pc_z, snh_876, snh_897, snh_1091, \
                         sog0_775, sog1_775, soh_1086, soh_1091 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_12 * snh_1091[k]
                    + f_3 * pc_x[k] * soh_1091[k];

        t_1449[k] = f_13 * snh_897[k]
                    + f_1 * sog0_775[k]
                    - f_2 * sog1_775[k]
                    + f_3 * pc_y[k] * soh_1086[k];

        t_1450[k] = f_19 * snh_876[k]
                    + f_3 * pc_z[k] * soh_1086[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pc_y, snh_899, snh_900, snh_901, sog0_777, \
                         sog0_778, sog0_779, sog1_777, sog1_778, sog1_779, soh_1088, soh_1089, \
                         soh_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_13 * snh_899[k]
                    + f_4 * sog0_777[k]
                    - f_5 * sog1_777[k]
                    + f_3 * pc_y[k] * soh_1088[k];

        t_1452[k] = f_13 * snh_900[k]
                    + f_6 * sog0_778[k]
                    - f_7 * sog1_778[k]
                    + f_3 * pc_y[k] * soh_1089[k];

        t_1453[k] = f_13 * snh_901[k]
                    + f_8 * sog0_779[k]
                    - f_9 * sog1_779[k]
                    + f_3 * pc_y[k] * soh_1090[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pc_x, pc_y, pc_z, snh_881, snh_902, snh_1092, \
                         sog0_779, sog0_780, sog1_779, sog1_780, soh_1091, \
                         soh_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_13 * snh_902[k]
                    + f_3 * pc_y[k] * soh_1091[k];

        t_1455[k] = f_19 * snh_881[k]
                    + f_1 * sog0_779[k]
                    - f_2 * sog1_779[k]
                    + f_3 * pc_z[k] * soh_1091[k];

        t_1456[k] = f_12 * snh_1092[k]
                    + f_1 * sog0_780[k]
                    - f_2 * sog1_780[k]
                    + f_3 * pc_x[k] * soh_1092[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t sog0, const size_t sog1,
                                                           const size_t soh, const size_t ncols,
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
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;

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
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sni0_1232 = buffer.data(sni0 + 1232);
    const auto *sni0_1235 = buffer.data(sni0 + 1235);
    const auto *sni0_1237 = buffer.data(sni0 + 1237);
    const auto *sni0_1238 = buffer.data(sni0 + 1238);
    const auto *sni0_1241 = buffer.data(sni0 + 1241);
    const auto *sni0_1242 = buffer.data(sni0 + 1242);
    const auto *sni0_1244 = buffer.data(sni0 + 1244);
    const auto *sni0_1246 = buffer.data(sni0 + 1246);
    const auto *sni0_1259 = buffer.data(sni0 + 1259);
    const auto *sni0_1260 = buffer.data(sni0 + 1260);
    const auto *sni0_1263 = buffer.data(sni0 + 1263);
    const auto *sni0_1266 = buffer.data(sni0 + 1266);
    const auto *sni0_1540 = buffer.data(sni0 + 1540);
    const auto *sni0_1543 = buffer.data(sni0 + 1543);
    const auto *sni0_1545 = buffer.data(sni0 + 1545);
    const auto *sni0_1546 = buffer.data(sni0 + 1546);
    const auto *sni0_1549 = buffer.data(sni0 + 1549);
    const auto *sni0_1550 = buffer.data(sni0 + 1550);
    const auto *sni0_1552 = buffer.data(sni0 + 1552);
    const auto *sni0_1554 = buffer.data(sni0 + 1554);
    const auto *sni0_1561 = buffer.data(sni0 + 1561);
    const auto *sni0_1563 = buffer.data(sni0 + 1563);
    const auto *sni0_1564 = buffer.data(sni0 + 1564);
    const auto *sni0_1565 = buffer.data(sni0 + 1565);
    const auto *sni0_1567 = buffer.data(sni0 + 1567);
    const auto *sni0_1573 = buffer.data(sni0 + 1573);

    const auto *snh_882 = buffer.data(snh + 882);
    const auto *snh_885 = buffer.data(snh + 885);
    const auto *snh_888 = buffer.data(snh + 888);
    const auto *snh_897 = buffer.data(snh + 897);
    const auto *snh_902 = buffer.data(snh + 902);
    const auto *snh_903 = buffer.data(snh + 903);
    const auto *snh_905 = buffer.data(snh + 905);
    const auto *snh_906 = buffer.data(snh + 906);
    const auto *snh_908 = buffer.data(snh + 908);
    const auto *snh_909 = buffer.data(snh + 909);
    const auto *snh_912 = buffer.data(snh + 912);
    const auto *snh_918 = buffer.data(snh + 918);
    const auto *snh_920 = buffer.data(snh + 920);
    const auto *snh_921 = buffer.data(snh + 921);
    const auto *snh_922 = buffer.data(snh + 922);
    const auto *snh_923 = buffer.data(snh + 923);
    const auto *snh_924 = buffer.data(snh + 924);
    const auto *snh_925 = buffer.data(snh + 925);
    const auto *snh_926 = buffer.data(snh + 926);
    const auto *snh_927 = buffer.data(snh + 927);
    const auto *snh_929 = buffer.data(snh + 929);
    const auto *snh_930 = buffer.data(snh + 930);
    const auto *snh_932 = buffer.data(snh + 932);
    const auto *snh_933 = buffer.data(snh + 933);
    const auto *snh_939 = buffer.data(snh + 939);
    const auto *snh_941 = buffer.data(snh + 941);
    const auto *snh_942 = buffer.data(snh + 942);
    const auto *snh_943 = buffer.data(snh + 943);
    const auto *snh_944 = buffer.data(snh + 944);
    const auto *snh_945 = buffer.data(snh + 945);
    const auto *snh_947 = buffer.data(snh + 947);
    const auto *snh_950 = buffer.data(snh + 950);
    const auto *snh_954 = buffer.data(snh + 954);
    const auto *snh_965 = buffer.data(snh + 965);
    const auto *snh_966 = buffer.data(snh + 966);
    const auto *snh_968 = buffer.data(snh + 968);
    const auto *snh_1095 = buffer.data(snh + 1095);
    const auto *snh_1097 = buffer.data(snh + 1097);
    const auto *snh_1098 = buffer.data(snh + 1098);
    const auto *snh_1101 = buffer.data(snh + 1101);
    const auto *snh_1102 = buffer.data(snh + 1102);
    const auto *snh_1104 = buffer.data(snh + 1104);
    const auto *snh_1106 = buffer.data(snh + 1106);
    const auto *snh_1107 = buffer.data(snh + 1107);
    const auto *snh_1108 = buffer.data(snh + 1108);
    const auto *snh_1109 = buffer.data(snh + 1109);
    const auto *snh_1110 = buffer.data(snh + 1110);
    const auto *snh_1111 = buffer.data(snh + 1111);
    const auto *snh_1112 = buffer.data(snh + 1112);
    const auto *snh_1128 = buffer.data(snh + 1128);
    const auto *snh_1129 = buffer.data(snh + 1129);
    const auto *snh_1130 = buffer.data(snh + 1130);
    const auto *snh_1131 = buffer.data(snh + 1131);
    const auto *snh_1132 = buffer.data(snh + 1132);
    const auto *snh_1133 = buffer.data(snh + 1133);
    const auto *snh_1134 = buffer.data(snh + 1134);
    const auto *snh_1137 = buffer.data(snh + 1137);
    const auto *snh_1139 = buffer.data(snh + 1139);
    const auto *snh_1140 = buffer.data(snh + 1140);
    const auto *snh_1143 = buffer.data(snh + 1143);
    const auto *snh_1144 = buffer.data(snh + 1144);
    const auto *snh_1146 = buffer.data(snh + 1146);
    const auto *snh_1148 = buffer.data(snh + 1148);
    const auto *snh_1149 = buffer.data(snh + 1149);
    const auto *snh_1150 = buffer.data(snh + 1150);
    const auto *snh_1151 = buffer.data(snh + 1151);
    const auto *snh_1152 = buffer.data(snh + 1152);
    const auto *snh_1153 = buffer.data(snh + 1153);
    const auto *snh_1154 = buffer.data(snh + 1154);
    const auto *snh_1155 = buffer.data(snh + 1155);
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
    const auto *snh_1181 = buffer.data(snh + 1181);

    const auto *sni1_1232 = buffer.data(sni1 + 1232);
    const auto *sni1_1235 = buffer.data(sni1 + 1235);
    const auto *sni1_1237 = buffer.data(sni1 + 1237);
    const auto *sni1_1238 = buffer.data(sni1 + 1238);
    const auto *sni1_1241 = buffer.data(sni1 + 1241);
    const auto *sni1_1242 = buffer.data(sni1 + 1242);
    const auto *sni1_1244 = buffer.data(sni1 + 1244);
    const auto *sni1_1246 = buffer.data(sni1 + 1246);
    const auto *sni1_1259 = buffer.data(sni1 + 1259);
    const auto *sni1_1260 = buffer.data(sni1 + 1260);
    const auto *sni1_1263 = buffer.data(sni1 + 1263);
    const auto *sni1_1266 = buffer.data(sni1 + 1266);
    const auto *sni1_1540 = buffer.data(sni1 + 1540);
    const auto *sni1_1543 = buffer.data(sni1 + 1543);
    const auto *sni1_1545 = buffer.data(sni1 + 1545);
    const auto *sni1_1546 = buffer.data(sni1 + 1546);
    const auto *sni1_1549 = buffer.data(sni1 + 1549);
    const auto *sni1_1550 = buffer.data(sni1 + 1550);
    const auto *sni1_1552 = buffer.data(sni1 + 1552);
    const auto *sni1_1554 = buffer.data(sni1 + 1554);
    const auto *sni1_1561 = buffer.data(sni1 + 1561);
    const auto *sni1_1563 = buffer.data(sni1 + 1563);
    const auto *sni1_1564 = buffer.data(sni1 + 1564);
    const auto *sni1_1565 = buffer.data(sni1 + 1565);
    const auto *sni1_1567 = buffer.data(sni1 + 1567);
    const auto *sni1_1573 = buffer.data(sni1 + 1573);

    const auto *sog0_783 = buffer.data(sog0 + 783);
    const auto *sog0_785 = buffer.data(sog0 + 785);
    const auto *sog0_786 = buffer.data(sog0 + 786);
    const auto *sog0_789 = buffer.data(sog0 + 789);
    const auto *sog0_790 = buffer.data(sog0 + 790);
    const auto *sog0_792 = buffer.data(sog0 + 792);
    const auto *sog0_793 = buffer.data(sog0 + 793);
    const auto *sog0_794 = buffer.data(sog0 + 794);
    const auto *sog0_805 = buffer.data(sog0 + 805);
    const auto *sog0_807 = buffer.data(sog0 + 807);
    const auto *sog0_808 = buffer.data(sog0 + 808);
    const auto *sog0_809 = buffer.data(sog0 + 809);
    const auto *sog0_810 = buffer.data(sog0 + 810);
    const auto *sog0_813 = buffer.data(sog0 + 813);
    const auto *sog0_815 = buffer.data(sog0 + 815);
    const auto *sog0_816 = buffer.data(sog0 + 816);
    const auto *sog0_819 = buffer.data(sog0 + 819);
    const auto *sog0_820 = buffer.data(sog0 + 820);
    const auto *sog0_822 = buffer.data(sog0 + 822);
    const auto *sog0_823 = buffer.data(sog0 + 823);
    const auto *sog0_824 = buffer.data(sog0 + 824);

    const auto *sog1_783 = buffer.data(sog1 + 783);
    const auto *sog1_785 = buffer.data(sog1 + 785);
    const auto *sog1_786 = buffer.data(sog1 + 786);
    const auto *sog1_789 = buffer.data(sog1 + 789);
    const auto *sog1_790 = buffer.data(sog1 + 790);
    const auto *sog1_792 = buffer.data(sog1 + 792);
    const auto *sog1_793 = buffer.data(sog1 + 793);
    const auto *sog1_794 = buffer.data(sog1 + 794);
    const auto *sog1_805 = buffer.data(sog1 + 805);
    const auto *sog1_807 = buffer.data(sog1 + 807);
    const auto *sog1_808 = buffer.data(sog1 + 808);
    const auto *sog1_809 = buffer.data(sog1 + 809);
    const auto *sog1_810 = buffer.data(sog1 + 810);
    const auto *sog1_813 = buffer.data(sog1 + 813);
    const auto *sog1_815 = buffer.data(sog1 + 815);
    const auto *sog1_816 = buffer.data(sog1 + 816);
    const auto *sog1_819 = buffer.data(sog1 + 819);
    const auto *sog1_820 = buffer.data(sog1 + 820);
    const auto *sog1_822 = buffer.data(sog1 + 822);
    const auto *sog1_823 = buffer.data(sog1 + 823);
    const auto *sog1_824 = buffer.data(sog1 + 824);

    const auto *soh_1092 = buffer.data(soh + 1092);
    const auto *soh_1094 = buffer.data(soh + 1094);
    const auto *soh_1095 = buffer.data(soh + 1095);
    const auto *soh_1097 = buffer.data(soh + 1097);
    const auto *soh_1098 = buffer.data(soh + 1098);
    const auto *soh_1101 = buffer.data(soh + 1101);
    const auto *soh_1102 = buffer.data(soh + 1102);
    const auto *soh_1104 = buffer.data(soh + 1104);
    const auto *soh_1106 = buffer.data(soh + 1106);
    const auto *soh_1107 = buffer.data(soh + 1107);
    const auto *soh_1108 = buffer.data(soh + 1108);
    const auto *soh_1109 = buffer.data(soh + 1109);
    const auto *soh_1110 = buffer.data(soh + 1110);
    const auto *soh_1111 = buffer.data(soh + 1111);
    const auto *soh_1112 = buffer.data(soh + 1112);
    const auto *soh_1113 = buffer.data(soh + 1113);
    const auto *soh_1115 = buffer.data(soh + 1115);
    const auto *soh_1116 = buffer.data(soh + 1116);
    const auto *soh_1118 = buffer.data(soh + 1118);
    const auto *soh_1119 = buffer.data(soh + 1119);
    const auto *soh_1122 = buffer.data(soh + 1122);
    const auto *soh_1128 = buffer.data(soh + 1128);
    const auto *soh_1129 = buffer.data(soh + 1129);
    const auto *soh_1130 = buffer.data(soh + 1130);
    const auto *soh_1131 = buffer.data(soh + 1131);
    const auto *soh_1132 = buffer.data(soh + 1132);
    const auto *soh_1133 = buffer.data(soh + 1133);
    const auto *soh_1134 = buffer.data(soh + 1134);
    const auto *soh_1136 = buffer.data(soh + 1136);
    const auto *soh_1137 = buffer.data(soh + 1137);
    const auto *soh_1139 = buffer.data(soh + 1139);
    const auto *soh_1140 = buffer.data(soh + 1140);
    const auto *soh_1143 = buffer.data(soh + 1143);
    const auto *soh_1144 = buffer.data(soh + 1144);
    const auto *soh_1146 = buffer.data(soh + 1146);
    const auto *soh_1148 = buffer.data(soh + 1148);
    const auto *soh_1149 = buffer.data(soh + 1149);
    const auto *soh_1150 = buffer.data(soh + 1150);
    const auto *soh_1151 = buffer.data(soh + 1151);
    const auto *soh_1152 = buffer.data(soh + 1152);
    const auto *soh_1153 = buffer.data(soh + 1153);
    const auto *soh_1154 = buffer.data(soh + 1154);
    const auto *soh_1155 = buffer.data(soh + 1155);
    const auto *soh_1157 = buffer.data(soh + 1157);
    const auto *soh_1158 = buffer.data(soh + 1158);
    const auto *soh_1160 = buffer.data(soh + 1160);
    const auto *soh_1161 = buffer.data(soh + 1161);
    const auto *soh_1164 = buffer.data(soh + 1164);
    const auto *soh_1170 = buffer.data(soh + 1170);
    const auto *soh_1171 = buffer.data(soh + 1171);
    const auto *soh_1172 = buffer.data(soh + 1172);
    const auto *soh_1173 = buffer.data(soh + 1173);
    const auto *soh_1174 = buffer.data(soh + 1174);
    const auto *soh_1175 = buffer.data(soh + 1175);
    const auto *soh_1176 = buffer.data(soh + 1176);
    const auto *soh_1178 = buffer.data(soh + 1178);

#pragma omp simd aligned(t_1457, t_1458, t_1459, t_1460, pc_x, pc_y, pc_z, snh_882, snh_903, \
                         snh_905, snh_1095, sog0_783, sog1_783, soh_1092, soh_1094, \
                         soh_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_12 * snh_903[k]
                    + f_3 * pc_y[k] * soh_1092[k];

        t_1458[k] = f_18 * snh_882[k]
                    + f_3 * pc_z[k] * soh_1092[k];

        t_1459[k] = f_12 * snh_1095[k]
                    + f_4 * sog0_783[k]
                    - f_5 * sog1_783[k]
                    + f_3 * pc_x[k] * soh_1095[k];

        t_1460[k] = f_12 * snh_905[k]
                    + f_3 * pc_y[k] * soh_1094[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pc_x, pc_z, snh_885, snh_1097, snh_1098, \
                         sog0_785, sog0_786, sog1_785, sog1_786, soh_1095, soh_1097, \
                         soh_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = f_12 * snh_1097[k]
                    + f_4 * sog0_785[k]
                    - f_5 * sog1_785[k]
                    + f_3 * pc_x[k] * soh_1097[k];

        t_1462[k] = f_12 * snh_1098[k]
                    + f_6 * sog0_786[k]
                    - f_7 * sog1_786[k]
                    + f_3 * pc_x[k] * soh_1098[k];

        t_1463[k] = f_18 * snh_885[k]
                    + f_3 * pc_z[k] * soh_1095[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pc_x, pc_y, snh_908, snh_1101, snh_1102, \
                         sog0_789, sog0_790, sog1_789, sog1_790, soh_1097, soh_1101, \
                         soh_1102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_12 * snh_908[k]
                    + f_3 * pc_y[k] * soh_1097[k];

        t_1465[k] = f_12 * snh_1101[k]
                    + f_6 * sog0_789[k]
                    - f_7 * sog1_789[k]
                    + f_3 * pc_x[k] * soh_1101[k];

        t_1466[k] = f_12 * snh_1102[k]
                    + f_8 * sog0_790[k]
                    - f_9 * sog1_790[k]
                    + f_3 * pc_x[k] * soh_1102[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, pc_x, pc_y, pc_z, snh_888, snh_912, snh_1104, \
                         sog0_792, sog1_792, soh_1098, soh_1101, \
                         soh_1104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_18 * snh_888[k]
                    + f_3 * pc_z[k] * soh_1098[k];

        t_1468[k] = f_12 * snh_1104[k]
                    + f_8 * sog0_792[k]
                    - f_9 * sog1_792[k]
                    + f_3 * pc_x[k] * soh_1104[k];

        t_1469[k] = f_12 * snh_912[k]
                    + f_3 * pc_y[k] * soh_1101[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, pc_x, snh_1106, snh_1107, snh_1108, \
                         snh_1109, sog0_794, sog1_794, soh_1106, soh_1107, soh_1108, \
                         soh_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_12 * snh_1106[k]
                    + f_8 * sog0_794[k]
                    - f_9 * sog1_794[k]
                    + f_3 * pc_x[k] * soh_1106[k];

        t_1471[k] = f_12 * snh_1107[k]
                    + f_3 * pc_x[k] * soh_1107[k];

        t_1472[k] = f_12 * snh_1108[k]
                    + f_3 * pc_x[k] * soh_1108[k];

        t_1473[k] = f_12 * snh_1109[k]
                    + f_3 * pc_x[k] * soh_1109[k];
    }

#pragma omp simd aligned(t_1474, t_1475, t_1476, t_1477, pc_x, pc_y, snh_918, snh_1110, \
                         snh_1111, snh_1112, sog0_790, sog1_790, soh_1107, soh_1110, soh_1111, \
                         soh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1474[k] = f_12 * snh_1110[k]
                    + f_3 * pc_x[k] * soh_1110[k];

        t_1475[k] = f_12 * snh_1111[k]
                    + f_3 * pc_x[k] * soh_1111[k];

        t_1476[k] = f_12 * snh_1112[k]
                    + f_3 * pc_x[k] * soh_1112[k];

        t_1477[k] = f_12 * snh_918[k]
                    + f_1 * sog0_790[k]
                    - f_2 * sog1_790[k]
                    + f_3 * pc_y[k] * soh_1107[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_y, pc_z, snh_897, snh_920, snh_921, \
                         sog0_792, sog0_793, sog1_792, sog1_793, soh_1107, soh_1109, \
                         soh_1110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_18 * snh_897[k]
                    + f_3 * pc_z[k] * soh_1107[k];

        t_1479[k] = f_12 * snh_920[k]
                    + f_4 * sog0_792[k]
                    - f_5 * sog1_792[k]
                    + f_3 * pc_y[k] * soh_1109[k];

        t_1480[k] = f_12 * snh_921[k]
                    + f_6 * sog0_793[k]
                    - f_7 * sog1_793[k]
                    + f_3 * pc_y[k] * soh_1110[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, t_1484, pb_y, pc_y, pc_z, sni0_1232, snh_902, \
                         snh_922, snh_923, sni1_1232, sog0_794, sog1_794, soh_1111, \
                         soh_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_12 * snh_922[k]
                    + f_8 * sog0_794[k]
                    - f_9 * sog1_794[k]
                    + f_3 * pc_y[k] * soh_1111[k];

        t_1482[k] = f_12 * snh_923[k]
                    + f_3 * pc_y[k] * soh_1112[k];

        t_1483[k] = f_18 * snh_902[k]
                    + f_1 * sog0_794[k]
                    - f_2 * sog1_794[k]
                    + f_3 * pc_z[k] * soh_1112[k];

        t_1484[k] = pb_y[k] * sni0_1232[k]
                    - f_10 * pc_y[k] * sni1_1232[k];
    }

#pragma omp simd aligned(t_1485, t_1486, t_1487, t_1488, pb_y, pc_y, pc_z, sni0_1235, snh_903, \
                         snh_924, snh_925, snh_926, sni1_1235, soh_1113, \
                         soh_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1485[k] = f_11 * snh_924[k]
                    + f_3 * pc_y[k] * soh_1113[k];

        t_1486[k] = f_17 * snh_903[k]
                    + f_3 * pc_z[k] * soh_1113[k];

        t_1487[k] = pb_y[k] * sni0_1235[k]
                    + f_12 * snh_925[k]
                    - f_10 * pc_y[k] * sni1_1235[k];

        t_1488[k] = f_11 * snh_926[k]
                    + f_3 * pc_y[k] * soh_1115[k];
    }

#pragma omp simd aligned(t_1489, t_1490, t_1491, t_1492, pb_y, pc_y, pc_z, sni0_1237, \
                         sni0_1238, snh_906, snh_927, snh_929, sni1_1237, sni1_1238, soh_1116, \
                         soh_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1489[k] = pb_y[k] * sni0_1237[k]
                    - f_10 * pc_y[k] * sni1_1237[k];

        t_1490[k] = pb_y[k] * sni0_1238[k]
                    + f_13 * snh_927[k]
                    - f_10 * pc_y[k] * sni1_1238[k];

        t_1491[k] = f_17 * snh_906[k]
                    + f_3 * pc_z[k] * soh_1116[k];

        t_1492[k] = f_11 * snh_929[k]
                    + f_3 * pc_y[k] * soh_1118[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pb_y, pc_y, pc_z, sni0_1241, sni0_1242, \
                         snh_909, snh_930, sni1_1241, sni1_1242, \
                         soh_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = pb_y[k] * sni0_1241[k]
                    - f_10 * pc_y[k] * sni1_1241[k];

        t_1494[k] = pb_y[k] * sni0_1242[k]
                    + f_14 * snh_930[k]
                    - f_10 * pc_y[k] * sni1_1242[k];

        t_1495[k] = f_17 * snh_909[k]
                    + f_3 * pc_z[k] * soh_1119[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, t_1499, pb_y, pc_x, pc_y, sni0_1244, \
                         sni0_1246, snh_932, snh_933, snh_1128, sni1_1244, sni1_1246, \
                         soh_1122, soh_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = pb_y[k] * sni0_1244[k]
                    + f_12 * snh_932[k]
                    - f_10 * pc_y[k] * sni1_1244[k];

        t_1497[k] = f_11 * snh_933[k]
                    + f_3 * pc_y[k] * soh_1122[k];

        t_1498[k] = pb_y[k] * sni0_1246[k]
                    - f_10 * pc_y[k] * sni1_1246[k];

        t_1499[k] = f_12 * snh_1128[k]
                    + f_3 * pc_x[k] * soh_1128[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, t_1503, t_1504, pc_x, snh_1129, snh_1130, \
                         snh_1131, snh_1132, snh_1133, soh_1129, soh_1130, soh_1131, soh_1132, \
                         soh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_12 * snh_1129[k]
                    + f_3 * pc_x[k] * soh_1129[k];

        t_1501[k] = f_12 * snh_1130[k]
                    + f_3 * pc_x[k] * soh_1130[k];

        t_1502[k] = f_12 * snh_1131[k]
                    + f_3 * pc_x[k] * soh_1131[k];

        t_1503[k] = f_12 * snh_1132[k]
                    + f_3 * pc_x[k] * soh_1132[k];

        t_1504[k] = f_12 * snh_1133[k]
                    + f_3 * pc_x[k] * soh_1133[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pc_y, pc_z, snh_918, snh_939, snh_941, \
                         sog0_805, sog0_807, sog1_805, sog1_807, soh_1128, \
                         soh_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_11 * snh_939[k]
                    + f_1 * sog0_805[k]
                    - f_2 * sog1_805[k]
                    + f_3 * pc_y[k] * soh_1128[k];

        t_1506[k] = f_17 * snh_918[k]
                    + f_3 * pc_z[k] * soh_1128[k];

        t_1507[k] = f_11 * snh_941[k]
                    + f_4 * sog0_807[k]
                    - f_5 * sog1_807[k]
                    + f_3 * pc_y[k] * soh_1130[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pc_y, snh_942, snh_943, snh_944, sog0_808, \
                         sog0_809, sog1_808, sog1_809, soh_1131, soh_1132, \
                         soh_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_11 * snh_942[k]
                    + f_6 * sog0_808[k]
                    - f_7 * sog1_808[k]
                    + f_3 * pc_y[k] * soh_1131[k];

        t_1509[k] = f_11 * snh_943[k]
                    + f_8 * sog0_809[k]
                    - f_9 * sog1_809[k]
                    + f_3 * pc_y[k] * soh_1132[k];

        t_1510[k] = f_11 * snh_944[k]
                    + f_3 * pc_y[k] * soh_1133[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pb_y, pc_x, pc_y, pc_z, sni0_1259, \
                         snh_924, snh_1134, sni1_1259, sog0_810, sog1_810, \
                         soh_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = pb_y[k] * sni0_1259[k]
                    - f_10 * pc_y[k] * sni1_1259[k];

        t_1512[k] = f_12 * snh_1134[k]
                    + f_1 * sog0_810[k]
                    - f_2 * sog1_810[k]
                    + f_3 * pc_x[k] * soh_1134[k];

        t_1513[k] = f_3 * pc_y[k] * soh_1134[k];

        t_1514[k] = f_16 * snh_924[k]
                    + f_3 * pc_z[k] * soh_1134[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, pc_x, pc_y, snh_1137, snh_1139, sog0_813, \
                         sog0_815, sog1_813, sog1_815, soh_1136, soh_1137, \
                         soh_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_12 * snh_1137[k]
                    + f_4 * sog0_813[k]
                    - f_5 * sog1_813[k]
                    + f_3 * pc_x[k] * soh_1137[k];

        t_1516[k] = f_3 * pc_y[k] * soh_1136[k];

        t_1517[k] = f_12 * snh_1139[k]
                    + f_4 * sog0_815[k]
                    - f_5 * sog1_815[k]
                    + f_3 * pc_x[k] * soh_1139[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, pc_x, pc_y, pc_z, snh_927, snh_1140, \
                         sog0_816, sog1_816, soh_1137, soh_1139, \
                         soh_1140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_12 * snh_1140[k]
                    + f_6 * sog0_816[k]
                    - f_7 * sog1_816[k]
                    + f_3 * pc_x[k] * soh_1140[k];

        t_1519[k] = f_16 * snh_927[k]
                    + f_3 * pc_z[k] * soh_1137[k];

        t_1520[k] = f_3 * pc_y[k] * soh_1139[k];
    }

#pragma omp simd aligned(t_1521, t_1522, t_1523, pc_x, pc_z, snh_930, snh_1143, snh_1144, \
                         sog0_819, sog0_820, sog1_819, sog1_820, soh_1140, soh_1143, \
                         soh_1144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1521[k] = f_12 * snh_1143[k]
                    + f_6 * sog0_819[k]
                    - f_7 * sog1_819[k]
                    + f_3 * pc_x[k] * soh_1143[k];

        t_1522[k] = f_12 * snh_1144[k]
                    + f_8 * sog0_820[k]
                    - f_9 * sog1_820[k]
                    + f_3 * pc_x[k] * soh_1144[k];

        t_1523[k] = f_16 * snh_930[k]
                    + f_3 * pc_z[k] * soh_1140[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, pc_x, pc_y, snh_1146, snh_1148, sog0_822, \
                         sog0_824, sog1_822, sog1_824, soh_1143, soh_1146, \
                         soh_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = f_12 * snh_1146[k]
                    + f_8 * sog0_822[k]
                    - f_9 * sog1_822[k]
                    + f_3 * pc_x[k] * soh_1146[k];

        t_1525[k] = f_3 * pc_y[k] * soh_1143[k];

        t_1526[k] = f_12 * snh_1148[k]
                    + f_8 * sog0_824[k]
                    - f_9 * sog1_824[k]
                    + f_3 * pc_x[k] * soh_1148[k];
    }

#pragma omp simd aligned(t_1527, t_1528, t_1529, t_1530, t_1531, pc_x, snh_1149, snh_1150, \
                         snh_1151, snh_1152, snh_1153, soh_1149, soh_1150, soh_1151, soh_1152, \
                         soh_1153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1527[k] = f_12 * snh_1149[k]
                    + f_3 * pc_x[k] * soh_1149[k];

        t_1528[k] = f_12 * snh_1150[k]
                    + f_3 * pc_x[k] * soh_1150[k];

        t_1529[k] = f_12 * snh_1151[k]
                    + f_3 * pc_x[k] * soh_1151[k];

        t_1530[k] = f_12 * snh_1152[k]
                    + f_3 * pc_x[k] * soh_1152[k];

        t_1531[k] = f_12 * snh_1153[k]
                    + f_3 * pc_x[k] * soh_1153[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, t_1535, pc_x, pc_y, pc_z, snh_939, snh_1154, \
                         sog0_820, sog0_822, sog1_820, sog1_822, soh_1149, soh_1151, \
                         soh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = f_12 * snh_1154[k]
                    + f_3 * pc_x[k] * soh_1154[k];

        t_1533[k] = f_1 * sog0_820[k]
                    - f_2 * sog1_820[k]
                    + f_3 * pc_y[k] * soh_1149[k];

        t_1534[k] = f_16 * snh_939[k]
                    + f_3 * pc_z[k] * soh_1149[k];

        t_1535[k] = f_4 * sog0_822[k]
                    - f_5 * sog1_822[k]
                    + f_3 * pc_y[k] * soh_1151[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, t_1539, pc_y, pc_z, snh_944, sog0_823, \
                         sog0_824, sog1_823, sog1_824, soh_1152, soh_1153, \
                         soh_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_6 * sog0_823[k]
                    - f_7 * sog1_823[k]
                    + f_3 * pc_y[k] * soh_1152[k];

        t_1537[k] = f_8 * sog0_824[k]
                    - f_9 * sog1_824[k]
                    + f_3 * pc_y[k] * soh_1153[k];

        t_1538[k] = f_3 * pc_y[k] * soh_1154[k];

        t_1539[k] = f_16 * snh_944[k]
                    + f_1 * sog0_824[k]
                    - f_2 * sog1_824[k]
                    + f_3 * pc_z[k] * soh_1154[k];
    }

#pragma omp simd aligned(t_1540, t_1541, t_1542, t_1543, pb_x, pc_x, pc_y, pc_z, sni0_1540, \
                         sni0_1543, snh_945, snh_1155, snh_1158, sni1_1540, sni1_1543, \
                         soh_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1540[k] = pb_x[k] * sni0_1540[k]
                    + f_19 * snh_1155[k]
                    - f_10 * pc_x[k] * sni1_1540[k];

        t_1541[k] = f_15 * snh_945[k]
                    + f_3 * pc_y[k] * soh_1155[k];

        t_1542[k] = f_3 * pc_z[k] * soh_1155[k];

        t_1543[k] = pb_x[k] * sni0_1543[k]
                    + f_14 * snh_1158[k]
                    - f_10 * pc_x[k] * sni1_1543[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, pb_x, pc_x, pc_y, sni0_1545, sni0_1546, \
                         snh_947, snh_1160, snh_1161, sni1_1545, sni1_1546, \
                         soh_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = f_15 * snh_947[k]
                    + f_3 * pc_y[k] * soh_1157[k];

        t_1545[k] = pb_x[k] * sni0_1545[k]
                    + f_14 * snh_1160[k]
                    - f_10 * pc_x[k] * sni1_1545[k];

        t_1546[k] = pb_x[k] * sni0_1546[k]
                    + f_13 * snh_1161[k]
                    - f_10 * pc_x[k] * sni1_1546[k];
    }

#pragma omp simd aligned(t_1547, t_1548, t_1549, pb_x, pc_x, pc_y, pc_z, sni0_1549, snh_950, \
                         snh_1164, sni1_1549, soh_1158, soh_1160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1547[k] = f_3 * pc_z[k] * soh_1158[k];

        t_1548[k] = f_15 * snh_950[k]
                    + f_3 * pc_y[k] * soh_1160[k];

        t_1549[k] = pb_x[k] * sni0_1549[k]
                    + f_13 * snh_1164[k]
                    - f_10 * pc_x[k] * sni1_1549[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pb_x, pc_x, pc_z, sni0_1550, sni0_1552, \
                         snh_1165, snh_1167, sni1_1550, sni1_1552, \
                         soh_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = pb_x[k] * sni0_1550[k]
                    + f_12 * snh_1165[k]
                    - f_10 * pc_x[k] * sni1_1550[k];

        t_1551[k] = f_3 * pc_z[k] * soh_1161[k];

        t_1552[k] = pb_x[k] * sni0_1552[k]
                    + f_12 * snh_1167[k]
                    - f_10 * pc_x[k] * sni1_1552[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, t_1556, pb_x, pc_x, pc_y, sni0_1554, snh_954, \
                         snh_1169, snh_1170, snh_1171, sni1_1554, soh_1164, soh_1170, \
                         soh_1171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = f_15 * snh_954[k]
                    + f_3 * pc_y[k] * soh_1164[k];

        t_1554[k] = pb_x[k] * sni0_1554[k]
                    + f_12 * snh_1169[k]
                    - f_10 * pc_x[k] * sni1_1554[k];

        t_1555[k] = f_11 * snh_1170[k]
                    + f_3 * pc_x[k] * soh_1170[k];

        t_1556[k] = f_11 * snh_1171[k]
                    + f_3 * pc_x[k] * soh_1171[k];
    }

#pragma omp simd aligned(t_1557, t_1558, t_1559, t_1560, pc_x, snh_1172, snh_1173, snh_1174, \
                         snh_1175, soh_1172, soh_1173, soh_1174, \
                         soh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1557[k] = f_11 * snh_1172[k]
                    + f_3 * pc_x[k] * soh_1172[k];

        t_1558[k] = f_11 * snh_1173[k]
                    + f_3 * pc_x[k] * soh_1173[k];

        t_1559[k] = f_11 * snh_1174[k]
                    + f_3 * pc_x[k] * soh_1174[k];

        t_1560[k] = f_11 * snh_1175[k]
                    + f_3 * pc_x[k] * soh_1175[k];
    }

#pragma omp simd aligned(t_1561, t_1562, t_1563, t_1564, pb_x, pc_x, pc_z, sni0_1561, \
                         sni0_1563, sni0_1564, sni1_1561, sni1_1563, sni1_1564, \
                         soh_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1561[k] = pb_x[k] * sni0_1561[k]
                    - f_10 * pc_x[k] * sni1_1561[k];

        t_1562[k] = f_3 * pc_z[k] * soh_1170[k];

        t_1563[k] = pb_x[k] * sni0_1563[k]
                    - f_10 * pc_x[k] * sni1_1563[k];

        t_1564[k] = pb_x[k] * sni0_1564[k]
                    - f_10 * pc_x[k] * sni1_1564[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, pb_x, pc_x, pc_y, sni0_1565, sni0_1567, \
                         snh_965, sni1_1565, sni1_1567, soh_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = pb_x[k] * sni0_1565[k]
                    - f_10 * pc_x[k] * sni1_1565[k];

        t_1566[k] = f_15 * snh_965[k]
                    + f_3 * pc_y[k] * soh_1175[k];

        t_1567[k] = pb_x[k] * sni0_1567[k]
                    - f_10 * pc_x[k] * sni1_1567[k];
    }

#pragma omp simd aligned(t_1568, t_1569, t_1570, t_1571, pb_z, pc_y, pc_z, sni0_1260, \
                         sni0_1263, snh_945, snh_966, sni1_1260, sni1_1263, \
                         soh_1176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1568[k] = pb_z[k] * sni0_1260[k]
                    - f_10 * pc_z[k] * sni1_1260[k];

        t_1569[k] = f_16 * snh_966[k]
                    + f_3 * pc_y[k] * soh_1176[k];

        t_1570[k] = f_11 * snh_945[k]
                    + f_3 * pc_z[k] * soh_1176[k];

        t_1571[k] = pb_z[k] * sni0_1263[k]
                    - f_10 * pc_z[k] * sni1_1263[k];
    }

#pragma omp simd aligned(t_1572, t_1573, t_1574, pb_x, pb_z, pc_x, pc_y, pc_z, sni0_1266, \
                         sni0_1573, snh_968, snh_1181, sni1_1266, sni1_1573, \
                         soh_1178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1572[k] = f_16 * snh_968[k]
                    + f_3 * pc_y[k] * soh_1178[k];

        t_1573[k] = pb_x[k] * sni0_1573[k]
                    + f_14 * snh_1181[k]
                    - f_10 * pc_x[k] * sni1_1573[k];

        t_1574[k] = pb_z[k] * sni0_1266[k]
                    - f_10 * pc_z[k] * sni1_1266[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t soh, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sni0_1270 = buffer.data(sni0 + 1270);
    const auto *sni0_1577 = buffer.data(sni0 + 1577);
    const auto *sni0_1580 = buffer.data(sni0 + 1580);
    const auto *sni0_1582 = buffer.data(sni0 + 1582);
    const auto *sni0_1589 = buffer.data(sni0 + 1589);
    const auto *sni0_1591 = buffer.data(sni0 + 1591);
    const auto *sni0_1592 = buffer.data(sni0 + 1592);
    const auto *sni0_1593 = buffer.data(sni0 + 1593);
    const auto *sni0_1595 = buffer.data(sni0 + 1595);
    const auto *sni0_1596 = buffer.data(sni0 + 1596);
    const auto *sni0_1599 = buffer.data(sni0 + 1599);
    const auto *sni0_1601 = buffer.data(sni0 + 1601);
    const auto *sni0_1602 = buffer.data(sni0 + 1602);
    const auto *sni0_1605 = buffer.data(sni0 + 1605);
    const auto *sni0_1606 = buffer.data(sni0 + 1606);
    const auto *sni0_1608 = buffer.data(sni0 + 1608);
    const auto *sni0_1610 = buffer.data(sni0 + 1610);
    const auto *sni0_1617 = buffer.data(sni0 + 1617);
    const auto *sni0_1619 = buffer.data(sni0 + 1619);
    const auto *sni0_1620 = buffer.data(sni0 + 1620);
    const auto *sni0_1621 = buffer.data(sni0 + 1621);
    const auto *sni0_1623 = buffer.data(sni0 + 1623);
    const auto *sni0_1624 = buffer.data(sni0 + 1624);
    const auto *sni0_1627 = buffer.data(sni0 + 1627);
    const auto *sni0_1629 = buffer.data(sni0 + 1629);
    const auto *sni0_1630 = buffer.data(sni0 + 1630);
    const auto *sni0_1633 = buffer.data(sni0 + 1633);
    const auto *sni0_1634 = buffer.data(sni0 + 1634);
    const auto *sni0_1636 = buffer.data(sni0 + 1636);
    const auto *sni0_1638 = buffer.data(sni0 + 1638);
    const auto *sni0_1645 = buffer.data(sni0 + 1645);
    const auto *sni0_1647 = buffer.data(sni0 + 1647);
    const auto *sni0_1648 = buffer.data(sni0 + 1648);
    const auto *sni0_1649 = buffer.data(sni0 + 1649);
    const auto *sni0_1651 = buffer.data(sni0 + 1651);
    const auto *sni0_1652 = buffer.data(sni0 + 1652);
    const auto *sni0_1655 = buffer.data(sni0 + 1655);
    const auto *sni0_1657 = buffer.data(sni0 + 1657);
    const auto *sni0_1658 = buffer.data(sni0 + 1658);
    const auto *sni0_1661 = buffer.data(sni0 + 1661);
    const auto *sni0_1662 = buffer.data(sni0 + 1662);
    const auto *sni0_1664 = buffer.data(sni0 + 1664);
    const auto *sni0_1666 = buffer.data(sni0 + 1666);
    const auto *sni0_1673 = buffer.data(sni0 + 1673);
    const auto *sni0_1675 = buffer.data(sni0 + 1675);
    const auto *sni0_1676 = buffer.data(sni0 + 1676);
    const auto *sni0_1677 = buffer.data(sni0 + 1677);
    const auto *sni0_1679 = buffer.data(sni0 + 1679);
    const auto *sni0_1680 = buffer.data(sni0 + 1680);
    const auto *sni0_1683 = buffer.data(sni0 + 1683);
    const auto *sni0_1685 = buffer.data(sni0 + 1685);
    const auto *sni0_1686 = buffer.data(sni0 + 1686);
    const auto *sni0_1689 = buffer.data(sni0 + 1689);
    const auto *sni0_1690 = buffer.data(sni0 + 1690);
    const auto *sni0_1692 = buffer.data(sni0 + 1692);

    const auto *snh_948 = buffer.data(snh + 948);
    const auto *snh_951 = buffer.data(snh + 951);
    const auto *snh_960 = buffer.data(snh + 960);
    const auto *snh_966 = buffer.data(snh + 966);
    const auto *snh_969 = buffer.data(snh + 969);
    const auto *snh_971 = buffer.data(snh + 971);
    const auto *snh_972 = buffer.data(snh + 972);
    const auto *snh_975 = buffer.data(snh + 975);
    const auto *snh_981 = buffer.data(snh + 981);
    const auto *snh_986 = buffer.data(snh + 986);
    const auto *snh_987 = buffer.data(snh + 987);
    const auto *snh_989 = buffer.data(snh + 989);
    const auto *snh_990 = buffer.data(snh + 990);
    const auto *snh_992 = buffer.data(snh + 992);
    const auto *snh_993 = buffer.data(snh + 993);
    const auto *snh_996 = buffer.data(snh + 996);
    const auto *snh_1002 = buffer.data(snh + 1002);
    const auto *snh_1007 = buffer.data(snh + 1007);
    const auto *snh_1008 = buffer.data(snh + 1008);
    const auto *snh_1010 = buffer.data(snh + 1010);
    const auto *snh_1011 = buffer.data(snh + 1011);
    const auto *snh_1013 = buffer.data(snh + 1013);
    const auto *snh_1014 = buffer.data(snh + 1014);
    const auto *snh_1017 = buffer.data(snh + 1017);
    const auto *snh_1023 = buffer.data(snh + 1023);
    const auto *snh_1028 = buffer.data(snh + 1028);
    const auto *snh_1029 = buffer.data(snh + 1029);
    const auto *snh_1031 = buffer.data(snh + 1031);
    const auto *snh_1032 = buffer.data(snh + 1032);
    const auto *snh_1034 = buffer.data(snh + 1034);
    const auto *snh_1035 = buffer.data(snh + 1035);
    const auto *snh_1038 = buffer.data(snh + 1038);
    const auto *snh_1049 = buffer.data(snh + 1049);
    const auto *snh_1050 = buffer.data(snh + 1050);
    const auto *snh_1052 = buffer.data(snh + 1052);
    const auto *snh_1055 = buffer.data(snh + 1055);
    const auto *snh_1059 = buffer.data(snh + 1059);
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
    const auto *snh_1200 = buffer.data(snh + 1200);
    const auto *snh_1202 = buffer.data(snh + 1202);
    const auto *snh_1203 = buffer.data(snh + 1203);
    const auto *snh_1206 = buffer.data(snh + 1206);
    const auto *snh_1207 = buffer.data(snh + 1207);
    const auto *snh_1209 = buffer.data(snh + 1209);
    const auto *snh_1211 = buffer.data(snh + 1211);
    const auto *snh_1212 = buffer.data(snh + 1212);
    const auto *snh_1213 = buffer.data(snh + 1213);
    const auto *snh_1214 = buffer.data(snh + 1214);
    const auto *snh_1215 = buffer.data(snh + 1215);
    const auto *snh_1216 = buffer.data(snh + 1216);
    const auto *snh_1217 = buffer.data(snh + 1217);
    const auto *snh_1218 = buffer.data(snh + 1218);
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
    const auto *snh_1263 = buffer.data(snh + 1263);
    const auto *snh_1265 = buffer.data(snh + 1265);
    const auto *snh_1266 = buffer.data(snh + 1266);
    const auto *snh_1269 = buffer.data(snh + 1269);
    const auto *snh_1270 = buffer.data(snh + 1270);
    const auto *snh_1272 = buffer.data(snh + 1272);

    const auto *sni1_1270 = buffer.data(sni1 + 1270);
    const auto *sni1_1577 = buffer.data(sni1 + 1577);
    const auto *sni1_1580 = buffer.data(sni1 + 1580);
    const auto *sni1_1582 = buffer.data(sni1 + 1582);
    const auto *sni1_1589 = buffer.data(sni1 + 1589);
    const auto *sni1_1591 = buffer.data(sni1 + 1591);
    const auto *sni1_1592 = buffer.data(sni1 + 1592);
    const auto *sni1_1593 = buffer.data(sni1 + 1593);
    const auto *sni1_1595 = buffer.data(sni1 + 1595);
    const auto *sni1_1596 = buffer.data(sni1 + 1596);
    const auto *sni1_1599 = buffer.data(sni1 + 1599);
    const auto *sni1_1601 = buffer.data(sni1 + 1601);
    const auto *sni1_1602 = buffer.data(sni1 + 1602);
    const auto *sni1_1605 = buffer.data(sni1 + 1605);
    const auto *sni1_1606 = buffer.data(sni1 + 1606);
    const auto *sni1_1608 = buffer.data(sni1 + 1608);
    const auto *sni1_1610 = buffer.data(sni1 + 1610);
    const auto *sni1_1617 = buffer.data(sni1 + 1617);
    const auto *sni1_1619 = buffer.data(sni1 + 1619);
    const auto *sni1_1620 = buffer.data(sni1 + 1620);
    const auto *sni1_1621 = buffer.data(sni1 + 1621);
    const auto *sni1_1623 = buffer.data(sni1 + 1623);
    const auto *sni1_1624 = buffer.data(sni1 + 1624);
    const auto *sni1_1627 = buffer.data(sni1 + 1627);
    const auto *sni1_1629 = buffer.data(sni1 + 1629);
    const auto *sni1_1630 = buffer.data(sni1 + 1630);
    const auto *sni1_1633 = buffer.data(sni1 + 1633);
    const auto *sni1_1634 = buffer.data(sni1 + 1634);
    const auto *sni1_1636 = buffer.data(sni1 + 1636);
    const auto *sni1_1638 = buffer.data(sni1 + 1638);
    const auto *sni1_1645 = buffer.data(sni1 + 1645);
    const auto *sni1_1647 = buffer.data(sni1 + 1647);
    const auto *sni1_1648 = buffer.data(sni1 + 1648);
    const auto *sni1_1649 = buffer.data(sni1 + 1649);
    const auto *sni1_1651 = buffer.data(sni1 + 1651);
    const auto *sni1_1652 = buffer.data(sni1 + 1652);
    const auto *sni1_1655 = buffer.data(sni1 + 1655);
    const auto *sni1_1657 = buffer.data(sni1 + 1657);
    const auto *sni1_1658 = buffer.data(sni1 + 1658);
    const auto *sni1_1661 = buffer.data(sni1 + 1661);
    const auto *sni1_1662 = buffer.data(sni1 + 1662);
    const auto *sni1_1664 = buffer.data(sni1 + 1664);
    const auto *sni1_1666 = buffer.data(sni1 + 1666);
    const auto *sni1_1673 = buffer.data(sni1 + 1673);
    const auto *sni1_1675 = buffer.data(sni1 + 1675);
    const auto *sni1_1676 = buffer.data(sni1 + 1676);
    const auto *sni1_1677 = buffer.data(sni1 + 1677);
    const auto *sni1_1679 = buffer.data(sni1 + 1679);
    const auto *sni1_1680 = buffer.data(sni1 + 1680);
    const auto *sni1_1683 = buffer.data(sni1 + 1683);
    const auto *sni1_1685 = buffer.data(sni1 + 1685);
    const auto *sni1_1686 = buffer.data(sni1 + 1686);
    const auto *sni1_1689 = buffer.data(sni1 + 1689);
    const auto *sni1_1690 = buffer.data(sni1 + 1690);
    const auto *sni1_1692 = buffer.data(sni1 + 1692);

    const auto *soh_1179 = buffer.data(soh + 1179);
    const auto *soh_1181 = buffer.data(soh + 1181);
    const auto *soh_1182 = buffer.data(soh + 1182);
    const auto *soh_1185 = buffer.data(soh + 1185);
    const auto *soh_1191 = buffer.data(soh + 1191);
    const auto *soh_1192 = buffer.data(soh + 1192);
    const auto *soh_1193 = buffer.data(soh + 1193);
    const auto *soh_1194 = buffer.data(soh + 1194);
    const auto *soh_1195 = buffer.data(soh + 1195);
    const auto *soh_1196 = buffer.data(soh + 1196);
    const auto *soh_1197 = buffer.data(soh + 1197);
    const auto *soh_1199 = buffer.data(soh + 1199);
    const auto *soh_1200 = buffer.data(soh + 1200);
    const auto *soh_1202 = buffer.data(soh + 1202);
    const auto *soh_1203 = buffer.data(soh + 1203);
    const auto *soh_1206 = buffer.data(soh + 1206);
    const auto *soh_1212 = buffer.data(soh + 1212);
    const auto *soh_1213 = buffer.data(soh + 1213);
    const auto *soh_1214 = buffer.data(soh + 1214);
    const auto *soh_1215 = buffer.data(soh + 1215);
    const auto *soh_1216 = buffer.data(soh + 1216);
    const auto *soh_1217 = buffer.data(soh + 1217);
    const auto *soh_1218 = buffer.data(soh + 1218);
    const auto *soh_1220 = buffer.data(soh + 1220);
    const auto *soh_1221 = buffer.data(soh + 1221);
    const auto *soh_1223 = buffer.data(soh + 1223);
    const auto *soh_1224 = buffer.data(soh + 1224);
    const auto *soh_1227 = buffer.data(soh + 1227);
    const auto *soh_1233 = buffer.data(soh + 1233);
    const auto *soh_1234 = buffer.data(soh + 1234);
    const auto *soh_1235 = buffer.data(soh + 1235);
    const auto *soh_1236 = buffer.data(soh + 1236);
    const auto *soh_1237 = buffer.data(soh + 1237);
    const auto *soh_1238 = buffer.data(soh + 1238);
    const auto *soh_1239 = buffer.data(soh + 1239);
    const auto *soh_1241 = buffer.data(soh + 1241);
    const auto *soh_1242 = buffer.data(soh + 1242);
    const auto *soh_1244 = buffer.data(soh + 1244);
    const auto *soh_1245 = buffer.data(soh + 1245);
    const auto *soh_1248 = buffer.data(soh + 1248);
    const auto *soh_1254 = buffer.data(soh + 1254);
    const auto *soh_1255 = buffer.data(soh + 1255);
    const auto *soh_1256 = buffer.data(soh + 1256);
    const auto *soh_1257 = buffer.data(soh + 1257);
    const auto *soh_1258 = buffer.data(soh + 1258);
    const auto *soh_1259 = buffer.data(soh + 1259);
    const auto *soh_1260 = buffer.data(soh + 1260);
    const auto *soh_1262 = buffer.data(soh + 1262);
    const auto *soh_1263 = buffer.data(soh + 1263);
    const auto *soh_1265 = buffer.data(soh + 1265);
    const auto *soh_1266 = buffer.data(soh + 1266);
    const auto *soh_1269 = buffer.data(soh + 1269);

#pragma omp simd aligned(t_1575, t_1576, t_1577, pb_x, pc_x, pc_y, pc_z, sni0_1577, snh_948, \
                         snh_971, snh_1185, sni1_1577, soh_1179, \
                         soh_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_11 * snh_948[k]
                    + f_3 * pc_z[k] * soh_1179[k];

        t_1576[k] = f_16 * snh_971[k]
                    + f_3 * pc_y[k] * soh_1181[k];

        t_1577[k] = pb_x[k] * sni0_1577[k]
                    + f_13 * snh_1185[k]
                    - f_10 * pc_x[k] * sni1_1577[k];
    }

#pragma omp simd aligned(t_1578, t_1579, t_1580, pb_x, pb_z, pc_x, pc_z, sni0_1270, sni0_1580, \
                         snh_951, snh_1188, sni1_1270, sni1_1580, \
                         soh_1182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1578[k] = pb_z[k] * sni0_1270[k]
                    - f_10 * pc_z[k] * sni1_1270[k];

        t_1579[k] = f_11 * snh_951[k]
                    + f_3 * pc_z[k] * soh_1182[k];

        t_1580[k] = pb_x[k] * sni0_1580[k]
                    + f_12 * snh_1188[k]
                    - f_10 * pc_x[k] * sni1_1580[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, t_1584, pb_x, pc_x, pc_y, sni0_1582, snh_975, \
                         snh_1190, snh_1191, snh_1192, sni1_1582, soh_1185, soh_1191, \
                         soh_1192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_16 * snh_975[k]
                    + f_3 * pc_y[k] * soh_1185[k];

        t_1582[k] = pb_x[k] * sni0_1582[k]
                    + f_12 * snh_1190[k]
                    - f_10 * pc_x[k] * sni1_1582[k];

        t_1583[k] = f_11 * snh_1191[k]
                    + f_3 * pc_x[k] * soh_1191[k];

        t_1584[k] = f_11 * snh_1192[k]
                    + f_3 * pc_x[k] * soh_1192[k];
    }

#pragma omp simd aligned(t_1585, t_1586, t_1587, t_1588, pc_x, snh_1193, snh_1194, snh_1195, \
                         snh_1196, soh_1193, soh_1194, soh_1195, \
                         soh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1585[k] = f_11 * snh_1193[k]
                    + f_3 * pc_x[k] * soh_1193[k];

        t_1586[k] = f_11 * snh_1194[k]
                    + f_3 * pc_x[k] * soh_1194[k];

        t_1587[k] = f_11 * snh_1195[k]
                    + f_3 * pc_x[k] * soh_1195[k];

        t_1588[k] = f_11 * snh_1196[k]
                    + f_3 * pc_x[k] * soh_1196[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, t_1592, pb_x, pc_x, pc_z, sni0_1589, \
                         sni0_1591, sni0_1592, snh_960, sni1_1589, sni1_1591, sni1_1592, \
                         soh_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = pb_x[k] * sni0_1589[k]
                    - f_10 * pc_x[k] * sni1_1589[k];

        t_1590[k] = f_11 * snh_960[k]
                    + f_3 * pc_z[k] * soh_1191[k];

        t_1591[k] = pb_x[k] * sni0_1591[k]
                    - f_10 * pc_x[k] * sni1_1591[k];

        t_1592[k] = pb_x[k] * sni0_1592[k]
                    - f_10 * pc_x[k] * sni1_1592[k];
    }

#pragma omp simd aligned(t_1593, t_1594, t_1595, t_1596, pb_x, pc_x, pc_y, sni0_1593, \
                         sni0_1595, sni0_1596, snh_986, snh_1197, sni1_1593, sni1_1595, \
                         sni1_1596, soh_1196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1593[k] = pb_x[k] * sni0_1593[k]
                    - f_10 * pc_x[k] * sni1_1593[k];

        t_1594[k] = f_16 * snh_986[k]
                    + f_3 * pc_y[k] * soh_1196[k];

        t_1595[k] = pb_x[k] * sni0_1595[k]
                    - f_10 * pc_x[k] * sni1_1595[k];

        t_1596[k] = pb_x[k] * sni0_1596[k]
                    + f_19 * snh_1197[k]
                    - f_10 * pc_x[k] * sni1_1596[k];
    }

#pragma omp simd aligned(t_1597, t_1598, t_1599, t_1600, pb_x, pc_x, pc_y, pc_z, sni0_1599, \
                         snh_966, snh_987, snh_989, snh_1200, sni1_1599, soh_1197, \
                         soh_1199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1597[k] = f_17 * snh_987[k]
                    + f_3 * pc_y[k] * soh_1197[k];

        t_1598[k] = f_12 * snh_966[k]
                    + f_3 * pc_z[k] * soh_1197[k];

        t_1599[k] = pb_x[k] * sni0_1599[k]
                    + f_14 * snh_1200[k]
                    - f_10 * pc_x[k] * sni1_1599[k];

        t_1600[k] = f_17 * snh_989[k]
                    + f_3 * pc_y[k] * soh_1199[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, pb_x, pc_x, pc_z, sni0_1601, sni0_1602, \
                         snh_969, snh_1202, snh_1203, sni1_1601, sni1_1602, \
                         soh_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = pb_x[k] * sni0_1601[k]
                    + f_14 * snh_1202[k]
                    - f_10 * pc_x[k] * sni1_1601[k];

        t_1602[k] = pb_x[k] * sni0_1602[k]
                    + f_13 * snh_1203[k]
                    - f_10 * pc_x[k] * sni1_1602[k];

        t_1603[k] = f_12 * snh_969[k]
                    + f_3 * pc_z[k] * soh_1200[k];
    }

#pragma omp simd aligned(t_1604, t_1605, t_1606, pb_x, pc_x, pc_y, sni0_1605, sni0_1606, \
                         snh_992, snh_1206, snh_1207, sni1_1605, sni1_1606, \
                         soh_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1604[k] = f_17 * snh_992[k]
                    + f_3 * pc_y[k] * soh_1202[k];

        t_1605[k] = pb_x[k] * sni0_1605[k]
                    + f_13 * snh_1206[k]
                    - f_10 * pc_x[k] * sni1_1605[k];

        t_1606[k] = pb_x[k] * sni0_1606[k]
                    + f_12 * snh_1207[k]
                    - f_10 * pc_x[k] * sni1_1606[k];
    }

#pragma omp simd aligned(t_1607, t_1608, t_1609, pb_x, pc_x, pc_y, pc_z, sni0_1608, snh_972, \
                         snh_996, snh_1209, sni1_1608, soh_1203, \
                         soh_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1607[k] = f_12 * snh_972[k]
                    + f_3 * pc_z[k] * soh_1203[k];

        t_1608[k] = pb_x[k] * sni0_1608[k]
                    + f_12 * snh_1209[k]
                    - f_10 * pc_x[k] * sni1_1608[k];

        t_1609[k] = f_17 * snh_996[k]
                    + f_3 * pc_y[k] * soh_1206[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, pb_x, pc_x, sni0_1610, snh_1211, \
                         snh_1212, snh_1213, snh_1214, sni1_1610, soh_1212, soh_1213, \
                         soh_1214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = pb_x[k] * sni0_1610[k]
                    + f_12 * snh_1211[k]
                    - f_10 * pc_x[k] * sni1_1610[k];

        t_1611[k] = f_11 * snh_1212[k]
                    + f_3 * pc_x[k] * soh_1212[k];

        t_1612[k] = f_11 * snh_1213[k]
                    + f_3 * pc_x[k] * soh_1213[k];

        t_1613[k] = f_11 * snh_1214[k]
                    + f_3 * pc_x[k] * soh_1214[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, t_1617, pb_x, pc_x, sni0_1617, snh_1215, \
                         snh_1216, snh_1217, sni1_1617, soh_1215, soh_1216, \
                         soh_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_11 * snh_1215[k]
                    + f_3 * pc_x[k] * soh_1215[k];

        t_1615[k] = f_11 * snh_1216[k]
                    + f_3 * pc_x[k] * soh_1216[k];

        t_1616[k] = f_11 * snh_1217[k]
                    + f_3 * pc_x[k] * soh_1217[k];

        t_1617[k] = pb_x[k] * sni0_1617[k]
                    - f_10 * pc_x[k] * sni1_1617[k];
    }

#pragma omp simd aligned(t_1618, t_1619, t_1620, t_1621, pb_x, pc_x, pc_z, sni0_1619, \
                         sni0_1620, sni0_1621, snh_981, sni1_1619, sni1_1620, sni1_1621, \
                         soh_1212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1618[k] = f_12 * snh_981[k]
                    + f_3 * pc_z[k] * soh_1212[k];

        t_1619[k] = pb_x[k] * sni0_1619[k]
                    - f_10 * pc_x[k] * sni1_1619[k];

        t_1620[k] = pb_x[k] * sni0_1620[k]
                    - f_10 * pc_x[k] * sni1_1620[k];

        t_1621[k] = pb_x[k] * sni0_1621[k]
                    - f_10 * pc_x[k] * sni1_1621[k];
    }

#pragma omp simd aligned(t_1622, t_1623, t_1624, t_1625, pb_x, pc_x, pc_y, sni0_1623, \
                         sni0_1624, snh_1007, snh_1008, snh_1218, sni1_1623, sni1_1624, \
                         soh_1217, soh_1218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1622[k] = f_17 * snh_1007[k]
                    + f_3 * pc_y[k] * soh_1217[k];

        t_1623[k] = pb_x[k] * sni0_1623[k]
                    - f_10 * pc_x[k] * sni1_1623[k];

        t_1624[k] = pb_x[k] * sni0_1624[k]
                    + f_19 * snh_1218[k]
                    - f_10 * pc_x[k] * sni1_1624[k];

        t_1625[k] = f_18 * snh_1008[k]
                    + f_3 * pc_y[k] * soh_1218[k];
    }

#pragma omp simd aligned(t_1626, t_1627, t_1628, pb_x, pc_x, pc_y, pc_z, sni0_1627, snh_987, \
                         snh_1010, snh_1221, sni1_1627, soh_1218, \
                         soh_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1626[k] = f_13 * snh_987[k]
                    + f_3 * pc_z[k] * soh_1218[k];

        t_1627[k] = pb_x[k] * sni0_1627[k]
                    + f_14 * snh_1221[k]
                    - f_10 * pc_x[k] * sni1_1627[k];

        t_1628[k] = f_18 * snh_1010[k]
                    + f_3 * pc_y[k] * soh_1220[k];
    }

#pragma omp simd aligned(t_1629, t_1630, t_1631, pb_x, pc_x, pc_z, sni0_1629, sni0_1630, \
                         snh_990, snh_1223, snh_1224, sni1_1629, sni1_1630, \
                         soh_1221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1629[k] = pb_x[k] * sni0_1629[k]
                    + f_14 * snh_1223[k]
                    - f_10 * pc_x[k] * sni1_1629[k];

        t_1630[k] = pb_x[k] * sni0_1630[k]
                    + f_13 * snh_1224[k]
                    - f_10 * pc_x[k] * sni1_1630[k];

        t_1631[k] = f_13 * snh_990[k]
                    + f_3 * pc_z[k] * soh_1221[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, pb_x, pc_x, pc_y, sni0_1633, sni0_1634, \
                         snh_1013, snh_1227, snh_1228, sni1_1633, sni1_1634, \
                         soh_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_18 * snh_1013[k]
                    + f_3 * pc_y[k] * soh_1223[k];

        t_1633[k] = pb_x[k] * sni0_1633[k]
                    + f_13 * snh_1227[k]
                    - f_10 * pc_x[k] * sni1_1633[k];

        t_1634[k] = pb_x[k] * sni0_1634[k]
                    + f_12 * snh_1228[k]
                    - f_10 * pc_x[k] * sni1_1634[k];
    }

#pragma omp simd aligned(t_1635, t_1636, t_1637, pb_x, pc_x, pc_y, pc_z, sni0_1636, snh_993, \
                         snh_1017, snh_1230, sni1_1636, soh_1224, \
                         soh_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1635[k] = f_13 * snh_993[k]
                    + f_3 * pc_z[k] * soh_1224[k];

        t_1636[k] = pb_x[k] * sni0_1636[k]
                    + f_12 * snh_1230[k]
                    - f_10 * pc_x[k] * sni1_1636[k];

        t_1637[k] = f_18 * snh_1017[k]
                    + f_3 * pc_y[k] * soh_1227[k];
    }

#pragma omp simd aligned(t_1638, t_1639, t_1640, t_1641, pb_x, pc_x, sni0_1638, snh_1232, \
                         snh_1233, snh_1234, snh_1235, sni1_1638, soh_1233, soh_1234, \
                         soh_1235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1638[k] = pb_x[k] * sni0_1638[k]
                    + f_12 * snh_1232[k]
                    - f_10 * pc_x[k] * sni1_1638[k];

        t_1639[k] = f_11 * snh_1233[k]
                    + f_3 * pc_x[k] * soh_1233[k];

        t_1640[k] = f_11 * snh_1234[k]
                    + f_3 * pc_x[k] * soh_1234[k];

        t_1641[k] = f_11 * snh_1235[k]
                    + f_3 * pc_x[k] * soh_1235[k];
    }

#pragma omp simd aligned(t_1642, t_1643, t_1644, t_1645, pb_x, pc_x, sni0_1645, snh_1236, \
                         snh_1237, snh_1238, sni1_1645, soh_1236, soh_1237, \
                         soh_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1642[k] = f_11 * snh_1236[k]
                    + f_3 * pc_x[k] * soh_1236[k];

        t_1643[k] = f_11 * snh_1237[k]
                    + f_3 * pc_x[k] * soh_1237[k];

        t_1644[k] = f_11 * snh_1238[k]
                    + f_3 * pc_x[k] * soh_1238[k];

        t_1645[k] = pb_x[k] * sni0_1645[k]
                    - f_10 * pc_x[k] * sni1_1645[k];
    }

#pragma omp simd aligned(t_1646, t_1647, t_1648, t_1649, pb_x, pc_x, pc_z, sni0_1647, \
                         sni0_1648, sni0_1649, snh_1002, sni1_1647, sni1_1648, sni1_1649, \
                         soh_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1646[k] = f_13 * snh_1002[k]
                    + f_3 * pc_z[k] * soh_1233[k];

        t_1647[k] = pb_x[k] * sni0_1647[k]
                    - f_10 * pc_x[k] * sni1_1647[k];

        t_1648[k] = pb_x[k] * sni0_1648[k]
                    - f_10 * pc_x[k] * sni1_1648[k];

        t_1649[k] = pb_x[k] * sni0_1649[k]
                    - f_10 * pc_x[k] * sni1_1649[k];
    }

#pragma omp simd aligned(t_1650, t_1651, t_1652, t_1653, pb_x, pc_x, pc_y, sni0_1651, \
                         sni0_1652, snh_1028, snh_1029, snh_1239, sni1_1651, sni1_1652, \
                         soh_1238, soh_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1650[k] = f_18 * snh_1028[k]
                    + f_3 * pc_y[k] * soh_1238[k];

        t_1651[k] = pb_x[k] * sni0_1651[k]
                    - f_10 * pc_x[k] * sni1_1651[k];

        t_1652[k] = pb_x[k] * sni0_1652[k]
                    + f_19 * snh_1239[k]
                    - f_10 * pc_x[k] * sni1_1652[k];

        t_1653[k] = f_19 * snh_1029[k]
                    + f_3 * pc_y[k] * soh_1239[k];
    }

#pragma omp simd aligned(t_1654, t_1655, t_1656, pb_x, pc_x, pc_y, pc_z, sni0_1655, snh_1008, \
                         snh_1031, snh_1242, sni1_1655, soh_1239, \
                         soh_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1654[k] = f_14 * snh_1008[k]
                    + f_3 * pc_z[k] * soh_1239[k];

        t_1655[k] = pb_x[k] * sni0_1655[k]
                    + f_14 * snh_1242[k]
                    - f_10 * pc_x[k] * sni1_1655[k];

        t_1656[k] = f_19 * snh_1031[k]
                    + f_3 * pc_y[k] * soh_1241[k];
    }

#pragma omp simd aligned(t_1657, t_1658, t_1659, pb_x, pc_x, pc_z, sni0_1657, sni0_1658, \
                         snh_1011, snh_1244, snh_1245, sni1_1657, sni1_1658, \
                         soh_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1657[k] = pb_x[k] * sni0_1657[k]
                    + f_14 * snh_1244[k]
                    - f_10 * pc_x[k] * sni1_1657[k];

        t_1658[k] = pb_x[k] * sni0_1658[k]
                    + f_13 * snh_1245[k]
                    - f_10 * pc_x[k] * sni1_1658[k];

        t_1659[k] = f_14 * snh_1011[k]
                    + f_3 * pc_z[k] * soh_1242[k];
    }

#pragma omp simd aligned(t_1660, t_1661, t_1662, pb_x, pc_x, pc_y, sni0_1661, sni0_1662, \
                         snh_1034, snh_1248, snh_1249, sni1_1661, sni1_1662, \
                         soh_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1660[k] = f_19 * snh_1034[k]
                    + f_3 * pc_y[k] * soh_1244[k];

        t_1661[k] = pb_x[k] * sni0_1661[k]
                    + f_13 * snh_1248[k]
                    - f_10 * pc_x[k] * sni1_1661[k];

        t_1662[k] = pb_x[k] * sni0_1662[k]
                    + f_12 * snh_1249[k]
                    - f_10 * pc_x[k] * sni1_1662[k];
    }

#pragma omp simd aligned(t_1663, t_1664, t_1665, pb_x, pc_x, pc_y, pc_z, sni0_1664, snh_1014, \
                         snh_1038, snh_1251, sni1_1664, soh_1245, \
                         soh_1248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1663[k] = f_14 * snh_1014[k]
                    + f_3 * pc_z[k] * soh_1245[k];

        t_1664[k] = pb_x[k] * sni0_1664[k]
                    + f_12 * snh_1251[k]
                    - f_10 * pc_x[k] * sni1_1664[k];

        t_1665[k] = f_19 * snh_1038[k]
                    + f_3 * pc_y[k] * soh_1248[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, t_1669, pb_x, pc_x, sni0_1666, snh_1253, \
                         snh_1254, snh_1255, snh_1256, sni1_1666, soh_1254, soh_1255, \
                         soh_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = pb_x[k] * sni0_1666[k]
                    + f_12 * snh_1253[k]
                    - f_10 * pc_x[k] * sni1_1666[k];

        t_1667[k] = f_11 * snh_1254[k]
                    + f_3 * pc_x[k] * soh_1254[k];

        t_1668[k] = f_11 * snh_1255[k]
                    + f_3 * pc_x[k] * soh_1255[k];

        t_1669[k] = f_11 * snh_1256[k]
                    + f_3 * pc_x[k] * soh_1256[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, t_1673, pb_x, pc_x, sni0_1673, snh_1257, \
                         snh_1258, snh_1259, sni1_1673, soh_1257, soh_1258, \
                         soh_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_11 * snh_1257[k]
                    + f_3 * pc_x[k] * soh_1257[k];

        t_1671[k] = f_11 * snh_1258[k]
                    + f_3 * pc_x[k] * soh_1258[k];

        t_1672[k] = f_11 * snh_1259[k]
                    + f_3 * pc_x[k] * soh_1259[k];

        t_1673[k] = pb_x[k] * sni0_1673[k]
                    - f_10 * pc_x[k] * sni1_1673[k];
    }

#pragma omp simd aligned(t_1674, t_1675, t_1676, t_1677, pb_x, pc_x, pc_z, sni0_1675, \
                         sni0_1676, sni0_1677, snh_1023, sni1_1675, sni1_1676, sni1_1677, \
                         soh_1254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1674[k] = f_14 * snh_1023[k]
                    + f_3 * pc_z[k] * soh_1254[k];

        t_1675[k] = pb_x[k] * sni0_1675[k]
                    - f_10 * pc_x[k] * sni1_1675[k];

        t_1676[k] = pb_x[k] * sni0_1676[k]
                    - f_10 * pc_x[k] * sni1_1676[k];

        t_1677[k] = pb_x[k] * sni0_1677[k]
                    - f_10 * pc_x[k] * sni1_1677[k];
    }

#pragma omp simd aligned(t_1678, t_1679, t_1680, t_1681, pb_x, pc_x, pc_y, sni0_1679, \
                         sni0_1680, snh_1049, snh_1050, snh_1260, sni1_1679, sni1_1680, \
                         soh_1259, soh_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1678[k] = f_19 * snh_1049[k]
                    + f_3 * pc_y[k] * soh_1259[k];

        t_1679[k] = pb_x[k] * sni0_1679[k]
                    - f_10 * pc_x[k] * sni1_1679[k];

        t_1680[k] = pb_x[k] * sni0_1680[k]
                    + f_19 * snh_1260[k]
                    - f_10 * pc_x[k] * sni1_1680[k];

        t_1681[k] = f_20 * snh_1050[k]
                    + f_3 * pc_y[k] * soh_1260[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, pb_x, pc_x, pc_y, pc_z, sni0_1683, snh_1029, \
                         snh_1052, snh_1263, sni1_1683, soh_1260, \
                         soh_1262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_20 * snh_1029[k]
                    + f_3 * pc_z[k] * soh_1260[k];

        t_1683[k] = pb_x[k] * sni0_1683[k]
                    + f_14 * snh_1263[k]
                    - f_10 * pc_x[k] * sni1_1683[k];

        t_1684[k] = f_20 * snh_1052[k]
                    + f_3 * pc_y[k] * soh_1262[k];
    }

#pragma omp simd aligned(t_1685, t_1686, t_1687, pb_x, pc_x, pc_z, sni0_1685, sni0_1686, \
                         snh_1032, snh_1265, snh_1266, sni1_1685, sni1_1686, \
                         soh_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1685[k] = pb_x[k] * sni0_1685[k]
                    + f_14 * snh_1265[k]
                    - f_10 * pc_x[k] * sni1_1685[k];

        t_1686[k] = pb_x[k] * sni0_1686[k]
                    + f_13 * snh_1266[k]
                    - f_10 * pc_x[k] * sni1_1686[k];

        t_1687[k] = f_20 * snh_1032[k]
                    + f_3 * pc_z[k] * soh_1263[k];
    }

#pragma omp simd aligned(t_1688, t_1689, t_1690, pb_x, pc_x, pc_y, sni0_1689, sni0_1690, \
                         snh_1055, snh_1269, snh_1270, sni1_1689, sni1_1690, \
                         soh_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1688[k] = f_20 * snh_1055[k]
                    + f_3 * pc_y[k] * soh_1265[k];

        t_1689[k] = pb_x[k] * sni0_1689[k]
                    + f_13 * snh_1269[k]
                    - f_10 * pc_x[k] * sni1_1689[k];

        t_1690[k] = pb_x[k] * sni0_1690[k]
                    + f_12 * snh_1270[k]
                    - f_10 * pc_x[k] * sni1_1690[k];
    }

#pragma omp simd aligned(t_1691, t_1692, t_1693, pb_x, pc_x, pc_y, pc_z, sni0_1692, snh_1035, \
                         snh_1059, snh_1272, sni1_1692, soh_1266, \
                         soh_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1691[k] = f_20 * snh_1035[k]
                    + f_3 * pc_z[k] * soh_1266[k];

        t_1692[k] = pb_x[k] * sni0_1692[k]
                    + f_12 * snh_1272[k]
                    - f_10 * pc_x[k] * sni1_1692[k];

        t_1693[k] = f_20 * snh_1059[k]
                    + f_3 * pc_y[k] * soh_1269[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t soh, const size_t ncols,
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
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sni0_1512 = buffer.data(sni0 + 1512);
    const auto *sni0_1517 = buffer.data(sni0 + 1517);
    const auto *sni0_1521 = buffer.data(sni0 + 1521);
    const auto *sni0_1526 = buffer.data(sni0 + 1526);
    const auto *sni0_1694 = buffer.data(sni0 + 1694);
    const auto *sni0_1701 = buffer.data(sni0 + 1701);
    const auto *sni0_1703 = buffer.data(sni0 + 1703);
    const auto *sni0_1704 = buffer.data(sni0 + 1704);
    const auto *sni0_1705 = buffer.data(sni0 + 1705);
    const auto *sni0_1707 = buffer.data(sni0 + 1707);
    const auto *sni0_1708 = buffer.data(sni0 + 1708);
    const auto *sni0_1711 = buffer.data(sni0 + 1711);
    const auto *sni0_1713 = buffer.data(sni0 + 1713);
    const auto *sni0_1714 = buffer.data(sni0 + 1714);
    const auto *sni0_1717 = buffer.data(sni0 + 1717);
    const auto *sni0_1718 = buffer.data(sni0 + 1718);
    const auto *sni0_1720 = buffer.data(sni0 + 1720);
    const auto *sni0_1722 = buffer.data(sni0 + 1722);
    const auto *sni0_1729 = buffer.data(sni0 + 1729);
    const auto *sni0_1731 = buffer.data(sni0 + 1731);
    const auto *sni0_1732 = buffer.data(sni0 + 1732);
    const auto *sni0_1733 = buffer.data(sni0 + 1733);
    const auto *sni0_1735 = buffer.data(sni0 + 1735);
    const auto *sni0_1736 = buffer.data(sni0 + 1736);
    const auto *sni0_1739 = buffer.data(sni0 + 1739);
    const auto *sni0_1741 = buffer.data(sni0 + 1741);
    const auto *sni0_1742 = buffer.data(sni0 + 1742);
    const auto *sni0_1745 = buffer.data(sni0 + 1745);
    const auto *sni0_1746 = buffer.data(sni0 + 1746);
    const auto *sni0_1748 = buffer.data(sni0 + 1748);
    const auto *sni0_1750 = buffer.data(sni0 + 1750);
    const auto *sni0_1757 = buffer.data(sni0 + 1757);
    const auto *sni0_1759 = buffer.data(sni0 + 1759);
    const auto *sni0_1760 = buffer.data(sni0 + 1760);
    const auto *sni0_1761 = buffer.data(sni0 + 1761);
    const auto *sni0_1763 = buffer.data(sni0 + 1763);
    const auto *sni0_1764 = buffer.data(sni0 + 1764);
    const auto *sni0_1767 = buffer.data(sni0 + 1767);
    const auto *sni0_1769 = buffer.data(sni0 + 1769);
    const auto *sni0_1770 = buffer.data(sni0 + 1770);
    const auto *sni0_1773 = buffer.data(sni0 + 1773);
    const auto *sni0_1774 = buffer.data(sni0 + 1774);
    const auto *sni0_1776 = buffer.data(sni0 + 1776);
    const auto *sni0_1778 = buffer.data(sni0 + 1778);
    const auto *sni0_1785 = buffer.data(sni0 + 1785);
    const auto *sni0_1787 = buffer.data(sni0 + 1787);
    const auto *sni0_1788 = buffer.data(sni0 + 1788);
    const auto *sni0_1789 = buffer.data(sni0 + 1789);
    const auto *sni0_1791 = buffer.data(sni0 + 1791);
    const auto *sni0_1795 = buffer.data(sni0 + 1795);
    const auto *sni0_1798 = buffer.data(sni0 + 1798);
    const auto *sni0_1802 = buffer.data(sni0 + 1802);
    const auto *sni0_1804 = buffer.data(sni0 + 1804);
    const auto *sni0_1813 = buffer.data(sni0 + 1813);

    const auto *snh_1044 = buffer.data(snh + 1044);
    const auto *snh_1050 = buffer.data(snh + 1050);
    const auto *snh_1053 = buffer.data(snh + 1053);
    const auto *snh_1056 = buffer.data(snh + 1056);
    const auto *snh_1065 = buffer.data(snh + 1065);
    const auto *snh_1070 = buffer.data(snh + 1070);
    const auto *snh_1071 = buffer.data(snh + 1071);
    const auto *snh_1073 = buffer.data(snh + 1073);
    const auto *snh_1074 = buffer.data(snh + 1074);
    const auto *snh_1076 = buffer.data(snh + 1076);
    const auto *snh_1077 = buffer.data(snh + 1077);
    const auto *snh_1080 = buffer.data(snh + 1080);
    const auto *snh_1086 = buffer.data(snh + 1086);
    const auto *snh_1091 = buffer.data(snh + 1091);
    const auto *snh_1092 = buffer.data(snh + 1092);
    const auto *snh_1094 = buffer.data(snh + 1094);
    const auto *snh_1095 = buffer.data(snh + 1095);
    const auto *snh_1097 = buffer.data(snh + 1097);
    const auto *snh_1098 = buffer.data(snh + 1098);
    const auto *snh_1101 = buffer.data(snh + 1101);
    const auto *snh_1107 = buffer.data(snh + 1107);
    const auto *snh_1112 = buffer.data(snh + 1112);
    const auto *snh_1113 = buffer.data(snh + 1113);
    const auto *snh_1115 = buffer.data(snh + 1115);
    const auto *snh_1116 = buffer.data(snh + 1116);
    const auto *snh_1118 = buffer.data(snh + 1118);
    const auto *snh_1119 = buffer.data(snh + 1119);
    const auto *snh_1122 = buffer.data(snh + 1122);
    const auto *snh_1133 = buffer.data(snh + 1133);
    const auto *snh_1134 = buffer.data(snh + 1134);
    const auto *snh_1136 = buffer.data(snh + 1136);
    const auto *snh_1139 = buffer.data(snh + 1139);
    const auto *snh_1143 = buffer.data(snh + 1143);
    const auto *snh_1274 = buffer.data(snh + 1274);
    const auto *snh_1275 = buffer.data(snh + 1275);
    const auto *snh_1276 = buffer.data(snh + 1276);
    const auto *snh_1277 = buffer.data(snh + 1277);
    const auto *snh_1278 = buffer.data(snh + 1278);
    const auto *snh_1279 = buffer.data(snh + 1279);
    const auto *snh_1280 = buffer.data(snh + 1280);
    const auto *snh_1281 = buffer.data(snh + 1281);
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
    const auto *snh_1302 = buffer.data(snh + 1302);
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
    const auto *snh_1347 = buffer.data(snh + 1347);
    const auto *snh_1350 = buffer.data(snh + 1350);
    const auto *snh_1354 = buffer.data(snh + 1354);
    const auto *snh_1356 = buffer.data(snh + 1356);
    const auto *snh_1359 = buffer.data(snh + 1359);
    const auto *snh_1360 = buffer.data(snh + 1360);
    const auto *snh_1361 = buffer.data(snh + 1361);
    const auto *snh_1362 = buffer.data(snh + 1362);
    const auto *snh_1363 = buffer.data(snh + 1363);
    const auto *snh_1364 = buffer.data(snh + 1364);

    const auto *sni1_1512 = buffer.data(sni1 + 1512);
    const auto *sni1_1517 = buffer.data(sni1 + 1517);
    const auto *sni1_1521 = buffer.data(sni1 + 1521);
    const auto *sni1_1526 = buffer.data(sni1 + 1526);
    const auto *sni1_1694 = buffer.data(sni1 + 1694);
    const auto *sni1_1701 = buffer.data(sni1 + 1701);
    const auto *sni1_1703 = buffer.data(sni1 + 1703);
    const auto *sni1_1704 = buffer.data(sni1 + 1704);
    const auto *sni1_1705 = buffer.data(sni1 + 1705);
    const auto *sni1_1707 = buffer.data(sni1 + 1707);
    const auto *sni1_1708 = buffer.data(sni1 + 1708);
    const auto *sni1_1711 = buffer.data(sni1 + 1711);
    const auto *sni1_1713 = buffer.data(sni1 + 1713);
    const auto *sni1_1714 = buffer.data(sni1 + 1714);
    const auto *sni1_1717 = buffer.data(sni1 + 1717);
    const auto *sni1_1718 = buffer.data(sni1 + 1718);
    const auto *sni1_1720 = buffer.data(sni1 + 1720);
    const auto *sni1_1722 = buffer.data(sni1 + 1722);
    const auto *sni1_1729 = buffer.data(sni1 + 1729);
    const auto *sni1_1731 = buffer.data(sni1 + 1731);
    const auto *sni1_1732 = buffer.data(sni1 + 1732);
    const auto *sni1_1733 = buffer.data(sni1 + 1733);
    const auto *sni1_1735 = buffer.data(sni1 + 1735);
    const auto *sni1_1736 = buffer.data(sni1 + 1736);
    const auto *sni1_1739 = buffer.data(sni1 + 1739);
    const auto *sni1_1741 = buffer.data(sni1 + 1741);
    const auto *sni1_1742 = buffer.data(sni1 + 1742);
    const auto *sni1_1745 = buffer.data(sni1 + 1745);
    const auto *sni1_1746 = buffer.data(sni1 + 1746);
    const auto *sni1_1748 = buffer.data(sni1 + 1748);
    const auto *sni1_1750 = buffer.data(sni1 + 1750);
    const auto *sni1_1757 = buffer.data(sni1 + 1757);
    const auto *sni1_1759 = buffer.data(sni1 + 1759);
    const auto *sni1_1760 = buffer.data(sni1 + 1760);
    const auto *sni1_1761 = buffer.data(sni1 + 1761);
    const auto *sni1_1763 = buffer.data(sni1 + 1763);
    const auto *sni1_1764 = buffer.data(sni1 + 1764);
    const auto *sni1_1767 = buffer.data(sni1 + 1767);
    const auto *sni1_1769 = buffer.data(sni1 + 1769);
    const auto *sni1_1770 = buffer.data(sni1 + 1770);
    const auto *sni1_1773 = buffer.data(sni1 + 1773);
    const auto *sni1_1774 = buffer.data(sni1 + 1774);
    const auto *sni1_1776 = buffer.data(sni1 + 1776);
    const auto *sni1_1778 = buffer.data(sni1 + 1778);
    const auto *sni1_1785 = buffer.data(sni1 + 1785);
    const auto *sni1_1787 = buffer.data(sni1 + 1787);
    const auto *sni1_1788 = buffer.data(sni1 + 1788);
    const auto *sni1_1789 = buffer.data(sni1 + 1789);
    const auto *sni1_1791 = buffer.data(sni1 + 1791);
    const auto *sni1_1795 = buffer.data(sni1 + 1795);
    const auto *sni1_1798 = buffer.data(sni1 + 1798);
    const auto *sni1_1802 = buffer.data(sni1 + 1802);
    const auto *sni1_1804 = buffer.data(sni1 + 1804);
    const auto *sni1_1813 = buffer.data(sni1 + 1813);

    const auto *soh_1275 = buffer.data(soh + 1275);
    const auto *soh_1276 = buffer.data(soh + 1276);
    const auto *soh_1277 = buffer.data(soh + 1277);
    const auto *soh_1278 = buffer.data(soh + 1278);
    const auto *soh_1279 = buffer.data(soh + 1279);
    const auto *soh_1280 = buffer.data(soh + 1280);
    const auto *soh_1281 = buffer.data(soh + 1281);
    const auto *soh_1283 = buffer.data(soh + 1283);
    const auto *soh_1284 = buffer.data(soh + 1284);
    const auto *soh_1286 = buffer.data(soh + 1286);
    const auto *soh_1287 = buffer.data(soh + 1287);
    const auto *soh_1290 = buffer.data(soh + 1290);
    const auto *soh_1296 = buffer.data(soh + 1296);
    const auto *soh_1297 = buffer.data(soh + 1297);
    const auto *soh_1298 = buffer.data(soh + 1298);
    const auto *soh_1299 = buffer.data(soh + 1299);
    const auto *soh_1300 = buffer.data(soh + 1300);
    const auto *soh_1301 = buffer.data(soh + 1301);
    const auto *soh_1302 = buffer.data(soh + 1302);
    const auto *soh_1304 = buffer.data(soh + 1304);
    const auto *soh_1305 = buffer.data(soh + 1305);
    const auto *soh_1307 = buffer.data(soh + 1307);
    const auto *soh_1308 = buffer.data(soh + 1308);
    const auto *soh_1311 = buffer.data(soh + 1311);
    const auto *soh_1317 = buffer.data(soh + 1317);
    const auto *soh_1318 = buffer.data(soh + 1318);
    const auto *soh_1319 = buffer.data(soh + 1319);
    const auto *soh_1320 = buffer.data(soh + 1320);
    const auto *soh_1321 = buffer.data(soh + 1321);
    const auto *soh_1322 = buffer.data(soh + 1322);
    const auto *soh_1323 = buffer.data(soh + 1323);
    const auto *soh_1325 = buffer.data(soh + 1325);
    const auto *soh_1326 = buffer.data(soh + 1326);
    const auto *soh_1328 = buffer.data(soh + 1328);
    const auto *soh_1329 = buffer.data(soh + 1329);
    const auto *soh_1332 = buffer.data(soh + 1332);
    const auto *soh_1338 = buffer.data(soh + 1338);
    const auto *soh_1339 = buffer.data(soh + 1339);
    const auto *soh_1340 = buffer.data(soh + 1340);
    const auto *soh_1341 = buffer.data(soh + 1341);
    const auto *soh_1342 = buffer.data(soh + 1342);
    const auto *soh_1343 = buffer.data(soh + 1343);
    const auto *soh_1344 = buffer.data(soh + 1344);
    const auto *soh_1346 = buffer.data(soh + 1346);
    const auto *soh_1347 = buffer.data(soh + 1347);
    const auto *soh_1349 = buffer.data(soh + 1349);
    const auto *soh_1350 = buffer.data(soh + 1350);
    const auto *soh_1353 = buffer.data(soh + 1353);
    const auto *soh_1359 = buffer.data(soh + 1359);
    const auto *soh_1360 = buffer.data(soh + 1360);
    const auto *soh_1361 = buffer.data(soh + 1361);
    const auto *soh_1362 = buffer.data(soh + 1362);
    const auto *soh_1363 = buffer.data(soh + 1363);
    const auto *soh_1364 = buffer.data(soh + 1364);

#pragma omp simd aligned(t_1694, t_1695, t_1696, t_1697, pb_x, pc_x, sni0_1694, snh_1274, \
                         snh_1275, snh_1276, snh_1277, sni1_1694, soh_1275, soh_1276, \
                         soh_1277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1694[k] = pb_x[k] * sni0_1694[k]
                    + f_12 * snh_1274[k]
                    - f_10 * pc_x[k] * sni1_1694[k];

        t_1695[k] = f_11 * snh_1275[k]
                    + f_3 * pc_x[k] * soh_1275[k];

        t_1696[k] = f_11 * snh_1276[k]
                    + f_3 * pc_x[k] * soh_1276[k];

        t_1697[k] = f_11 * snh_1277[k]
                    + f_3 * pc_x[k] * soh_1277[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, t_1701, pb_x, pc_x, sni0_1701, snh_1278, \
                         snh_1279, snh_1280, sni1_1701, soh_1278, soh_1279, \
                         soh_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_11 * snh_1278[k]
                    + f_3 * pc_x[k] * soh_1278[k];

        t_1699[k] = f_11 * snh_1279[k]
                    + f_3 * pc_x[k] * soh_1279[k];

        t_1700[k] = f_11 * snh_1280[k]
                    + f_3 * pc_x[k] * soh_1280[k];

        t_1701[k] = pb_x[k] * sni0_1701[k]
                    - f_10 * pc_x[k] * sni1_1701[k];
    }

#pragma omp simd aligned(t_1702, t_1703, t_1704, t_1705, pb_x, pc_x, pc_z, sni0_1703, \
                         sni0_1704, sni0_1705, snh_1044, sni1_1703, sni1_1704, sni1_1705, \
                         soh_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1702[k] = f_20 * snh_1044[k]
                    + f_3 * pc_z[k] * soh_1275[k];

        t_1703[k] = pb_x[k] * sni0_1703[k]
                    - f_10 * pc_x[k] * sni1_1703[k];

        t_1704[k] = pb_x[k] * sni0_1704[k]
                    - f_10 * pc_x[k] * sni1_1704[k];

        t_1705[k] = pb_x[k] * sni0_1705[k]
                    - f_10 * pc_x[k] * sni1_1705[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, t_1709, pb_x, pc_x, pc_y, sni0_1707, \
                         sni0_1708, snh_1070, snh_1071, snh_1281, sni1_1707, sni1_1708, \
                         soh_1280, soh_1281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = f_20 * snh_1070[k]
                    + f_3 * pc_y[k] * soh_1280[k];

        t_1707[k] = pb_x[k] * sni0_1707[k]
                    - f_10 * pc_x[k] * sni1_1707[k];

        t_1708[k] = pb_x[k] * sni0_1708[k]
                    + f_19 * snh_1281[k]
                    - f_10 * pc_x[k] * sni1_1708[k];

        t_1709[k] = f_14 * snh_1071[k]
                    + f_3 * pc_y[k] * soh_1281[k];
    }

#pragma omp simd aligned(t_1710, t_1711, t_1712, pb_x, pc_x, pc_y, pc_z, sni0_1711, snh_1050, \
                         snh_1073, snh_1284, sni1_1711, soh_1281, \
                         soh_1283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1710[k] = f_19 * snh_1050[k]
                    + f_3 * pc_z[k] * soh_1281[k];

        t_1711[k] = pb_x[k] * sni0_1711[k]
                    + f_14 * snh_1284[k]
                    - f_10 * pc_x[k] * sni1_1711[k];

        t_1712[k] = f_14 * snh_1073[k]
                    + f_3 * pc_y[k] * soh_1283[k];
    }

#pragma omp simd aligned(t_1713, t_1714, t_1715, pb_x, pc_x, pc_z, sni0_1713, sni0_1714, \
                         snh_1053, snh_1286, snh_1287, sni1_1713, sni1_1714, \
                         soh_1284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1713[k] = pb_x[k] * sni0_1713[k]
                    + f_14 * snh_1286[k]
                    - f_10 * pc_x[k] * sni1_1713[k];

        t_1714[k] = pb_x[k] * sni0_1714[k]
                    + f_13 * snh_1287[k]
                    - f_10 * pc_x[k] * sni1_1714[k];

        t_1715[k] = f_19 * snh_1053[k]
                    + f_3 * pc_z[k] * soh_1284[k];
    }

#pragma omp simd aligned(t_1716, t_1717, t_1718, pb_x, pc_x, pc_y, sni0_1717, sni0_1718, \
                         snh_1076, snh_1290, snh_1291, sni1_1717, sni1_1718, \
                         soh_1286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1716[k] = f_14 * snh_1076[k]
                    + f_3 * pc_y[k] * soh_1286[k];

        t_1717[k] = pb_x[k] * sni0_1717[k]
                    + f_13 * snh_1290[k]
                    - f_10 * pc_x[k] * sni1_1717[k];

        t_1718[k] = pb_x[k] * sni0_1718[k]
                    + f_12 * snh_1291[k]
                    - f_10 * pc_x[k] * sni1_1718[k];
    }

#pragma omp simd aligned(t_1719, t_1720, t_1721, pb_x, pc_x, pc_y, pc_z, sni0_1720, snh_1056, \
                         snh_1080, snh_1293, sni1_1720, soh_1287, \
                         soh_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1719[k] = f_19 * snh_1056[k]
                    + f_3 * pc_z[k] * soh_1287[k];

        t_1720[k] = pb_x[k] * sni0_1720[k]
                    + f_12 * snh_1293[k]
                    - f_10 * pc_x[k] * sni1_1720[k];

        t_1721[k] = f_14 * snh_1080[k]
                    + f_3 * pc_y[k] * soh_1290[k];
    }

#pragma omp simd aligned(t_1722, t_1723, t_1724, t_1725, pb_x, pc_x, sni0_1722, snh_1295, \
                         snh_1296, snh_1297, snh_1298, sni1_1722, soh_1296, soh_1297, \
                         soh_1298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1722[k] = pb_x[k] * sni0_1722[k]
                    + f_12 * snh_1295[k]
                    - f_10 * pc_x[k] * sni1_1722[k];

        t_1723[k] = f_11 * snh_1296[k]
                    + f_3 * pc_x[k] * soh_1296[k];

        t_1724[k] = f_11 * snh_1297[k]
                    + f_3 * pc_x[k] * soh_1297[k];

        t_1725[k] = f_11 * snh_1298[k]
                    + f_3 * pc_x[k] * soh_1298[k];
    }

#pragma omp simd aligned(t_1726, t_1727, t_1728, t_1729, pb_x, pc_x, sni0_1729, snh_1299, \
                         snh_1300, snh_1301, sni1_1729, soh_1299, soh_1300, \
                         soh_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1726[k] = f_11 * snh_1299[k]
                    + f_3 * pc_x[k] * soh_1299[k];

        t_1727[k] = f_11 * snh_1300[k]
                    + f_3 * pc_x[k] * soh_1300[k];

        t_1728[k] = f_11 * snh_1301[k]
                    + f_3 * pc_x[k] * soh_1301[k];

        t_1729[k] = pb_x[k] * sni0_1729[k]
                    - f_10 * pc_x[k] * sni1_1729[k];
    }

#pragma omp simd aligned(t_1730, t_1731, t_1732, t_1733, pb_x, pc_x, pc_z, sni0_1731, \
                         sni0_1732, sni0_1733, snh_1065, sni1_1731, sni1_1732, sni1_1733, \
                         soh_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1730[k] = f_19 * snh_1065[k]
                    + f_3 * pc_z[k] * soh_1296[k];

        t_1731[k] = pb_x[k] * sni0_1731[k]
                    - f_10 * pc_x[k] * sni1_1731[k];

        t_1732[k] = pb_x[k] * sni0_1732[k]
                    - f_10 * pc_x[k] * sni1_1732[k];

        t_1733[k] = pb_x[k] * sni0_1733[k]
                    - f_10 * pc_x[k] * sni1_1733[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, t_1737, pb_x, pc_x, pc_y, sni0_1735, \
                         sni0_1736, snh_1091, snh_1092, snh_1302, sni1_1735, sni1_1736, \
                         soh_1301, soh_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = f_14 * snh_1091[k]
                    + f_3 * pc_y[k] * soh_1301[k];

        t_1735[k] = pb_x[k] * sni0_1735[k]
                    - f_10 * pc_x[k] * sni1_1735[k];

        t_1736[k] = pb_x[k] * sni0_1736[k]
                    + f_19 * snh_1302[k]
                    - f_10 * pc_x[k] * sni1_1736[k];

        t_1737[k] = f_13 * snh_1092[k]
                    + f_3 * pc_y[k] * soh_1302[k];
    }

#pragma omp simd aligned(t_1738, t_1739, t_1740, pb_x, pc_x, pc_y, pc_z, sni0_1739, snh_1071, \
                         snh_1094, snh_1305, sni1_1739, soh_1302, \
                         soh_1304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1738[k] = f_18 * snh_1071[k]
                    + f_3 * pc_z[k] * soh_1302[k];

        t_1739[k] = pb_x[k] * sni0_1739[k]
                    + f_14 * snh_1305[k]
                    - f_10 * pc_x[k] * sni1_1739[k];

        t_1740[k] = f_13 * snh_1094[k]
                    + f_3 * pc_y[k] * soh_1304[k];
    }

#pragma omp simd aligned(t_1741, t_1742, t_1743, pb_x, pc_x, pc_z, sni0_1741, sni0_1742, \
                         snh_1074, snh_1307, snh_1308, sni1_1741, sni1_1742, \
                         soh_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1741[k] = pb_x[k] * sni0_1741[k]
                    + f_14 * snh_1307[k]
                    - f_10 * pc_x[k] * sni1_1741[k];

        t_1742[k] = pb_x[k] * sni0_1742[k]
                    + f_13 * snh_1308[k]
                    - f_10 * pc_x[k] * sni1_1742[k];

        t_1743[k] = f_18 * snh_1074[k]
                    + f_3 * pc_z[k] * soh_1305[k];
    }

#pragma omp simd aligned(t_1744, t_1745, t_1746, pb_x, pc_x, pc_y, sni0_1745, sni0_1746, \
                         snh_1097, snh_1311, snh_1312, sni1_1745, sni1_1746, \
                         soh_1307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1744[k] = f_13 * snh_1097[k]
                    + f_3 * pc_y[k] * soh_1307[k];

        t_1745[k] = pb_x[k] * sni0_1745[k]
                    + f_13 * snh_1311[k]
                    - f_10 * pc_x[k] * sni1_1745[k];

        t_1746[k] = pb_x[k] * sni0_1746[k]
                    + f_12 * snh_1312[k]
                    - f_10 * pc_x[k] * sni1_1746[k];
    }

#pragma omp simd aligned(t_1747, t_1748, t_1749, pb_x, pc_x, pc_y, pc_z, sni0_1748, snh_1077, \
                         snh_1101, snh_1314, sni1_1748, soh_1308, \
                         soh_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1747[k] = f_18 * snh_1077[k]
                    + f_3 * pc_z[k] * soh_1308[k];

        t_1748[k] = pb_x[k] * sni0_1748[k]
                    + f_12 * snh_1314[k]
                    - f_10 * pc_x[k] * sni1_1748[k];

        t_1749[k] = f_13 * snh_1101[k]
                    + f_3 * pc_y[k] * soh_1311[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, pb_x, pc_x, sni0_1750, snh_1316, \
                         snh_1317, snh_1318, snh_1319, sni1_1750, soh_1317, soh_1318, \
                         soh_1319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = pb_x[k] * sni0_1750[k]
                    + f_12 * snh_1316[k]
                    - f_10 * pc_x[k] * sni1_1750[k];

        t_1751[k] = f_11 * snh_1317[k]
                    + f_3 * pc_x[k] * soh_1317[k];

        t_1752[k] = f_11 * snh_1318[k]
                    + f_3 * pc_x[k] * soh_1318[k];

        t_1753[k] = f_11 * snh_1319[k]
                    + f_3 * pc_x[k] * soh_1319[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, pb_x, pc_x, sni0_1757, snh_1320, \
                         snh_1321, snh_1322, sni1_1757, soh_1320, soh_1321, \
                         soh_1322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = f_11 * snh_1320[k]
                    + f_3 * pc_x[k] * soh_1320[k];

        t_1755[k] = f_11 * snh_1321[k]
                    + f_3 * pc_x[k] * soh_1321[k];

        t_1756[k] = f_11 * snh_1322[k]
                    + f_3 * pc_x[k] * soh_1322[k];

        t_1757[k] = pb_x[k] * sni0_1757[k]
                    - f_10 * pc_x[k] * sni1_1757[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, t_1761, pb_x, pc_x, pc_z, sni0_1759, \
                         sni0_1760, sni0_1761, snh_1086, sni1_1759, sni1_1760, sni1_1761, \
                         soh_1317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = f_18 * snh_1086[k]
                    + f_3 * pc_z[k] * soh_1317[k];

        t_1759[k] = pb_x[k] * sni0_1759[k]
                    - f_10 * pc_x[k] * sni1_1759[k];

        t_1760[k] = pb_x[k] * sni0_1760[k]
                    - f_10 * pc_x[k] * sni1_1760[k];

        t_1761[k] = pb_x[k] * sni0_1761[k]
                    - f_10 * pc_x[k] * sni1_1761[k];
    }

#pragma omp simd aligned(t_1762, t_1763, t_1764, t_1765, pb_x, pc_x, pc_y, sni0_1763, \
                         sni0_1764, snh_1112, snh_1113, snh_1323, sni1_1763, sni1_1764, \
                         soh_1322, soh_1323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1762[k] = f_13 * snh_1112[k]
                    + f_3 * pc_y[k] * soh_1322[k];

        t_1763[k] = pb_x[k] * sni0_1763[k]
                    - f_10 * pc_x[k] * sni1_1763[k];

        t_1764[k] = pb_x[k] * sni0_1764[k]
                    + f_19 * snh_1323[k]
                    - f_10 * pc_x[k] * sni1_1764[k];

        t_1765[k] = f_12 * snh_1113[k]
                    + f_3 * pc_y[k] * soh_1323[k];
    }

#pragma omp simd aligned(t_1766, t_1767, t_1768, pb_x, pc_x, pc_y, pc_z, sni0_1767, snh_1092, \
                         snh_1115, snh_1326, sni1_1767, soh_1323, \
                         soh_1325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1766[k] = f_17 * snh_1092[k]
                    + f_3 * pc_z[k] * soh_1323[k];

        t_1767[k] = pb_x[k] * sni0_1767[k]
                    + f_14 * snh_1326[k]
                    - f_10 * pc_x[k] * sni1_1767[k];

        t_1768[k] = f_12 * snh_1115[k]
                    + f_3 * pc_y[k] * soh_1325[k];
    }

#pragma omp simd aligned(t_1769, t_1770, t_1771, pb_x, pc_x, pc_z, sni0_1769, sni0_1770, \
                         snh_1095, snh_1328, snh_1329, sni1_1769, sni1_1770, \
                         soh_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1769[k] = pb_x[k] * sni0_1769[k]
                    + f_14 * snh_1328[k]
                    - f_10 * pc_x[k] * sni1_1769[k];

        t_1770[k] = pb_x[k] * sni0_1770[k]
                    + f_13 * snh_1329[k]
                    - f_10 * pc_x[k] * sni1_1770[k];

        t_1771[k] = f_17 * snh_1095[k]
                    + f_3 * pc_z[k] * soh_1326[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, pb_x, pc_x, pc_y, sni0_1773, sni0_1774, \
                         snh_1118, snh_1332, snh_1333, sni1_1773, sni1_1774, \
                         soh_1328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = f_12 * snh_1118[k]
                    + f_3 * pc_y[k] * soh_1328[k];

        t_1773[k] = pb_x[k] * sni0_1773[k]
                    + f_13 * snh_1332[k]
                    - f_10 * pc_x[k] * sni1_1773[k];

        t_1774[k] = pb_x[k] * sni0_1774[k]
                    + f_12 * snh_1333[k]
                    - f_10 * pc_x[k] * sni1_1774[k];
    }

#pragma omp simd aligned(t_1775, t_1776, t_1777, pb_x, pc_x, pc_y, pc_z, sni0_1776, snh_1098, \
                         snh_1122, snh_1335, sni1_1776, soh_1329, \
                         soh_1332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1775[k] = f_17 * snh_1098[k]
                    + f_3 * pc_z[k] * soh_1329[k];

        t_1776[k] = pb_x[k] * sni0_1776[k]
                    + f_12 * snh_1335[k]
                    - f_10 * pc_x[k] * sni1_1776[k];

        t_1777[k] = f_12 * snh_1122[k]
                    + f_3 * pc_y[k] * soh_1332[k];
    }

#pragma omp simd aligned(t_1778, t_1779, t_1780, t_1781, pb_x, pc_x, sni0_1778, snh_1337, \
                         snh_1338, snh_1339, snh_1340, sni1_1778, soh_1338, soh_1339, \
                         soh_1340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1778[k] = pb_x[k] * sni0_1778[k]
                    + f_12 * snh_1337[k]
                    - f_10 * pc_x[k] * sni1_1778[k];

        t_1779[k] = f_11 * snh_1338[k]
                    + f_3 * pc_x[k] * soh_1338[k];

        t_1780[k] = f_11 * snh_1339[k]
                    + f_3 * pc_x[k] * soh_1339[k];

        t_1781[k] = f_11 * snh_1340[k]
                    + f_3 * pc_x[k] * soh_1340[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, pb_x, pc_x, sni0_1785, snh_1341, \
                         snh_1342, snh_1343, sni1_1785, soh_1341, soh_1342, \
                         soh_1343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_11 * snh_1341[k]
                    + f_3 * pc_x[k] * soh_1341[k];

        t_1783[k] = f_11 * snh_1342[k]
                    + f_3 * pc_x[k] * soh_1342[k];

        t_1784[k] = f_11 * snh_1343[k]
                    + f_3 * pc_x[k] * soh_1343[k];

        t_1785[k] = pb_x[k] * sni0_1785[k]
                    - f_10 * pc_x[k] * sni1_1785[k];
    }

#pragma omp simd aligned(t_1786, t_1787, t_1788, t_1789, pb_x, pc_x, pc_z, sni0_1787, \
                         sni0_1788, sni0_1789, snh_1107, sni1_1787, sni1_1788, sni1_1789, \
                         soh_1338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1786[k] = f_17 * snh_1107[k]
                    + f_3 * pc_z[k] * soh_1338[k];

        t_1787[k] = pb_x[k] * sni0_1787[k]
                    - f_10 * pc_x[k] * sni1_1787[k];

        t_1788[k] = pb_x[k] * sni0_1788[k]
                    - f_10 * pc_x[k] * sni1_1788[k];

        t_1789[k] = pb_x[k] * sni0_1789[k]
                    - f_10 * pc_x[k] * sni1_1789[k];
    }

#pragma omp simd aligned(t_1790, t_1791, t_1792, t_1793, pb_x, pb_y, pc_x, pc_y, sni0_1512, \
                         sni0_1791, snh_1133, snh_1134, sni1_1512, sni1_1791, soh_1343, \
                         soh_1344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1790[k] = f_12 * snh_1133[k]
                    + f_3 * pc_y[k] * soh_1343[k];

        t_1791[k] = pb_x[k] * sni0_1791[k]
                    - f_10 * pc_x[k] * sni1_1791[k];

        t_1792[k] = pb_y[k] * sni0_1512[k]
                    - f_10 * pc_y[k] * sni1_1512[k];

        t_1793[k] = f_11 * snh_1134[k]
                    + f_3 * pc_y[k] * soh_1344[k];
    }

#pragma omp simd aligned(t_1794, t_1795, t_1796, pb_x, pc_x, pc_y, pc_z, sni0_1795, snh_1113, \
                         snh_1136, snh_1347, sni1_1795, soh_1344, \
                         soh_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1794[k] = f_16 * snh_1113[k]
                    + f_3 * pc_z[k] * soh_1344[k];

        t_1795[k] = pb_x[k] * sni0_1795[k]
                    + f_14 * snh_1347[k]
                    - f_10 * pc_x[k] * sni1_1795[k];

        t_1796[k] = f_11 * snh_1136[k]
                    + f_3 * pc_y[k] * soh_1346[k];
    }

#pragma omp simd aligned(t_1797, t_1798, t_1799, pb_x, pb_y, pc_x, pc_y, pc_z, sni0_1517, \
                         sni0_1798, snh_1116, snh_1350, sni1_1517, sni1_1798, \
                         soh_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1797[k] = pb_y[k] * sni0_1517[k]
                    - f_10 * pc_y[k] * sni1_1517[k];

        t_1798[k] = pb_x[k] * sni0_1798[k]
                    + f_13 * snh_1350[k]
                    - f_10 * pc_x[k] * sni1_1798[k];

        t_1799[k] = f_16 * snh_1116[k]
                    + f_3 * pc_z[k] * soh_1347[k];
    }

#pragma omp simd aligned(t_1800, t_1801, t_1802, pb_x, pb_y, pc_x, pc_y, sni0_1521, sni0_1802, \
                         snh_1139, snh_1354, sni1_1521, sni1_1802, \
                         soh_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1800[k] = f_11 * snh_1139[k]
                    + f_3 * pc_y[k] * soh_1349[k];

        t_1801[k] = pb_y[k] * sni0_1521[k]
                    - f_10 * pc_y[k] * sni1_1521[k];

        t_1802[k] = pb_x[k] * sni0_1802[k]
                    + f_12 * snh_1354[k]
                    - f_10 * pc_x[k] * sni1_1802[k];
    }

#pragma omp simd aligned(t_1803, t_1804, t_1805, pb_x, pc_x, pc_y, pc_z, sni0_1804, snh_1119, \
                         snh_1143, snh_1356, sni1_1804, soh_1350, \
                         soh_1353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1803[k] = f_16 * snh_1119[k]
                    + f_3 * pc_z[k] * soh_1350[k];

        t_1804[k] = pb_x[k] * sni0_1804[k]
                    + f_12 * snh_1356[k]
                    - f_10 * pc_x[k] * sni1_1804[k];

        t_1805[k] = f_11 * snh_1143[k]
                    + f_3 * pc_y[k] * soh_1353[k];
    }

#pragma omp simd aligned(t_1806, t_1807, t_1808, t_1809, pb_y, pc_x, pc_y, sni0_1526, \
                         snh_1359, snh_1360, snh_1361, sni1_1526, soh_1359, soh_1360, \
                         soh_1361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1806[k] = pb_y[k] * sni0_1526[k]
                    - f_10 * pc_y[k] * sni1_1526[k];

        t_1807[k] = f_11 * snh_1359[k]
                    + f_3 * pc_x[k] * soh_1359[k];

        t_1808[k] = f_11 * snh_1360[k]
                    + f_3 * pc_x[k] * soh_1360[k];

        t_1809[k] = f_11 * snh_1361[k]
                    + f_3 * pc_x[k] * soh_1361[k];
    }

#pragma omp simd aligned(t_1810, t_1811, t_1812, t_1813, pb_x, pc_x, sni0_1813, snh_1362, \
                         snh_1363, snh_1364, sni1_1813, soh_1362, soh_1363, \
                         soh_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1810[k] = f_11 * snh_1362[k]
                    + f_3 * pc_x[k] * soh_1362[k];

        t_1811[k] = f_11 * snh_1363[k]
                    + f_3 * pc_x[k] * soh_1363[k];

        t_1812[k] = f_11 * snh_1364[k]
                    + f_3 * pc_x[k] * soh_1364[k];

        t_1813[k] = pb_x[k] * sni0_1813[k]
                    - f_10 * pc_x[k] * sni1_1813[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t sog0, const size_t sog1,
                                                           const size_t soh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sni0_1540 = buffer.data(sni0 + 1540);
    const auto *sni0_1543 = buffer.data(sni0 + 1543);
    const auto *sni0_1546 = buffer.data(sni0 + 1546);
    const auto *sni0_1550 = buffer.data(sni0 + 1550);
    const auto *sni0_1561 = buffer.data(sni0 + 1561);
    const auto *sni0_1563 = buffer.data(sni0 + 1563);
    const auto *sni0_1564 = buffer.data(sni0 + 1564);
    const auto *sni0_1565 = buffer.data(sni0 + 1565);
    const auto *sni0_1815 = buffer.data(sni0 + 1815);
    const auto *sni0_1816 = buffer.data(sni0 + 1816);
    const auto *sni0_1817 = buffer.data(sni0 + 1817);
    const auto *sni0_1819 = buffer.data(sni0 + 1819);
    const auto *sni0_1820 = buffer.data(sni0 + 1820);
    const auto *sni0_1823 = buffer.data(sni0 + 1823);
    const auto *sni0_1825 = buffer.data(sni0 + 1825);
    const auto *sni0_1826 = buffer.data(sni0 + 1826);
    const auto *sni0_1829 = buffer.data(sni0 + 1829);
    const auto *sni0_1830 = buffer.data(sni0 + 1830);
    const auto *sni0_1832 = buffer.data(sni0 + 1832);
    const auto *sni0_1834 = buffer.data(sni0 + 1834);
    const auto *sni0_1841 = buffer.data(sni0 + 1841);
    const auto *sni0_1843 = buffer.data(sni0 + 1843);
    const auto *sni0_1844 = buffer.data(sni0 + 1844);
    const auto *sni0_1845 = buffer.data(sni0 + 1845);
    const auto *sni0_1847 = buffer.data(sni0 + 1847);

    const auto *snh_1128 = buffer.data(snh + 1128);
    const auto *snh_1134 = buffer.data(snh + 1134);
    const auto *snh_1137 = buffer.data(snh + 1137);
    const auto *snh_1140 = buffer.data(snh + 1140);
    const auto *snh_1149 = buffer.data(snh + 1149);
    const auto *snh_1154 = buffer.data(snh + 1154);
    const auto *snh_1155 = buffer.data(snh + 1155);
    const auto *snh_1157 = buffer.data(snh + 1157);
    const auto *snh_1158 = buffer.data(snh + 1158);
    const auto *snh_1160 = buffer.data(snh + 1160);
    const auto *snh_1161 = buffer.data(snh + 1161);
    const auto *snh_1164 = buffer.data(snh + 1164);
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
    const auto *snh_1191 = buffer.data(snh + 1191);
    const auto *snh_1196 = buffer.data(snh + 1196);
    const auto *snh_1197 = buffer.data(snh + 1197);
    const auto *snh_1199 = buffer.data(snh + 1199);
    const auto *snh_1200 = buffer.data(snh + 1200);
    const auto *snh_1202 = buffer.data(snh + 1202);
    const auto *snh_1206 = buffer.data(snh + 1206);
    const auto *snh_1212 = buffer.data(snh + 1212);
    const auto *snh_1214 = buffer.data(snh + 1214);
    const auto *snh_1215 = buffer.data(snh + 1215);
    const auto *snh_1216 = buffer.data(snh + 1216);
    const auto *snh_1217 = buffer.data(snh + 1217);
    const auto *snh_1218 = buffer.data(snh + 1218);
    const auto *snh_1220 = buffer.data(snh + 1220);
    const auto *snh_1223 = buffer.data(snh + 1223);
    const auto *snh_1365 = buffer.data(snh + 1365);
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

    const auto *sni1_1540 = buffer.data(sni1 + 1540);
    const auto *sni1_1543 = buffer.data(sni1 + 1543);
    const auto *sni1_1546 = buffer.data(sni1 + 1546);
    const auto *sni1_1550 = buffer.data(sni1 + 1550);
    const auto *sni1_1561 = buffer.data(sni1 + 1561);
    const auto *sni1_1563 = buffer.data(sni1 + 1563);
    const auto *sni1_1564 = buffer.data(sni1 + 1564);
    const auto *sni1_1565 = buffer.data(sni1 + 1565);
    const auto *sni1_1815 = buffer.data(sni1 + 1815);
    const auto *sni1_1816 = buffer.data(sni1 + 1816);
    const auto *sni1_1817 = buffer.data(sni1 + 1817);
    const auto *sni1_1819 = buffer.data(sni1 + 1819);
    const auto *sni1_1820 = buffer.data(sni1 + 1820);
    const auto *sni1_1823 = buffer.data(sni1 + 1823);
    const auto *sni1_1825 = buffer.data(sni1 + 1825);
    const auto *sni1_1826 = buffer.data(sni1 + 1826);
    const auto *sni1_1829 = buffer.data(sni1 + 1829);
    const auto *sni1_1830 = buffer.data(sni1 + 1830);
    const auto *sni1_1832 = buffer.data(sni1 + 1832);
    const auto *sni1_1834 = buffer.data(sni1 + 1834);
    const auto *sni1_1841 = buffer.data(sni1 + 1841);
    const auto *sni1_1843 = buffer.data(sni1 + 1843);
    const auto *sni1_1844 = buffer.data(sni1 + 1844);
    const auto *sni1_1845 = buffer.data(sni1 + 1845);
    const auto *sni1_1847 = buffer.data(sni1 + 1847);

    const auto *sog0_990 = buffer.data(sog0 + 990);
    const auto *sog0_993 = buffer.data(sog0 + 993);
    const auto *sog0_995 = buffer.data(sog0 + 995);
    const auto *sog0_996 = buffer.data(sog0 + 996);
    const auto *sog0_999 = buffer.data(sog0 + 999);
    const auto *sog0_1000 = buffer.data(sog0 + 1000);
    const auto *sog0_1002 = buffer.data(sog0 + 1002);
    const auto *sog0_1003 = buffer.data(sog0 + 1003);
    const auto *sog0_1004 = buffer.data(sog0 + 1004);
    const auto *sog0_1010 = buffer.data(sog0 + 1010);
    const auto *sog0_1014 = buffer.data(sog0 + 1014);
    const auto *sog0_1017 = buffer.data(sog0 + 1017);
    const auto *sog0_1019 = buffer.data(sog0 + 1019);
    const auto *sog0_1020 = buffer.data(sog0 + 1020);
    const auto *sog0_1023 = buffer.data(sog0 + 1023);
    const auto *sog0_1025 = buffer.data(sog0 + 1025);
    const auto *sog0_1026 = buffer.data(sog0 + 1026);
    const auto *sog0_1029 = buffer.data(sog0 + 1029);
    const auto *sog0_1030 = buffer.data(sog0 + 1030);
    const auto *sog0_1032 = buffer.data(sog0 + 1032);
    const auto *sog0_1033 = buffer.data(sog0 + 1033);
    const auto *sog0_1034 = buffer.data(sog0 + 1034);
    const auto *sog0_1035 = buffer.data(sog0 + 1035);
    const auto *sog0_1038 = buffer.data(sog0 + 1038);
    const auto *sog0_1040 = buffer.data(sog0 + 1040);
    const auto *sog0_1041 = buffer.data(sog0 + 1041);

    const auto *sog1_990 = buffer.data(sog1 + 990);
    const auto *sog1_993 = buffer.data(sog1 + 993);
    const auto *sog1_995 = buffer.data(sog1 + 995);
    const auto *sog1_996 = buffer.data(sog1 + 996);
    const auto *sog1_999 = buffer.data(sog1 + 999);
    const auto *sog1_1000 = buffer.data(sog1 + 1000);
    const auto *sog1_1002 = buffer.data(sog1 + 1002);
    const auto *sog1_1003 = buffer.data(sog1 + 1003);
    const auto *sog1_1004 = buffer.data(sog1 + 1004);
    const auto *sog1_1010 = buffer.data(sog1 + 1010);
    const auto *sog1_1014 = buffer.data(sog1 + 1014);
    const auto *sog1_1017 = buffer.data(sog1 + 1017);
    const auto *sog1_1019 = buffer.data(sog1 + 1019);
    const auto *sog1_1020 = buffer.data(sog1 + 1020);
    const auto *sog1_1023 = buffer.data(sog1 + 1023);
    const auto *sog1_1025 = buffer.data(sog1 + 1025);
    const auto *sog1_1026 = buffer.data(sog1 + 1026);
    const auto *sog1_1029 = buffer.data(sog1 + 1029);
    const auto *sog1_1030 = buffer.data(sog1 + 1030);
    const auto *sog1_1032 = buffer.data(sog1 + 1032);
    const auto *sog1_1033 = buffer.data(sog1 + 1033);
    const auto *sog1_1034 = buffer.data(sog1 + 1034);
    const auto *sog1_1035 = buffer.data(sog1 + 1035);
    const auto *sog1_1038 = buffer.data(sog1 + 1038);
    const auto *sog1_1040 = buffer.data(sog1 + 1040);
    const auto *sog1_1041 = buffer.data(sog1 + 1041);

    const auto *soh_1359 = buffer.data(soh + 1359);
    const auto *soh_1364 = buffer.data(soh + 1364);
    const auto *soh_1365 = buffer.data(soh + 1365);
    const auto *soh_1367 = buffer.data(soh + 1367);
    const auto *soh_1368 = buffer.data(soh + 1368);
    const auto *soh_1370 = buffer.data(soh + 1370);
    const auto *soh_1371 = buffer.data(soh + 1371);
    const auto *soh_1374 = buffer.data(soh + 1374);
    const auto *soh_1380 = buffer.data(soh + 1380);
    const auto *soh_1381 = buffer.data(soh + 1381);
    const auto *soh_1382 = buffer.data(soh + 1382);
    const auto *soh_1383 = buffer.data(soh + 1383);
    const auto *soh_1384 = buffer.data(soh + 1384);
    const auto *soh_1385 = buffer.data(soh + 1385);
    const auto *soh_1386 = buffer.data(soh + 1386);
    const auto *soh_1388 = buffer.data(soh + 1388);
    const auto *soh_1389 = buffer.data(soh + 1389);
    const auto *soh_1391 = buffer.data(soh + 1391);
    const auto *soh_1392 = buffer.data(soh + 1392);
    const auto *soh_1395 = buffer.data(soh + 1395);
    const auto *soh_1396 = buffer.data(soh + 1396);
    const auto *soh_1398 = buffer.data(soh + 1398);
    const auto *soh_1400 = buffer.data(soh + 1400);
    const auto *soh_1401 = buffer.data(soh + 1401);
    const auto *soh_1402 = buffer.data(soh + 1402);
    const auto *soh_1403 = buffer.data(soh + 1403);
    const auto *soh_1404 = buffer.data(soh + 1404);
    const auto *soh_1405 = buffer.data(soh + 1405);
    const auto *soh_1406 = buffer.data(soh + 1406);
    const auto *soh_1407 = buffer.data(soh + 1407);
    const auto *soh_1409 = buffer.data(soh + 1409);
    const auto *soh_1410 = buffer.data(soh + 1410);
    const auto *soh_1412 = buffer.data(soh + 1412);
    const auto *soh_1413 = buffer.data(soh + 1413);
    const auto *soh_1416 = buffer.data(soh + 1416);
    const auto *soh_1419 = buffer.data(soh + 1419);
    const auto *soh_1421 = buffer.data(soh + 1421);
    const auto *soh_1422 = buffer.data(soh + 1422);
    const auto *soh_1423 = buffer.data(soh + 1423);
    const auto *soh_1424 = buffer.data(soh + 1424);
    const auto *soh_1425 = buffer.data(soh + 1425);
    const auto *soh_1426 = buffer.data(soh + 1426);
    const auto *soh_1427 = buffer.data(soh + 1427);
    const auto *soh_1428 = buffer.data(soh + 1428);
    const auto *soh_1430 = buffer.data(soh + 1430);
    const auto *soh_1431 = buffer.data(soh + 1431);
    const auto *soh_1433 = buffer.data(soh + 1433);
    const auto *soh_1434 = buffer.data(soh + 1434);
    const auto *soh_1437 = buffer.data(soh + 1437);
    const auto *soh_1438 = buffer.data(soh + 1438);
    const auto *soh_1440 = buffer.data(soh + 1440);
    const auto *soh_1442 = buffer.data(soh + 1442);
    const auto *soh_1443 = buffer.data(soh + 1443);
    const auto *soh_1444 = buffer.data(soh + 1444);
    const auto *soh_1445 = buffer.data(soh + 1445);
    const auto *soh_1446 = buffer.data(soh + 1446);
    const auto *soh_1447 = buffer.data(soh + 1447);
    const auto *soh_1448 = buffer.data(soh + 1448);
    const auto *soh_1449 = buffer.data(soh + 1449);
    const auto *soh_1451 = buffer.data(soh + 1451);
    const auto *soh_1452 = buffer.data(soh + 1452);
    const auto *soh_1454 = buffer.data(soh + 1454);
    const auto *soh_1455 = buffer.data(soh + 1455);

#pragma omp simd aligned(t_1814, t_1815, t_1816, t_1817, pb_x, pc_x, pc_z, sni0_1815, \
                         sni0_1816, sni0_1817, snh_1128, sni1_1815, sni1_1816, sni1_1817, \
                         soh_1359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1814[k] = f_16 * snh_1128[k]
                    + f_3 * pc_z[k] * soh_1359[k];

        t_1815[k] = pb_x[k] * sni0_1815[k]
                    - f_10 * pc_x[k] * sni1_1815[k];

        t_1816[k] = pb_x[k] * sni0_1816[k]
                    - f_10 * pc_x[k] * sni1_1816[k];

        t_1817[k] = pb_x[k] * sni0_1817[k]
                    - f_10 * pc_x[k] * sni1_1817[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, t_1821, pb_x, pc_x, pc_y, sni0_1819, \
                         sni0_1820, snh_1154, snh_1365, sni1_1819, sni1_1820, soh_1364, \
                         soh_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = f_11 * snh_1154[k]
                    + f_3 * pc_y[k] * soh_1364[k];

        t_1819[k] = pb_x[k] * sni0_1819[k]
                    - f_10 * pc_x[k] * sni1_1819[k];

        t_1820[k] = pb_x[k] * sni0_1820[k]
                    + f_19 * snh_1365[k]
                    - f_10 * pc_x[k] * sni1_1820[k];

        t_1821[k] = f_3 * pc_y[k] * soh_1365[k];
    }

#pragma omp simd aligned(t_1822, t_1823, t_1824, pb_x, pc_x, pc_y, pc_z, sni0_1823, snh_1134, \
                         snh_1368, sni1_1823, soh_1365, soh_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1822[k] = f_15 * snh_1134[k]
                    + f_3 * pc_z[k] * soh_1365[k];

        t_1823[k] = pb_x[k] * sni0_1823[k]
                    + f_14 * snh_1368[k]
                    - f_10 * pc_x[k] * sni1_1823[k];

        t_1824[k] = f_3 * pc_y[k] * soh_1367[k];
    }

#pragma omp simd aligned(t_1825, t_1826, t_1827, pb_x, pc_x, pc_z, sni0_1825, sni0_1826, \
                         snh_1137, snh_1370, snh_1371, sni1_1825, sni1_1826, \
                         soh_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1825[k] = pb_x[k] * sni0_1825[k]
                    + f_14 * snh_1370[k]
                    - f_10 * pc_x[k] * sni1_1825[k];

        t_1826[k] = pb_x[k] * sni0_1826[k]
                    + f_13 * snh_1371[k]
                    - f_10 * pc_x[k] * sni1_1826[k];

        t_1827[k] = f_15 * snh_1137[k]
                    + f_3 * pc_z[k] * soh_1368[k];
    }

#pragma omp simd aligned(t_1828, t_1829, t_1830, pb_x, pc_x, pc_y, sni0_1829, sni0_1830, \
                         snh_1374, snh_1375, sni1_1829, sni1_1830, \
                         soh_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1828[k] = f_3 * pc_y[k] * soh_1370[k];

        t_1829[k] = pb_x[k] * sni0_1829[k]
                    + f_13 * snh_1374[k]
                    - f_10 * pc_x[k] * sni1_1829[k];

        t_1830[k] = pb_x[k] * sni0_1830[k]
                    + f_12 * snh_1375[k]
                    - f_10 * pc_x[k] * sni1_1830[k];
    }

#pragma omp simd aligned(t_1831, t_1832, t_1833, pb_x, pc_x, pc_y, pc_z, sni0_1832, snh_1140, \
                         snh_1377, sni1_1832, soh_1371, soh_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1831[k] = f_15 * snh_1140[k]
                    + f_3 * pc_z[k] * soh_1371[k];

        t_1832[k] = pb_x[k] * sni0_1832[k]
                    + f_12 * snh_1377[k]
                    - f_10 * pc_x[k] * sni1_1832[k];

        t_1833[k] = f_3 * pc_y[k] * soh_1374[k];
    }

#pragma omp simd aligned(t_1834, t_1835, t_1836, t_1837, pb_x, pc_x, sni0_1834, snh_1379, \
                         snh_1380, snh_1381, snh_1382, sni1_1834, soh_1380, soh_1381, \
                         soh_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = pb_x[k] * sni0_1834[k]
                    + f_12 * snh_1379[k]
                    - f_10 * pc_x[k] * sni1_1834[k];

        t_1835[k] = f_11 * snh_1380[k]
                    + f_3 * pc_x[k] * soh_1380[k];

        t_1836[k] = f_11 * snh_1381[k]
                    + f_3 * pc_x[k] * soh_1381[k];

        t_1837[k] = f_11 * snh_1382[k]
                    + f_3 * pc_x[k] * soh_1382[k];
    }

#pragma omp simd aligned(t_1838, t_1839, t_1840, t_1841, pb_x, pc_x, sni0_1841, snh_1383, \
                         snh_1384, snh_1385, sni1_1841, soh_1383, soh_1384, \
                         soh_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1838[k] = f_11 * snh_1383[k]
                    + f_3 * pc_x[k] * soh_1383[k];

        t_1839[k] = f_11 * snh_1384[k]
                    + f_3 * pc_x[k] * soh_1384[k];

        t_1840[k] = f_11 * snh_1385[k]
                    + f_3 * pc_x[k] * soh_1385[k];

        t_1841[k] = pb_x[k] * sni0_1841[k]
                    - f_10 * pc_x[k] * sni1_1841[k];
    }

#pragma omp simd aligned(t_1842, t_1843, t_1844, t_1845, pb_x, pc_x, pc_z, sni0_1843, \
                         sni0_1844, sni0_1845, snh_1149, sni1_1843, sni1_1844, sni1_1845, \
                         soh_1380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1842[k] = f_15 * snh_1149[k]
                    + f_3 * pc_z[k] * soh_1380[k];

        t_1843[k] = pb_x[k] * sni0_1843[k]
                    - f_10 * pc_x[k] * sni1_1843[k];

        t_1844[k] = pb_x[k] * sni0_1844[k]
                    - f_10 * pc_x[k] * sni1_1844[k];

        t_1845[k] = pb_x[k] * sni0_1845[k]
                    - f_10 * pc_x[k] * sni1_1845[k];
    }

#pragma omp simd aligned(t_1846, t_1847, t_1848, t_1849, t_1850, pb_x, pc_x, pc_y, pc_z, \
                         sni0_1847, snh_1155, sni1_1847, sog0_990, sog1_990, soh_1385, \
                         soh_1386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1846[k] = f_3 * pc_y[k] * soh_1385[k];

        t_1847[k] = pb_x[k] * sni0_1847[k]
                    - f_10 * pc_x[k] * sni1_1847[k];

        t_1848[k] = f_1 * sog0_990[k]
                    - f_2 * sog1_990[k]
                    + f_3 * pc_x[k] * soh_1386[k];

        t_1849[k] = f_0 * snh_1155[k]
                    + f_3 * pc_y[k] * soh_1386[k];

        t_1850[k] = f_3 * pc_z[k] * soh_1386[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pc_x, pc_y, snh_1157, sog0_993, sog0_995, \
                         sog1_993, sog1_995, soh_1388, soh_1389, \
                         soh_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = f_4 * sog0_993[k]
                    - f_5 * sog1_993[k]
                    + f_3 * pc_x[k] * soh_1389[k];

        t_1852[k] = f_0 * snh_1157[k]
                    + f_3 * pc_y[k] * soh_1388[k];

        t_1853[k] = f_4 * sog0_995[k]
                    - f_5 * sog1_995[k]
                    + f_3 * pc_x[k] * soh_1391[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, t_1857, pc_x, pc_y, pc_z, snh_1160, sog0_996, \
                         sog0_999, sog1_996, sog1_999, soh_1389, soh_1391, soh_1392, \
                         soh_1395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_6 * sog0_996[k]
                    - f_7 * sog1_996[k]
                    + f_3 * pc_x[k] * soh_1392[k];

        t_1855[k] = f_3 * pc_z[k] * soh_1389[k];

        t_1856[k] = f_0 * snh_1160[k]
                    + f_3 * pc_y[k] * soh_1391[k];

        t_1857[k] = f_6 * sog0_999[k]
                    - f_7 * sog1_999[k]
                    + f_3 * pc_x[k] * soh_1395[k];
    }

#pragma omp simd aligned(t_1858, t_1859, t_1860, t_1861, pc_x, pc_y, pc_z, snh_1164, \
                         sog0_1000, sog0_1002, sog1_1000, sog1_1002, soh_1392, soh_1395, \
                         soh_1396, soh_1398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1858[k] = f_8 * sog0_1000[k]
                    - f_9 * sog1_1000[k]
                    + f_3 * pc_x[k] * soh_1396[k];

        t_1859[k] = f_3 * pc_z[k] * soh_1392[k];

        t_1860[k] = f_8 * sog0_1002[k]
                    - f_9 * sog1_1002[k]
                    + f_3 * pc_x[k] * soh_1398[k];

        t_1861[k] = f_0 * snh_1164[k]
                    + f_3 * pc_y[k] * soh_1395[k];
    }

#pragma omp simd aligned(t_1862, t_1863, t_1864, t_1865, t_1866, t_1867, pc_x, sog0_1004, \
                         sog1_1004, soh_1400, soh_1401, soh_1402, soh_1403, soh_1404, \
                         soh_1405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1862[k] = f_8 * sog0_1004[k]
                    - f_9 * sog1_1004[k]
                    + f_3 * pc_x[k] * soh_1400[k];

        t_1863[k] = f_3 * pc_x[k] * soh_1401[k];

        t_1864[k] = f_3 * pc_x[k] * soh_1402[k];

        t_1865[k] = f_3 * pc_x[k] * soh_1403[k];

        t_1866[k] = f_3 * pc_x[k] * soh_1404[k];

        t_1867[k] = f_3 * pc_x[k] * soh_1405[k];
    }

#pragma omp simd aligned(t_1868, t_1869, t_1870, t_1871, pc_x, pc_y, pc_z, snh_1170, snh_1172, \
                         sog0_1000, sog0_1002, sog1_1000, sog1_1002, soh_1401, soh_1403, \
                         soh_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1868[k] = f_3 * pc_x[k] * soh_1406[k];

        t_1869[k] = f_0 * snh_1170[k]
                    + f_1 * sog0_1000[k]
                    - f_2 * sog1_1000[k]
                    + f_3 * pc_y[k] * soh_1401[k];

        t_1870[k] = f_3 * pc_z[k] * soh_1401[k];

        t_1871[k] = f_0 * snh_1172[k]
                    + f_4 * sog0_1002[k]
                    - f_5 * sog1_1002[k]
                    + f_3 * pc_y[k] * soh_1403[k];
    }

#pragma omp simd aligned(t_1872, t_1873, t_1874, t_1875, pc_y, pc_z, snh_1173, snh_1174, \
                         snh_1175, sog0_1003, sog0_1004, sog1_1003, sog1_1004, soh_1404, \
                         soh_1405, soh_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1872[k] = f_0 * snh_1173[k]
                    + f_6 * sog0_1003[k]
                    - f_7 * sog1_1003[k]
                    + f_3 * pc_y[k] * soh_1404[k];

        t_1873[k] = f_0 * snh_1174[k]
                    + f_8 * sog0_1004[k]
                    - f_9 * sog1_1004[k]
                    + f_3 * pc_y[k] * soh_1405[k];

        t_1874[k] = f_0 * snh_1175[k]
                    + f_3 * pc_y[k] * soh_1406[k];

        t_1875[k] = f_1 * sog0_1004[k]
                    - f_2 * sog1_1004[k]
                    + f_3 * pc_z[k] * soh_1406[k];
    }

#pragma omp simd aligned(t_1876, t_1877, t_1878, t_1879, pb_z, pc_y, pc_z, sni0_1540, \
                         sni0_1543, snh_1155, snh_1176, sni1_1540, sni1_1543, \
                         soh_1407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1876[k] = pb_z[k] * sni0_1540[k]
                    - f_10 * pc_z[k] * sni1_1540[k];

        t_1877[k] = f_15 * snh_1176[k]
                    + f_3 * pc_y[k] * soh_1407[k];

        t_1878[k] = f_11 * snh_1155[k]
                    + f_3 * pc_z[k] * soh_1407[k];

        t_1879[k] = pb_z[k] * sni0_1543[k]
                    - f_10 * pc_z[k] * sni1_1543[k];
    }

#pragma omp simd aligned(t_1880, t_1881, t_1882, pb_z, pc_x, pc_y, pc_z, sni0_1546, snh_1178, \
                         sni1_1546, sog0_1010, sog1_1010, soh_1409, \
                         soh_1412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1880[k] = f_15 * snh_1178[k]
                    + f_3 * pc_y[k] * soh_1409[k];

        t_1881[k] = f_4 * sog0_1010[k]
                    - f_5 * sog1_1010[k]
                    + f_3 * pc_x[k] * soh_1412[k];

        t_1882[k] = pb_z[k] * sni0_1546[k]
                    - f_10 * pc_z[k] * sni1_1546[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pc_x, pc_y, pc_z, snh_1158, snh_1181, \
                         sog0_1014, sog1_1014, soh_1410, soh_1412, \
                         soh_1416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_11 * snh_1158[k]
                    + f_3 * pc_z[k] * soh_1410[k];

        t_1884[k] = f_15 * snh_1181[k]
                    + f_3 * pc_y[k] * soh_1412[k];

        t_1885[k] = f_6 * sog0_1014[k]
                    - f_7 * sog1_1014[k]
                    + f_3 * pc_x[k] * soh_1416[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pb_z, pc_x, pc_z, sni0_1550, snh_1161, \
                         sni1_1550, sog0_1017, sog1_1017, soh_1413, \
                         soh_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = pb_z[k] * sni0_1550[k]
                    - f_10 * pc_z[k] * sni1_1550[k];

        t_1887[k] = f_11 * snh_1161[k]
                    + f_3 * pc_z[k] * soh_1413[k];

        t_1888[k] = f_8 * sog0_1017[k]
                    - f_9 * sog1_1017[k]
                    + f_3 * pc_x[k] * soh_1419[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, t_1892, t_1893, pc_x, pc_y, snh_1185, \
                         sog0_1019, sog1_1019, soh_1416, soh_1421, soh_1422, soh_1423, \
                         soh_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = f_15 * snh_1185[k]
                    + f_3 * pc_y[k] * soh_1416[k];

        t_1890[k] = f_8 * sog0_1019[k]
                    - f_9 * sog1_1019[k]
                    + f_3 * pc_x[k] * soh_1421[k];

        t_1891[k] = f_3 * pc_x[k] * soh_1422[k];

        t_1892[k] = f_3 * pc_x[k] * soh_1423[k];

        t_1893[k] = f_3 * pc_x[k] * soh_1424[k];
    }

#pragma omp simd aligned(t_1894, t_1895, t_1896, t_1897, t_1898, pb_z, pc_x, pc_z, sni0_1561, \
                         snh_1170, sni1_1561, soh_1422, soh_1425, soh_1426, \
                         soh_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1894[k] = f_3 * pc_x[k] * soh_1425[k];

        t_1895[k] = f_3 * pc_x[k] * soh_1426[k];

        t_1896[k] = f_3 * pc_x[k] * soh_1427[k];

        t_1897[k] = pb_z[k] * sni0_1561[k]
                    - f_10 * pc_z[k] * sni1_1561[k];

        t_1898[k] = f_11 * snh_1170[k]
                    + f_3 * pc_z[k] * soh_1422[k];
    }

#pragma omp simd aligned(t_1899, t_1900, t_1901, pb_z, pc_z, sni0_1563, sni0_1564, sni0_1565, \
                         snh_1171, snh_1172, snh_1173, sni1_1563, sni1_1564, \
                         sni1_1565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1899[k] = pb_z[k] * sni0_1563[k]
                    + f_12 * snh_1171[k]
                    - f_10 * pc_z[k] * sni1_1563[k];

        t_1900[k] = pb_z[k] * sni0_1564[k]
                    + f_13 * snh_1172[k]
                    - f_10 * pc_z[k] * sni1_1564[k];

        t_1901[k] = pb_z[k] * sni0_1565[k]
                    + f_14 * snh_1173[k]
                    - f_10 * pc_z[k] * sni1_1565[k];
    }

#pragma omp simd aligned(t_1902, t_1903, t_1904, t_1905, pc_x, pc_y, pc_z, snh_1175, snh_1196, \
                         snh_1197, sog0_1019, sog0_1020, sog1_1019, sog1_1020, soh_1427, \
                         soh_1428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1902[k] = f_15 * snh_1196[k]
                    + f_3 * pc_y[k] * soh_1427[k];

        t_1903[k] = f_11 * snh_1175[k]
                    + f_1 * sog0_1019[k]
                    - f_2 * sog1_1019[k]
                    + f_3 * pc_z[k] * soh_1427[k];

        t_1904[k] = f_1 * sog0_1020[k]
                    - f_2 * sog1_1020[k]
                    + f_3 * pc_x[k] * soh_1428[k];

        t_1905[k] = f_16 * snh_1197[k]
                    + f_3 * pc_y[k] * soh_1428[k];
    }

#pragma omp simd aligned(t_1906, t_1907, t_1908, pc_x, pc_y, pc_z, snh_1176, snh_1199, \
                         sog0_1023, sog1_1023, soh_1428, soh_1430, \
                         soh_1431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1906[k] = f_12 * snh_1176[k]
                    + f_3 * pc_z[k] * soh_1428[k];

        t_1907[k] = f_4 * sog0_1023[k]
                    - f_5 * sog1_1023[k]
                    + f_3 * pc_x[k] * soh_1431[k];

        t_1908[k] = f_16 * snh_1199[k]
                    + f_3 * pc_y[k] * soh_1430[k];
    }

#pragma omp simd aligned(t_1909, t_1910, t_1911, t_1912, pc_x, pc_y, pc_z, snh_1179, snh_1202, \
                         sog0_1025, sog0_1026, sog1_1025, sog1_1026, soh_1431, soh_1433, \
                         soh_1434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1909[k] = f_4 * sog0_1025[k]
                    - f_5 * sog1_1025[k]
                    + f_3 * pc_x[k] * soh_1433[k];

        t_1910[k] = f_6 * sog0_1026[k]
                    - f_7 * sog1_1026[k]
                    + f_3 * pc_x[k] * soh_1434[k];

        t_1911[k] = f_12 * snh_1179[k]
                    + f_3 * pc_z[k] * soh_1431[k];

        t_1912[k] = f_16 * snh_1202[k]
                    + f_3 * pc_y[k] * soh_1433[k];
    }

#pragma omp simd aligned(t_1913, t_1914, t_1915, pc_x, pc_z, snh_1182, sog0_1029, sog0_1030, \
                         sog1_1029, sog1_1030, soh_1434, soh_1437, \
                         soh_1438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1913[k] = f_6 * sog0_1029[k]
                    - f_7 * sog1_1029[k]
                    + f_3 * pc_x[k] * soh_1437[k];

        t_1914[k] = f_8 * sog0_1030[k]
                    - f_9 * sog1_1030[k]
                    + f_3 * pc_x[k] * soh_1438[k];

        t_1915[k] = f_12 * snh_1182[k]
                    + f_3 * pc_z[k] * soh_1434[k];
    }

#pragma omp simd aligned(t_1916, t_1917, t_1918, t_1919, pc_x, pc_y, snh_1206, sog0_1032, \
                         sog0_1034, sog1_1032, sog1_1034, soh_1437, soh_1440, soh_1442, \
                         soh_1443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1916[k] = f_8 * sog0_1032[k]
                    - f_9 * sog1_1032[k]
                    + f_3 * pc_x[k] * soh_1440[k];

        t_1917[k] = f_16 * snh_1206[k]
                    + f_3 * pc_y[k] * soh_1437[k];

        t_1918[k] = f_8 * sog0_1034[k]
                    - f_9 * sog1_1034[k]
                    + f_3 * pc_x[k] * soh_1442[k];

        t_1919[k] = f_3 * pc_x[k] * soh_1443[k];
    }

#pragma omp simd aligned(t_1920, t_1921, t_1922, t_1923, t_1924, pc_x, soh_1444, soh_1445, \
                         soh_1446, soh_1447, soh_1448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1920[k] = f_3 * pc_x[k] * soh_1444[k];

        t_1921[k] = f_3 * pc_x[k] * soh_1445[k];

        t_1922[k] = f_3 * pc_x[k] * soh_1446[k];

        t_1923[k] = f_3 * pc_x[k] * soh_1447[k];

        t_1924[k] = f_3 * pc_x[k] * soh_1448[k];
    }

#pragma omp simd aligned(t_1925, t_1926, t_1927, pc_y, pc_z, snh_1191, snh_1212, snh_1214, \
                         sog0_1030, sog0_1032, sog1_1030, sog1_1032, soh_1443, \
                         soh_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1925[k] = f_16 * snh_1212[k]
                    + f_1 * sog0_1030[k]
                    - f_2 * sog1_1030[k]
                    + f_3 * pc_y[k] * soh_1443[k];

        t_1926[k] = f_12 * snh_1191[k]
                    + f_3 * pc_z[k] * soh_1443[k];

        t_1927[k] = f_16 * snh_1214[k]
                    + f_4 * sog0_1032[k]
                    - f_5 * sog1_1032[k]
                    + f_3 * pc_y[k] * soh_1445[k];
    }

#pragma omp simd aligned(t_1928, t_1929, t_1930, pc_y, snh_1215, snh_1216, snh_1217, \
                         sog0_1033, sog0_1034, sog1_1033, sog1_1034, soh_1446, soh_1447, \
                         soh_1448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1928[k] = f_16 * snh_1215[k]
                    + f_6 * sog0_1033[k]
                    - f_7 * sog1_1033[k]
                    + f_3 * pc_y[k] * soh_1446[k];

        t_1929[k] = f_16 * snh_1216[k]
                    + f_8 * sog0_1034[k]
                    - f_9 * sog1_1034[k]
                    + f_3 * pc_y[k] * soh_1447[k];

        t_1930[k] = f_16 * snh_1217[k]
                    + f_3 * pc_y[k] * soh_1448[k];
    }

#pragma omp simd aligned(t_1931, t_1932, t_1933, t_1934, pc_x, pc_y, pc_z, snh_1196, snh_1197, \
                         snh_1218, sog0_1034, sog0_1035, sog1_1034, sog1_1035, soh_1448, \
                         soh_1449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1931[k] = f_12 * snh_1196[k]
                    + f_1 * sog0_1034[k]
                    - f_2 * sog1_1034[k]
                    + f_3 * pc_z[k] * soh_1448[k];

        t_1932[k] = f_1 * sog0_1035[k]
                    - f_2 * sog1_1035[k]
                    + f_3 * pc_x[k] * soh_1449[k];

        t_1933[k] = f_17 * snh_1218[k]
                    + f_3 * pc_y[k] * soh_1449[k];

        t_1934[k] = f_13 * snh_1197[k]
                    + f_3 * pc_z[k] * soh_1449[k];
    }

#pragma omp simd aligned(t_1935, t_1936, t_1937, pc_x, pc_y, snh_1220, sog0_1038, sog0_1040, \
                         sog1_1038, sog1_1040, soh_1451, soh_1452, \
                         soh_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1935[k] = f_4 * sog0_1038[k]
                    - f_5 * sog1_1038[k]
                    + f_3 * pc_x[k] * soh_1452[k];

        t_1936[k] = f_17 * snh_1220[k]
                    + f_3 * pc_y[k] * soh_1451[k];

        t_1937[k] = f_4 * sog0_1040[k]
                    - f_5 * sog1_1040[k]
                    + f_3 * pc_x[k] * soh_1454[k];
    }

#pragma omp simd aligned(t_1938, t_1939, t_1940, pc_x, pc_y, pc_z, snh_1200, snh_1223, \
                         sog0_1041, sog1_1041, soh_1452, soh_1454, \
                         soh_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1938[k] = f_6 * sog0_1041[k]
                    - f_7 * sog1_1041[k]
                    + f_3 * pc_x[k] * soh_1455[k];

        t_1939[k] = f_13 * snh_1200[k]
                    + f_3 * pc_z[k] * soh_1452[k];

        t_1940[k] = f_17 * snh_1223[k]
                    + f_3 * pc_y[k] * soh_1454[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t snh, const size_t sog0,
                                                           const size_t sog1, const size_t soh,
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
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;
    const auto f_20 = 2.5 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snh_1203 = buffer.data(snh + 1203);
    const auto *snh_1212 = buffer.data(snh + 1212);
    const auto *snh_1217 = buffer.data(snh + 1217);
    const auto *snh_1218 = buffer.data(snh + 1218);
    const auto *snh_1221 = buffer.data(snh + 1221);
    const auto *snh_1224 = buffer.data(snh + 1224);
    const auto *snh_1227 = buffer.data(snh + 1227);
    const auto *snh_1233 = buffer.data(snh + 1233);
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
    const auto *snh_1254 = buffer.data(snh + 1254);
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
    const auto *snh_1275 = buffer.data(snh + 1275);
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
    const auto *snh_1296 = buffer.data(snh + 1296);
    const auto *snh_1298 = buffer.data(snh + 1298);
    const auto *snh_1299 = buffer.data(snh + 1299);
    const auto *snh_1300 = buffer.data(snh + 1300);
    const auto *snh_1301 = buffer.data(snh + 1301);
    const auto *snh_1302 = buffer.data(snh + 1302);
    const auto *snh_1304 = buffer.data(snh + 1304);
    const auto *snh_1307 = buffer.data(snh + 1307);
    const auto *snh_1311 = buffer.data(snh + 1311);

    const auto *sog0_1044 = buffer.data(sog0 + 1044);
    const auto *sog0_1045 = buffer.data(sog0 + 1045);
    const auto *sog0_1047 = buffer.data(sog0 + 1047);
    const auto *sog0_1048 = buffer.data(sog0 + 1048);
    const auto *sog0_1049 = buffer.data(sog0 + 1049);
    const auto *sog0_1050 = buffer.data(sog0 + 1050);
    const auto *sog0_1053 = buffer.data(sog0 + 1053);
    const auto *sog0_1055 = buffer.data(sog0 + 1055);
    const auto *sog0_1056 = buffer.data(sog0 + 1056);
    const auto *sog0_1059 = buffer.data(sog0 + 1059);
    const auto *sog0_1060 = buffer.data(sog0 + 1060);
    const auto *sog0_1062 = buffer.data(sog0 + 1062);
    const auto *sog0_1063 = buffer.data(sog0 + 1063);
    const auto *sog0_1064 = buffer.data(sog0 + 1064);
    const auto *sog0_1065 = buffer.data(sog0 + 1065);
    const auto *sog0_1068 = buffer.data(sog0 + 1068);
    const auto *sog0_1070 = buffer.data(sog0 + 1070);
    const auto *sog0_1071 = buffer.data(sog0 + 1071);
    const auto *sog0_1074 = buffer.data(sog0 + 1074);
    const auto *sog0_1075 = buffer.data(sog0 + 1075);
    const auto *sog0_1077 = buffer.data(sog0 + 1077);
    const auto *sog0_1078 = buffer.data(sog0 + 1078);
    const auto *sog0_1079 = buffer.data(sog0 + 1079);
    const auto *sog0_1080 = buffer.data(sog0 + 1080);
    const auto *sog0_1083 = buffer.data(sog0 + 1083);
    const auto *sog0_1085 = buffer.data(sog0 + 1085);
    const auto *sog0_1086 = buffer.data(sog0 + 1086);
    const auto *sog0_1089 = buffer.data(sog0 + 1089);
    const auto *sog0_1090 = buffer.data(sog0 + 1090);
    const auto *sog0_1092 = buffer.data(sog0 + 1092);
    const auto *sog0_1093 = buffer.data(sog0 + 1093);
    const auto *sog0_1094 = buffer.data(sog0 + 1094);
    const auto *sog0_1095 = buffer.data(sog0 + 1095);
    const auto *sog0_1098 = buffer.data(sog0 + 1098);
    const auto *sog0_1100 = buffer.data(sog0 + 1100);
    const auto *sog0_1101 = buffer.data(sog0 + 1101);
    const auto *sog0_1104 = buffer.data(sog0 + 1104);
    const auto *sog0_1105 = buffer.data(sog0 + 1105);
    const auto *sog0_1107 = buffer.data(sog0 + 1107);
    const auto *sog0_1109 = buffer.data(sog0 + 1109);

    const auto *sog1_1044 = buffer.data(sog1 + 1044);
    const auto *sog1_1045 = buffer.data(sog1 + 1045);
    const auto *sog1_1047 = buffer.data(sog1 + 1047);
    const auto *sog1_1048 = buffer.data(sog1 + 1048);
    const auto *sog1_1049 = buffer.data(sog1 + 1049);
    const auto *sog1_1050 = buffer.data(sog1 + 1050);
    const auto *sog1_1053 = buffer.data(sog1 + 1053);
    const auto *sog1_1055 = buffer.data(sog1 + 1055);
    const auto *sog1_1056 = buffer.data(sog1 + 1056);
    const auto *sog1_1059 = buffer.data(sog1 + 1059);
    const auto *sog1_1060 = buffer.data(sog1 + 1060);
    const auto *sog1_1062 = buffer.data(sog1 + 1062);
    const auto *sog1_1063 = buffer.data(sog1 + 1063);
    const auto *sog1_1064 = buffer.data(sog1 + 1064);
    const auto *sog1_1065 = buffer.data(sog1 + 1065);
    const auto *sog1_1068 = buffer.data(sog1 + 1068);
    const auto *sog1_1070 = buffer.data(sog1 + 1070);
    const auto *sog1_1071 = buffer.data(sog1 + 1071);
    const auto *sog1_1074 = buffer.data(sog1 + 1074);
    const auto *sog1_1075 = buffer.data(sog1 + 1075);
    const auto *sog1_1077 = buffer.data(sog1 + 1077);
    const auto *sog1_1078 = buffer.data(sog1 + 1078);
    const auto *sog1_1079 = buffer.data(sog1 + 1079);
    const auto *sog1_1080 = buffer.data(sog1 + 1080);
    const auto *sog1_1083 = buffer.data(sog1 + 1083);
    const auto *sog1_1085 = buffer.data(sog1 + 1085);
    const auto *sog1_1086 = buffer.data(sog1 + 1086);
    const auto *sog1_1089 = buffer.data(sog1 + 1089);
    const auto *sog1_1090 = buffer.data(sog1 + 1090);
    const auto *sog1_1092 = buffer.data(sog1 + 1092);
    const auto *sog1_1093 = buffer.data(sog1 + 1093);
    const auto *sog1_1094 = buffer.data(sog1 + 1094);
    const auto *sog1_1095 = buffer.data(sog1 + 1095);
    const auto *sog1_1098 = buffer.data(sog1 + 1098);
    const auto *sog1_1100 = buffer.data(sog1 + 1100);
    const auto *sog1_1101 = buffer.data(sog1 + 1101);
    const auto *sog1_1104 = buffer.data(sog1 + 1104);
    const auto *sog1_1105 = buffer.data(sog1 + 1105);
    const auto *sog1_1107 = buffer.data(sog1 + 1107);
    const auto *sog1_1109 = buffer.data(sog1 + 1109);

    const auto *soh_1455 = buffer.data(soh + 1455);
    const auto *soh_1458 = buffer.data(soh + 1458);
    const auto *soh_1459 = buffer.data(soh + 1459);
    const auto *soh_1461 = buffer.data(soh + 1461);
    const auto *soh_1463 = buffer.data(soh + 1463);
    const auto *soh_1464 = buffer.data(soh + 1464);
    const auto *soh_1465 = buffer.data(soh + 1465);
    const auto *soh_1466 = buffer.data(soh + 1466);
    const auto *soh_1467 = buffer.data(soh + 1467);
    const auto *soh_1468 = buffer.data(soh + 1468);
    const auto *soh_1469 = buffer.data(soh + 1469);
    const auto *soh_1470 = buffer.data(soh + 1470);
    const auto *soh_1472 = buffer.data(soh + 1472);
    const auto *soh_1473 = buffer.data(soh + 1473);
    const auto *soh_1475 = buffer.data(soh + 1475);
    const auto *soh_1476 = buffer.data(soh + 1476);
    const auto *soh_1479 = buffer.data(soh + 1479);
    const auto *soh_1480 = buffer.data(soh + 1480);
    const auto *soh_1482 = buffer.data(soh + 1482);
    const auto *soh_1484 = buffer.data(soh + 1484);
    const auto *soh_1485 = buffer.data(soh + 1485);
    const auto *soh_1486 = buffer.data(soh + 1486);
    const auto *soh_1487 = buffer.data(soh + 1487);
    const auto *soh_1488 = buffer.data(soh + 1488);
    const auto *soh_1489 = buffer.data(soh + 1489);
    const auto *soh_1490 = buffer.data(soh + 1490);
    const auto *soh_1491 = buffer.data(soh + 1491);
    const auto *soh_1493 = buffer.data(soh + 1493);
    const auto *soh_1494 = buffer.data(soh + 1494);
    const auto *soh_1496 = buffer.data(soh + 1496);
    const auto *soh_1497 = buffer.data(soh + 1497);
    const auto *soh_1500 = buffer.data(soh + 1500);
    const auto *soh_1501 = buffer.data(soh + 1501);
    const auto *soh_1503 = buffer.data(soh + 1503);
    const auto *soh_1505 = buffer.data(soh + 1505);
    const auto *soh_1506 = buffer.data(soh + 1506);
    const auto *soh_1507 = buffer.data(soh + 1507);
    const auto *soh_1508 = buffer.data(soh + 1508);
    const auto *soh_1509 = buffer.data(soh + 1509);
    const auto *soh_1510 = buffer.data(soh + 1510);
    const auto *soh_1511 = buffer.data(soh + 1511);
    const auto *soh_1512 = buffer.data(soh + 1512);
    const auto *soh_1514 = buffer.data(soh + 1514);
    const auto *soh_1515 = buffer.data(soh + 1515);
    const auto *soh_1517 = buffer.data(soh + 1517);
    const auto *soh_1518 = buffer.data(soh + 1518);
    const auto *soh_1521 = buffer.data(soh + 1521);
    const auto *soh_1522 = buffer.data(soh + 1522);
    const auto *soh_1524 = buffer.data(soh + 1524);
    const auto *soh_1526 = buffer.data(soh + 1526);
    const auto *soh_1527 = buffer.data(soh + 1527);
    const auto *soh_1528 = buffer.data(soh + 1528);
    const auto *soh_1529 = buffer.data(soh + 1529);
    const auto *soh_1530 = buffer.data(soh + 1530);
    const auto *soh_1531 = buffer.data(soh + 1531);
    const auto *soh_1532 = buffer.data(soh + 1532);
    const auto *soh_1533 = buffer.data(soh + 1533);
    const auto *soh_1535 = buffer.data(soh + 1535);
    const auto *soh_1536 = buffer.data(soh + 1536);
    const auto *soh_1538 = buffer.data(soh + 1538);
    const auto *soh_1539 = buffer.data(soh + 1539);
    const auto *soh_1542 = buffer.data(soh + 1542);
    const auto *soh_1543 = buffer.data(soh + 1543);
    const auto *soh_1545 = buffer.data(soh + 1545);
    const auto *soh_1547 = buffer.data(soh + 1547);
    const auto *soh_1548 = buffer.data(soh + 1548);
    const auto *soh_1549 = buffer.data(soh + 1549);
    const auto *soh_1550 = buffer.data(soh + 1550);
    const auto *soh_1551 = buffer.data(soh + 1551);
    const auto *soh_1552 = buffer.data(soh + 1552);
    const auto *soh_1553 = buffer.data(soh + 1553);

#pragma omp simd aligned(t_1941, t_1942, t_1943, pc_x, pc_z, snh_1203, sog0_1044, sog0_1045, \
                         sog1_1044, sog1_1045, soh_1455, soh_1458, \
                         soh_1459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1941[k] = f_6 * sog0_1044[k]
                    - f_7 * sog1_1044[k]
                    + f_3 * pc_x[k] * soh_1458[k];

        t_1942[k] = f_8 * sog0_1045[k]
                    - f_9 * sog1_1045[k]
                    + f_3 * pc_x[k] * soh_1459[k];

        t_1943[k] = f_13 * snh_1203[k]
                    + f_3 * pc_z[k] * soh_1455[k];
    }

#pragma omp simd aligned(t_1944, t_1945, t_1946, t_1947, pc_x, pc_y, snh_1227, sog0_1047, \
                         sog0_1049, sog1_1047, sog1_1049, soh_1458, soh_1461, soh_1463, \
                         soh_1464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1944[k] = f_8 * sog0_1047[k]
                    - f_9 * sog1_1047[k]
                    + f_3 * pc_x[k] * soh_1461[k];

        t_1945[k] = f_17 * snh_1227[k]
                    + f_3 * pc_y[k] * soh_1458[k];

        t_1946[k] = f_8 * sog0_1049[k]
                    - f_9 * sog1_1049[k]
                    + f_3 * pc_x[k] * soh_1463[k];

        t_1947[k] = f_3 * pc_x[k] * soh_1464[k];
    }

#pragma omp simd aligned(t_1948, t_1949, t_1950, t_1951, t_1952, pc_x, soh_1465, soh_1466, \
                         soh_1467, soh_1468, soh_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1948[k] = f_3 * pc_x[k] * soh_1465[k];

        t_1949[k] = f_3 * pc_x[k] * soh_1466[k];

        t_1950[k] = f_3 * pc_x[k] * soh_1467[k];

        t_1951[k] = f_3 * pc_x[k] * soh_1468[k];

        t_1952[k] = f_3 * pc_x[k] * soh_1469[k];
    }

#pragma omp simd aligned(t_1953, t_1954, t_1955, pc_y, pc_z, snh_1212, snh_1233, snh_1235, \
                         sog0_1045, sog0_1047, sog1_1045, sog1_1047, soh_1464, \
                         soh_1466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1953[k] = f_17 * snh_1233[k]
                    + f_1 * sog0_1045[k]
                    - f_2 * sog1_1045[k]
                    + f_3 * pc_y[k] * soh_1464[k];

        t_1954[k] = f_13 * snh_1212[k]
                    + f_3 * pc_z[k] * soh_1464[k];

        t_1955[k] = f_17 * snh_1235[k]
                    + f_4 * sog0_1047[k]
                    - f_5 * sog1_1047[k]
                    + f_3 * pc_y[k] * soh_1466[k];
    }

#pragma omp simd aligned(t_1956, t_1957, t_1958, pc_y, snh_1236, snh_1237, snh_1238, \
                         sog0_1048, sog0_1049, sog1_1048, sog1_1049, soh_1467, soh_1468, \
                         soh_1469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1956[k] = f_17 * snh_1236[k]
                    + f_6 * sog0_1048[k]
                    - f_7 * sog1_1048[k]
                    + f_3 * pc_y[k] * soh_1467[k];

        t_1957[k] = f_17 * snh_1237[k]
                    + f_8 * sog0_1049[k]
                    - f_9 * sog1_1049[k]
                    + f_3 * pc_y[k] * soh_1468[k];

        t_1958[k] = f_17 * snh_1238[k]
                    + f_3 * pc_y[k] * soh_1469[k];
    }

#pragma omp simd aligned(t_1959, t_1960, t_1961, t_1962, pc_x, pc_y, pc_z, snh_1217, snh_1218, \
                         snh_1239, sog0_1049, sog0_1050, sog1_1049, sog1_1050, soh_1469, \
                         soh_1470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1959[k] = f_13 * snh_1217[k]
                    + f_1 * sog0_1049[k]
                    - f_2 * sog1_1049[k]
                    + f_3 * pc_z[k] * soh_1469[k];

        t_1960[k] = f_1 * sog0_1050[k]
                    - f_2 * sog1_1050[k]
                    + f_3 * pc_x[k] * soh_1470[k];

        t_1961[k] = f_18 * snh_1239[k]
                    + f_3 * pc_y[k] * soh_1470[k];

        t_1962[k] = f_14 * snh_1218[k]
                    + f_3 * pc_z[k] * soh_1470[k];
    }

#pragma omp simd aligned(t_1963, t_1964, t_1965, pc_x, pc_y, snh_1241, sog0_1053, sog0_1055, \
                         sog1_1053, sog1_1055, soh_1472, soh_1473, \
                         soh_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1963[k] = f_4 * sog0_1053[k]
                    - f_5 * sog1_1053[k]
                    + f_3 * pc_x[k] * soh_1473[k];

        t_1964[k] = f_18 * snh_1241[k]
                    + f_3 * pc_y[k] * soh_1472[k];

        t_1965[k] = f_4 * sog0_1055[k]
                    - f_5 * sog1_1055[k]
                    + f_3 * pc_x[k] * soh_1475[k];
    }

#pragma omp simd aligned(t_1966, t_1967, t_1968, pc_x, pc_y, pc_z, snh_1221, snh_1244, \
                         sog0_1056, sog1_1056, soh_1473, soh_1475, \
                         soh_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1966[k] = f_6 * sog0_1056[k]
                    - f_7 * sog1_1056[k]
                    + f_3 * pc_x[k] * soh_1476[k];

        t_1967[k] = f_14 * snh_1221[k]
                    + f_3 * pc_z[k] * soh_1473[k];

        t_1968[k] = f_18 * snh_1244[k]
                    + f_3 * pc_y[k] * soh_1475[k];
    }

#pragma omp simd aligned(t_1969, t_1970, t_1971, pc_x, pc_z, snh_1224, sog0_1059, sog0_1060, \
                         sog1_1059, sog1_1060, soh_1476, soh_1479, \
                         soh_1480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1969[k] = f_6 * sog0_1059[k]
                    - f_7 * sog1_1059[k]
                    + f_3 * pc_x[k] * soh_1479[k];

        t_1970[k] = f_8 * sog0_1060[k]
                    - f_9 * sog1_1060[k]
                    + f_3 * pc_x[k] * soh_1480[k];

        t_1971[k] = f_14 * snh_1224[k]
                    + f_3 * pc_z[k] * soh_1476[k];
    }

#pragma omp simd aligned(t_1972, t_1973, t_1974, t_1975, pc_x, pc_y, snh_1248, sog0_1062, \
                         sog0_1064, sog1_1062, sog1_1064, soh_1479, soh_1482, soh_1484, \
                         soh_1485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1972[k] = f_8 * sog0_1062[k]
                    - f_9 * sog1_1062[k]
                    + f_3 * pc_x[k] * soh_1482[k];

        t_1973[k] = f_18 * snh_1248[k]
                    + f_3 * pc_y[k] * soh_1479[k];

        t_1974[k] = f_8 * sog0_1064[k]
                    - f_9 * sog1_1064[k]
                    + f_3 * pc_x[k] * soh_1484[k];

        t_1975[k] = f_3 * pc_x[k] * soh_1485[k];
    }

#pragma omp simd aligned(t_1976, t_1977, t_1978, t_1979, t_1980, pc_x, soh_1486, soh_1487, \
                         soh_1488, soh_1489, soh_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1976[k] = f_3 * pc_x[k] * soh_1486[k];

        t_1977[k] = f_3 * pc_x[k] * soh_1487[k];

        t_1978[k] = f_3 * pc_x[k] * soh_1488[k];

        t_1979[k] = f_3 * pc_x[k] * soh_1489[k];

        t_1980[k] = f_3 * pc_x[k] * soh_1490[k];
    }

#pragma omp simd aligned(t_1981, t_1982, t_1983, pc_y, pc_z, snh_1233, snh_1254, snh_1256, \
                         sog0_1060, sog0_1062, sog1_1060, sog1_1062, soh_1485, \
                         soh_1487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1981[k] = f_18 * snh_1254[k]
                    + f_1 * sog0_1060[k]
                    - f_2 * sog1_1060[k]
                    + f_3 * pc_y[k] * soh_1485[k];

        t_1982[k] = f_14 * snh_1233[k]
                    + f_3 * pc_z[k] * soh_1485[k];

        t_1983[k] = f_18 * snh_1256[k]
                    + f_4 * sog0_1062[k]
                    - f_5 * sog1_1062[k]
                    + f_3 * pc_y[k] * soh_1487[k];
    }

#pragma omp simd aligned(t_1984, t_1985, t_1986, pc_y, snh_1257, snh_1258, snh_1259, \
                         sog0_1063, sog0_1064, sog1_1063, sog1_1064, soh_1488, soh_1489, \
                         soh_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1984[k] = f_18 * snh_1257[k]
                    + f_6 * sog0_1063[k]
                    - f_7 * sog1_1063[k]
                    + f_3 * pc_y[k] * soh_1488[k];

        t_1985[k] = f_18 * snh_1258[k]
                    + f_8 * sog0_1064[k]
                    - f_9 * sog1_1064[k]
                    + f_3 * pc_y[k] * soh_1489[k];

        t_1986[k] = f_18 * snh_1259[k]
                    + f_3 * pc_y[k] * soh_1490[k];
    }

#pragma omp simd aligned(t_1987, t_1988, t_1989, t_1990, pc_x, pc_y, pc_z, snh_1238, snh_1239, \
                         snh_1260, sog0_1064, sog0_1065, sog1_1064, sog1_1065, soh_1490, \
                         soh_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1987[k] = f_14 * snh_1238[k]
                    + f_1 * sog0_1064[k]
                    - f_2 * sog1_1064[k]
                    + f_3 * pc_z[k] * soh_1490[k];

        t_1988[k] = f_1 * sog0_1065[k]
                    - f_2 * sog1_1065[k]
                    + f_3 * pc_x[k] * soh_1491[k];

        t_1989[k] = f_19 * snh_1260[k]
                    + f_3 * pc_y[k] * soh_1491[k];

        t_1990[k] = f_20 * snh_1239[k]
                    + f_3 * pc_z[k] * soh_1491[k];
    }

#pragma omp simd aligned(t_1991, t_1992, t_1993, pc_x, pc_y, snh_1262, sog0_1068, sog0_1070, \
                         sog1_1068, sog1_1070, soh_1493, soh_1494, \
                         soh_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1991[k] = f_4 * sog0_1068[k]
                    - f_5 * sog1_1068[k]
                    + f_3 * pc_x[k] * soh_1494[k];

        t_1992[k] = f_19 * snh_1262[k]
                    + f_3 * pc_y[k] * soh_1493[k];

        t_1993[k] = f_4 * sog0_1070[k]
                    - f_5 * sog1_1070[k]
                    + f_3 * pc_x[k] * soh_1496[k];
    }

#pragma omp simd aligned(t_1994, t_1995, t_1996, pc_x, pc_y, pc_z, snh_1242, snh_1265, \
                         sog0_1071, sog1_1071, soh_1494, soh_1496, \
                         soh_1497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1994[k] = f_6 * sog0_1071[k]
                    - f_7 * sog1_1071[k]
                    + f_3 * pc_x[k] * soh_1497[k];

        t_1995[k] = f_20 * snh_1242[k]
                    + f_3 * pc_z[k] * soh_1494[k];

        t_1996[k] = f_19 * snh_1265[k]
                    + f_3 * pc_y[k] * soh_1496[k];
    }

#pragma omp simd aligned(t_1997, t_1998, t_1999, pc_x, pc_z, snh_1245, sog0_1074, sog0_1075, \
                         sog1_1074, sog1_1075, soh_1497, soh_1500, \
                         soh_1501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1997[k] = f_6 * sog0_1074[k]
                    - f_7 * sog1_1074[k]
                    + f_3 * pc_x[k] * soh_1500[k];

        t_1998[k] = f_8 * sog0_1075[k]
                    - f_9 * sog1_1075[k]
                    + f_3 * pc_x[k] * soh_1501[k];

        t_1999[k] = f_20 * snh_1245[k]
                    + f_3 * pc_z[k] * soh_1497[k];
    }

#pragma omp simd aligned(t_2000, t_2001, t_2002, t_2003, pc_x, pc_y, snh_1269, sog0_1077, \
                         sog0_1079, sog1_1077, sog1_1079, soh_1500, soh_1503, soh_1505, \
                         soh_1506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2000[k] = f_8 * sog0_1077[k]
                    - f_9 * sog1_1077[k]
                    + f_3 * pc_x[k] * soh_1503[k];

        t_2001[k] = f_19 * snh_1269[k]
                    + f_3 * pc_y[k] * soh_1500[k];

        t_2002[k] = f_8 * sog0_1079[k]
                    - f_9 * sog1_1079[k]
                    + f_3 * pc_x[k] * soh_1505[k];

        t_2003[k] = f_3 * pc_x[k] * soh_1506[k];
    }

#pragma omp simd aligned(t_2004, t_2005, t_2006, t_2007, t_2008, pc_x, soh_1507, soh_1508, \
                         soh_1509, soh_1510, soh_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2004[k] = f_3 * pc_x[k] * soh_1507[k];

        t_2005[k] = f_3 * pc_x[k] * soh_1508[k];

        t_2006[k] = f_3 * pc_x[k] * soh_1509[k];

        t_2007[k] = f_3 * pc_x[k] * soh_1510[k];

        t_2008[k] = f_3 * pc_x[k] * soh_1511[k];
    }

#pragma omp simd aligned(t_2009, t_2010, t_2011, pc_y, pc_z, snh_1254, snh_1275, snh_1277, \
                         sog0_1075, sog0_1077, sog1_1075, sog1_1077, soh_1506, \
                         soh_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2009[k] = f_19 * snh_1275[k]
                    + f_1 * sog0_1075[k]
                    - f_2 * sog1_1075[k]
                    + f_3 * pc_y[k] * soh_1506[k];

        t_2010[k] = f_20 * snh_1254[k]
                    + f_3 * pc_z[k] * soh_1506[k];

        t_2011[k] = f_19 * snh_1277[k]
                    + f_4 * sog0_1077[k]
                    - f_5 * sog1_1077[k]
                    + f_3 * pc_y[k] * soh_1508[k];
    }

#pragma omp simd aligned(t_2012, t_2013, t_2014, pc_y, snh_1278, snh_1279, snh_1280, \
                         sog0_1078, sog0_1079, sog1_1078, sog1_1079, soh_1509, soh_1510, \
                         soh_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2012[k] = f_19 * snh_1278[k]
                    + f_6 * sog0_1078[k]
                    - f_7 * sog1_1078[k]
                    + f_3 * pc_y[k] * soh_1509[k];

        t_2013[k] = f_19 * snh_1279[k]
                    + f_8 * sog0_1079[k]
                    - f_9 * sog1_1079[k]
                    + f_3 * pc_y[k] * soh_1510[k];

        t_2014[k] = f_19 * snh_1280[k]
                    + f_3 * pc_y[k] * soh_1511[k];
    }

#pragma omp simd aligned(t_2015, t_2016, t_2017, t_2018, pc_x, pc_y, pc_z, snh_1259, snh_1260, \
                         snh_1281, sog0_1079, sog0_1080, sog1_1079, sog1_1080, soh_1511, \
                         soh_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2015[k] = f_20 * snh_1259[k]
                    + f_1 * sog0_1079[k]
                    - f_2 * sog1_1079[k]
                    + f_3 * pc_z[k] * soh_1511[k];

        t_2016[k] = f_1 * sog0_1080[k]
                    - f_2 * sog1_1080[k]
                    + f_3 * pc_x[k] * soh_1512[k];

        t_2017[k] = f_20 * snh_1281[k]
                    + f_3 * pc_y[k] * soh_1512[k];

        t_2018[k] = f_19 * snh_1260[k]
                    + f_3 * pc_z[k] * soh_1512[k];
    }

#pragma omp simd aligned(t_2019, t_2020, t_2021, pc_x, pc_y, snh_1283, sog0_1083, sog0_1085, \
                         sog1_1083, sog1_1085, soh_1514, soh_1515, \
                         soh_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2019[k] = f_4 * sog0_1083[k]
                    - f_5 * sog1_1083[k]
                    + f_3 * pc_x[k] * soh_1515[k];

        t_2020[k] = f_20 * snh_1283[k]
                    + f_3 * pc_y[k] * soh_1514[k];

        t_2021[k] = f_4 * sog0_1085[k]
                    - f_5 * sog1_1085[k]
                    + f_3 * pc_x[k] * soh_1517[k];
    }

#pragma omp simd aligned(t_2022, t_2023, t_2024, pc_x, pc_y, pc_z, snh_1263, snh_1286, \
                         sog0_1086, sog1_1086, soh_1515, soh_1517, \
                         soh_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2022[k] = f_6 * sog0_1086[k]
                    - f_7 * sog1_1086[k]
                    + f_3 * pc_x[k] * soh_1518[k];

        t_2023[k] = f_19 * snh_1263[k]
                    + f_3 * pc_z[k] * soh_1515[k];

        t_2024[k] = f_20 * snh_1286[k]
                    + f_3 * pc_y[k] * soh_1517[k];
    }

#pragma omp simd aligned(t_2025, t_2026, t_2027, pc_x, pc_z, snh_1266, sog0_1089, sog0_1090, \
                         sog1_1089, sog1_1090, soh_1518, soh_1521, \
                         soh_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2025[k] = f_6 * sog0_1089[k]
                    - f_7 * sog1_1089[k]
                    + f_3 * pc_x[k] * soh_1521[k];

        t_2026[k] = f_8 * sog0_1090[k]
                    - f_9 * sog1_1090[k]
                    + f_3 * pc_x[k] * soh_1522[k];

        t_2027[k] = f_19 * snh_1266[k]
                    + f_3 * pc_z[k] * soh_1518[k];
    }

#pragma omp simd aligned(t_2028, t_2029, t_2030, t_2031, pc_x, pc_y, snh_1290, sog0_1092, \
                         sog0_1094, sog1_1092, sog1_1094, soh_1521, soh_1524, soh_1526, \
                         soh_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2028[k] = f_8 * sog0_1092[k]
                    - f_9 * sog1_1092[k]
                    + f_3 * pc_x[k] * soh_1524[k];

        t_2029[k] = f_20 * snh_1290[k]
                    + f_3 * pc_y[k] * soh_1521[k];

        t_2030[k] = f_8 * sog0_1094[k]
                    - f_9 * sog1_1094[k]
                    + f_3 * pc_x[k] * soh_1526[k];

        t_2031[k] = f_3 * pc_x[k] * soh_1527[k];
    }

#pragma omp simd aligned(t_2032, t_2033, t_2034, t_2035, t_2036, pc_x, soh_1528, soh_1529, \
                         soh_1530, soh_1531, soh_1532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2032[k] = f_3 * pc_x[k] * soh_1528[k];

        t_2033[k] = f_3 * pc_x[k] * soh_1529[k];

        t_2034[k] = f_3 * pc_x[k] * soh_1530[k];

        t_2035[k] = f_3 * pc_x[k] * soh_1531[k];

        t_2036[k] = f_3 * pc_x[k] * soh_1532[k];
    }

#pragma omp simd aligned(t_2037, t_2038, t_2039, pc_y, pc_z, snh_1275, snh_1296, snh_1298, \
                         sog0_1090, sog0_1092, sog1_1090, sog1_1092, soh_1527, \
                         soh_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2037[k] = f_20 * snh_1296[k]
                    + f_1 * sog0_1090[k]
                    - f_2 * sog1_1090[k]
                    + f_3 * pc_y[k] * soh_1527[k];

        t_2038[k] = f_19 * snh_1275[k]
                    + f_3 * pc_z[k] * soh_1527[k];

        t_2039[k] = f_20 * snh_1298[k]
                    + f_4 * sog0_1092[k]
                    - f_5 * sog1_1092[k]
                    + f_3 * pc_y[k] * soh_1529[k];
    }

#pragma omp simd aligned(t_2040, t_2041, t_2042, pc_y, snh_1299, snh_1300, snh_1301, \
                         sog0_1093, sog0_1094, sog1_1093, sog1_1094, soh_1530, soh_1531, \
                         soh_1532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2040[k] = f_20 * snh_1299[k]
                    + f_6 * sog0_1093[k]
                    - f_7 * sog1_1093[k]
                    + f_3 * pc_y[k] * soh_1530[k];

        t_2041[k] = f_20 * snh_1300[k]
                    + f_8 * sog0_1094[k]
                    - f_9 * sog1_1094[k]
                    + f_3 * pc_y[k] * soh_1531[k];

        t_2042[k] = f_20 * snh_1301[k]
                    + f_3 * pc_y[k] * soh_1532[k];
    }

#pragma omp simd aligned(t_2043, t_2044, t_2045, t_2046, pc_x, pc_y, pc_z, snh_1280, snh_1281, \
                         snh_1302, sog0_1094, sog0_1095, sog1_1094, sog1_1095, soh_1532, \
                         soh_1533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2043[k] = f_19 * snh_1280[k]
                    + f_1 * sog0_1094[k]
                    - f_2 * sog1_1094[k]
                    + f_3 * pc_z[k] * soh_1532[k];

        t_2044[k] = f_1 * sog0_1095[k]
                    - f_2 * sog1_1095[k]
                    + f_3 * pc_x[k] * soh_1533[k];

        t_2045[k] = f_14 * snh_1302[k]
                    + f_3 * pc_y[k] * soh_1533[k];

        t_2046[k] = f_18 * snh_1281[k]
                    + f_3 * pc_z[k] * soh_1533[k];
    }

#pragma omp simd aligned(t_2047, t_2048, t_2049, pc_x, pc_y, snh_1304, sog0_1098, sog0_1100, \
                         sog1_1098, sog1_1100, soh_1535, soh_1536, \
                         soh_1538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2047[k] = f_4 * sog0_1098[k]
                    - f_5 * sog1_1098[k]
                    + f_3 * pc_x[k] * soh_1536[k];

        t_2048[k] = f_14 * snh_1304[k]
                    + f_3 * pc_y[k] * soh_1535[k];

        t_2049[k] = f_4 * sog0_1100[k]
                    - f_5 * sog1_1100[k]
                    + f_3 * pc_x[k] * soh_1538[k];
    }

#pragma omp simd aligned(t_2050, t_2051, t_2052, pc_x, pc_y, pc_z, snh_1284, snh_1307, \
                         sog0_1101, sog1_1101, soh_1536, soh_1538, \
                         soh_1539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2050[k] = f_6 * sog0_1101[k]
                    - f_7 * sog1_1101[k]
                    + f_3 * pc_x[k] * soh_1539[k];

        t_2051[k] = f_18 * snh_1284[k]
                    + f_3 * pc_z[k] * soh_1536[k];

        t_2052[k] = f_14 * snh_1307[k]
                    + f_3 * pc_y[k] * soh_1538[k];
    }

#pragma omp simd aligned(t_2053, t_2054, t_2055, pc_x, pc_z, snh_1287, sog0_1104, sog0_1105, \
                         sog1_1104, sog1_1105, soh_1539, soh_1542, \
                         soh_1543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2053[k] = f_6 * sog0_1104[k]
                    - f_7 * sog1_1104[k]
                    + f_3 * pc_x[k] * soh_1542[k];

        t_2054[k] = f_8 * sog0_1105[k]
                    - f_9 * sog1_1105[k]
                    + f_3 * pc_x[k] * soh_1543[k];

        t_2055[k] = f_18 * snh_1287[k]
                    + f_3 * pc_z[k] * soh_1539[k];
    }

#pragma omp simd aligned(t_2056, t_2057, t_2058, t_2059, pc_x, pc_y, snh_1311, sog0_1107, \
                         sog0_1109, sog1_1107, sog1_1109, soh_1542, soh_1545, soh_1547, \
                         soh_1548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2056[k] = f_8 * sog0_1107[k]
                    - f_9 * sog1_1107[k]
                    + f_3 * pc_x[k] * soh_1545[k];

        t_2057[k] = f_14 * snh_1311[k]
                    + f_3 * pc_y[k] * soh_1542[k];

        t_2058[k] = f_8 * sog0_1109[k]
                    - f_9 * sog1_1109[k]
                    + f_3 * pc_x[k] * soh_1547[k];

        t_2059[k] = f_3 * pc_x[k] * soh_1548[k];
    }

#pragma omp simd aligned(t_2060, t_2061, t_2062, t_2063, t_2064, pc_x, soh_1549, soh_1550, \
                         soh_1551, soh_1552, soh_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2060[k] = f_3 * pc_x[k] * soh_1549[k];

        t_2061[k] = f_3 * pc_x[k] * soh_1550[k];

        t_2062[k] = f_3 * pc_x[k] * soh_1551[k];

        t_2063[k] = f_3 * pc_x[k] * soh_1552[k];

        t_2064[k] = f_3 * pc_x[k] * soh_1553[k];
    }
}

static auto
compute_prim_soi_three_center_electron_repulsion_0_piece18(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sni0,
                                                           const size_t snh, const size_t sni1,
                                                           const size_t sog0, const size_t sog1,
                                                           const size_t soh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
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
    const auto f_15 = 5.0 / q;
    const auto f_16 = 4.5 / q;
    const auto f_17 = 4.0 / q;
    const auto f_18 = 3.5 / q;
    const auto f_19 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sni0_1820 = buffer.data(sni0 + 1820);
    const auto *sni0_1825 = buffer.data(sni0 + 1825);
    const auto *sni0_1829 = buffer.data(sni0 + 1829);
    const auto *sni0_1834 = buffer.data(sni0 + 1834);
    const auto *sni0_1841 = buffer.data(sni0 + 1841);
    const auto *sni0_1843 = buffer.data(sni0 + 1843);
    const auto *sni0_1844 = buffer.data(sni0 + 1844);
    const auto *sni0_1845 = buffer.data(sni0 + 1845);
    const auto *sni0_1847 = buffer.data(sni0 + 1847);

    const auto *snh_1296 = buffer.data(snh + 1296);
    const auto *snh_1301 = buffer.data(snh + 1301);
    const auto *snh_1302 = buffer.data(snh + 1302);
    const auto *snh_1305 = buffer.data(snh + 1305);
    const auto *snh_1308 = buffer.data(snh + 1308);
    const auto *snh_1317 = buffer.data(snh + 1317);
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
    const auto *snh_1338 = buffer.data(snh + 1338);
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
    const auto *snh_1359 = buffer.data(snh + 1359);
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
    const auto *snh_1380 = buffer.data(snh + 1380);
    const auto *snh_1382 = buffer.data(snh + 1382);
    const auto *snh_1383 = buffer.data(snh + 1383);
    const auto *snh_1384 = buffer.data(snh + 1384);
    const auto *snh_1385 = buffer.data(snh + 1385);

    const auto *sni1_1820 = buffer.data(sni1 + 1820);
    const auto *sni1_1825 = buffer.data(sni1 + 1825);
    const auto *sni1_1829 = buffer.data(sni1 + 1829);
    const auto *sni1_1834 = buffer.data(sni1 + 1834);
    const auto *sni1_1841 = buffer.data(sni1 + 1841);
    const auto *sni1_1843 = buffer.data(sni1 + 1843);
    const auto *sni1_1844 = buffer.data(sni1 + 1844);
    const auto *sni1_1845 = buffer.data(sni1 + 1845);
    const auto *sni1_1847 = buffer.data(sni1 + 1847);

    const auto *sog0_1105 = buffer.data(sog0 + 1105);
    const auto *sog0_1107 = buffer.data(sog0 + 1107);
    const auto *sog0_1108 = buffer.data(sog0 + 1108);
    const auto *sog0_1109 = buffer.data(sog0 + 1109);
    const auto *sog0_1110 = buffer.data(sog0 + 1110);
    const auto *sog0_1113 = buffer.data(sog0 + 1113);
    const auto *sog0_1115 = buffer.data(sog0 + 1115);
    const auto *sog0_1116 = buffer.data(sog0 + 1116);
    const auto *sog0_1119 = buffer.data(sog0 + 1119);
    const auto *sog0_1120 = buffer.data(sog0 + 1120);
    const auto *sog0_1122 = buffer.data(sog0 + 1122);
    const auto *sog0_1123 = buffer.data(sog0 + 1123);
    const auto *sog0_1124 = buffer.data(sog0 + 1124);
    const auto *sog0_1125 = buffer.data(sog0 + 1125);
    const auto *sog0_1128 = buffer.data(sog0 + 1128);
    const auto *sog0_1130 = buffer.data(sog0 + 1130);
    const auto *sog0_1131 = buffer.data(sog0 + 1131);
    const auto *sog0_1134 = buffer.data(sog0 + 1134);
    const auto *sog0_1135 = buffer.data(sog0 + 1135);
    const auto *sog0_1137 = buffer.data(sog0 + 1137);
    const auto *sog0_1138 = buffer.data(sog0 + 1138);
    const auto *sog0_1139 = buffer.data(sog0 + 1139);
    const auto *sog0_1143 = buffer.data(sog0 + 1143);
    const auto *sog0_1146 = buffer.data(sog0 + 1146);
    const auto *sog0_1150 = buffer.data(sog0 + 1150);
    const auto *sog0_1152 = buffer.data(sog0 + 1152);
    const auto *sog0_1155 = buffer.data(sog0 + 1155);
    const auto *sog0_1158 = buffer.data(sog0 + 1158);
    const auto *sog0_1160 = buffer.data(sog0 + 1160);
    const auto *sog0_1161 = buffer.data(sog0 + 1161);
    const auto *sog0_1164 = buffer.data(sog0 + 1164);
    const auto *sog0_1165 = buffer.data(sog0 + 1165);
    const auto *sog0_1167 = buffer.data(sog0 + 1167);
    const auto *sog0_1168 = buffer.data(sog0 + 1168);
    const auto *sog0_1169 = buffer.data(sog0 + 1169);

    const auto *sog1_1105 = buffer.data(sog1 + 1105);
    const auto *sog1_1107 = buffer.data(sog1 + 1107);
    const auto *sog1_1108 = buffer.data(sog1 + 1108);
    const auto *sog1_1109 = buffer.data(sog1 + 1109);
    const auto *sog1_1110 = buffer.data(sog1 + 1110);
    const auto *sog1_1113 = buffer.data(sog1 + 1113);
    const auto *sog1_1115 = buffer.data(sog1 + 1115);
    const auto *sog1_1116 = buffer.data(sog1 + 1116);
    const auto *sog1_1119 = buffer.data(sog1 + 1119);
    const auto *sog1_1120 = buffer.data(sog1 + 1120);
    const auto *sog1_1122 = buffer.data(sog1 + 1122);
    const auto *sog1_1123 = buffer.data(sog1 + 1123);
    const auto *sog1_1124 = buffer.data(sog1 + 1124);
    const auto *sog1_1125 = buffer.data(sog1 + 1125);
    const auto *sog1_1128 = buffer.data(sog1 + 1128);
    const auto *sog1_1130 = buffer.data(sog1 + 1130);
    const auto *sog1_1131 = buffer.data(sog1 + 1131);
    const auto *sog1_1134 = buffer.data(sog1 + 1134);
    const auto *sog1_1135 = buffer.data(sog1 + 1135);
    const auto *sog1_1137 = buffer.data(sog1 + 1137);
    const auto *sog1_1138 = buffer.data(sog1 + 1138);
    const auto *sog1_1139 = buffer.data(sog1 + 1139);
    const auto *sog1_1143 = buffer.data(sog1 + 1143);
    const auto *sog1_1146 = buffer.data(sog1 + 1146);
    const auto *sog1_1150 = buffer.data(sog1 + 1150);
    const auto *sog1_1152 = buffer.data(sog1 + 1152);
    const auto *sog1_1155 = buffer.data(sog1 + 1155);
    const auto *sog1_1158 = buffer.data(sog1 + 1158);
    const auto *sog1_1160 = buffer.data(sog1 + 1160);
    const auto *sog1_1161 = buffer.data(sog1 + 1161);
    const auto *sog1_1164 = buffer.data(sog1 + 1164);
    const auto *sog1_1165 = buffer.data(sog1 + 1165);
    const auto *sog1_1167 = buffer.data(sog1 + 1167);
    const auto *sog1_1168 = buffer.data(sog1 + 1168);
    const auto *sog1_1169 = buffer.data(sog1 + 1169);

    const auto *soh_1548 = buffer.data(soh + 1548);
    const auto *soh_1550 = buffer.data(soh + 1550);
    const auto *soh_1551 = buffer.data(soh + 1551);
    const auto *soh_1552 = buffer.data(soh + 1552);
    const auto *soh_1553 = buffer.data(soh + 1553);
    const auto *soh_1554 = buffer.data(soh + 1554);
    const auto *soh_1556 = buffer.data(soh + 1556);
    const auto *soh_1557 = buffer.data(soh + 1557);
    const auto *soh_1559 = buffer.data(soh + 1559);
    const auto *soh_1560 = buffer.data(soh + 1560);
    const auto *soh_1563 = buffer.data(soh + 1563);
    const auto *soh_1564 = buffer.data(soh + 1564);
    const auto *soh_1566 = buffer.data(soh + 1566);
    const auto *soh_1568 = buffer.data(soh + 1568);
    const auto *soh_1569 = buffer.data(soh + 1569);
    const auto *soh_1570 = buffer.data(soh + 1570);
    const auto *soh_1571 = buffer.data(soh + 1571);
    const auto *soh_1572 = buffer.data(soh + 1572);
    const auto *soh_1573 = buffer.data(soh + 1573);
    const auto *soh_1574 = buffer.data(soh + 1574);
    const auto *soh_1575 = buffer.data(soh + 1575);
    const auto *soh_1577 = buffer.data(soh + 1577);
    const auto *soh_1578 = buffer.data(soh + 1578);
    const auto *soh_1580 = buffer.data(soh + 1580);
    const auto *soh_1581 = buffer.data(soh + 1581);
    const auto *soh_1584 = buffer.data(soh + 1584);
    const auto *soh_1585 = buffer.data(soh + 1585);
    const auto *soh_1587 = buffer.data(soh + 1587);
    const auto *soh_1589 = buffer.data(soh + 1589);
    const auto *soh_1590 = buffer.data(soh + 1590);
    const auto *soh_1591 = buffer.data(soh + 1591);
    const auto *soh_1592 = buffer.data(soh + 1592);
    const auto *soh_1593 = buffer.data(soh + 1593);
    const auto *soh_1594 = buffer.data(soh + 1594);
    const auto *soh_1595 = buffer.data(soh + 1595);
    const auto *soh_1596 = buffer.data(soh + 1596);
    const auto *soh_1598 = buffer.data(soh + 1598);
    const auto *soh_1599 = buffer.data(soh + 1599);
    const auto *soh_1601 = buffer.data(soh + 1601);
    const auto *soh_1602 = buffer.data(soh + 1602);
    const auto *soh_1605 = buffer.data(soh + 1605);
    const auto *soh_1606 = buffer.data(soh + 1606);
    const auto *soh_1608 = buffer.data(soh + 1608);
    const auto *soh_1611 = buffer.data(soh + 1611);
    const auto *soh_1612 = buffer.data(soh + 1612);
    const auto *soh_1613 = buffer.data(soh + 1613);
    const auto *soh_1614 = buffer.data(soh + 1614);
    const auto *soh_1615 = buffer.data(soh + 1615);
    const auto *soh_1616 = buffer.data(soh + 1616);
    const auto *soh_1617 = buffer.data(soh + 1617);
    const auto *soh_1619 = buffer.data(soh + 1619);
    const auto *soh_1620 = buffer.data(soh + 1620);
    const auto *soh_1622 = buffer.data(soh + 1622);
    const auto *soh_1623 = buffer.data(soh + 1623);
    const auto *soh_1626 = buffer.data(soh + 1626);
    const auto *soh_1627 = buffer.data(soh + 1627);
    const auto *soh_1629 = buffer.data(soh + 1629);
    const auto *soh_1631 = buffer.data(soh + 1631);
    const auto *soh_1632 = buffer.data(soh + 1632);
    const auto *soh_1633 = buffer.data(soh + 1633);
    const auto *soh_1634 = buffer.data(soh + 1634);
    const auto *soh_1635 = buffer.data(soh + 1635);
    const auto *soh_1636 = buffer.data(soh + 1636);
    const auto *soh_1637 = buffer.data(soh + 1637);

#pragma omp simd aligned(t_2065, t_2066, t_2067, pc_y, pc_z, snh_1296, snh_1317, snh_1319, \
                         sog0_1105, sog0_1107, sog1_1105, sog1_1107, soh_1548, \
                         soh_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2065[k] = f_14 * snh_1317[k]
                    + f_1 * sog0_1105[k]
                    - f_2 * sog1_1105[k]
                    + f_3 * pc_y[k] * soh_1548[k];

        t_2066[k] = f_18 * snh_1296[k]
                    + f_3 * pc_z[k] * soh_1548[k];

        t_2067[k] = f_14 * snh_1319[k]
                    + f_4 * sog0_1107[k]
                    - f_5 * sog1_1107[k]
                    + f_3 * pc_y[k] * soh_1550[k];
    }

#pragma omp simd aligned(t_2068, t_2069, t_2070, pc_y, snh_1320, snh_1321, snh_1322, \
                         sog0_1108, sog0_1109, sog1_1108, sog1_1109, soh_1551, soh_1552, \
                         soh_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2068[k] = f_14 * snh_1320[k]
                    + f_6 * sog0_1108[k]
                    - f_7 * sog1_1108[k]
                    + f_3 * pc_y[k] * soh_1551[k];

        t_2069[k] = f_14 * snh_1321[k]
                    + f_8 * sog0_1109[k]
                    - f_9 * sog1_1109[k]
                    + f_3 * pc_y[k] * soh_1552[k];

        t_2070[k] = f_14 * snh_1322[k]
                    + f_3 * pc_y[k] * soh_1553[k];
    }

#pragma omp simd aligned(t_2071, t_2072, t_2073, t_2074, pc_x, pc_y, pc_z, snh_1301, snh_1302, \
                         snh_1323, sog0_1109, sog0_1110, sog1_1109, sog1_1110, soh_1553, \
                         soh_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2071[k] = f_18 * snh_1301[k]
                    + f_1 * sog0_1109[k]
                    - f_2 * sog1_1109[k]
                    + f_3 * pc_z[k] * soh_1553[k];

        t_2072[k] = f_1 * sog0_1110[k]
                    - f_2 * sog1_1110[k]
                    + f_3 * pc_x[k] * soh_1554[k];

        t_2073[k] = f_13 * snh_1323[k]
                    + f_3 * pc_y[k] * soh_1554[k];

        t_2074[k] = f_17 * snh_1302[k]
                    + f_3 * pc_z[k] * soh_1554[k];
    }

#pragma omp simd aligned(t_2075, t_2076, t_2077, pc_x, pc_y, snh_1325, sog0_1113, sog0_1115, \
                         sog1_1113, sog1_1115, soh_1556, soh_1557, \
                         soh_1559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2075[k] = f_4 * sog0_1113[k]
                    - f_5 * sog1_1113[k]
                    + f_3 * pc_x[k] * soh_1557[k];

        t_2076[k] = f_13 * snh_1325[k]
                    + f_3 * pc_y[k] * soh_1556[k];

        t_2077[k] = f_4 * sog0_1115[k]
                    - f_5 * sog1_1115[k]
                    + f_3 * pc_x[k] * soh_1559[k];
    }

#pragma omp simd aligned(t_2078, t_2079, t_2080, pc_x, pc_y, pc_z, snh_1305, snh_1328, \
                         sog0_1116, sog1_1116, soh_1557, soh_1559, \
                         soh_1560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2078[k] = f_6 * sog0_1116[k]
                    - f_7 * sog1_1116[k]
                    + f_3 * pc_x[k] * soh_1560[k];

        t_2079[k] = f_17 * snh_1305[k]
                    + f_3 * pc_z[k] * soh_1557[k];

        t_2080[k] = f_13 * snh_1328[k]
                    + f_3 * pc_y[k] * soh_1559[k];
    }

#pragma omp simd aligned(t_2081, t_2082, t_2083, pc_x, pc_z, snh_1308, sog0_1119, sog0_1120, \
                         sog1_1119, sog1_1120, soh_1560, soh_1563, \
                         soh_1564 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2081[k] = f_6 * sog0_1119[k]
                    - f_7 * sog1_1119[k]
                    + f_3 * pc_x[k] * soh_1563[k];

        t_2082[k] = f_8 * sog0_1120[k]
                    - f_9 * sog1_1120[k]
                    + f_3 * pc_x[k] * soh_1564[k];

        t_2083[k] = f_17 * snh_1308[k]
                    + f_3 * pc_z[k] * soh_1560[k];
    }

#pragma omp simd aligned(t_2084, t_2085, t_2086, t_2087, pc_x, pc_y, snh_1332, sog0_1122, \
                         sog0_1124, sog1_1122, sog1_1124, soh_1563, soh_1566, soh_1568, \
                         soh_1569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2084[k] = f_8 * sog0_1122[k]
                    - f_9 * sog1_1122[k]
                    + f_3 * pc_x[k] * soh_1566[k];

        t_2085[k] = f_13 * snh_1332[k]
                    + f_3 * pc_y[k] * soh_1563[k];

        t_2086[k] = f_8 * sog0_1124[k]
                    - f_9 * sog1_1124[k]
                    + f_3 * pc_x[k] * soh_1568[k];

        t_2087[k] = f_3 * pc_x[k] * soh_1569[k];
    }

#pragma omp simd aligned(t_2088, t_2089, t_2090, t_2091, t_2092, pc_x, soh_1570, soh_1571, \
                         soh_1572, soh_1573, soh_1574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2088[k] = f_3 * pc_x[k] * soh_1570[k];

        t_2089[k] = f_3 * pc_x[k] * soh_1571[k];

        t_2090[k] = f_3 * pc_x[k] * soh_1572[k];

        t_2091[k] = f_3 * pc_x[k] * soh_1573[k];

        t_2092[k] = f_3 * pc_x[k] * soh_1574[k];
    }

#pragma omp simd aligned(t_2093, t_2094, t_2095, pc_y, pc_z, snh_1317, snh_1338, snh_1340, \
                         sog0_1120, sog0_1122, sog1_1120, sog1_1122, soh_1569, \
                         soh_1571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2093[k] = f_13 * snh_1338[k]
                    + f_1 * sog0_1120[k]
                    - f_2 * sog1_1120[k]
                    + f_3 * pc_y[k] * soh_1569[k];

        t_2094[k] = f_17 * snh_1317[k]
                    + f_3 * pc_z[k] * soh_1569[k];

        t_2095[k] = f_13 * snh_1340[k]
                    + f_4 * sog0_1122[k]
                    - f_5 * sog1_1122[k]
                    + f_3 * pc_y[k] * soh_1571[k];
    }

#pragma omp simd aligned(t_2096, t_2097, t_2098, pc_y, snh_1341, snh_1342, snh_1343, \
                         sog0_1123, sog0_1124, sog1_1123, sog1_1124, soh_1572, soh_1573, \
                         soh_1574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2096[k] = f_13 * snh_1341[k]
                    + f_6 * sog0_1123[k]
                    - f_7 * sog1_1123[k]
                    + f_3 * pc_y[k] * soh_1572[k];

        t_2097[k] = f_13 * snh_1342[k]
                    + f_8 * sog0_1124[k]
                    - f_9 * sog1_1124[k]
                    + f_3 * pc_y[k] * soh_1573[k];

        t_2098[k] = f_13 * snh_1343[k]
                    + f_3 * pc_y[k] * soh_1574[k];
    }

#pragma omp simd aligned(t_2099, t_2100, t_2101, t_2102, pc_x, pc_y, pc_z, snh_1322, snh_1323, \
                         snh_1344, sog0_1124, sog0_1125, sog1_1124, sog1_1125, soh_1574, \
                         soh_1575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2099[k] = f_17 * snh_1322[k]
                    + f_1 * sog0_1124[k]
                    - f_2 * sog1_1124[k]
                    + f_3 * pc_z[k] * soh_1574[k];

        t_2100[k] = f_1 * sog0_1125[k]
                    - f_2 * sog1_1125[k]
                    + f_3 * pc_x[k] * soh_1575[k];

        t_2101[k] = f_12 * snh_1344[k]
                    + f_3 * pc_y[k] * soh_1575[k];

        t_2102[k] = f_16 * snh_1323[k]
                    + f_3 * pc_z[k] * soh_1575[k];
    }

#pragma omp simd aligned(t_2103, t_2104, t_2105, pc_x, pc_y, snh_1346, sog0_1128, sog0_1130, \
                         sog1_1128, sog1_1130, soh_1577, soh_1578, \
                         soh_1580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2103[k] = f_4 * sog0_1128[k]
                    - f_5 * sog1_1128[k]
                    + f_3 * pc_x[k] * soh_1578[k];

        t_2104[k] = f_12 * snh_1346[k]
                    + f_3 * pc_y[k] * soh_1577[k];

        t_2105[k] = f_4 * sog0_1130[k]
                    - f_5 * sog1_1130[k]
                    + f_3 * pc_x[k] * soh_1580[k];
    }

#pragma omp simd aligned(t_2106, t_2107, t_2108, pc_x, pc_y, pc_z, snh_1326, snh_1349, \
                         sog0_1131, sog1_1131, soh_1578, soh_1580, \
                         soh_1581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2106[k] = f_6 * sog0_1131[k]
                    - f_7 * sog1_1131[k]
                    + f_3 * pc_x[k] * soh_1581[k];

        t_2107[k] = f_16 * snh_1326[k]
                    + f_3 * pc_z[k] * soh_1578[k];

        t_2108[k] = f_12 * snh_1349[k]
                    + f_3 * pc_y[k] * soh_1580[k];
    }

#pragma omp simd aligned(t_2109, t_2110, t_2111, pc_x, pc_z, snh_1329, sog0_1134, sog0_1135, \
                         sog1_1134, sog1_1135, soh_1581, soh_1584, \
                         soh_1585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2109[k] = f_6 * sog0_1134[k]
                    - f_7 * sog1_1134[k]
                    + f_3 * pc_x[k] * soh_1584[k];

        t_2110[k] = f_8 * sog0_1135[k]
                    - f_9 * sog1_1135[k]
                    + f_3 * pc_x[k] * soh_1585[k];

        t_2111[k] = f_16 * snh_1329[k]
                    + f_3 * pc_z[k] * soh_1581[k];
    }

#pragma omp simd aligned(t_2112, t_2113, t_2114, t_2115, pc_x, pc_y, snh_1353, sog0_1137, \
                         sog0_1139, sog1_1137, sog1_1139, soh_1584, soh_1587, soh_1589, \
                         soh_1590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2112[k] = f_8 * sog0_1137[k]
                    - f_9 * sog1_1137[k]
                    + f_3 * pc_x[k] * soh_1587[k];

        t_2113[k] = f_12 * snh_1353[k]
                    + f_3 * pc_y[k] * soh_1584[k];

        t_2114[k] = f_8 * sog0_1139[k]
                    - f_9 * sog1_1139[k]
                    + f_3 * pc_x[k] * soh_1589[k];

        t_2115[k] = f_3 * pc_x[k] * soh_1590[k];
    }

#pragma omp simd aligned(t_2116, t_2117, t_2118, t_2119, t_2120, pc_x, soh_1591, soh_1592, \
                         soh_1593, soh_1594, soh_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2116[k] = f_3 * pc_x[k] * soh_1591[k];

        t_2117[k] = f_3 * pc_x[k] * soh_1592[k];

        t_2118[k] = f_3 * pc_x[k] * soh_1593[k];

        t_2119[k] = f_3 * pc_x[k] * soh_1594[k];

        t_2120[k] = f_3 * pc_x[k] * soh_1595[k];
    }

#pragma omp simd aligned(t_2121, t_2122, t_2123, pc_y, pc_z, snh_1338, snh_1359, snh_1361, \
                         sog0_1135, sog0_1137, sog1_1135, sog1_1137, soh_1590, \
                         soh_1592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2121[k] = f_12 * snh_1359[k]
                    + f_1 * sog0_1135[k]
                    - f_2 * sog1_1135[k]
                    + f_3 * pc_y[k] * soh_1590[k];

        t_2122[k] = f_16 * snh_1338[k]
                    + f_3 * pc_z[k] * soh_1590[k];

        t_2123[k] = f_12 * snh_1361[k]
                    + f_4 * sog0_1137[k]
                    - f_5 * sog1_1137[k]
                    + f_3 * pc_y[k] * soh_1592[k];
    }

#pragma omp simd aligned(t_2124, t_2125, t_2126, pc_y, snh_1362, snh_1363, snh_1364, \
                         sog0_1138, sog0_1139, sog1_1138, sog1_1139, soh_1593, soh_1594, \
                         soh_1595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2124[k] = f_12 * snh_1362[k]
                    + f_6 * sog0_1138[k]
                    - f_7 * sog1_1138[k]
                    + f_3 * pc_y[k] * soh_1593[k];

        t_2125[k] = f_12 * snh_1363[k]
                    + f_8 * sog0_1139[k]
                    - f_9 * sog1_1139[k]
                    + f_3 * pc_y[k] * soh_1594[k];

        t_2126[k] = f_12 * snh_1364[k]
                    + f_3 * pc_y[k] * soh_1595[k];
    }

#pragma omp simd aligned(t_2127, t_2128, t_2129, t_2130, pb_y, pc_y, pc_z, sni0_1820, \
                         snh_1343, snh_1344, snh_1365, sni1_1820, sog0_1139, sog1_1139, \
                         soh_1595, soh_1596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2127[k] = f_16 * snh_1343[k]
                    + f_1 * sog0_1139[k]
                    - f_2 * sog1_1139[k]
                    + f_3 * pc_z[k] * soh_1595[k];

        t_2128[k] = pb_y[k] * sni0_1820[k]
                    - f_10 * pc_y[k] * sni1_1820[k];

        t_2129[k] = f_11 * snh_1365[k]
                    + f_3 * pc_y[k] * soh_1596[k];

        t_2130[k] = f_15 * snh_1344[k]
                    + f_3 * pc_z[k] * soh_1596[k];
    }

#pragma omp simd aligned(t_2131, t_2132, t_2133, pb_y, pc_x, pc_y, sni0_1825, snh_1367, \
                         sni1_1825, sog0_1143, sog1_1143, soh_1598, \
                         soh_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2131[k] = f_4 * sog0_1143[k]
                    - f_5 * sog1_1143[k]
                    + f_3 * pc_x[k] * soh_1599[k];

        t_2132[k] = f_11 * snh_1367[k]
                    + f_3 * pc_y[k] * soh_1598[k];

        t_2133[k] = pb_y[k] * sni0_1825[k]
                    - f_10 * pc_y[k] * sni1_1825[k];
    }

#pragma omp simd aligned(t_2134, t_2135, t_2136, pc_x, pc_y, pc_z, snh_1347, snh_1370, \
                         sog0_1146, sog1_1146, soh_1599, soh_1601, \
                         soh_1602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2134[k] = f_6 * sog0_1146[k]
                    - f_7 * sog1_1146[k]
                    + f_3 * pc_x[k] * soh_1602[k];

        t_2135[k] = f_15 * snh_1347[k]
                    + f_3 * pc_z[k] * soh_1599[k];

        t_2136[k] = f_11 * snh_1370[k]
                    + f_3 * pc_y[k] * soh_1601[k];
    }

#pragma omp simd aligned(t_2137, t_2138, t_2139, pb_y, pc_x, pc_y, pc_z, sni0_1829, snh_1350, \
                         sni1_1829, sog0_1150, sog1_1150, soh_1602, \
                         soh_1606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2137[k] = pb_y[k] * sni0_1829[k]
                    - f_10 * pc_y[k] * sni1_1829[k];

        t_2138[k] = f_8 * sog0_1150[k]
                    - f_9 * sog1_1150[k]
                    + f_3 * pc_x[k] * soh_1606[k];

        t_2139[k] = f_15 * snh_1350[k]
                    + f_3 * pc_z[k] * soh_1602[k];
    }

#pragma omp simd aligned(t_2140, t_2141, t_2142, t_2143, pb_y, pc_x, pc_y, sni0_1834, \
                         snh_1374, sni1_1834, sog0_1152, sog1_1152, soh_1605, soh_1608, \
                         soh_1611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2140[k] = f_8 * sog0_1152[k]
                    - f_9 * sog1_1152[k]
                    + f_3 * pc_x[k] * soh_1608[k];

        t_2141[k] = f_11 * snh_1374[k]
                    + f_3 * pc_y[k] * soh_1605[k];

        t_2142[k] = pb_y[k] * sni0_1834[k]
                    - f_10 * pc_y[k] * sni1_1834[k];

        t_2143[k] = f_3 * pc_x[k] * soh_1611[k];
    }

#pragma omp simd aligned(t_2144, t_2145, t_2146, t_2147, t_2148, pc_x, soh_1612, soh_1613, \
                         soh_1614, soh_1615, soh_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2144[k] = f_3 * pc_x[k] * soh_1612[k];

        t_2145[k] = f_3 * pc_x[k] * soh_1613[k];

        t_2146[k] = f_3 * pc_x[k] * soh_1614[k];

        t_2147[k] = f_3 * pc_x[k] * soh_1615[k];

        t_2148[k] = f_3 * pc_x[k] * soh_1616[k];
    }

#pragma omp simd aligned(t_2149, t_2150, t_2151, pb_y, pc_y, pc_z, sni0_1841, sni0_1843, \
                         snh_1359, snh_1380, snh_1382, sni1_1841, sni1_1843, \
                         soh_1611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2149[k] = pb_y[k] * sni0_1841[k]
                    + f_19 * snh_1380[k]
                    - f_10 * pc_y[k] * sni1_1841[k];

        t_2150[k] = f_15 * snh_1359[k]
                    + f_3 * pc_z[k] * soh_1611[k];

        t_2151[k] = pb_y[k] * sni0_1843[k]
                    + f_14 * snh_1382[k]
                    - f_10 * pc_y[k] * sni1_1843[k];
    }

#pragma omp simd aligned(t_2152, t_2153, t_2154, t_2155, pb_y, pc_y, sni0_1844, sni0_1845, \
                         sni0_1847, snh_1383, snh_1384, snh_1385, sni1_1844, sni1_1845, \
                         sni1_1847, soh_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2152[k] = pb_y[k] * sni0_1844[k]
                    + f_13 * snh_1383[k]
                    - f_10 * pc_y[k] * sni1_1844[k];

        t_2153[k] = pb_y[k] * sni0_1845[k]
                    + f_12 * snh_1384[k]
                    - f_10 * pc_y[k] * sni1_1845[k];

        t_2154[k] = f_11 * snh_1385[k]
                    + f_3 * pc_y[k] * soh_1616[k];

        t_2155[k] = pb_y[k] * sni0_1847[k]
                    - f_10 * pc_y[k] * sni1_1847[k];
    }

#pragma omp simd aligned(t_2156, t_2157, t_2158, t_2159, t_2160, pc_x, pc_y, pc_z, snh_1365, \
                         sog0_1155, sog0_1158, sog1_1155, sog1_1158, soh_1617, soh_1619, \
                         soh_1620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2156[k] = f_1 * sog0_1155[k]
                    - f_2 * sog1_1155[k]
                    + f_3 * pc_x[k] * soh_1617[k];

        t_2157[k] = f_3 * pc_y[k] * soh_1617[k];

        t_2158[k] = f_0 * snh_1365[k]
                    + f_3 * pc_z[k] * soh_1617[k];

        t_2159[k] = f_4 * sog0_1158[k]
                    - f_5 * sog1_1158[k]
                    + f_3 * pc_x[k] * soh_1620[k];

        t_2160[k] = f_3 * pc_y[k] * soh_1619[k];
    }

#pragma omp simd aligned(t_2161, t_2162, t_2163, t_2164, pc_x, pc_y, pc_z, snh_1368, \
                         sog0_1160, sog0_1161, sog1_1160, sog1_1161, soh_1620, soh_1622, \
                         soh_1623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2161[k] = f_4 * sog0_1160[k]
                    - f_5 * sog1_1160[k]
                    + f_3 * pc_x[k] * soh_1622[k];

        t_2162[k] = f_6 * sog0_1161[k]
                    - f_7 * sog1_1161[k]
                    + f_3 * pc_x[k] * soh_1623[k];

        t_2163[k] = f_0 * snh_1368[k]
                    + f_3 * pc_z[k] * soh_1620[k];

        t_2164[k] = f_3 * pc_y[k] * soh_1622[k];
    }

#pragma omp simd aligned(t_2165, t_2166, t_2167, pc_x, pc_z, snh_1371, sog0_1164, sog0_1165, \
                         sog1_1164, sog1_1165, soh_1623, soh_1626, \
                         soh_1627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2165[k] = f_6 * sog0_1164[k]
                    - f_7 * sog1_1164[k]
                    + f_3 * pc_x[k] * soh_1626[k];

        t_2166[k] = f_8 * sog0_1165[k]
                    - f_9 * sog1_1165[k]
                    + f_3 * pc_x[k] * soh_1627[k];

        t_2167[k] = f_0 * snh_1371[k]
                    + f_3 * pc_z[k] * soh_1623[k];
    }

#pragma omp simd aligned(t_2168, t_2169, t_2170, t_2171, t_2172, pc_x, pc_y, sog0_1167, \
                         sog0_1169, sog1_1167, sog1_1169, soh_1626, soh_1629, soh_1631, \
                         soh_1632, soh_1633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2168[k] = f_8 * sog0_1167[k]
                    - f_9 * sog1_1167[k]
                    + f_3 * pc_x[k] * soh_1629[k];

        t_2169[k] = f_3 * pc_y[k] * soh_1626[k];

        t_2170[k] = f_8 * sog0_1169[k]
                    - f_9 * sog1_1169[k]
                    + f_3 * pc_x[k] * soh_1631[k];

        t_2171[k] = f_3 * pc_x[k] * soh_1632[k];

        t_2172[k] = f_3 * pc_x[k] * soh_1633[k];
    }

#pragma omp simd aligned(t_2173, t_2174, t_2175, t_2176, t_2177, pc_x, pc_y, sog0_1165, \
                         sog1_1165, soh_1632, soh_1634, soh_1635, soh_1636, \
                         soh_1637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2173[k] = f_3 * pc_x[k] * soh_1634[k];

        t_2174[k] = f_3 * pc_x[k] * soh_1635[k];

        t_2175[k] = f_3 * pc_x[k] * soh_1636[k];

        t_2176[k] = f_3 * pc_x[k] * soh_1637[k];

        t_2177[k] = f_1 * sog0_1165[k]
                    - f_2 * sog1_1165[k]
                    + f_3 * pc_y[k] * soh_1632[k];
    }

#pragma omp simd aligned(t_2178, t_2179, t_2180, pc_y, pc_z, snh_1380, sog0_1167, sog0_1168, \
                         sog1_1167, sog1_1168, soh_1632, soh_1634, \
                         soh_1635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2178[k] = f_0 * snh_1380[k]
                    + f_3 * pc_z[k] * soh_1632[k];

        t_2179[k] = f_4 * sog0_1167[k]
                    - f_5 * sog1_1167[k]
                    + f_3 * pc_y[k] * soh_1634[k];

        t_2180[k] = f_6 * sog0_1168[k]
                    - f_7 * sog1_1168[k]
                    + f_3 * pc_y[k] * soh_1635[k];
    }

#pragma omp simd aligned(t_2181, t_2182, t_2183, pc_y, pc_z, snh_1385, sog0_1169, sog1_1169, \
                         soh_1636, soh_1637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2181[k] = f_8 * sog0_1169[k]
                    - f_9 * sog1_1169[k]
                    + f_3 * pc_y[k] * soh_1636[k];

        t_2182[k] = f_3 * pc_y[k] * soh_1637[k];

        t_2183[k] = f_0 * snh_1385[k]
                    + f_1 * sog0_1169[k]
                    - f_2 * sog1_1169[k]
                    + f_3 * pc_z[k] * soh_1637[k];
    }
}

auto
compute_prim_soi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sni0, const size_t snh,
                                                   const size_t sni1, const size_t sog0,
                                                   const size_t sog1, const size_t soh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_soi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, sni0, snh,
                                                              sni1, sog0, sog1, soh, ncols,
                                                              gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, sni0,
                                                               snh, sni1, sog0, sog1, soh,
                                                               ncols, gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, sni0,
                                                               snh, sni1, sog0, sog1, soh,
                                                               ncols, gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece12(buffer, target, pc, snh, sog0,
                                                               sog1, soh, ncols, gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, sni0,
                                                               snh, sni1, sog0, sog1, soh,
                                                               ncols, gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece14(buffer, target, pb, pc, sni0,
                                                               snh, sni1, soh, ncols, gamma, p,
                                                               q);

    compute_prim_soi_three_center_electron_repulsion_0_piece15(buffer, target, pb, pc, sni0,
                                                               snh, sni1, soh, ncols, gamma, p,
                                                               q);

    compute_prim_soi_three_center_electron_repulsion_0_piece16(buffer, target, pb, pc, sni0,
                                                               snh, sni1, sog0, sog1, soh,
                                                               ncols, gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece17(buffer, target, pc, snh, sog0,
                                                               sog1, soh, ncols, gamma, p, q);

    compute_prim_soi_three_center_electron_repulsion_0_piece18(buffer, target, pb, pc, sni0,
                                                               snh, sni1, sog0, sog1, soh,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
