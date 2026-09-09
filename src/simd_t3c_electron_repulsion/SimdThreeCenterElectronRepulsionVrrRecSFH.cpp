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


#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sfh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdh0,
                                                          const size_t sdg, const size_t sdh1,
                                                          const size_t sff0, const size_t sff1,
                                                          const size_t sfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;

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
    auto *t_128 = buffer.data(target + 128);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh0_0 = buffer.data(sdh0 + 0);
    const auto *sdh0_3 = buffer.data(sdh0 + 3);
    const auto *sdh0_5 = buffer.data(sdh0 + 5);
    const auto *sdh0_6 = buffer.data(sdh0 + 6);
    const auto *sdh0_9 = buffer.data(sdh0 + 9);
    const auto *sdh0_15 = buffer.data(sdh0 + 15);
    const auto *sdh0_20 = buffer.data(sdh0 + 20);
    const auto *sdh0_24 = buffer.data(sdh0 + 24);
    const auto *sdh0_27 = buffer.data(sdh0 + 27);
    const auto *sdh0_42 = buffer.data(sdh0 + 42);
    const auto *sdh0_47 = buffer.data(sdh0 + 47);
    const auto *sdh0_51 = buffer.data(sdh0 + 51);
    const auto *sdh0_63 = buffer.data(sdh0 + 63);
    const auto *sdh0_66 = buffer.data(sdh0 + 66);
    const auto *sdh0_68 = buffer.data(sdh0 + 68);
    const auto *sdh0_69 = buffer.data(sdh0 + 69);
    const auto *sdh0_72 = buffer.data(sdh0 + 72);
    const auto *sdh0_78 = buffer.data(sdh0 + 78);
    const auto *sdh0_80 = buffer.data(sdh0 + 80);
    const auto *sdh0_81 = buffer.data(sdh0 + 81);
    const auto *sdh0_83 = buffer.data(sdh0 + 83);
    const auto *sdh0_99 = buffer.data(sdh0 + 99);
    const auto *sdh0_101 = buffer.data(sdh0 + 101);
    const auto *sdh0_102 = buffer.data(sdh0 + 102);
    const auto *sdh0_104 = buffer.data(sdh0 + 104);
    const auto *sdh0_105 = buffer.data(sdh0 + 105);
    const auto *sdh0_108 = buffer.data(sdh0 + 108);
    const auto *sdh0_110 = buffer.data(sdh0 + 110);
    const auto *sdh0_111 = buffer.data(sdh0 + 111);
    const auto *sdh0_114 = buffer.data(sdh0 + 114);
    const auto *sdh0_120 = buffer.data(sdh0 + 120);
    const auto *sdh0_122 = buffer.data(sdh0 + 122);
    const auto *sdh0_123 = buffer.data(sdh0 + 123);
    const auto *sdh0_125 = buffer.data(sdh0 + 125);

    const auto *sdg_0 = buffer.data(sdg + 0);
    const auto *sdg_1 = buffer.data(sdg + 1);
    const auto *sdg_2 = buffer.data(sdg + 2);
    const auto *sdg_3 = buffer.data(sdg + 3);
    const auto *sdg_5 = buffer.data(sdg + 5);
    const auto *sdg_6 = buffer.data(sdg + 6);
    const auto *sdg_9 = buffer.data(sdg + 9);
    const auto *sdg_10 = buffer.data(sdg + 10);
    const auto *sdg_11 = buffer.data(sdg + 11);
    const auto *sdg_12 = buffer.data(sdg + 12);
    const auto *sdg_13 = buffer.data(sdg + 13);
    const auto *sdg_14 = buffer.data(sdg + 14);
    const auto *sdg_15 = buffer.data(sdg + 15);
    const auto *sdg_17 = buffer.data(sdg + 17);
    const auto *sdg_18 = buffer.data(sdg + 18);
    const auto *sdg_20 = buffer.data(sdg + 20);
    const auto *sdg_25 = buffer.data(sdg + 25);
    const auto *sdg_26 = buffer.data(sdg + 26);
    const auto *sdg_27 = buffer.data(sdg + 27);
    const auto *sdg_28 = buffer.data(sdg + 28);
    const auto *sdg_29 = buffer.data(sdg + 29);
    const auto *sdg_30 = buffer.data(sdg + 30);
    const auto *sdg_32 = buffer.data(sdg + 32);
    const auto *sdg_33 = buffer.data(sdg + 33);
    const auto *sdg_35 = buffer.data(sdg + 35);
    const auto *sdg_40 = buffer.data(sdg + 40);
    const auto *sdg_41 = buffer.data(sdg + 41);
    const auto *sdg_42 = buffer.data(sdg + 42);
    const auto *sdg_43 = buffer.data(sdg + 43);
    const auto *sdg_44 = buffer.data(sdg + 44);
    const auto *sdg_45 = buffer.data(sdg + 45);
    const auto *sdg_48 = buffer.data(sdg + 48);
    const auto *sdg_50 = buffer.data(sdg + 50);
    const auto *sdg_51 = buffer.data(sdg + 51);
    const auto *sdg_54 = buffer.data(sdg + 54);
    const auto *sdg_55 = buffer.data(sdg + 55);
    const auto *sdg_56 = buffer.data(sdg + 56);
    const auto *sdg_57 = buffer.data(sdg + 57);
    const auto *sdg_58 = buffer.data(sdg + 58);
    const auto *sdg_59 = buffer.data(sdg + 59);
    const auto *sdg_70 = buffer.data(sdg + 70);
    const auto *sdg_71 = buffer.data(sdg + 71);
    const auto *sdg_72 = buffer.data(sdg + 72);
    const auto *sdg_73 = buffer.data(sdg + 73);
    const auto *sdg_74 = buffer.data(sdg + 74);
    const auto *sdg_75 = buffer.data(sdg + 75);
    const auto *sdg_78 = buffer.data(sdg + 78);
    const auto *sdg_80 = buffer.data(sdg + 80);
    const auto *sdg_81 = buffer.data(sdg + 81);
    const auto *sdg_84 = buffer.data(sdg + 84);
    const auto *sdg_85 = buffer.data(sdg + 85);
    const auto *sdg_86 = buffer.data(sdg + 86);
    const auto *sdg_87 = buffer.data(sdg + 87);
    const auto *sdg_88 = buffer.data(sdg + 88);
    const auto *sdg_89 = buffer.data(sdg + 89);

    const auto *sdh1_0 = buffer.data(sdh1 + 0);
    const auto *sdh1_3 = buffer.data(sdh1 + 3);
    const auto *sdh1_5 = buffer.data(sdh1 + 5);
    const auto *sdh1_6 = buffer.data(sdh1 + 6);
    const auto *sdh1_9 = buffer.data(sdh1 + 9);
    const auto *sdh1_15 = buffer.data(sdh1 + 15);
    const auto *sdh1_20 = buffer.data(sdh1 + 20);
    const auto *sdh1_24 = buffer.data(sdh1 + 24);
    const auto *sdh1_27 = buffer.data(sdh1 + 27);
    const auto *sdh1_42 = buffer.data(sdh1 + 42);
    const auto *sdh1_47 = buffer.data(sdh1 + 47);
    const auto *sdh1_51 = buffer.data(sdh1 + 51);
    const auto *sdh1_63 = buffer.data(sdh1 + 63);
    const auto *sdh1_66 = buffer.data(sdh1 + 66);
    const auto *sdh1_68 = buffer.data(sdh1 + 68);
    const auto *sdh1_69 = buffer.data(sdh1 + 69);
    const auto *sdh1_72 = buffer.data(sdh1 + 72);
    const auto *sdh1_78 = buffer.data(sdh1 + 78);
    const auto *sdh1_80 = buffer.data(sdh1 + 80);
    const auto *sdh1_81 = buffer.data(sdh1 + 81);
    const auto *sdh1_83 = buffer.data(sdh1 + 83);
    const auto *sdh1_99 = buffer.data(sdh1 + 99);
    const auto *sdh1_101 = buffer.data(sdh1 + 101);
    const auto *sdh1_102 = buffer.data(sdh1 + 102);
    const auto *sdh1_104 = buffer.data(sdh1 + 104);
    const auto *sdh1_105 = buffer.data(sdh1 + 105);
    const auto *sdh1_108 = buffer.data(sdh1 + 108);
    const auto *sdh1_110 = buffer.data(sdh1 + 110);
    const auto *sdh1_111 = buffer.data(sdh1 + 111);
    const auto *sdh1_114 = buffer.data(sdh1 + 114);
    const auto *sdh1_120 = buffer.data(sdh1 + 120);
    const auto *sdh1_122 = buffer.data(sdh1 + 122);
    const auto *sdh1_123 = buffer.data(sdh1 + 123);
    const auto *sdh1_125 = buffer.data(sdh1 + 125);

    const auto *sff0_0 = buffer.data(sff0 + 0);
    const auto *sff0_3 = buffer.data(sff0 + 3);
    const auto *sff0_5 = buffer.data(sff0 + 5);
    const auto *sff0_6 = buffer.data(sff0 + 6);
    const auto *sff0_8 = buffer.data(sff0 + 8);
    const auto *sff0_9 = buffer.data(sff0 + 9);
    const auto *sff0_16 = buffer.data(sff0 + 16);
    const auto *sff0_18 = buffer.data(sff0 + 18);
    const auto *sff0_19 = buffer.data(sff0 + 19);
    const auto *sff0_28 = buffer.data(sff0 + 28);
    const auto *sff0_29 = buffer.data(sff0 + 29);
    const auto *sff0_60 = buffer.data(sff0 + 60);

    const auto *sff1_0 = buffer.data(sff1 + 0);
    const auto *sff1_3 = buffer.data(sff1 + 3);
    const auto *sff1_5 = buffer.data(sff1 + 5);
    const auto *sff1_6 = buffer.data(sff1 + 6);
    const auto *sff1_8 = buffer.data(sff1 + 8);
    const auto *sff1_9 = buffer.data(sff1 + 9);
    const auto *sff1_16 = buffer.data(sff1 + 16);
    const auto *sff1_18 = buffer.data(sff1 + 18);
    const auto *sff1_19 = buffer.data(sff1 + 19);
    const auto *sff1_28 = buffer.data(sff1 + 28);
    const auto *sff1_29 = buffer.data(sff1 + 29);
    const auto *sff1_60 = buffer.data(sff1 + 60);

    const auto *sfg_0 = buffer.data(sfg + 0);
    const auto *sfg_2 = buffer.data(sfg + 2);
    const auto *sfg_3 = buffer.data(sfg + 3);
    const auto *sfg_5 = buffer.data(sfg + 5);
    const auto *sfg_6 = buffer.data(sfg + 6);
    const auto *sfg_9 = buffer.data(sfg + 9);
    const auto *sfg_10 = buffer.data(sfg + 10);
    const auto *sfg_11 = buffer.data(sfg + 11);
    const auto *sfg_12 = buffer.data(sfg + 12);
    const auto *sfg_13 = buffer.data(sfg + 13);
    const auto *sfg_14 = buffer.data(sfg + 14);
    const auto *sfg_15 = buffer.data(sfg + 15);
    const auto *sfg_17 = buffer.data(sfg + 17);
    const auto *sfg_18 = buffer.data(sfg + 18);
    const auto *sfg_20 = buffer.data(sfg + 20);
    const auto *sfg_25 = buffer.data(sfg + 25);
    const auto *sfg_26 = buffer.data(sfg + 26);
    const auto *sfg_27 = buffer.data(sfg + 27);
    const auto *sfg_28 = buffer.data(sfg + 28);
    const auto *sfg_29 = buffer.data(sfg + 29);
    const auto *sfg_30 = buffer.data(sfg + 30);
    const auto *sfg_32 = buffer.data(sfg + 32);
    const auto *sfg_33 = buffer.data(sfg + 33);
    const auto *sfg_35 = buffer.data(sfg + 35);
    const auto *sfg_40 = buffer.data(sfg + 40);
    const auto *sfg_41 = buffer.data(sfg + 41);
    const auto *sfg_42 = buffer.data(sfg + 42);
    const auto *sfg_43 = buffer.data(sfg + 43);
    const auto *sfg_44 = buffer.data(sfg + 44);
    const auto *sfg_45 = buffer.data(sfg + 45);
    const auto *sfg_47 = buffer.data(sfg + 47);
    const auto *sfg_48 = buffer.data(sfg + 48);
    const auto *sfg_50 = buffer.data(sfg + 50);
    const auto *sfg_55 = buffer.data(sfg + 55);
    const auto *sfg_56 = buffer.data(sfg + 56);
    const auto *sfg_57 = buffer.data(sfg + 57);
    const auto *sfg_58 = buffer.data(sfg + 58);
    const auto *sfg_59 = buffer.data(sfg + 59);
    const auto *sfg_60 = buffer.data(sfg + 60);
    const auto *sfg_62 = buffer.data(sfg + 62);
    const auto *sfg_63 = buffer.data(sfg + 63);
    const auto *sfg_65 = buffer.data(sfg + 65);
    const auto *sfg_70 = buffer.data(sfg + 70);
    const auto *sfg_71 = buffer.data(sfg + 71);
    const auto *sfg_72 = buffer.data(sfg + 72);
    const auto *sfg_73 = buffer.data(sfg + 73);
    const auto *sfg_74 = buffer.data(sfg + 74);
    const auto *sfg_75 = buffer.data(sfg + 75);
    const auto *sfg_77 = buffer.data(sfg + 77);
    const auto *sfg_78 = buffer.data(sfg + 78);
    const auto *sfg_80 = buffer.data(sfg + 80);
    const auto *sfg_85 = buffer.data(sfg + 85);
    const auto *sfg_86 = buffer.data(sfg + 86);
    const auto *sfg_87 = buffer.data(sfg + 87);
    const auto *sfg_88 = buffer.data(sfg + 88);
    const auto *sfg_89 = buffer.data(sfg + 89);
    const auto *sfg_90 = buffer.data(sfg + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sdg_0, sdg_3, sff0_0, sff0_3, \
                         sff1_0, sff1_3, sfg_0, sfg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdg_0[k]
                 + f_1 * sff0_0[k]
                 - f_2 * sff1_0[k]
                 + f_3 * pc_x[k] * sfg_0[k];

        t_1[k] = f_3 * pc_y[k] * sfg_0[k];

        t_2[k] = f_3 * pc_z[k] * sfg_0[k];

        t_3[k] = f_0 * sdg_3[k]
                 + f_4 * sff0_3[k]
                 - f_5 * sff1_3[k]
                 + f_3 * pc_x[k] * sfg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sdg_5, sdg_6, sff0_5, sff0_6, sff1_5, \
                         sff1_6, sfg_2, sfg_5, sfg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sfg_2[k];

        t_5[k] = f_0 * sdg_5[k]
                 + f_4 * sff0_5[k]
                 - f_5 * sff1_5[k]
                 + f_3 * pc_x[k] * sfg_5[k];

        t_6[k] = f_0 * sdg_6[k]
                 + f_6 * sff0_6[k]
                 - f_7 * sff1_6[k]
                 + f_3 * pc_x[k] * sfg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, sdg_9, sdg_10, sff0_9, sff1_9, \
                         sfg_3, sfg_5, sfg_9, sfg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sfg_3[k];

        t_8[k] = f_3 * pc_y[k] * sfg_5[k];

        t_9[k] = f_0 * sdg_9[k]
                 + f_6 * sff0_9[k]
                 - f_7 * sff1_9[k]
                 + f_3 * pc_x[k] * sfg_9[k];

        t_10[k] = f_0 * sdg_10[k]
                  + f_3 * pc_x[k] * sfg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, sdg_11, sdg_12, sdg_13, sdg_14, sfg_11, \
                         sfg_12, sfg_13, sfg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sdg_11[k]
                  + f_3 * pc_x[k] * sfg_11[k];

        t_12[k] = f_0 * sdg_12[k]
                  + f_3 * pc_x[k] * sfg_12[k];

        t_13[k] = f_0 * sdg_13[k]
                  + f_3 * pc_x[k] * sfg_13[k];

        t_14[k] = f_0 * sdg_14[k]
                  + f_3 * pc_x[k] * sfg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, sff0_6, sff0_8, sff0_9, sff1_6, \
                         sff1_8, sff1_9, sfg_10, sfg_12, sfg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * sff0_6[k]
                  - f_2 * sff1_6[k]
                  + f_3 * pc_y[k] * sfg_10[k];

        t_16[k] = f_3 * pc_z[k] * sfg_10[k];

        t_17[k] = f_4 * sff0_8[k]
                  - f_5 * sff1_8[k]
                  + f_3 * pc_y[k] * sfg_12[k];

        t_18[k] = f_6 * sff0_9[k]
                  - f_7 * sff1_9[k]
                  + f_3 * pc_y[k] * sfg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, sdh0_0, sdg_0, \
                         sdh1_0, sff0_9, sff1_9, sfg_14, sfg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sfg_14[k];

        t_20[k] = f_1 * sff0_9[k]
                  - f_2 * sff1_9[k]
                  + f_3 * pc_z[k] * sfg_14[k];

        t_21[k] = pb_y[k] * sdh0_0[k]
                  - f_8 * pc_y[k] * sdh1_0[k];

        t_22[k] = f_9 * sdg_0[k]
                  + f_3 * pc_y[k] * sfg_15[k];

        t_23[k] = f_3 * pc_z[k] * sfg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, sdh0_3, sdh0_5, sdh0_6, sdg_1, \
                         sdg_2, sdg_3, sdh1_3, sdh1_5, sdh1_6, sfg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sdh0_3[k]
                  + f_10 * sdg_1[k]
                  - f_8 * pc_y[k] * sdh1_3[k];

        t_25[k] = f_9 * sdg_2[k]
                  + f_3 * pc_y[k] * sfg_17[k];

        t_26[k] = pb_y[k] * sdh0_5[k]
                  - f_8 * pc_y[k] * sdh1_5[k];

        t_27[k] = pb_y[k] * sdh0_6[k]
                  + f_0 * sdg_3[k]
                  - f_8 * pc_y[k] * sdh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, sdh0_9, sdg_5, \
                         sdg_25, sdh1_9, sfg_18, sfg_20, sfg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * sfg_18[k];

        t_29[k] = f_9 * sdg_5[k]
                  + f_3 * pc_y[k] * sfg_20[k];

        t_30[k] = pb_y[k] * sdh0_9[k]
                  - f_8 * pc_y[k] * sdh1_9[k];

        t_31[k] = f_10 * sdg_25[k]
                  + f_3 * pc_x[k] * sfg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, sdg_26, sdg_27, sdg_28, sdg_29, sfg_26, \
                         sfg_27, sfg_28, sfg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_10 * sdg_26[k]
                  + f_3 * pc_x[k] * sfg_26[k];

        t_33[k] = f_10 * sdg_27[k]
                  + f_3 * pc_x[k] * sfg_27[k];

        t_34[k] = f_10 * sdg_28[k]
                  + f_3 * pc_x[k] * sfg_28[k];

        t_35[k] = f_10 * sdg_29[k]
                  + f_3 * pc_x[k] * sfg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, sdg_10, sdg_12, sff0_16, sff0_18, \
                         sff1_16, sff1_18, sfg_25, sfg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * sdg_10[k]
                  + f_1 * sff0_16[k]
                  - f_2 * sff1_16[k]
                  + f_3 * pc_y[k] * sfg_25[k];

        t_37[k] = f_3 * pc_z[k] * sfg_25[k];

        t_38[k] = f_9 * sdg_12[k]
                  + f_4 * sff0_18[k]
                  - f_5 * sff1_18[k]
                  + f_3 * pc_y[k] * sfg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, sdh0_20, sdg_13, sdg_14, sdh1_20, \
                         sff0_19, sff1_19, sfg_28, sfg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * sdg_13[k]
                  + f_6 * sff0_19[k]
                  - f_7 * sff1_19[k]
                  + f_3 * pc_y[k] * sfg_28[k];

        t_40[k] = f_9 * sdg_14[k]
                  + f_3 * pc_y[k] * sfg_29[k];

        t_41[k] = pb_y[k] * sdh0_20[k]
                  - f_8 * pc_y[k] * sdh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, sdh0_0, sdh0_3, \
                         sdg_0, sdh1_0, sdh1_3, sfg_30, sfg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sdh0_0[k]
                  - f_8 * pc_z[k] * sdh1_0[k];

        t_43[k] = f_3 * pc_y[k] * sfg_30[k];

        t_44[k] = f_9 * sdg_0[k]
                  + f_3 * pc_z[k] * sfg_30[k];

        t_45[k] = pb_z[k] * sdh0_3[k]
                  - f_8 * pc_z[k] * sdh1_3[k];

        t_46[k] = f_3 * pc_y[k] * sfg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, sdh0_5, sdh0_6, sdg_2, \
                         sdg_3, sdh1_5, sdh1_6, sfg_33, sfg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * sdh0_5[k]
                  + f_10 * sdg_2[k]
                  - f_8 * pc_z[k] * sdh1_5[k];

        t_48[k] = pb_z[k] * sdh0_6[k]
                  - f_8 * pc_z[k] * sdh1_6[k];

        t_49[k] = f_9 * sdg_3[k]
                  + f_3 * pc_z[k] * sfg_33[k];

        t_50[k] = f_3 * pc_y[k] * sfg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, sdh0_9, sdg_5, sdg_40, \
                         sdg_41, sdg_42, sdh1_9, sfg_40, sfg_41, \
                         sfg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * sdh0_9[k]
                  + f_0 * sdg_5[k]
                  - f_8 * pc_z[k] * sdh1_9[k];

        t_52[k] = f_10 * sdg_40[k]
                  + f_3 * pc_x[k] * sfg_40[k];

        t_53[k] = f_10 * sdg_41[k]
                  + f_3 * pc_x[k] * sfg_41[k];

        t_54[k] = f_10 * sdg_42[k]
                  + f_3 * pc_x[k] * sfg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, sdh0_15, sdg_10, sdg_43, \
                         sdg_44, sdh1_15, sfg_40, sfg_43, sfg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_10 * sdg_43[k]
                  + f_3 * pc_x[k] * sfg_43[k];

        t_56[k] = f_10 * sdg_44[k]
                  + f_3 * pc_x[k] * sfg_44[k];

        t_57[k] = pb_z[k] * sdh0_15[k]
                  - f_8 * pc_z[k] * sdh1_15[k];

        t_58[k] = f_9 * sdg_10[k]
                  + f_3 * pc_z[k] * sfg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, sdg_14, sff0_28, sff0_29, \
                         sff1_28, sff1_29, sfg_42, sfg_43, sfg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * sff0_28[k]
                  - f_5 * sff1_28[k]
                  + f_3 * pc_y[k] * sfg_42[k];

        t_60[k] = f_6 * sff0_29[k]
                  - f_7 * sff1_29[k]
                  + f_3 * pc_y[k] * sfg_43[k];

        t_61[k] = f_3 * pc_y[k] * sfg_44[k];

        t_62[k] = f_9 * sdg_14[k]
                  + f_1 * sff0_29[k]
                  - f_2 * sff1_29[k]
                  + f_3 * pc_z[k] * sfg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_x, pc_x, pc_y, pc_z, sdh0_63, sdh0_66, \
                         sdg_15, sdg_45, sdg_48, sdh1_63, sdh1_66, \
                         sfg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_x[k] * sdh0_63[k]
                  + f_11 * sdg_45[k]
                  - f_8 * pc_x[k] * sdh1_63[k];

        t_64[k] = f_10 * sdg_15[k]
                  + f_3 * pc_y[k] * sfg_45[k];

        t_65[k] = f_3 * pc_z[k] * sfg_45[k];

        t_66[k] = pb_x[k] * sdh0_66[k]
                  + f_0 * sdg_48[k]
                  - f_8 * pc_x[k] * sdh1_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pc_x, pc_y, sdh0_68, sdh0_69, sdg_17, sdg_50, \
                         sdg_51, sdh1_68, sdh1_69, sfg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * sdg_17[k]
                  + f_3 * pc_y[k] * sfg_47[k];

        t_68[k] = pb_x[k] * sdh0_68[k]
                  + f_0 * sdg_50[k]
                  - f_8 * pc_x[k] * sdh1_68[k];

        t_69[k] = pb_x[k] * sdh0_69[k]
                  + f_10 * sdg_51[k]
                  - f_8 * pc_x[k] * sdh1_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pc_x, pc_y, pc_z, sdh0_72, sdg_20, \
                         sdg_54, sdg_55, sdh1_72, sfg_48, sfg_50, \
                         sfg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * sfg_48[k];

        t_71[k] = f_10 * sdg_20[k]
                  + f_3 * pc_y[k] * sfg_50[k];

        t_72[k] = pb_x[k] * sdh0_72[k]
                  + f_10 * sdg_54[k]
                  - f_8 * pc_x[k] * sdh1_72[k];

        t_73[k] = f_9 * sdg_55[k]
                  + f_3 * pc_x[k] * sfg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, sdg_56, sdg_57, sdg_58, sdg_59, sfg_56, \
                         sfg_57, sfg_58, sfg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_9 * sdg_56[k]
                  + f_3 * pc_x[k] * sfg_56[k];

        t_75[k] = f_9 * sdg_57[k]
                  + f_3 * pc_x[k] * sfg_57[k];

        t_76[k] = f_9 * sdg_58[k]
                  + f_3 * pc_x[k] * sfg_58[k];

        t_77[k] = f_9 * sdg_59[k]
                  + f_3 * pc_x[k] * sfg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pc_x, pc_z, sdh0_78, sdh0_80, sdh0_81, \
                         sdh1_78, sdh1_80, sdh1_81, sfg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_x[k] * sdh0_78[k]
                  - f_8 * pc_x[k] * sdh1_78[k];

        t_79[k] = f_3 * pc_z[k] * sfg_55[k];

        t_80[k] = pb_x[k] * sdh0_80[k]
                  - f_8 * pc_x[k] * sdh1_80[k];

        t_81[k] = pb_x[k] * sdh0_81[k]
                  - f_8 * pc_x[k] * sdh1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_x, pb_y, pc_x, pc_y, sdh0_42, sdh0_83, \
                         sdg_29, sdg_30, sdh1_42, sdh1_83, sfg_59, \
                         sfg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_10 * sdg_29[k]
                  + f_3 * pc_y[k] * sfg_59[k];

        t_83[k] = pb_x[k] * sdh0_83[k]
                  - f_8 * pc_x[k] * sdh1_83[k];

        t_84[k] = pb_y[k] * sdh0_42[k]
                  - f_8 * pc_y[k] * sdh1_42[k];

        t_85[k] = f_9 * sdg_30[k]
                  + f_3 * pc_y[k] * sfg_60[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_y, pb_z, pc_y, pc_z, sdh0_24, sdh0_47, \
                         sdg_15, sdg_32, sdh1_24, sdh1_47, sfg_60, \
                         sfg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_9 * sdg_15[k]
                  + f_3 * pc_z[k] * sfg_60[k];

        t_87[k] = pb_z[k] * sdh0_24[k]
                  - f_8 * pc_z[k] * sdh1_24[k];

        t_88[k] = f_9 * sdg_32[k]
                  + f_3 * pc_y[k] * sfg_62[k];

        t_89[k] = pb_y[k] * sdh0_47[k]
                  - f_8 * pc_y[k] * sdh1_47[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_y, pb_z, pc_y, pc_z, sdh0_27, sdh0_51, \
                         sdg_18, sdg_35, sdh1_27, sdh1_51, sfg_63, \
                         sfg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * sdh0_27[k]
                  - f_8 * pc_z[k] * sdh1_27[k];

        t_91[k] = f_9 * sdg_18[k]
                  + f_3 * pc_z[k] * sfg_63[k];

        t_92[k] = f_9 * sdg_35[k]
                  + f_3 * pc_y[k] * sfg_65[k];

        t_93[k] = pb_y[k] * sdh0_51[k]
                  - f_8 * pc_y[k] * sdh1_51[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, sdg_70, sdg_71, sdg_72, sdg_73, \
                         sdg_74, sfg_70, sfg_71, sfg_72, sfg_73, \
                         sfg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * sdg_70[k]
                  + f_3 * pc_x[k] * sfg_70[k];

        t_95[k] = f_9 * sdg_71[k]
                  + f_3 * pc_x[k] * sfg_71[k];

        t_96[k] = f_9 * sdg_72[k]
                  + f_3 * pc_x[k] * sfg_72[k];

        t_97[k] = f_9 * sdg_73[k]
                  + f_3 * pc_x[k] * sfg_73[k];

        t_98[k] = f_9 * sdg_74[k]
                  + f_3 * pc_x[k] * sfg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pb_x, pc_x, pc_z, sdh0_99, sdh0_101, \
                         sdh0_102, sdg_25, sdh1_99, sdh1_101, sdh1_102, \
                         sfg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pb_x[k] * sdh0_99[k]
                  - f_8 * pc_x[k] * sdh1_99[k];

        t_100[k] = f_9 * sdg_25[k]
                   + f_3 * pc_z[k] * sfg_70[k];

        t_101[k] = pb_x[k] * sdh0_101[k]
                   - f_8 * pc_x[k] * sdh1_101[k];

        t_102[k] = pb_x[k] * sdh0_102[k]
                   - f_8 * pc_x[k] * sdh1_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pb_x, pc_x, pc_y, sdh0_104, sdh0_105, \
                         sdg_44, sdg_75, sdh1_104, sdh1_105, sfg_74, \
                         sfg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_9 * sdg_44[k]
                   + f_3 * pc_y[k] * sfg_74[k];

        t_104[k] = pb_x[k] * sdh0_104[k]
                   - f_8 * pc_x[k] * sdh1_104[k];

        t_105[k] = pb_x[k] * sdh0_105[k]
                   + f_11 * sdg_75[k]
                   - f_8 * pc_x[k] * sdh1_105[k];

        t_106[k] = f_3 * pc_y[k] * sfg_75[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, pc_x, pc_y, pc_z, sdh0_108, sdg_30, \
                         sdg_78, sdh1_108, sfg_75, sfg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_10 * sdg_30[k]
                   + f_3 * pc_z[k] * sfg_75[k];

        t_108[k] = pb_x[k] * sdh0_108[k]
                   + f_0 * sdg_78[k]
                   - f_8 * pc_x[k] * sdh1_108[k];

        t_109[k] = f_3 * pc_y[k] * sfg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, pc_x, pc_z, sdh0_110, sdh0_111, sdg_33, \
                         sdg_80, sdg_81, sdh1_110, sdh1_111, sfg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * sdh0_110[k]
                   + f_0 * sdg_80[k]
                   - f_8 * pc_x[k] * sdh1_110[k];

        t_111[k] = pb_x[k] * sdh0_111[k]
                   + f_10 * sdg_81[k]
                   - f_8 * pc_x[k] * sdh1_111[k];

        t_112[k] = f_10 * sdg_33[k]
                   + f_3 * pc_z[k] * sfg_78[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pc_x, pc_y, sdh0_114, sdg_84, \
                         sdg_85, sdg_86, sdh1_114, sfg_80, sfg_85, \
                         sfg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * pc_y[k] * sfg_80[k];

        t_114[k] = pb_x[k] * sdh0_114[k]
                   + f_10 * sdg_84[k]
                   - f_8 * pc_x[k] * sdh1_114[k];

        t_115[k] = f_9 * sdg_85[k]
                   + f_3 * pc_x[k] * sfg_85[k];

        t_116[k] = f_9 * sdg_86[k]
                   + f_3 * pc_x[k] * sfg_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pc_x, sdh0_120, sdg_87, sdg_88, \
                         sdg_89, sdh1_120, sfg_87, sfg_88, sfg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_9 * sdg_87[k]
                   + f_3 * pc_x[k] * sfg_87[k];

        t_118[k] = f_9 * sdg_88[k]
                   + f_3 * pc_x[k] * sfg_88[k];

        t_119[k] = f_9 * sdg_89[k]
                   + f_3 * pc_x[k] * sfg_89[k];

        t_120[k] = pb_x[k] * sdh0_120[k]
                   - f_8 * pc_x[k] * sdh1_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pc_x, pc_y, pc_z, sdh0_122, \
                         sdh0_123, sdg_40, sdh1_122, sdh1_123, sfg_85, \
                         sfg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * sdg_40[k]
                   + f_3 * pc_z[k] * sfg_85[k];

        t_122[k] = pb_x[k] * sdh0_122[k]
                   - f_8 * pc_x[k] * sdh1_122[k];

        t_123[k] = pb_x[k] * sdh0_123[k]
                   - f_8 * pc_x[k] * sdh1_123[k];

        t_124[k] = f_3 * pc_y[k] * sfg_89[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_x, pc_x, pc_y, pc_z, sdh0_125, sdg_45, \
                         sdh1_125, sff0_60, sff1_60, sfg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_x[k] * sdh0_125[k]
                   - f_8 * pc_x[k] * sdh1_125[k];

        t_126[k] = f_1 * sff0_60[k]
                   - f_2 * sff1_60[k]
                   + f_3 * pc_x[k] * sfg_90[k];

        t_127[k] = f_0 * sdg_45[k]
                   + f_3 * pc_y[k] * sfg_90[k];

        t_128[k] = f_3 * pc_z[k] * sfg_90[k];
    }
}

static auto
compute_prim_sfh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sdh0,
                                                          const size_t sdg, const size_t sdh1,
                                                          const size_t sff0, const size_t sff1,
                                                          const size_t sfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh0_63 = buffer.data(sdh0 + 63);
    const auto *sdh0_66 = buffer.data(sdh0 + 66);
    const auto *sdh0_69 = buffer.data(sdh0 + 69);
    const auto *sdh0_78 = buffer.data(sdh0 + 78);
    const auto *sdh0_80 = buffer.data(sdh0 + 80);
    const auto *sdh0_81 = buffer.data(sdh0 + 81);
    const auto *sdh0_105 = buffer.data(sdh0 + 105);
    const auto *sdh0_110 = buffer.data(sdh0 + 110);
    const auto *sdh0_114 = buffer.data(sdh0 + 114);
    const auto *sdh0_120 = buffer.data(sdh0 + 120);
    const auto *sdh0_122 = buffer.data(sdh0 + 122);
    const auto *sdh0_123 = buffer.data(sdh0 + 123);
    const auto *sdh0_125 = buffer.data(sdh0 + 125);

    const auto *sdg_45 = buffer.data(sdg + 45);
    const auto *sdg_47 = buffer.data(sdg + 47);
    const auto *sdg_48 = buffer.data(sdg + 48);
    const auto *sdg_50 = buffer.data(sdg + 50);
    const auto *sdg_55 = buffer.data(sdg + 55);
    const auto *sdg_56 = buffer.data(sdg + 56);
    const auto *sdg_57 = buffer.data(sdg + 57);
    const auto *sdg_58 = buffer.data(sdg + 58);
    const auto *sdg_59 = buffer.data(sdg + 59);
    const auto *sdg_60 = buffer.data(sdg + 60);
    const auto *sdg_62 = buffer.data(sdg + 62);
    const auto *sdg_63 = buffer.data(sdg + 63);
    const auto *sdg_65 = buffer.data(sdg + 65);
    const auto *sdg_70 = buffer.data(sdg + 70);
    const auto *sdg_74 = buffer.data(sdg + 74);
    const auto *sdg_75 = buffer.data(sdg + 75);
    const auto *sdg_77 = buffer.data(sdg + 77);
    const auto *sdg_78 = buffer.data(sdg + 78);
    const auto *sdg_80 = buffer.data(sdg + 80);
    const auto *sdg_85 = buffer.data(sdg + 85);
    const auto *sdg_87 = buffer.data(sdg + 87);
    const auto *sdg_88 = buffer.data(sdg + 88);
    const auto *sdg_89 = buffer.data(sdg + 89);

    const auto *sdh1_63 = buffer.data(sdh1 + 63);
    const auto *sdh1_66 = buffer.data(sdh1 + 66);
    const auto *sdh1_69 = buffer.data(sdh1 + 69);
    const auto *sdh1_78 = buffer.data(sdh1 + 78);
    const auto *sdh1_80 = buffer.data(sdh1 + 80);
    const auto *sdh1_81 = buffer.data(sdh1 + 81);
    const auto *sdh1_105 = buffer.data(sdh1 + 105);
    const auto *sdh1_110 = buffer.data(sdh1 + 110);
    const auto *sdh1_114 = buffer.data(sdh1 + 114);
    const auto *sdh1_120 = buffer.data(sdh1 + 120);
    const auto *sdh1_122 = buffer.data(sdh1 + 122);
    const auto *sdh1_123 = buffer.data(sdh1 + 123);
    const auto *sdh1_125 = buffer.data(sdh1 + 125);

    const auto *sff0_63 = buffer.data(sff0 + 63);
    const auto *sff0_65 = buffer.data(sff0 + 65);
    const auto *sff0_66 = buffer.data(sff0 + 66);
    const auto *sff0_68 = buffer.data(sff0 + 68);
    const auto *sff0_69 = buffer.data(sff0 + 69);
    const auto *sff0_75 = buffer.data(sff0 + 75);
    const auto *sff0_79 = buffer.data(sff0 + 79);
    const auto *sff0_83 = buffer.data(sff0 + 83);
    const auto *sff0_86 = buffer.data(sff0 + 86);
    const auto *sff0_90 = buffer.data(sff0 + 90);
    const auto *sff0_93 = buffer.data(sff0 + 93);
    const auto *sff0_95 = buffer.data(sff0 + 95);
    const auto *sff0_96 = buffer.data(sff0 + 96);
    const auto *sff0_98 = buffer.data(sff0 + 98);
    const auto *sff0_99 = buffer.data(sff0 + 99);

    const auto *sff1_63 = buffer.data(sff1 + 63);
    const auto *sff1_65 = buffer.data(sff1 + 65);
    const auto *sff1_66 = buffer.data(sff1 + 66);
    const auto *sff1_68 = buffer.data(sff1 + 68);
    const auto *sff1_69 = buffer.data(sff1 + 69);
    const auto *sff1_75 = buffer.data(sff1 + 75);
    const auto *sff1_79 = buffer.data(sff1 + 79);
    const auto *sff1_83 = buffer.data(sff1 + 83);
    const auto *sff1_86 = buffer.data(sff1 + 86);
    const auto *sff1_90 = buffer.data(sff1 + 90);
    const auto *sff1_93 = buffer.data(sff1 + 93);
    const auto *sff1_95 = buffer.data(sff1 + 95);
    const auto *sff1_96 = buffer.data(sff1 + 96);
    const auto *sff1_98 = buffer.data(sff1 + 98);
    const auto *sff1_99 = buffer.data(sff1 + 99);

    const auto *sfg_92 = buffer.data(sfg + 92);
    const auto *sfg_93 = buffer.data(sfg + 93);
    const auto *sfg_95 = buffer.data(sfg + 95);
    const auto *sfg_96 = buffer.data(sfg + 96);
    const auto *sfg_99 = buffer.data(sfg + 99);
    const auto *sfg_100 = buffer.data(sfg + 100);
    const auto *sfg_101 = buffer.data(sfg + 101);
    const auto *sfg_102 = buffer.data(sfg + 102);
    const auto *sfg_103 = buffer.data(sfg + 103);
    const auto *sfg_104 = buffer.data(sfg + 104);
    const auto *sfg_105 = buffer.data(sfg + 105);
    const auto *sfg_107 = buffer.data(sfg + 107);
    const auto *sfg_108 = buffer.data(sfg + 108);
    const auto *sfg_110 = buffer.data(sfg + 110);
    const auto *sfg_114 = buffer.data(sfg + 114);
    const auto *sfg_115 = buffer.data(sfg + 115);
    const auto *sfg_116 = buffer.data(sfg + 116);
    const auto *sfg_117 = buffer.data(sfg + 117);
    const auto *sfg_118 = buffer.data(sfg + 118);
    const auto *sfg_119 = buffer.data(sfg + 119);
    const auto *sfg_120 = buffer.data(sfg + 120);
    const auto *sfg_122 = buffer.data(sfg + 122);
    const auto *sfg_123 = buffer.data(sfg + 123);
    const auto *sfg_125 = buffer.data(sfg + 125);
    const auto *sfg_126 = buffer.data(sfg + 126);
    const auto *sfg_130 = buffer.data(sfg + 130);
    const auto *sfg_131 = buffer.data(sfg + 131);
    const auto *sfg_132 = buffer.data(sfg + 132);
    const auto *sfg_133 = buffer.data(sfg + 133);
    const auto *sfg_134 = buffer.data(sfg + 134);
    const auto *sfg_135 = buffer.data(sfg + 135);
    const auto *sfg_137 = buffer.data(sfg + 137);
    const auto *sfg_138 = buffer.data(sfg + 138);
    const auto *sfg_140 = buffer.data(sfg + 140);
    const auto *sfg_141 = buffer.data(sfg + 141);
    const auto *sfg_144 = buffer.data(sfg + 144);
    const auto *sfg_145 = buffer.data(sfg + 145);
    const auto *sfg_146 = buffer.data(sfg + 146);
    const auto *sfg_147 = buffer.data(sfg + 147);
    const auto *sfg_148 = buffer.data(sfg + 148);
    const auto *sfg_149 = buffer.data(sfg + 149);

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, sdg_47, sff0_63, sff0_65, sff1_63, \
                         sff1_65, sfg_92, sfg_93, sfg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_4 * sff0_63[k]
                   - f_5 * sff1_63[k]
                   + f_3 * pc_x[k] * sfg_93[k];

        t_130[k] = f_0 * sdg_47[k]
                   + f_3 * pc_y[k] * sfg_92[k];

        t_131[k] = f_4 * sff0_65[k]
                   - f_5 * sff1_65[k]
                   + f_3 * pc_x[k] * sfg_95[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_y, pc_z, sdg_50, sff0_66, \
                         sff0_69, sff1_66, sff1_69, sfg_93, sfg_95, sfg_96, \
                         sfg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_6 * sff0_66[k]
                   - f_7 * sff1_66[k]
                   + f_3 * pc_x[k] * sfg_96[k];

        t_133[k] = f_3 * pc_z[k] * sfg_93[k];

        t_134[k] = f_0 * sdg_50[k]
                   + f_3 * pc_y[k] * sfg_95[k];

        t_135[k] = f_6 * sff0_69[k]
                   - f_7 * sff1_69[k]
                   + f_3 * pc_x[k] * sfg_99[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, t_141, pc_x, pc_y, sdg_55, \
                         sff0_66, sff1_66, sfg_100, sfg_101, sfg_102, sfg_103, \
                         sfg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_3 * pc_x[k] * sfg_100[k];

        t_137[k] = f_3 * pc_x[k] * sfg_101[k];

        t_138[k] = f_3 * pc_x[k] * sfg_102[k];

        t_139[k] = f_3 * pc_x[k] * sfg_103[k];

        t_140[k] = f_3 * pc_x[k] * sfg_104[k];

        t_141[k] = f_0 * sdg_55[k]
                   + f_1 * sff0_66[k]
                   - f_2 * sff1_66[k]
                   + f_3 * pc_y[k] * sfg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pc_y, pc_z, sdg_57, sdg_58, sff0_68, sff0_69, \
                         sff1_68, sff1_69, sfg_100, sfg_102, sfg_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * sfg_100[k];

        t_143[k] = f_0 * sdg_57[k]
                   + f_4 * sff0_68[k]
                   - f_5 * sff1_68[k]
                   + f_3 * pc_y[k] * sfg_102[k];

        t_144[k] = f_0 * sdg_58[k]
                   + f_6 * sff0_69[k]
                   - f_7 * sff1_69[k]
                   + f_3 * pc_y[k] * sfg_103[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pb_z, pc_y, pc_z, sdh0_63, sdg_59, \
                         sdg_60, sdh1_63, sff0_69, sff1_69, sfg_104, \
                         sfg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_0 * sdg_59[k]
                   + f_3 * pc_y[k] * sfg_104[k];

        t_146[k] = f_1 * sff0_69[k]
                   - f_2 * sff1_69[k]
                   + f_3 * pc_z[k] * sfg_104[k];

        t_147[k] = pb_z[k] * sdh0_63[k]
                   - f_8 * pc_z[k] * sdh1_63[k];

        t_148[k] = f_10 * sdg_60[k]
                   + f_3 * pc_y[k] * sfg_105[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_z, pc_y, pc_z, sdh0_66, sdg_45, sdg_62, \
                         sdh1_66, sfg_105, sfg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_9 * sdg_45[k]
                   + f_3 * pc_z[k] * sfg_105[k];

        t_150[k] = pb_z[k] * sdh0_66[k]
                   - f_8 * pc_z[k] * sdh1_66[k];

        t_151[k] = f_10 * sdg_62[k]
                   + f_3 * pc_y[k] * sfg_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_z, pc_x, pc_y, pc_z, sdh0_69, sdg_48, \
                         sdg_65, sdh1_69, sff0_75, sff1_75, sfg_108, \
                         sfg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * sff0_75[k]
                   - f_5 * sff1_75[k]
                   + f_3 * pc_x[k] * sfg_110[k];

        t_153[k] = pb_z[k] * sdh0_69[k]
                   - f_8 * pc_z[k] * sdh1_69[k];

        t_154[k] = f_9 * sdg_48[k]
                   + f_3 * pc_z[k] * sfg_108[k];

        t_155[k] = f_10 * sdg_65[k]
                   + f_3 * pc_y[k] * sfg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, pc_x, sff0_79, sff1_79, \
                         sfg_114, sfg_115, sfg_116, sfg_117, sfg_118, \
                         sfg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_6 * sff0_79[k]
                   - f_7 * sff1_79[k]
                   + f_3 * pc_x[k] * sfg_114[k];

        t_157[k] = f_3 * pc_x[k] * sfg_115[k];

        t_158[k] = f_3 * pc_x[k] * sfg_116[k];

        t_159[k] = f_3 * pc_x[k] * sfg_117[k];

        t_160[k] = f_3 * pc_x[k] * sfg_118[k];

        t_161[k] = f_3 * pc_x[k] * sfg_119[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_z, sdh0_78, sdh0_80, sdh0_81, \
                         sdg_55, sdg_56, sdg_57, sdh1_78, sdh1_80, sdh1_81, \
                         sfg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pb_z[k] * sdh0_78[k]
                   - f_8 * pc_z[k] * sdh1_78[k];

        t_163[k] = f_9 * sdg_55[k]
                   + f_3 * pc_z[k] * sfg_115[k];

        t_164[k] = pb_z[k] * sdh0_80[k]
                   + f_10 * sdg_56[k]
                   - f_8 * pc_z[k] * sdh1_80[k];

        t_165[k] = pb_z[k] * sdh0_81[k]
                   + f_0 * sdg_57[k]
                   - f_8 * pc_z[k] * sdh1_81[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, sdh0_105, sdg_59, \
                         sdg_74, sdg_75, sdh1_105, sff0_79, sff1_79, sfg_119, \
                         sfg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * sdg_74[k]
                   + f_3 * pc_y[k] * sfg_119[k];

        t_167[k] = f_9 * sdg_59[k]
                   + f_1 * sff0_79[k]
                   - f_2 * sff1_79[k]
                   + f_3 * pc_z[k] * sfg_119[k];

        t_168[k] = pb_y[k] * sdh0_105[k]
                   - f_8 * pc_y[k] * sdh1_105[k];

        t_169[k] = f_9 * sdg_75[k]
                   + f_3 * pc_y[k] * sfg_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pc_x, pc_y, pc_z, sdg_60, sdg_77, sff0_83, \
                         sff1_83, sfg_120, sfg_122, sfg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * sdg_60[k]
                   + f_3 * pc_z[k] * sfg_120[k];

        t_171[k] = f_4 * sff0_83[k]
                   - f_5 * sff1_83[k]
                   + f_3 * pc_x[k] * sfg_123[k];

        t_172[k] = f_9 * sdg_77[k]
                   + f_3 * pc_y[k] * sfg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pb_y, pc_x, pc_y, pc_z, sdh0_110, sdg_63, \
                         sdh1_110, sff0_86, sff1_86, sfg_123, sfg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pb_y[k] * sdh0_110[k]
                   - f_8 * pc_y[k] * sdh1_110[k];

        t_174[k] = f_6 * sff0_86[k]
                   - f_7 * sff1_86[k]
                   + f_3 * pc_x[k] * sfg_126[k];

        t_175[k] = f_10 * sdg_63[k]
                   + f_3 * pc_z[k] * sfg_123[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, pb_y, pc_x, pc_y, sdh0_114, \
                         sdg_80, sdh1_114, sfg_125, sfg_130, sfg_131, \
                         sfg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * sdg_80[k]
                   + f_3 * pc_y[k] * sfg_125[k];

        t_177[k] = pb_y[k] * sdh0_114[k]
                   - f_8 * pc_y[k] * sdh1_114[k];

        t_178[k] = f_3 * pc_x[k] * sfg_130[k];

        t_179[k] = f_3 * pc_x[k] * sfg_131[k];

        t_180[k] = f_3 * pc_x[k] * sfg_132[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_y, pc_x, pc_y, pc_z, sdh0_120, sdg_70, \
                         sdg_85, sdh1_120, sfg_130, sfg_133, sfg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_3 * pc_x[k] * sfg_133[k];

        t_182[k] = f_3 * pc_x[k] * sfg_134[k];

        t_183[k] = pb_y[k] * sdh0_120[k]
                   + f_11 * sdg_85[k]
                   - f_8 * pc_y[k] * sdh1_120[k];

        t_184[k] = f_10 * sdg_70[k]
                   + f_3 * pc_z[k] * sfg_130[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pc_y, sdh0_122, sdh0_123, sdh0_125, \
                         sdg_87, sdg_88, sdg_89, sdh1_122, sdh1_123, sdh1_125, \
                         sfg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * sdh0_122[k]
                   + f_0 * sdg_87[k]
                   - f_8 * pc_y[k] * sdh1_122[k];

        t_186[k] = pb_y[k] * sdh0_123[k]
                   + f_10 * sdg_88[k]
                   - f_8 * pc_y[k] * sdh1_123[k];

        t_187[k] = f_9 * sdg_89[k]
                   + f_3 * pc_y[k] * sfg_134[k];

        t_188[k] = pb_y[k] * sdh0_125[k]
                   - f_8 * pc_y[k] * sdh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, sdg_75, sff0_90, \
                         sff0_93, sff1_90, sff1_93, sfg_135, sfg_137, \
                         sfg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_1 * sff0_90[k]
                   - f_2 * sff1_90[k]
                   + f_3 * pc_x[k] * sfg_135[k];

        t_190[k] = f_3 * pc_y[k] * sfg_135[k];

        t_191[k] = f_0 * sdg_75[k]
                   + f_3 * pc_z[k] * sfg_135[k];

        t_192[k] = f_4 * sff0_93[k]
                   - f_5 * sff1_93[k]
                   + f_3 * pc_x[k] * sfg_138[k];

        t_193[k] = f_3 * pc_y[k] * sfg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, sdg_78, sff0_95, \
                         sff0_96, sff1_95, sff1_96, sfg_138, sfg_140, \
                         sfg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_4 * sff0_95[k]
                   - f_5 * sff1_95[k]
                   + f_3 * pc_x[k] * sfg_140[k];

        t_195[k] = f_6 * sff0_96[k]
                   - f_7 * sff1_96[k]
                   + f_3 * pc_x[k] * sfg_141[k];

        t_196[k] = f_0 * sdg_78[k]
                   + f_3 * pc_z[k] * sfg_138[k];

        t_197[k] = f_3 * pc_y[k] * sfg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, pc_x, sff0_99, sff1_99, \
                         sfg_144, sfg_145, sfg_146, sfg_147, sfg_148, \
                         sfg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_6 * sff0_99[k]
                   - f_7 * sff1_99[k]
                   + f_3 * pc_x[k] * sfg_144[k];

        t_199[k] = f_3 * pc_x[k] * sfg_145[k];

        t_200[k] = f_3 * pc_x[k] * sfg_146[k];

        t_201[k] = f_3 * pc_x[k] * sfg_147[k];

        t_202[k] = f_3 * pc_x[k] * sfg_148[k];

        t_203[k] = f_3 * pc_x[k] * sfg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, sdg_85, sff0_96, sff0_98, \
                         sff0_99, sff1_96, sff1_98, sff1_99, sfg_145, sfg_147, \
                         sfg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * sff0_96[k]
                   - f_2 * sff1_96[k]
                   + f_3 * pc_y[k] * sfg_145[k];

        t_205[k] = f_0 * sdg_85[k]
                   + f_3 * pc_z[k] * sfg_145[k];

        t_206[k] = f_4 * sff0_98[k]
                   - f_5 * sff1_98[k]
                   + f_3 * pc_y[k] * sfg_147[k];

        t_207[k] = f_6 * sff0_99[k]
                   - f_7 * sff1_99[k]
                   + f_3 * pc_y[k] * sfg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, pc_y, pc_z, sdg_89, sff0_99, sff1_99, \
                         sfg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * sfg_149[k];

        t_209[k] = f_0 * sdg_89[k]
                   + f_1 * sff0_99[k]
                   - f_2 * sff1_99[k]
                   + f_3 * pc_z[k] * sfg_149[k];
    }
}

auto
compute_prim_sfh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdh0, const size_t sdg,
                                                   const size_t sdh1, const size_t sff0,
                                                   const size_t sff1, const size_t sfg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sfh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sdh0, sdg,
                                                              sdh1, sff0, sff1, sfg, ncols,
                                                              gamma, p, q);

    compute_prim_sfh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sdh0, sdg,
                                                              sdh1, sff0, sff1, sfg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
