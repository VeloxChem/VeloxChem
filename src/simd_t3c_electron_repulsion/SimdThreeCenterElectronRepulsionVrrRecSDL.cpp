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


#include "SimdThreeCenterElectronRepulsionVrrRecSDL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sdl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t spl0,
                                                          const size_t spk, const size_t spl1,
                                                          const size_t sdi0, const size_t sdi1,
                                                          const size_t sdk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spl0_0 = buffer.data(spl0 + 0);
    const auto *spl0_3 = buffer.data(spl0 + 3);
    const auto *spl0_5 = buffer.data(spl0 + 5);
    const auto *spl0_6 = buffer.data(spl0 + 6);
    const auto *spl0_9 = buffer.data(spl0 + 9);
    const auto *spl0_10 = buffer.data(spl0 + 10);
    const auto *spl0_14 = buffer.data(spl0 + 14);
    const auto *spl0_15 = buffer.data(spl0 + 15);
    const auto *spl0_20 = buffer.data(spl0 + 20);
    const auto *spl0_21 = buffer.data(spl0 + 21);
    const auto *spl0_27 = buffer.data(spl0 + 27);
    const auto *spl0_48 = buffer.data(spl0 + 48);
    const auto *spl0_51 = buffer.data(spl0 + 51);
    const auto *spl0_55 = buffer.data(spl0 + 55);
    const auto *spl0_57 = buffer.data(spl0 + 57);
    const auto *spl0_60 = buffer.data(spl0 + 60);
    const auto *spl0_62 = buffer.data(spl0 + 62);
    const auto *spl0_63 = buffer.data(spl0 + 63);
    const auto *spl0_66 = buffer.data(spl0 + 66);
    const auto *spl0_68 = buffer.data(spl0 + 68);
    const auto *spl0_69 = buffer.data(spl0 + 69);
    const auto *spl0_70 = buffer.data(spl0 + 70);
    const auto *spl0_81 = buffer.data(spl0 + 81);
    const auto *spl0_83 = buffer.data(spl0 + 83);
    const auto *spl0_84 = buffer.data(spl0 + 84);
    const auto *spl0_85 = buffer.data(spl0 + 85);
    const auto *spl0_86 = buffer.data(spl0 + 86);
    const auto *spl0_87 = buffer.data(spl0 + 87);
    const auto *spl0_89 = buffer.data(spl0 + 89);
    const auto *spl0_95 = buffer.data(spl0 + 95);
    const auto *spl0_99 = buffer.data(spl0 + 99);
    const auto *spl0_102 = buffer.data(spl0 + 102);
    const auto *spl0_104 = buffer.data(spl0 + 104);
    const auto *spl0_107 = buffer.data(spl0 + 107);
    const auto *spl0_108 = buffer.data(spl0 + 108);
    const auto *spl0_110 = buffer.data(spl0 + 110);
    const auto *spl0_113 = buffer.data(spl0 + 113);
    const auto *spl0_114 = buffer.data(spl0 + 114);
    const auto *spl0_115 = buffer.data(spl0 + 115);
    const auto *spl0_117 = buffer.data(spl0 + 117);

    const auto *spk_0 = buffer.data(spk + 0);
    const auto *spk_2 = buffer.data(spk + 2);
    const auto *spk_3 = buffer.data(spk + 3);
    const auto *spk_5 = buffer.data(spk + 5);
    const auto *spk_6 = buffer.data(spk + 6);
    const auto *spk_9 = buffer.data(spk + 9);
    const auto *spk_10 = buffer.data(spk + 10);
    const auto *spk_12 = buffer.data(spk + 12);
    const auto *spk_14 = buffer.data(spk + 14);
    const auto *spk_15 = buffer.data(spk + 15);
    const auto *spk_17 = buffer.data(spk + 17);
    const auto *spk_18 = buffer.data(spk + 18);
    const auto *spk_20 = buffer.data(spk + 20);
    const auto *spk_21 = buffer.data(spk + 21);
    const auto *spk_23 = buffer.data(spk + 23);
    const auto *spk_24 = buffer.data(spk + 24);
    const auto *spk_25 = buffer.data(spk + 25);
    const auto *spk_27 = buffer.data(spk + 27);
    const auto *spk_28 = buffer.data(spk + 28);
    const auto *spk_29 = buffer.data(spk + 29);
    const auto *spk_30 = buffer.data(spk + 30);
    const auto *spk_31 = buffer.data(spk + 31);
    const auto *spk_32 = buffer.data(spk + 32);
    const auto *spk_33 = buffer.data(spk + 33);
    const auto *spk_34 = buffer.data(spk + 34);
    const auto *spk_35 = buffer.data(spk + 35);
    const auto *spk_39 = buffer.data(spk + 39);
    const auto *spk_42 = buffer.data(spk + 42);
    const auto *spk_46 = buffer.data(spk + 46);
    const auto *spk_48 = buffer.data(spk + 48);
    const auto *spk_51 = buffer.data(spk + 51);
    const auto *spk_53 = buffer.data(spk + 53);
    const auto *spk_54 = buffer.data(spk + 54);
    const auto *spk_57 = buffer.data(spk + 57);
    const auto *spk_59 = buffer.data(spk + 59);
    const auto *spk_60 = buffer.data(spk + 60);
    const auto *spk_61 = buffer.data(spk + 61);
    const auto *spk_64 = buffer.data(spk + 64);
    const auto *spk_65 = buffer.data(spk + 65);
    const auto *spk_66 = buffer.data(spk + 66);
    const auto *spk_67 = buffer.data(spk + 67);
    const auto *spk_68 = buffer.data(spk + 68);
    const auto *spk_69 = buffer.data(spk + 69);
    const auto *spk_70 = buffer.data(spk + 70);
    const auto *spk_71 = buffer.data(spk + 71);
    const auto *spk_77 = buffer.data(spk + 77);
    const auto *spk_81 = buffer.data(spk + 81);
    const auto *spk_84 = buffer.data(spk + 84);
    const auto *spk_86 = buffer.data(spk + 86);
    const auto *spk_89 = buffer.data(spk + 89);
    const auto *spk_90 = buffer.data(spk + 90);
    const auto *spk_92 = buffer.data(spk + 92);
    const auto *spk_95 = buffer.data(spk + 95);
    const auto *spk_96 = buffer.data(spk + 96);
    const auto *spk_97 = buffer.data(spk + 97);
    const auto *spk_99 = buffer.data(spk + 99);
    const auto *spk_100 = buffer.data(spk + 100);
    const auto *spk_101 = buffer.data(spk + 101);

    const auto *spl1_0 = buffer.data(spl1 + 0);
    const auto *spl1_3 = buffer.data(spl1 + 3);
    const auto *spl1_5 = buffer.data(spl1 + 5);
    const auto *spl1_6 = buffer.data(spl1 + 6);
    const auto *spl1_9 = buffer.data(spl1 + 9);
    const auto *spl1_10 = buffer.data(spl1 + 10);
    const auto *spl1_14 = buffer.data(spl1 + 14);
    const auto *spl1_15 = buffer.data(spl1 + 15);
    const auto *spl1_20 = buffer.data(spl1 + 20);
    const auto *spl1_21 = buffer.data(spl1 + 21);
    const auto *spl1_27 = buffer.data(spl1 + 27);
    const auto *spl1_48 = buffer.data(spl1 + 48);
    const auto *spl1_51 = buffer.data(spl1 + 51);
    const auto *spl1_55 = buffer.data(spl1 + 55);
    const auto *spl1_57 = buffer.data(spl1 + 57);
    const auto *spl1_60 = buffer.data(spl1 + 60);
    const auto *spl1_62 = buffer.data(spl1 + 62);
    const auto *spl1_63 = buffer.data(spl1 + 63);
    const auto *spl1_66 = buffer.data(spl1 + 66);
    const auto *spl1_68 = buffer.data(spl1 + 68);
    const auto *spl1_69 = buffer.data(spl1 + 69);
    const auto *spl1_70 = buffer.data(spl1 + 70);
    const auto *spl1_81 = buffer.data(spl1 + 81);
    const auto *spl1_83 = buffer.data(spl1 + 83);
    const auto *spl1_84 = buffer.data(spl1 + 84);
    const auto *spl1_85 = buffer.data(spl1 + 85);
    const auto *spl1_86 = buffer.data(spl1 + 86);
    const auto *spl1_87 = buffer.data(spl1 + 87);
    const auto *spl1_89 = buffer.data(spl1 + 89);
    const auto *spl1_95 = buffer.data(spl1 + 95);
    const auto *spl1_99 = buffer.data(spl1 + 99);
    const auto *spl1_102 = buffer.data(spl1 + 102);
    const auto *spl1_104 = buffer.data(spl1 + 104);
    const auto *spl1_107 = buffer.data(spl1 + 107);
    const auto *spl1_108 = buffer.data(spl1 + 108);
    const auto *spl1_110 = buffer.data(spl1 + 110);
    const auto *spl1_113 = buffer.data(spl1 + 113);
    const auto *spl1_114 = buffer.data(spl1 + 114);
    const auto *spl1_115 = buffer.data(spl1 + 115);
    const auto *spl1_117 = buffer.data(spl1 + 117);

    const auto *sdi0_0 = buffer.data(sdi0 + 0);
    const auto *sdi0_3 = buffer.data(sdi0 + 3);
    const auto *sdi0_5 = buffer.data(sdi0 + 5);
    const auto *sdi0_6 = buffer.data(sdi0 + 6);
    const auto *sdi0_9 = buffer.data(sdi0 + 9);
    const auto *sdi0_10 = buffer.data(sdi0 + 10);
    const auto *sdi0_12 = buffer.data(sdi0 + 12);
    const auto *sdi0_14 = buffer.data(sdi0 + 14);
    const auto *sdi0_15 = buffer.data(sdi0 + 15);
    const auto *sdi0_17 = buffer.data(sdi0 + 17);
    const auto *sdi0_18 = buffer.data(sdi0 + 18);
    const auto *sdi0_20 = buffer.data(sdi0 + 20);
    const auto *sdi0_21 = buffer.data(sdi0 + 21);
    const auto *sdi0_23 = buffer.data(sdi0 + 23);
    const auto *sdi0_24 = buffer.data(sdi0 + 24);
    const auto *sdi0_25 = buffer.data(sdi0 + 25);
    const auto *sdi0_26 = buffer.data(sdi0 + 26);
    const auto *sdi0_27 = buffer.data(sdi0 + 27);

    const auto *sdi1_0 = buffer.data(sdi1 + 0);
    const auto *sdi1_3 = buffer.data(sdi1 + 3);
    const auto *sdi1_5 = buffer.data(sdi1 + 5);
    const auto *sdi1_6 = buffer.data(sdi1 + 6);
    const auto *sdi1_9 = buffer.data(sdi1 + 9);
    const auto *sdi1_10 = buffer.data(sdi1 + 10);
    const auto *sdi1_12 = buffer.data(sdi1 + 12);
    const auto *sdi1_14 = buffer.data(sdi1 + 14);
    const auto *sdi1_15 = buffer.data(sdi1 + 15);
    const auto *sdi1_17 = buffer.data(sdi1 + 17);
    const auto *sdi1_18 = buffer.data(sdi1 + 18);
    const auto *sdi1_20 = buffer.data(sdi1 + 20);
    const auto *sdi1_21 = buffer.data(sdi1 + 21);
    const auto *sdi1_23 = buffer.data(sdi1 + 23);
    const auto *sdi1_24 = buffer.data(sdi1 + 24);
    const auto *sdi1_25 = buffer.data(sdi1 + 25);
    const auto *sdi1_26 = buffer.data(sdi1 + 26);
    const auto *sdi1_27 = buffer.data(sdi1 + 27);

    const auto *sdk_0 = buffer.data(sdk + 0);
    const auto *sdk_2 = buffer.data(sdk + 2);
    const auto *sdk_3 = buffer.data(sdk + 3);
    const auto *sdk_5 = buffer.data(sdk + 5);
    const auto *sdk_6 = buffer.data(sdk + 6);
    const auto *sdk_9 = buffer.data(sdk + 9);
    const auto *sdk_10 = buffer.data(sdk + 10);
    const auto *sdk_12 = buffer.data(sdk + 12);
    const auto *sdk_14 = buffer.data(sdk + 14);
    const auto *sdk_15 = buffer.data(sdk + 15);
    const auto *sdk_17 = buffer.data(sdk + 17);
    const auto *sdk_18 = buffer.data(sdk + 18);
    const auto *sdk_20 = buffer.data(sdk + 20);
    const auto *sdk_21 = buffer.data(sdk + 21);
    const auto *sdk_23 = buffer.data(sdk + 23);
    const auto *sdk_24 = buffer.data(sdk + 24);
    const auto *sdk_25 = buffer.data(sdk + 25);
    const auto *sdk_27 = buffer.data(sdk + 27);
    const auto *sdk_28 = buffer.data(sdk + 28);
    const auto *sdk_29 = buffer.data(sdk + 29);
    const auto *sdk_30 = buffer.data(sdk + 30);
    const auto *sdk_31 = buffer.data(sdk + 31);
    const auto *sdk_32 = buffer.data(sdk + 32);
    const auto *sdk_33 = buffer.data(sdk + 33);
    const auto *sdk_34 = buffer.data(sdk + 34);
    const auto *sdk_35 = buffer.data(sdk + 35);
    const auto *sdk_36 = buffer.data(sdk + 36);
    const auto *sdk_38 = buffer.data(sdk + 38);
    const auto *sdk_39 = buffer.data(sdk + 39);
    const auto *sdk_41 = buffer.data(sdk + 41);
    const auto *sdk_42 = buffer.data(sdk + 42);
    const auto *sdk_45 = buffer.data(sdk + 45);
    const auto *sdk_46 = buffer.data(sdk + 46);
    const auto *sdk_50 = buffer.data(sdk + 50);
    const auto *sdk_51 = buffer.data(sdk + 51);
    const auto *sdk_56 = buffer.data(sdk + 56);
    const auto *sdk_64 = buffer.data(sdk + 64);
    const auto *sdk_65 = buffer.data(sdk + 65);
    const auto *sdk_66 = buffer.data(sdk + 66);
    const auto *sdk_67 = buffer.data(sdk + 67);
    const auto *sdk_68 = buffer.data(sdk + 68);
    const auto *sdk_69 = buffer.data(sdk + 69);
    const auto *sdk_70 = buffer.data(sdk + 70);
    const auto *sdk_71 = buffer.data(sdk + 71);
    const auto *sdk_72 = buffer.data(sdk + 72);
    const auto *sdk_74 = buffer.data(sdk + 74);
    const auto *sdk_75 = buffer.data(sdk + 75);
    const auto *sdk_77 = buffer.data(sdk + 77);
    const auto *sdk_78 = buffer.data(sdk + 78);
    const auto *sdk_81 = buffer.data(sdk + 81);
    const auto *sdk_82 = buffer.data(sdk + 82);
    const auto *sdk_86 = buffer.data(sdk + 86);
    const auto *sdk_87 = buffer.data(sdk + 87);
    const auto *sdk_92 = buffer.data(sdk + 92);
    const auto *sdk_100 = buffer.data(sdk + 100);
    const auto *sdk_101 = buffer.data(sdk + 101);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, spk_0, spk_3, sdi0_0, sdi0_3, \
                         sdi1_0, sdi1_3, sdk_0, sdk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spk_0[k]
                 + f_1 * sdi0_0[k]
                 - f_2 * sdi1_0[k]
                 + f_3 * pc_x[k] * sdk_0[k];

        t_1[k] = f_3 * pc_y[k] * sdk_0[k];

        t_2[k] = f_3 * pc_z[k] * sdk_0[k];

        t_3[k] = f_0 * spk_3[k]
                 + f_4 * sdi0_3[k]
                 - f_5 * sdi1_3[k]
                 + f_3 * pc_x[k] * sdk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, spk_5, spk_6, sdi0_5, sdi0_6, sdi1_5, \
                         sdi1_6, sdk_2, sdk_5, sdk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sdk_2[k];

        t_5[k] = f_0 * spk_5[k]
                 + f_4 * sdi0_5[k]
                 - f_5 * sdi1_5[k]
                 + f_3 * pc_x[k] * sdk_5[k];

        t_6[k] = f_0 * spk_6[k]
                 + f_6 * sdi0_6[k]
                 - f_7 * sdi1_6[k]
                 + f_3 * pc_x[k] * sdk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, spk_9, sdi0_9, sdi1_9, sdk_3, sdk_5, \
                         sdk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sdk_3[k];

        t_8[k] = f_3 * pc_y[k] * sdk_5[k];

        t_9[k] = f_0 * spk_9[k]
                 + f_6 * sdi0_9[k]
                 - f_7 * sdi1_9[k]
                 + f_3 * pc_x[k] * sdk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, spk_10, spk_12, sdi0_10, sdi0_12, \
                         sdi1_10, sdi1_12, sdk_6, sdk_10, sdk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * spk_10[k]
                  + f_8 * sdi0_10[k]
                  - f_9 * sdi1_10[k]
                  + f_3 * pc_x[k] * sdk_10[k];

        t_11[k] = f_3 * pc_z[k] * sdk_6[k];

        t_12[k] = f_0 * spk_12[k]
                  + f_8 * sdi0_12[k]
                  - f_9 * sdi1_12[k]
                  + f_3 * pc_x[k] * sdk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, spk_14, spk_15, sdi0_14, sdi0_15, \
                         sdi1_14, sdi1_15, sdk_9, sdk_14, sdk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sdk_9[k];

        t_14[k] = f_0 * spk_14[k]
                  + f_8 * sdi0_14[k]
                  - f_9 * sdi1_14[k]
                  + f_3 * pc_x[k] * sdk_14[k];

        t_15[k] = f_0 * spk_15[k]
                  + f_10 * sdi0_15[k]
                  - f_11 * sdi1_15[k]
                  + f_3 * pc_x[k] * sdk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, spk_17, spk_18, sdi0_17, sdi0_18, \
                         sdi1_17, sdi1_18, sdk_10, sdk_17, sdk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sdk_10[k];

        t_17[k] = f_0 * spk_17[k]
                  + f_10 * sdi0_17[k]
                  - f_11 * sdi1_17[k]
                  + f_3 * pc_x[k] * sdk_17[k];

        t_18[k] = f_0 * spk_18[k]
                  + f_10 * sdi0_18[k]
                  - f_11 * sdi1_18[k]
                  + f_3 * pc_x[k] * sdk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, spk_20, spk_21, sdi0_20, sdi0_21, \
                         sdi1_20, sdi1_21, sdk_14, sdk_20, sdk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sdk_14[k];

        t_20[k] = f_0 * spk_20[k]
                  + f_10 * sdi0_20[k]
                  - f_11 * sdi1_20[k]
                  + f_3 * pc_x[k] * sdk_20[k];

        t_21[k] = f_0 * spk_21[k]
                  + f_12 * sdi0_21[k]
                  - f_13 * sdi1_21[k]
                  + f_3 * pc_x[k] * sdk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, spk_23, spk_24, sdi0_23, sdi0_24, \
                         sdi1_23, sdi1_24, sdk_15, sdk_23, sdk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * sdk_15[k];

        t_23[k] = f_0 * spk_23[k]
                  + f_12 * sdi0_23[k]
                  - f_13 * sdi1_23[k]
                  + f_3 * pc_x[k] * sdk_23[k];

        t_24[k] = f_0 * spk_24[k]
                  + f_12 * sdi0_24[k]
                  - f_13 * sdi1_24[k]
                  + f_3 * pc_x[k] * sdk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, spk_25, spk_27, sdi0_25, sdi0_27, \
                         sdi1_25, sdi1_27, sdk_20, sdk_25, sdk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * spk_25[k]
                  + f_12 * sdi0_25[k]
                  - f_13 * sdi1_25[k]
                  + f_3 * pc_x[k] * sdk_25[k];

        t_26[k] = f_3 * pc_y[k] * sdk_20[k];

        t_27[k] = f_0 * spk_27[k]
                  + f_12 * sdi0_27[k]
                  - f_13 * sdi1_27[k]
                  + f_3 * pc_x[k] * sdk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, spk_28, spk_29, spk_30, spk_31, \
                         spk_32, sdk_28, sdk_29, sdk_30, sdk_31, \
                         sdk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * spk_28[k]
                  + f_3 * pc_x[k] * sdk_28[k];

        t_29[k] = f_0 * spk_29[k]
                  + f_3 * pc_x[k] * sdk_29[k];

        t_30[k] = f_0 * spk_30[k]
                  + f_3 * pc_x[k] * sdk_30[k];

        t_31[k] = f_0 * spk_31[k]
                  + f_3 * pc_x[k] * sdk_31[k];

        t_32[k] = f_0 * spk_32[k]
                  + f_3 * pc_x[k] * sdk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, spk_33, spk_34, spk_35, sdi0_21, \
                         sdi1_21, sdk_28, sdk_33, sdk_34, sdk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * spk_33[k]
                  + f_3 * pc_x[k] * sdk_33[k];

        t_34[k] = f_0 * spk_34[k]
                  + f_3 * pc_x[k] * sdk_34[k];

        t_35[k] = f_0 * spk_35[k]
                  + f_3 * pc_x[k] * sdk_35[k];

        t_36[k] = f_1 * sdi0_21[k]
                  - f_2 * sdi1_21[k]
                  + f_3 * pc_y[k] * sdk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, sdi0_23, sdi0_24, sdi0_25, \
                         sdi1_23, sdi1_24, sdi1_25, sdk_28, sdk_30, sdk_31, \
                         sdk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * sdk_28[k];

        t_38[k] = f_4 * sdi0_23[k]
                  - f_5 * sdi1_23[k]
                  + f_3 * pc_y[k] * sdk_30[k];

        t_39[k] = f_6 * sdi0_24[k]
                  - f_7 * sdi1_24[k]
                  + f_3 * pc_y[k] * sdk_31[k];

        t_40[k] = f_8 * sdi0_25[k]
                  - f_9 * sdi1_25[k]
                  + f_3 * pc_y[k] * sdk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sdi0_26, sdi0_27, sdi1_26, \
                         sdi1_27, sdk_33, sdk_34, sdk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * sdi0_26[k]
                  - f_11 * sdi1_26[k]
                  + f_3 * pc_y[k] * sdk_33[k];

        t_42[k] = f_12 * sdi0_27[k]
                  - f_13 * sdi1_27[k]
                  + f_3 * pc_y[k] * sdk_34[k];

        t_43[k] = f_3 * pc_y[k] * sdk_35[k];

        t_44[k] = f_1 * sdi0_27[k]
                  - f_2 * sdi1_27[k]
                  + f_3 * pc_z[k] * sdk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_x, pb_y, pc_x, pc_y, pc_z, spl0_0, \
                         spl0_48, spk_0, spk_39, spl1_0, spl1_48, \
                         sdk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * spl0_0[k]
                  - f_14 * pc_y[k] * spl1_0[k];

        t_46[k] = f_15 * spk_0[k]
                  + f_3 * pc_y[k] * sdk_36[k];

        t_47[k] = f_3 * pc_z[k] * sdk_36[k];

        t_48[k] = pb_x[k] * spl0_48[k]
                  + f_16 * spk_39[k]
                  - f_14 * pc_x[k] * spl1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_y, pc_x, pc_y, spl0_5, spl0_51, spk_2, \
                         spk_42, spl1_5, spl1_51, sdk_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * spk_2[k]
                  + f_3 * pc_y[k] * sdk_38[k];

        t_50[k] = pb_y[k] * spl0_5[k]
                  - f_14 * pc_y[k] * spl1_5[k];

        t_51[k] = pb_x[k] * spl0_51[k]
                  + f_17 * spk_42[k]
                  - f_14 * pc_x[k] * spl1_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_y, pc_y, pc_z, spl0_9, spk_5, spl1_9, sdk_39, \
                         sdk_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sdk_39[k];

        t_53[k] = f_15 * spk_5[k]
                  + f_3 * pc_y[k] * sdk_41[k];

        t_54[k] = pb_y[k] * spl0_9[k]
                  - f_14 * pc_y[k] * spl1_9[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pc_x, pc_z, spl0_55, spl0_57, spk_46, spk_48, \
                         spl1_55, spl1_57, sdk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * spl0_55[k]
                  + f_18 * spk_46[k]
                  - f_14 * pc_x[k] * spl1_55[k];

        t_56[k] = f_3 * pc_z[k] * sdk_42[k];

        t_57[k] = pb_x[k] * spl0_57[k]
                  + f_18 * spk_48[k]
                  - f_14 * pc_x[k] * spl1_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, pc_x, pc_y, spl0_14, spl0_60, spk_9, \
                         spk_51, spl1_14, spl1_60, sdk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_15 * spk_9[k]
                  + f_3 * pc_y[k] * sdk_45[k];

        t_59[k] = pb_y[k] * spl0_14[k]
                  - f_14 * pc_y[k] * spl1_14[k];

        t_60[k] = pb_x[k] * spl0_60[k]
                  + f_19 * spk_51[k]
                  - f_14 * pc_x[k] * spl1_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, pc_x, pc_z, spl0_62, spl0_63, spk_53, spk_54, \
                         spl1_62, spl1_63, sdk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * sdk_46[k];

        t_62[k] = pb_x[k] * spl0_62[k]
                  + f_19 * spk_53[k]
                  - f_14 * pc_x[k] * spl1_62[k];

        t_63[k] = pb_x[k] * spl0_63[k]
                  + f_19 * spk_54[k]
                  - f_14 * pc_x[k] * spl1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, pc_x, pc_y, spl0_20, spl0_66, spk_14, \
                         spk_57, spl1_20, spl1_66, sdk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_15 * spk_14[k]
                  + f_3 * pc_y[k] * sdk_50[k];

        t_65[k] = pb_y[k] * spl0_20[k]
                  - f_14 * pc_y[k] * spl1_20[k];

        t_66[k] = pb_x[k] * spl0_66[k]
                  + f_0 * spk_57[k]
                  - f_14 * pc_x[k] * spl1_66[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pc_x, pc_z, spl0_68, spl0_69, spk_59, spk_60, \
                         spl1_68, spl1_69, sdk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * sdk_51[k];

        t_68[k] = pb_x[k] * spl0_68[k]
                  + f_0 * spk_59[k]
                  - f_14 * pc_x[k] * spl1_68[k];

        t_69[k] = pb_x[k] * spl0_69[k]
                  + f_0 * spk_60[k]
                  - f_14 * pc_x[k] * spl1_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_x, pb_y, pc_x, pc_y, spl0_27, spl0_70, spk_20, \
                         spk_61, spl1_27, spl1_70, sdk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_x[k] * spl0_70[k]
                  + f_0 * spk_61[k]
                  - f_14 * pc_x[k] * spl1_70[k];

        t_71[k] = f_15 * spk_20[k]
                  + f_3 * pc_y[k] * sdk_56[k];

        t_72[k] = pb_y[k] * spl0_27[k]
                  - f_14 * pc_y[k] * spl1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, spk_64, spk_65, spk_66, spk_67, \
                         spk_68, sdk_64, sdk_65, sdk_66, sdk_67, \
                         sdk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_15 * spk_64[k]
                  + f_3 * pc_x[k] * sdk_64[k];

        t_74[k] = f_15 * spk_65[k]
                  + f_3 * pc_x[k] * sdk_65[k];

        t_75[k] = f_15 * spk_66[k]
                  + f_3 * pc_x[k] * sdk_66[k];

        t_76[k] = f_15 * spk_67[k]
                  + f_3 * pc_x[k] * sdk_67[k];

        t_77[k] = f_15 * spk_68[k]
                  + f_3 * pc_x[k] * sdk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pc_x, spl0_81, spk_69, spk_70, spk_71, \
                         spl1_81, sdk_69, sdk_70, sdk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_15 * spk_69[k]
                  + f_3 * pc_x[k] * sdk_69[k];

        t_79[k] = f_15 * spk_70[k]
                  + f_3 * pc_x[k] * sdk_70[k];

        t_80[k] = f_15 * spk_71[k]
                  + f_3 * pc_x[k] * sdk_71[k];

        t_81[k] = pb_x[k] * spl0_81[k]
                  - f_14 * pc_x[k] * spl1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_x, pc_x, pc_z, spl0_83, spl0_84, spl0_85, \
                         spl1_83, spl1_84, spl1_85, sdk_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * sdk_64[k];

        t_83[k] = pb_x[k] * spl0_83[k]
                  - f_14 * pc_x[k] * spl1_83[k];

        t_84[k] = pb_x[k] * spl0_84[k]
                  - f_14 * pc_x[k] * spl1_84[k];

        t_85[k] = pb_x[k] * spl0_85[k]
                  - f_14 * pc_x[k] * spl1_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pc_x, pc_y, spl0_86, spl0_87, spl0_89, \
                         spk_35, spl1_86, spl1_87, spl1_89, sdk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * spl0_86[k]
                  - f_14 * pc_x[k] * spl1_86[k];

        t_87[k] = pb_x[k] * spl0_87[k]
                  - f_14 * pc_x[k] * spl1_87[k];

        t_88[k] = f_15 * spk_35[k]
                  + f_3 * pc_y[k] * sdk_71[k];

        t_89[k] = pb_x[k] * spl0_89[k]
                  - f_14 * pc_x[k] * spl1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pb_z, pc_y, pc_z, spl0_0, spl0_3, \
                         spk_0, spl1_0, spl1_3, sdk_72, sdk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * spl0_0[k]
                  - f_14 * pc_z[k] * spl1_0[k];

        t_91[k] = f_3 * pc_y[k] * sdk_72[k];

        t_92[k] = f_15 * spk_0[k]
                  + f_3 * pc_z[k] * sdk_72[k];

        t_93[k] = pb_z[k] * spl0_3[k]
                  - f_14 * pc_z[k] * spl1_3[k];

        t_94[k] = f_3 * pc_y[k] * sdk_74[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_x, pb_z, pc_x, pc_z, spl0_6, spl0_95, spk_3, \
                         spk_77, spl1_6, spl1_95, sdk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_x[k] * spl0_95[k]
                  + f_16 * spk_77[k]
                  - f_14 * pc_x[k] * spl1_95[k];

        t_96[k] = pb_z[k] * spl0_6[k]
                  - f_14 * pc_z[k] * spl1_6[k];

        t_97[k] = f_15 * spk_3[k]
                  + f_3 * pc_z[k] * sdk_75[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pb_x, pb_z, pc_x, pc_y, pc_z, spl0_10, spl0_99, \
                         spk_81, spl1_10, spl1_99, sdk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_3 * pc_y[k] * sdk_77[k];

        t_99[k] = pb_x[k] * spl0_99[k]
                  + f_17 * spk_81[k]
                  - f_14 * pc_x[k] * spl1_99[k];

        t_100[k] = pb_z[k] * spl0_10[k]
                   - f_14 * pc_z[k] * spl1_10[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pb_x, pc_x, pc_y, pc_z, spl0_102, spk_6, spk_84, \
                         spl1_102, sdk_78, sdk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_15 * spk_6[k]
                   + f_3 * pc_z[k] * sdk_78[k];

        t_102[k] = pb_x[k] * spl0_102[k]
                   + f_18 * spk_84[k]
                   - f_14 * pc_x[k] * spl1_102[k];

        t_103[k] = f_3 * pc_y[k] * sdk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_x, pb_z, pc_x, pc_z, spl0_15, spl0_104, \
                         spk_10, spk_86, spl1_15, spl1_104, sdk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * spl0_104[k]
                   + f_18 * spk_86[k]
                   - f_14 * pc_x[k] * spl1_104[k];

        t_105[k] = pb_z[k] * spl0_15[k]
                   - f_14 * pc_z[k] * spl1_15[k];

        t_106[k] = f_15 * spk_10[k]
                   + f_3 * pc_z[k] * sdk_82[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, pc_x, pc_y, spl0_107, spl0_108, spk_89, \
                         spk_90, spl1_107, spl1_108, sdk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pb_x[k] * spl0_107[k]
                   + f_19 * spk_89[k]
                   - f_14 * pc_x[k] * spl1_107[k];

        t_108[k] = pb_x[k] * spl0_108[k]
                   + f_19 * spk_90[k]
                   - f_14 * pc_x[k] * spl1_108[k];

        t_109[k] = f_3 * pc_y[k] * sdk_86[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, pb_z, pc_x, pc_z, spl0_21, spl0_110, \
                         spk_15, spk_92, spl1_21, spl1_110, sdk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * spl0_110[k]
                   + f_19 * spk_92[k]
                   - f_14 * pc_x[k] * spl1_110[k];

        t_111[k] = pb_z[k] * spl0_21[k]
                   - f_14 * pc_z[k] * spl1_21[k];

        t_112[k] = f_15 * spk_15[k]
                   + f_3 * pc_z[k] * sdk_87[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_x, pc_x, spl0_113, spl0_114, spl0_115, \
                         spk_95, spk_96, spk_97, spl1_113, spl1_114, \
                         spl1_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pb_x[k] * spl0_113[k]
                   + f_0 * spk_95[k]
                   - f_14 * pc_x[k] * spl1_113[k];

        t_114[k] = pb_x[k] * spl0_114[k]
                   + f_0 * spk_96[k]
                   - f_14 * pc_x[k] * spl1_114[k];

        t_115[k] = pb_x[k] * spl0_115[k]
                   + f_0 * spk_97[k]
                   - f_14 * pc_x[k] * spl1_115[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_x, pc_x, pc_y, spl0_117, spk_99, \
                         spk_100, spk_101, spl1_117, sdk_92, sdk_100, \
                         sdk_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_3 * pc_y[k] * sdk_92[k];

        t_117[k] = pb_x[k] * spl0_117[k]
                   + f_0 * spk_99[k]
                   - f_14 * pc_x[k] * spl1_117[k];

        t_118[k] = f_15 * spk_100[k]
                   + f_3 * pc_x[k] * sdk_100[k];

        t_119[k] = f_15 * spk_101[k]
                   + f_3 * pc_x[k] * sdk_101[k];
    }
}

static auto
compute_prim_sdl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t spl0,
                                                          const size_t spk, const size_t spl1,
                                                          const size_t sdi0, const size_t sdi1,
                                                          const size_t sdk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_16 = 3.0 / q;
    const auto f_17 = 2.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spl0_48 = buffer.data(spl0 + 48);
    const auto *spl0_51 = buffer.data(spl0 + 51);
    const auto *spl0_55 = buffer.data(spl0 + 55);
    const auto *spl0_60 = buffer.data(spl0 + 60);
    const auto *spl0_66 = buffer.data(spl0 + 66);
    const auto *spl0_81 = buffer.data(spl0 + 81);
    const auto *spl0_90 = buffer.data(spl0 + 90);
    const auto *spl0_95 = buffer.data(spl0 + 95);
    const auto *spl0_99 = buffer.data(spl0 + 99);
    const auto *spl0_104 = buffer.data(spl0 + 104);
    const auto *spl0_110 = buffer.data(spl0 + 110);
    const auto *spl0_117 = buffer.data(spl0 + 117);
    const auto *spl0_126 = buffer.data(spl0 + 126);
    const auto *spl0_128 = buffer.data(spl0 + 128);
    const auto *spl0_129 = buffer.data(spl0 + 129);
    const auto *spl0_130 = buffer.data(spl0 + 130);
    const auto *spl0_131 = buffer.data(spl0 + 131);
    const auto *spl0_132 = buffer.data(spl0 + 132);
    const auto *spl0_134 = buffer.data(spl0 + 134);

    const auto *spk_28 = buffer.data(spk + 28);
    const auto *spk_36 = buffer.data(spk + 36);
    const auto *spk_38 = buffer.data(spk + 38);
    const auto *spk_39 = buffer.data(spk + 39);
    const auto *spk_41 = buffer.data(spk + 41);
    const auto *spk_42 = buffer.data(spk + 42);
    const auto *spk_45 = buffer.data(spk + 45);
    const auto *spk_46 = buffer.data(spk + 46);
    const auto *spk_50 = buffer.data(spk + 50);
    const auto *spk_51 = buffer.data(spk + 51);
    const auto *spk_56 = buffer.data(spk + 56);
    const auto *spk_64 = buffer.data(spk + 64);
    const auto *spk_66 = buffer.data(spk + 66);
    const auto *spk_67 = buffer.data(spk + 67);
    const auto *spk_68 = buffer.data(spk + 68);
    const auto *spk_69 = buffer.data(spk + 69);
    const auto *spk_70 = buffer.data(spk + 70);
    const auto *spk_71 = buffer.data(spk + 71);
    const auto *spk_72 = buffer.data(spk + 72);
    const auto *spk_74 = buffer.data(spk + 74);
    const auto *spk_75 = buffer.data(spk + 75);
    const auto *spk_77 = buffer.data(spk + 77);
    const auto *spk_78 = buffer.data(spk + 78);
    const auto *spk_81 = buffer.data(spk + 81);
    const auto *spk_82 = buffer.data(spk + 82);
    const auto *spk_86 = buffer.data(spk + 86);
    const auto *spk_87 = buffer.data(spk + 87);
    const auto *spk_92 = buffer.data(spk + 92);
    const auto *spk_102 = buffer.data(spk + 102);
    const auto *spk_103 = buffer.data(spk + 103);
    const auto *spk_104 = buffer.data(spk + 104);
    const auto *spk_105 = buffer.data(spk + 105);
    const auto *spk_106 = buffer.data(spk + 106);
    const auto *spk_107 = buffer.data(spk + 107);

    const auto *spl1_48 = buffer.data(spl1 + 48);
    const auto *spl1_51 = buffer.data(spl1 + 51);
    const auto *spl1_55 = buffer.data(spl1 + 55);
    const auto *spl1_60 = buffer.data(spl1 + 60);
    const auto *spl1_66 = buffer.data(spl1 + 66);
    const auto *spl1_81 = buffer.data(spl1 + 81);
    const auto *spl1_90 = buffer.data(spl1 + 90);
    const auto *spl1_95 = buffer.data(spl1 + 95);
    const auto *spl1_99 = buffer.data(spl1 + 99);
    const auto *spl1_104 = buffer.data(spl1 + 104);
    const auto *spl1_110 = buffer.data(spl1 + 110);
    const auto *spl1_117 = buffer.data(spl1 + 117);
    const auto *spl1_126 = buffer.data(spl1 + 126);
    const auto *spl1_128 = buffer.data(spl1 + 128);
    const auto *spl1_129 = buffer.data(spl1 + 129);
    const auto *spl1_130 = buffer.data(spl1 + 130);
    const auto *spl1_131 = buffer.data(spl1 + 131);
    const auto *spl1_132 = buffer.data(spl1 + 132);
    const auto *spl1_134 = buffer.data(spl1 + 134);

    const auto *sdi0_84 = buffer.data(sdi0 + 84);
    const auto *sdi0_87 = buffer.data(sdi0 + 87);
    const auto *sdi0_89 = buffer.data(sdi0 + 89);
    const auto *sdi0_90 = buffer.data(sdi0 + 90);
    const auto *sdi0_93 = buffer.data(sdi0 + 93);
    const auto *sdi0_94 = buffer.data(sdi0 + 94);
    const auto *sdi0_96 = buffer.data(sdi0 + 96);
    const auto *sdi0_98 = buffer.data(sdi0 + 98);
    const auto *sdi0_99 = buffer.data(sdi0 + 99);
    const auto *sdi0_101 = buffer.data(sdi0 + 101);
    const auto *sdi0_102 = buffer.data(sdi0 + 102);
    const auto *sdi0_104 = buffer.data(sdi0 + 104);
    const auto *sdi0_105 = buffer.data(sdi0 + 105);
    const auto *sdi0_107 = buffer.data(sdi0 + 107);
    const auto *sdi0_108 = buffer.data(sdi0 + 108);
    const auto *sdi0_109 = buffer.data(sdi0 + 109);
    const auto *sdi0_110 = buffer.data(sdi0 + 110);
    const auto *sdi0_111 = buffer.data(sdi0 + 111);
    const auto *sdi0_124 = buffer.data(sdi0 + 124);
    const auto *sdi0_129 = buffer.data(sdi0 + 129);
    const auto *sdi0_130 = buffer.data(sdi0 + 130);
    const auto *sdi0_135 = buffer.data(sdi0 + 135);
    const auto *sdi0_136 = buffer.data(sdi0 + 136);
    const auto *sdi0_137 = buffer.data(sdi0 + 137);
    const auto *sdi0_140 = buffer.data(sdi0 + 140);
    const auto *sdi0_143 = buffer.data(sdi0 + 143);
    const auto *sdi0_145 = buffer.data(sdi0 + 145);
    const auto *sdi0_146 = buffer.data(sdi0 + 146);
    const auto *sdi0_149 = buffer.data(sdi0 + 149);
    const auto *sdi0_150 = buffer.data(sdi0 + 150);
    const auto *sdi0_152 = buffer.data(sdi0 + 152);
    const auto *sdi0_154 = buffer.data(sdi0 + 154);
    const auto *sdi0_155 = buffer.data(sdi0 + 155);
    const auto *sdi0_157 = buffer.data(sdi0 + 157);
    const auto *sdi0_158 = buffer.data(sdi0 + 158);
    const auto *sdi0_160 = buffer.data(sdi0 + 160);
    const auto *sdi0_161 = buffer.data(sdi0 + 161);

    const auto *sdi1_84 = buffer.data(sdi1 + 84);
    const auto *sdi1_87 = buffer.data(sdi1 + 87);
    const auto *sdi1_89 = buffer.data(sdi1 + 89);
    const auto *sdi1_90 = buffer.data(sdi1 + 90);
    const auto *sdi1_93 = buffer.data(sdi1 + 93);
    const auto *sdi1_94 = buffer.data(sdi1 + 94);
    const auto *sdi1_96 = buffer.data(sdi1 + 96);
    const auto *sdi1_98 = buffer.data(sdi1 + 98);
    const auto *sdi1_99 = buffer.data(sdi1 + 99);
    const auto *sdi1_101 = buffer.data(sdi1 + 101);
    const auto *sdi1_102 = buffer.data(sdi1 + 102);
    const auto *sdi1_104 = buffer.data(sdi1 + 104);
    const auto *sdi1_105 = buffer.data(sdi1 + 105);
    const auto *sdi1_107 = buffer.data(sdi1 + 107);
    const auto *sdi1_108 = buffer.data(sdi1 + 108);
    const auto *sdi1_109 = buffer.data(sdi1 + 109);
    const auto *sdi1_110 = buffer.data(sdi1 + 110);
    const auto *sdi1_111 = buffer.data(sdi1 + 111);
    const auto *sdi1_124 = buffer.data(sdi1 + 124);
    const auto *sdi1_129 = buffer.data(sdi1 + 129);
    const auto *sdi1_130 = buffer.data(sdi1 + 130);
    const auto *sdi1_135 = buffer.data(sdi1 + 135);
    const auto *sdi1_136 = buffer.data(sdi1 + 136);
    const auto *sdi1_137 = buffer.data(sdi1 + 137);
    const auto *sdi1_140 = buffer.data(sdi1 + 140);
    const auto *sdi1_143 = buffer.data(sdi1 + 143);
    const auto *sdi1_145 = buffer.data(sdi1 + 145);
    const auto *sdi1_146 = buffer.data(sdi1 + 146);
    const auto *sdi1_149 = buffer.data(sdi1 + 149);
    const auto *sdi1_150 = buffer.data(sdi1 + 150);
    const auto *sdi1_152 = buffer.data(sdi1 + 152);
    const auto *sdi1_154 = buffer.data(sdi1 + 154);
    const auto *sdi1_155 = buffer.data(sdi1 + 155);
    const auto *sdi1_157 = buffer.data(sdi1 + 157);
    const auto *sdi1_158 = buffer.data(sdi1 + 158);
    const auto *sdi1_160 = buffer.data(sdi1 + 160);
    const auto *sdi1_161 = buffer.data(sdi1 + 161);

    const auto *sdk_100 = buffer.data(sdk + 100);
    const auto *sdk_102 = buffer.data(sdk + 102);
    const auto *sdk_103 = buffer.data(sdk + 103);
    const auto *sdk_104 = buffer.data(sdk + 104);
    const auto *sdk_105 = buffer.data(sdk + 105);
    const auto *sdk_106 = buffer.data(sdk + 106);
    const auto *sdk_107 = buffer.data(sdk + 107);
    const auto *sdk_108 = buffer.data(sdk + 108);
    const auto *sdk_110 = buffer.data(sdk + 110);
    const auto *sdk_111 = buffer.data(sdk + 111);
    const auto *sdk_113 = buffer.data(sdk + 113);
    const auto *sdk_114 = buffer.data(sdk + 114);
    const auto *sdk_117 = buffer.data(sdk + 117);
    const auto *sdk_118 = buffer.data(sdk + 118);
    const auto *sdk_120 = buffer.data(sdk + 120);
    const auto *sdk_122 = buffer.data(sdk + 122);
    const auto *sdk_123 = buffer.data(sdk + 123);
    const auto *sdk_125 = buffer.data(sdk + 125);
    const auto *sdk_126 = buffer.data(sdk + 126);
    const auto *sdk_128 = buffer.data(sdk + 128);
    const auto *sdk_129 = buffer.data(sdk + 129);
    const auto *sdk_131 = buffer.data(sdk + 131);
    const auto *sdk_132 = buffer.data(sdk + 132);
    const auto *sdk_133 = buffer.data(sdk + 133);
    const auto *sdk_135 = buffer.data(sdk + 135);
    const auto *sdk_136 = buffer.data(sdk + 136);
    const auto *sdk_137 = buffer.data(sdk + 137);
    const auto *sdk_138 = buffer.data(sdk + 138);
    const auto *sdk_139 = buffer.data(sdk + 139);
    const auto *sdk_140 = buffer.data(sdk + 140);
    const auto *sdk_141 = buffer.data(sdk + 141);
    const auto *sdk_142 = buffer.data(sdk + 142);
    const auto *sdk_143 = buffer.data(sdk + 143);
    const auto *sdk_144 = buffer.data(sdk + 144);
    const auto *sdk_146 = buffer.data(sdk + 146);
    const auto *sdk_147 = buffer.data(sdk + 147);
    const auto *sdk_149 = buffer.data(sdk + 149);
    const auto *sdk_150 = buffer.data(sdk + 150);
    const auto *sdk_153 = buffer.data(sdk + 153);
    const auto *sdk_154 = buffer.data(sdk + 154);
    const auto *sdk_156 = buffer.data(sdk + 156);
    const auto *sdk_158 = buffer.data(sdk + 158);
    const auto *sdk_159 = buffer.data(sdk + 159);
    const auto *sdk_161 = buffer.data(sdk + 161);
    const auto *sdk_162 = buffer.data(sdk + 162);
    const auto *sdk_164 = buffer.data(sdk + 164);
    const auto *sdk_167 = buffer.data(sdk + 167);
    const auto *sdk_168 = buffer.data(sdk + 168);
    const auto *sdk_169 = buffer.data(sdk + 169);
    const auto *sdk_172 = buffer.data(sdk + 172);
    const auto *sdk_173 = buffer.data(sdk + 173);
    const auto *sdk_174 = buffer.data(sdk + 174);
    const auto *sdk_175 = buffer.data(sdk + 175);
    const auto *sdk_176 = buffer.data(sdk + 176);
    const auto *sdk_177 = buffer.data(sdk + 177);
    const auto *sdk_178 = buffer.data(sdk + 178);
    const auto *sdk_179 = buffer.data(sdk + 179);
    const auto *sdk_180 = buffer.data(sdk + 180);
    const auto *sdk_182 = buffer.data(sdk + 182);
    const auto *sdk_183 = buffer.data(sdk + 183);
    const auto *sdk_185 = buffer.data(sdk + 185);
    const auto *sdk_186 = buffer.data(sdk + 186);
    const auto *sdk_189 = buffer.data(sdk + 189);
    const auto *sdk_190 = buffer.data(sdk + 190);
    const auto *sdk_192 = buffer.data(sdk + 192);
    const auto *sdk_194 = buffer.data(sdk + 194);
    const auto *sdk_195 = buffer.data(sdk + 195);
    const auto *sdk_197 = buffer.data(sdk + 197);
    const auto *sdk_198 = buffer.data(sdk + 198);
    const auto *sdk_200 = buffer.data(sdk + 200);
    const auto *sdk_201 = buffer.data(sdk + 201);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, spk_102, spk_103, spk_104, \
                         spk_105, spk_106, sdk_102, sdk_103, sdk_104, sdk_105, \
                         sdk_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_15 * spk_102[k]
                   + f_3 * pc_x[k] * sdk_102[k];

        t_121[k] = f_15 * spk_103[k]
                   + f_3 * pc_x[k] * sdk_103[k];

        t_122[k] = f_15 * spk_104[k]
                   + f_3 * pc_x[k] * sdk_104[k];

        t_123[k] = f_15 * spk_105[k]
                   + f_3 * pc_x[k] * sdk_105[k];

        t_124[k] = f_15 * spk_106[k]
                   + f_3 * pc_x[k] * sdk_106[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_x, pc_x, pc_z, spl0_126, spl0_128, \
                         spk_28, spk_107, spl1_126, spl1_128, sdk_100, \
                         sdk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_15 * spk_107[k]
                   + f_3 * pc_x[k] * sdk_107[k];

        t_126[k] = pb_x[k] * spl0_126[k]
                   - f_14 * pc_x[k] * spl1_126[k];

        t_127[k] = f_15 * spk_28[k]
                   + f_3 * pc_z[k] * sdk_100[k];

        t_128[k] = pb_x[k] * spl0_128[k]
                   - f_14 * pc_x[k] * spl1_128[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pb_x, pc_x, spl0_129, spl0_130, spl0_131, \
                         spl0_132, spl1_129, spl1_130, spl1_131, \
                         spl1_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pb_x[k] * spl0_129[k]
                   - f_14 * pc_x[k] * spl1_129[k];

        t_130[k] = pb_x[k] * spl0_130[k]
                   - f_14 * pc_x[k] * spl1_130[k];

        t_131[k] = pb_x[k] * spl0_131[k]
                   - f_14 * pc_x[k] * spl1_131[k];

        t_132[k] = pb_x[k] * spl0_132[k]
                   - f_14 * pc_x[k] * spl1_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pb_x, pc_x, pc_y, pc_z, spl0_134, \
                         spk_36, spl1_134, sdi0_84, sdi1_84, sdk_107, \
                         sdk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_y[k] * sdk_107[k];

        t_134[k] = pb_x[k] * spl0_134[k]
                   - f_14 * pc_x[k] * spl1_134[k];

        t_135[k] = f_1 * sdi0_84[k]
                   - f_2 * sdi1_84[k]
                   + f_3 * pc_x[k] * sdk_108[k];

        t_136[k] = f_0 * spk_36[k]
                   + f_3 * pc_y[k] * sdk_108[k];

        t_137[k] = f_3 * pc_z[k] * sdk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, spk_38, sdi0_87, sdi0_89, sdi1_87, \
                         sdi1_89, sdk_110, sdk_111, sdk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_4 * sdi0_87[k]
                   - f_5 * sdi1_87[k]
                   + f_3 * pc_x[k] * sdk_111[k];

        t_139[k] = f_0 * spk_38[k]
                   + f_3 * pc_y[k] * sdk_110[k];

        t_140[k] = f_4 * sdi0_89[k]
                   - f_5 * sdi1_89[k]
                   + f_3 * pc_x[k] * sdk_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, spk_41, sdi0_90, \
                         sdi0_93, sdi1_90, sdi1_93, sdk_111, sdk_113, sdk_114, \
                         sdk_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_6 * sdi0_90[k]
                   - f_7 * sdi1_90[k]
                   + f_3 * pc_x[k] * sdk_114[k];

        t_142[k] = f_3 * pc_z[k] * sdk_111[k];

        t_143[k] = f_0 * spk_41[k]
                   + f_3 * pc_y[k] * sdk_113[k];

        t_144[k] = f_6 * sdi0_93[k]
                   - f_7 * sdi1_93[k]
                   + f_3 * pc_x[k] * sdk_117[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, pc_z, spk_45, sdi0_94, \
                         sdi0_96, sdi1_94, sdi1_96, sdk_114, sdk_117, sdk_118, \
                         sdk_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * sdi0_94[k]
                   - f_9 * sdi1_94[k]
                   + f_3 * pc_x[k] * sdk_118[k];

        t_146[k] = f_3 * pc_z[k] * sdk_114[k];

        t_147[k] = f_8 * sdi0_96[k]
                   - f_9 * sdi1_96[k]
                   + f_3 * pc_x[k] * sdk_120[k];

        t_148[k] = f_0 * spk_45[k]
                   + f_3 * pc_y[k] * sdk_117[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_z, sdi0_98, sdi0_99, sdi0_101, \
                         sdi1_98, sdi1_99, sdi1_101, sdk_118, sdk_122, sdk_123, \
                         sdk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_8 * sdi0_98[k]
                   - f_9 * sdi1_98[k]
                   + f_3 * pc_x[k] * sdk_122[k];

        t_150[k] = f_10 * sdi0_99[k]
                   - f_11 * sdi1_99[k]
                   + f_3 * pc_x[k] * sdk_123[k];

        t_151[k] = f_3 * pc_z[k] * sdk_118[k];

        t_152[k] = f_10 * sdi0_101[k]
                   - f_11 * sdi1_101[k]
                   + f_3 * pc_x[k] * sdk_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, spk_50, sdi0_102, sdi0_104, \
                         sdi1_102, sdi1_104, sdk_122, sdk_126, \
                         sdk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_10 * sdi0_102[k]
                   - f_11 * sdi1_102[k]
                   + f_3 * pc_x[k] * sdk_126[k];

        t_154[k] = f_0 * spk_50[k]
                   + f_3 * pc_y[k] * sdk_122[k];

        t_155[k] = f_10 * sdi0_104[k]
                   - f_11 * sdi1_104[k]
                   + f_3 * pc_x[k] * sdk_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, pc_z, sdi0_105, sdi0_107, sdi0_108, \
                         sdi1_105, sdi1_107, sdi1_108, sdk_123, sdk_129, sdk_131, \
                         sdk_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_12 * sdi0_105[k]
                   - f_13 * sdi1_105[k]
                   + f_3 * pc_x[k] * sdk_129[k];

        t_157[k] = f_3 * pc_z[k] * sdk_123[k];

        t_158[k] = f_12 * sdi0_107[k]
                   - f_13 * sdi1_107[k]
                   + f_3 * pc_x[k] * sdk_131[k];

        t_159[k] = f_12 * sdi0_108[k]
                   - f_13 * sdi1_108[k]
                   + f_3 * pc_x[k] * sdk_132[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pc_x, pc_y, spk_56, sdi0_109, sdi0_111, \
                         sdi1_109, sdi1_111, sdk_128, sdk_133, sdk_135, \
                         sdk_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_12 * sdi0_109[k]
                   - f_13 * sdi1_109[k]
                   + f_3 * pc_x[k] * sdk_133[k];

        t_161[k] = f_0 * spk_56[k]
                   + f_3 * pc_y[k] * sdk_128[k];

        t_162[k] = f_12 * sdi0_111[k]
                   - f_13 * sdi1_111[k]
                   + f_3 * pc_x[k] * sdk_135[k];

        t_163[k] = f_3 * pc_x[k] * sdk_136[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, t_169, t_170, pc_x, sdk_137, \
                         sdk_138, sdk_139, sdk_140, sdk_141, sdk_142, \
                         sdk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_3 * pc_x[k] * sdk_137[k];

        t_165[k] = f_3 * pc_x[k] * sdk_138[k];

        t_166[k] = f_3 * pc_x[k] * sdk_139[k];

        t_167[k] = f_3 * pc_x[k] * sdk_140[k];

        t_168[k] = f_3 * pc_x[k] * sdk_141[k];

        t_169[k] = f_3 * pc_x[k] * sdk_142[k];

        t_170[k] = f_3 * pc_x[k] * sdk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, spk_64, spk_66, sdi0_105, sdi0_107, \
                         sdi1_105, sdi1_107, sdk_136, sdk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_0 * spk_64[k]
                   + f_1 * sdi0_105[k]
                   - f_2 * sdi1_105[k]
                   + f_3 * pc_y[k] * sdk_136[k];

        t_172[k] = f_3 * pc_z[k] * sdk_136[k];

        t_173[k] = f_0 * spk_66[k]
                   + f_4 * sdi0_107[k]
                   - f_5 * sdi1_107[k]
                   + f_3 * pc_y[k] * sdk_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, spk_67, spk_68, spk_69, sdi0_108, \
                         sdi0_109, sdi0_110, sdi1_108, sdi1_109, sdi1_110, sdk_139, sdk_140, \
                         sdk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_0 * spk_67[k]
                   + f_6 * sdi0_108[k]
                   - f_7 * sdi1_108[k]
                   + f_3 * pc_y[k] * sdk_139[k];

        t_175[k] = f_0 * spk_68[k]
                   + f_8 * sdi0_109[k]
                   - f_9 * sdi1_109[k]
                   + f_3 * pc_y[k] * sdk_140[k];

        t_176[k] = f_0 * spk_69[k]
                   + f_10 * sdi0_110[k]
                   - f_11 * sdi1_110[k]
                   + f_3 * pc_y[k] * sdk_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, spl0_90, spk_70, \
                         spk_71, spl1_90, sdi0_111, sdi1_111, sdk_142, \
                         sdk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_0 * spk_70[k]
                   + f_12 * sdi0_111[k]
                   - f_13 * sdi1_111[k]
                   + f_3 * pc_y[k] * sdk_142[k];

        t_178[k] = f_0 * spk_71[k]
                   + f_3 * pc_y[k] * sdk_143[k];

        t_179[k] = f_1 * sdi0_111[k]
                   - f_2 * sdi1_111[k]
                   + f_3 * pc_z[k] * sdk_143[k];

        t_180[k] = pb_y[k] * spl0_90[k]
                   - f_14 * pc_y[k] * spl1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, spl0_48, spk_36, \
                         spk_72, spk_74, spl1_48, sdk_144, sdk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * spk_72[k]
                   + f_3 * pc_y[k] * sdk_144[k];

        t_182[k] = f_15 * spk_36[k]
                   + f_3 * pc_z[k] * sdk_144[k];

        t_183[k] = pb_z[k] * spl0_48[k]
                   - f_14 * pc_z[k] * spl1_48[k];

        t_184[k] = f_15 * spk_74[k]
                   + f_3 * pc_y[k] * sdk_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, spl0_51, spl0_95, \
                         spk_39, spk_77, spl1_51, spl1_95, sdk_147, \
                         sdk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * spl0_95[k]
                   - f_14 * pc_y[k] * spl1_95[k];

        t_186[k] = pb_z[k] * spl0_51[k]
                   - f_14 * pc_z[k] * spl1_51[k];

        t_187[k] = f_15 * spk_39[k]
                   + f_3 * pc_z[k] * sdk_147[k];

        t_188[k] = f_15 * spk_77[k]
                   + f_3 * pc_y[k] * sdk_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, spl0_55, spl0_99, \
                         spk_42, spl1_55, spl1_99, sdk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * spl0_99[k]
                   - f_14 * pc_y[k] * spl1_99[k];

        t_190[k] = pb_z[k] * spl0_55[k]
                   - f_14 * pc_z[k] * spl1_55[k];

        t_191[k] = f_15 * spk_42[k]
                   + f_3 * pc_z[k] * sdk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_x, pc_y, spl0_104, spk_81, spl1_104, \
                         sdi0_124, sdi1_124, sdk_153, sdk_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * sdi0_124[k]
                   - f_9 * sdi1_124[k]
                   + f_3 * pc_x[k] * sdk_156[k];

        t_193[k] = f_15 * spk_81[k]
                   + f_3 * pc_y[k] * sdk_153[k];

        t_194[k] = pb_y[k] * spl0_104[k]
                   - f_14 * pc_y[k] * spl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_z, pc_x, pc_z, spl0_60, spk_46, spl1_60, \
                         sdi0_129, sdi1_129, sdk_154, sdk_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * spl0_60[k]
                   - f_14 * pc_z[k] * spl1_60[k];

        t_196[k] = f_15 * spk_46[k]
                   + f_3 * pc_z[k] * sdk_154[k];

        t_197[k] = f_10 * sdi0_129[k]
                   - f_11 * sdi1_129[k]
                   + f_3 * pc_x[k] * sdk_161[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_x, pc_y, spl0_110, spk_86, spl1_110, \
                         sdi0_130, sdi1_130, sdk_158, sdk_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_10 * sdi0_130[k]
                   - f_11 * sdi1_130[k]
                   + f_3 * pc_x[k] * sdk_162[k];

        t_199[k] = f_15 * spk_86[k]
                   + f_3 * pc_y[k] * sdk_158[k];

        t_200[k] = pb_y[k] * spl0_110[k]
                   - f_14 * pc_y[k] * spl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_z, pc_x, pc_z, spl0_66, spk_51, spl1_66, \
                         sdi0_135, sdi1_135, sdk_159, sdk_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * spl0_66[k]
                   - f_14 * pc_z[k] * spl1_66[k];

        t_202[k] = f_15 * spk_51[k]
                   + f_3 * pc_z[k] * sdk_159[k];

        t_203[k] = f_12 * sdi0_135[k]
                   - f_13 * sdi1_135[k]
                   + f_3 * pc_x[k] * sdk_167[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, spk_92, sdi0_136, sdi0_137, \
                         sdi1_136, sdi1_137, sdk_164, sdk_168, \
                         sdk_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * sdi0_136[k]
                   - f_13 * sdi1_136[k]
                   + f_3 * pc_x[k] * sdk_168[k];

        t_205[k] = f_12 * sdi0_137[k]
                   - f_13 * sdi1_137[k]
                   + f_3 * pc_x[k] * sdk_169[k];

        t_206[k] = f_15 * spk_92[k]
                   + f_3 * pc_y[k] * sdk_164[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, pb_y, pc_x, pc_y, spl0_117, \
                         spl1_117, sdk_172, sdk_173, sdk_174, sdk_175, \
                         sdk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pb_y[k] * spl0_117[k]
                   - f_14 * pc_y[k] * spl1_117[k];

        t_208[k] = f_3 * pc_x[k] * sdk_172[k];

        t_209[k] = f_3 * pc_x[k] * sdk_173[k];

        t_210[k] = f_3 * pc_x[k] * sdk_174[k];

        t_211[k] = f_3 * pc_x[k] * sdk_175[k];

        t_212[k] = f_3 * pc_x[k] * sdk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pb_z, pc_x, pc_z, spl0_81, spk_64, \
                         spl1_81, sdk_172, sdk_177, sdk_178, sdk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_3 * pc_x[k] * sdk_177[k];

        t_214[k] = f_3 * pc_x[k] * sdk_178[k];

        t_215[k] = f_3 * pc_x[k] * sdk_179[k];

        t_216[k] = pb_z[k] * spl0_81[k]
                   - f_14 * pc_z[k] * spl1_81[k];

        t_217[k] = f_15 * spk_64[k]
                   + f_3 * pc_z[k] * sdk_172[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pb_y, pc_y, spl0_128, spl0_129, spl0_130, \
                         spk_102, spk_103, spk_104, spl1_128, spl1_129, \
                         spl1_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = pb_y[k] * spl0_128[k]
                   + f_16 * spk_102[k]
                   - f_14 * pc_y[k] * spl1_128[k];

        t_219[k] = pb_y[k] * spl0_129[k]
                   + f_17 * spk_103[k]
                   - f_14 * pc_y[k] * spl1_129[k];

        t_220[k] = pb_y[k] * spl0_130[k]
                   + f_18 * spk_104[k]
                   - f_14 * pc_y[k] * spl1_130[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pb_y, pc_y, spl0_131, spl0_132, spl0_134, \
                         spk_105, spk_106, spk_107, spl1_131, spl1_132, spl1_134, \
                         sdk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pb_y[k] * spl0_131[k]
                   + f_19 * spk_105[k]
                   - f_14 * pc_y[k] * spl1_131[k];

        t_222[k] = pb_y[k] * spl0_132[k]
                   + f_0 * spk_106[k]
                   - f_14 * pc_y[k] * spl1_132[k];

        t_223[k] = f_15 * spk_107[k]
                   + f_3 * pc_y[k] * sdk_179[k];

        t_224[k] = pb_y[k] * spl0_134[k]
                   - f_14 * pc_y[k] * spl1_134[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, spk_72, \
                         sdi0_140, sdi0_143, sdi1_140, sdi1_143, sdk_180, sdk_182, \
                         sdk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_1 * sdi0_140[k]
                   - f_2 * sdi1_140[k]
                   + f_3 * pc_x[k] * sdk_180[k];

        t_226[k] = f_3 * pc_y[k] * sdk_180[k];

        t_227[k] = f_0 * spk_72[k]
                   + f_3 * pc_z[k] * sdk_180[k];

        t_228[k] = f_4 * sdi0_143[k]
                   - f_5 * sdi1_143[k]
                   + f_3 * pc_x[k] * sdk_183[k];

        t_229[k] = f_3 * pc_y[k] * sdk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, spk_75, sdi0_145, \
                         sdi0_146, sdi1_145, sdi1_146, sdk_183, sdk_185, \
                         sdk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * sdi0_145[k]
                   - f_5 * sdi1_145[k]
                   + f_3 * pc_x[k] * sdk_185[k];

        t_231[k] = f_6 * sdi0_146[k]
                   - f_7 * sdi1_146[k]
                   + f_3 * pc_x[k] * sdk_186[k];

        t_232[k] = f_0 * spk_75[k]
                   + f_3 * pc_z[k] * sdk_183[k];

        t_233[k] = f_3 * pc_y[k] * sdk_185[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_z, spk_78, sdi0_149, sdi0_150, \
                         sdi1_149, sdi1_150, sdk_186, sdk_189, \
                         sdk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_6 * sdi0_149[k]
                   - f_7 * sdi1_149[k]
                   + f_3 * pc_x[k] * sdk_189[k];

        t_235[k] = f_8 * sdi0_150[k]
                   - f_9 * sdi1_150[k]
                   + f_3 * pc_x[k] * sdk_190[k];

        t_236[k] = f_0 * spk_78[k]
                   + f_3 * pc_z[k] * sdk_186[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pc_x, pc_y, sdi0_152, sdi0_154, sdi0_155, \
                         sdi1_152, sdi1_154, sdi1_155, sdk_189, sdk_192, sdk_194, \
                         sdk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_8 * sdi0_152[k]
                   - f_9 * sdi1_152[k]
                   + f_3 * pc_x[k] * sdk_192[k];

        t_238[k] = f_3 * pc_y[k] * sdk_189[k];

        t_239[k] = f_8 * sdi0_154[k]
                   - f_9 * sdi1_154[k]
                   + f_3 * pc_x[k] * sdk_194[k];

        t_240[k] = f_10 * sdi0_155[k]
                   - f_11 * sdi1_155[k]
                   + f_3 * pc_x[k] * sdk_195[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, spk_82, sdi0_157, \
                         sdi0_158, sdi1_157, sdi1_158, sdk_190, sdk_194, sdk_197, \
                         sdk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_0 * spk_82[k]
                   + f_3 * pc_z[k] * sdk_190[k];

        t_242[k] = f_10 * sdi0_157[k]
                   - f_11 * sdi1_157[k]
                   + f_3 * pc_x[k] * sdk_197[k];

        t_243[k] = f_10 * sdi0_158[k]
                   - f_11 * sdi1_158[k]
                   + f_3 * pc_x[k] * sdk_198[k];

        t_244[k] = f_3 * pc_y[k] * sdk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, spk_87, sdi0_160, sdi0_161, \
                         sdi1_160, sdi1_161, sdk_195, sdk_200, \
                         sdk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_10 * sdi0_160[k]
                   - f_11 * sdi1_160[k]
                   + f_3 * pc_x[k] * sdk_200[k];

        t_246[k] = f_12 * sdi0_161[k]
                   - f_13 * sdi1_161[k]
                   + f_3 * pc_x[k] * sdk_201[k];

        t_247[k] = f_0 * spk_87[k]
                   + f_3 * pc_z[k] * sdk_195[k];
    }
}

static auto
compute_prim_sdl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t spk, const size_t sdi0,
                                                          const size_t sdi1, const size_t sdk,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spk_100 = buffer.data(spk + 100);
    const auto *spk_107 = buffer.data(spk + 107);

    const auto *sdi0_161 = buffer.data(sdi0 + 161);
    const auto *sdi0_163 = buffer.data(sdi0 + 163);
    const auto *sdi0_164 = buffer.data(sdi0 + 164);
    const auto *sdi0_165 = buffer.data(sdi0 + 165);
    const auto *sdi0_166 = buffer.data(sdi0 + 166);
    const auto *sdi0_167 = buffer.data(sdi0 + 167);

    const auto *sdi1_161 = buffer.data(sdi1 + 161);
    const auto *sdi1_163 = buffer.data(sdi1 + 163);
    const auto *sdi1_164 = buffer.data(sdi1 + 164);
    const auto *sdi1_165 = buffer.data(sdi1 + 165);
    const auto *sdi1_166 = buffer.data(sdi1 + 166);
    const auto *sdi1_167 = buffer.data(sdi1 + 167);

    const auto *sdk_200 = buffer.data(sdk + 200);
    const auto *sdk_203 = buffer.data(sdk + 203);
    const auto *sdk_204 = buffer.data(sdk + 204);
    const auto *sdk_205 = buffer.data(sdk + 205);
    const auto *sdk_207 = buffer.data(sdk + 207);
    const auto *sdk_208 = buffer.data(sdk + 208);
    const auto *sdk_209 = buffer.data(sdk + 209);
    const auto *sdk_210 = buffer.data(sdk + 210);
    const auto *sdk_211 = buffer.data(sdk + 211);
    const auto *sdk_212 = buffer.data(sdk + 212);
    const auto *sdk_213 = buffer.data(sdk + 213);
    const auto *sdk_214 = buffer.data(sdk + 214);
    const auto *sdk_215 = buffer.data(sdk + 215);

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pc_x, pc_y, sdi0_163, sdi0_164, sdi0_165, \
                         sdi1_163, sdi1_164, sdi1_165, sdk_200, sdk_203, sdk_204, \
                         sdk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_12 * sdi0_163[k]
                   - f_13 * sdi1_163[k]
                   + f_3 * pc_x[k] * sdk_203[k];

        t_249[k] = f_12 * sdi0_164[k]
                   - f_13 * sdi1_164[k]
                   + f_3 * pc_x[k] * sdk_204[k];

        t_250[k] = f_12 * sdi0_165[k]
                   - f_13 * sdi1_165[k]
                   + f_3 * pc_x[k] * sdk_205[k];

        t_251[k] = f_3 * pc_y[k] * sdk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, t_257, pc_x, sdi0_167, sdi1_167, \
                         sdk_207, sdk_208, sdk_209, sdk_210, sdk_211, \
                         sdk_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_12 * sdi0_167[k]
                   - f_13 * sdi1_167[k]
                   + f_3 * pc_x[k] * sdk_207[k];

        t_253[k] = f_3 * pc_x[k] * sdk_208[k];

        t_254[k] = f_3 * pc_x[k] * sdk_209[k];

        t_255[k] = f_3 * pc_x[k] * sdk_210[k];

        t_256[k] = f_3 * pc_x[k] * sdk_211[k];

        t_257[k] = f_3 * pc_x[k] * sdk_212[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pc_x, pc_y, pc_z, spk_100, \
                         sdi0_161, sdi1_161, sdk_208, sdk_213, sdk_214, \
                         sdk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_3 * pc_x[k] * sdk_213[k];

        t_259[k] = f_3 * pc_x[k] * sdk_214[k];

        t_260[k] = f_3 * pc_x[k] * sdk_215[k];

        t_261[k] = f_1 * sdi0_161[k]
                   - f_2 * sdi1_161[k]
                   + f_3 * pc_y[k] * sdk_208[k];

        t_262[k] = f_0 * spk_100[k]
                   + f_3 * pc_z[k] * sdk_208[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pc_y, sdi0_163, sdi0_164, sdi0_165, sdi1_163, \
                         sdi1_164, sdi1_165, sdk_210, sdk_211, \
                         sdk_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_4 * sdi0_163[k]
                   - f_5 * sdi1_163[k]
                   + f_3 * pc_y[k] * sdk_210[k];

        t_264[k] = f_6 * sdi0_164[k]
                   - f_7 * sdi1_164[k]
                   + f_3 * pc_y[k] * sdk_211[k];

        t_265[k] = f_8 * sdi0_165[k]
                   - f_9 * sdi1_165[k]
                   + f_3 * pc_y[k] * sdk_212[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pc_y, pc_z, spk_107, sdi0_166, sdi0_167, \
                         sdi1_166, sdi1_167, sdk_213, sdk_214, \
                         sdk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_10 * sdi0_166[k]
                   - f_11 * sdi1_166[k]
                   + f_3 * pc_y[k] * sdk_213[k];

        t_267[k] = f_12 * sdi0_167[k]
                   - f_13 * sdi1_167[k]
                   + f_3 * pc_y[k] * sdk_214[k];

        t_268[k] = f_3 * pc_y[k] * sdk_215[k];

        t_269[k] = f_0 * spk_107[k]
                   + f_1 * sdi0_167[k]
                   - f_2 * sdi1_167[k]
                   + f_3 * pc_z[k] * sdk_215[k];
    }
}

auto
compute_prim_sdl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t spl0, const size_t spk,
                                                   const size_t spl1, const size_t sdi0,
                                                   const size_t sdi1, const size_t sdk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sdl_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, spl0, spk,
                                                              spl1, sdi0, sdi1, sdk, ncols,
                                                              gamma, p, q);

    compute_prim_sdl_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, spl0, spk,
                                                              spl1, sdi0, sdi1, sdk, ncols,
                                                              gamma, p, q);

    compute_prim_sdl_three_center_electron_repulsion_0_piece2(buffer, target, pc, spk, sdi0,
                                                              sdi1, sdk, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
