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


#include "SimdThreeCenterElectronRepulsionVrrRecPPI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ppi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t spi0, const size_t sph,
                                                          const size_t spi1, const size_t psi0,
                                                          const size_t psh, const size_t psi1,
                                                          const size_t ppg0, const size_t ppg1,
                                                          const size_t pph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
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
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 1.0 / q;
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);

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
    auto *t_129 = buffer.data(target + 129);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spi0_0 = buffer.data(spi0 + 0);
    const auto *spi0_1 = buffer.data(spi0 + 1);
    const auto *spi0_3 = buffer.data(spi0 + 3);
    const auto *spi0_5 = buffer.data(spi0 + 5);
    const auto *spi0_6 = buffer.data(spi0 + 6);
    const auto *spi0_8 = buffer.data(spi0 + 8);
    const auto *spi0_9 = buffer.data(spi0 + 9);
    const auto *spi0_10 = buffer.data(spi0 + 10);
    const auto *spi0_12 = buffer.data(spi0 + 12);
    const auto *spi0_13 = buffer.data(spi0 + 13);
    const auto *spi0_14 = buffer.data(spi0 + 14);
    const auto *spi0_27 = buffer.data(spi0 + 27);
    const auto *spi0_31 = buffer.data(spi0 + 31);
    const auto *spi0_34 = buffer.data(spi0 + 34);
    const auto *spi0_38 = buffer.data(spi0 + 38);
    const auto *spi0_40 = buffer.data(spi0 + 40);
    const auto *spi0_49 = buffer.data(spi0 + 49);
    const auto *spi0_51 = buffer.data(spi0 + 51);
    const auto *spi0_52 = buffer.data(spi0 + 52);
    const auto *spi0_53 = buffer.data(spi0 + 53);
    const auto *spi0_55 = buffer.data(spi0 + 55);
    const auto *spi0_61 = buffer.data(spi0 + 61);
    const auto *spi0_65 = buffer.data(spi0 + 65);
    const auto *spi0_70 = buffer.data(spi0 + 70);
    const auto *spi0_77 = buffer.data(spi0 + 77);
    const auto *spi0_79 = buffer.data(spi0 + 79);
    const auto *spi0_80 = buffer.data(spi0 + 80);
    const auto *spi0_81 = buffer.data(spi0 + 81);
    const auto *spi0_83 = buffer.data(spi0 + 83);

    const auto *sph_0 = buffer.data(sph + 0);
    const auto *sph_1 = buffer.data(sph + 1);
    const auto *sph_3 = buffer.data(sph + 3);
    const auto *sph_5 = buffer.data(sph + 5);
    const auto *sph_6 = buffer.data(sph + 6);
    const auto *sph_8 = buffer.data(sph + 8);
    const auto *sph_9 = buffer.data(sph + 9);
    const auto *sph_15 = buffer.data(sph + 15);
    const auto *sph_17 = buffer.data(sph + 17);
    const auto *sph_18 = buffer.data(sph + 18);
    const auto *sph_20 = buffer.data(sph + 20);
    const auto *sph_24 = buffer.data(sph + 24);
    const auto *sph_27 = buffer.data(sph + 27);
    const auto *sph_31 = buffer.data(sph + 31);
    const auto *sph_33 = buffer.data(sph + 33);
    const auto *sph_36 = buffer.data(sph + 36);
    const auto *sph_38 = buffer.data(sph + 38);
    const auto *sph_39 = buffer.data(sph + 39);
    const auto *sph_41 = buffer.data(sph + 41);
    const auto *sph_47 = buffer.data(sph + 47);
    const auto *sph_51 = buffer.data(sph + 51);
    const auto *sph_56 = buffer.data(sph + 56);
    const auto *sph_57 = buffer.data(sph + 57);
    const auto *sph_59 = buffer.data(sph + 59);
    const auto *sph_60 = buffer.data(sph + 60);
    const auto *sph_62 = buffer.data(sph + 62);

    const auto *spi1_0 = buffer.data(spi1 + 0);
    const auto *spi1_1 = buffer.data(spi1 + 1);
    const auto *spi1_3 = buffer.data(spi1 + 3);
    const auto *spi1_5 = buffer.data(spi1 + 5);
    const auto *spi1_6 = buffer.data(spi1 + 6);
    const auto *spi1_8 = buffer.data(spi1 + 8);
    const auto *spi1_9 = buffer.data(spi1 + 9);
    const auto *spi1_10 = buffer.data(spi1 + 10);
    const auto *spi1_12 = buffer.data(spi1 + 12);
    const auto *spi1_13 = buffer.data(spi1 + 13);
    const auto *spi1_14 = buffer.data(spi1 + 14);
    const auto *spi1_27 = buffer.data(spi1 + 27);
    const auto *spi1_31 = buffer.data(spi1 + 31);
    const auto *spi1_34 = buffer.data(spi1 + 34);
    const auto *spi1_38 = buffer.data(spi1 + 38);
    const auto *spi1_40 = buffer.data(spi1 + 40);
    const auto *spi1_49 = buffer.data(spi1 + 49);
    const auto *spi1_51 = buffer.data(spi1 + 51);
    const auto *spi1_52 = buffer.data(spi1 + 52);
    const auto *spi1_53 = buffer.data(spi1 + 53);
    const auto *spi1_55 = buffer.data(spi1 + 55);
    const auto *spi1_61 = buffer.data(spi1 + 61);
    const auto *spi1_65 = buffer.data(spi1 + 65);
    const auto *spi1_70 = buffer.data(spi1 + 70);
    const auto *spi1_77 = buffer.data(spi1 + 77);
    const auto *spi1_79 = buffer.data(spi1 + 79);
    const auto *spi1_80 = buffer.data(spi1 + 80);
    const auto *spi1_81 = buffer.data(spi1 + 81);
    const auto *spi1_83 = buffer.data(spi1 + 83);

    const auto *psi0_0 = buffer.data(psi0 + 0);
    const auto *psi0_3 = buffer.data(psi0 + 3);
    const auto *psi0_5 = buffer.data(psi0 + 5);
    const auto *psi0_6 = buffer.data(psi0 + 6);
    const auto *psi0_9 = buffer.data(psi0 + 9);
    const auto *psi0_10 = buffer.data(psi0 + 10);
    const auto *psi0_14 = buffer.data(psi0 + 14);
    const auto *psi0_49 = buffer.data(psi0 + 49);
    const auto *psi0_51 = buffer.data(psi0 + 51);
    const auto *psi0_52 = buffer.data(psi0 + 52);
    const auto *psi0_53 = buffer.data(psi0 + 53);

    const auto *psh_0 = buffer.data(psh + 0);
    const auto *psh_2 = buffer.data(psh + 2);
    const auto *psh_3 = buffer.data(psh + 3);
    const auto *psh_5 = buffer.data(psh + 5);
    const auto *psh_6 = buffer.data(psh + 6);
    const auto *psh_9 = buffer.data(psh + 9);
    const auto *psh_10 = buffer.data(psh + 10);
    const auto *psh_14 = buffer.data(psh + 14);
    const auto *psh_15 = buffer.data(psh + 15);
    const auto *psh_17 = buffer.data(psh + 17);
    const auto *psh_18 = buffer.data(psh + 18);
    const auto *psh_20 = buffer.data(psh + 20);
    const auto *psh_36 = buffer.data(psh + 36);
    const auto *psh_37 = buffer.data(psh + 37);
    const auto *psh_38 = buffer.data(psh + 38);
    const auto *psh_39 = buffer.data(psh + 39);
    const auto *psh_40 = buffer.data(psh + 40);
    const auto *psh_41 = buffer.data(psh + 41);

    const auto *psi1_0 = buffer.data(psi1 + 0);
    const auto *psi1_3 = buffer.data(psi1 + 3);
    const auto *psi1_5 = buffer.data(psi1 + 5);
    const auto *psi1_6 = buffer.data(psi1 + 6);
    const auto *psi1_9 = buffer.data(psi1 + 9);
    const auto *psi1_10 = buffer.data(psi1 + 10);
    const auto *psi1_14 = buffer.data(psi1 + 14);
    const auto *psi1_49 = buffer.data(psi1 + 49);
    const auto *psi1_51 = buffer.data(psi1 + 51);
    const auto *psi1_52 = buffer.data(psi1 + 52);
    const auto *psi1_53 = buffer.data(psi1 + 53);

    const auto *ppg0_0 = buffer.data(ppg0 + 0);
    const auto *ppg0_1 = buffer.data(ppg0 + 1);
    const auto *ppg0_2 = buffer.data(ppg0 + 2);
    const auto *ppg0_3 = buffer.data(ppg0 + 3);
    const auto *ppg0_5 = buffer.data(ppg0 + 5);
    const auto *ppg0_10 = buffer.data(ppg0 + 10);
    const auto *ppg0_12 = buffer.data(ppg0 + 12);
    const auto *ppg0_13 = buffer.data(ppg0 + 13);
    const auto *ppg0_14 = buffer.data(ppg0 + 14);
    const auto *ppg0_35 = buffer.data(ppg0 + 35);
    const auto *ppg0_60 = buffer.data(ppg0 + 60);
    const auto *ppg0_61 = buffer.data(ppg0 + 61);
    const auto *ppg0_63 = buffer.data(ppg0 + 63);
    const auto *ppg0_65 = buffer.data(ppg0 + 65);
    const auto *ppg0_66 = buffer.data(ppg0 + 66);
    const auto *ppg0_68 = buffer.data(ppg0 + 68);
    const auto *ppg0_69 = buffer.data(ppg0 + 69);
    const auto *ppg0_70 = buffer.data(ppg0 + 70);
    const auto *ppg0_72 = buffer.data(ppg0 + 72);
    const auto *ppg0_73 = buffer.data(ppg0 + 73);
    const auto *ppg0_74 = buffer.data(ppg0 + 74);

    const auto *ppg1_0 = buffer.data(ppg1 + 0);
    const auto *ppg1_1 = buffer.data(ppg1 + 1);
    const auto *ppg1_2 = buffer.data(ppg1 + 2);
    const auto *ppg1_3 = buffer.data(ppg1 + 3);
    const auto *ppg1_5 = buffer.data(ppg1 + 5);
    const auto *ppg1_10 = buffer.data(ppg1 + 10);
    const auto *ppg1_12 = buffer.data(ppg1 + 12);
    const auto *ppg1_13 = buffer.data(ppg1 + 13);
    const auto *ppg1_14 = buffer.data(ppg1 + 14);
    const auto *ppg1_35 = buffer.data(ppg1 + 35);
    const auto *ppg1_60 = buffer.data(ppg1 + 60);
    const auto *ppg1_61 = buffer.data(ppg1 + 61);
    const auto *ppg1_63 = buffer.data(ppg1 + 63);
    const auto *ppg1_65 = buffer.data(ppg1 + 65);
    const auto *ppg1_66 = buffer.data(ppg1 + 66);
    const auto *ppg1_68 = buffer.data(ppg1 + 68);
    const auto *ppg1_69 = buffer.data(ppg1 + 69);
    const auto *ppg1_70 = buffer.data(ppg1 + 70);
    const auto *ppg1_72 = buffer.data(ppg1 + 72);
    const auto *ppg1_73 = buffer.data(ppg1 + 73);
    const auto *ppg1_74 = buffer.data(ppg1 + 74);

    const auto *pph_0 = buffer.data(pph + 0);
    const auto *pph_1 = buffer.data(pph + 1);
    const auto *pph_2 = buffer.data(pph + 2);
    const auto *pph_3 = buffer.data(pph + 3);
    const auto *pph_5 = buffer.data(pph + 5);
    const auto *pph_6 = buffer.data(pph + 6);
    const auto *pph_8 = buffer.data(pph + 8);
    const auto *pph_9 = buffer.data(pph + 9);
    const auto *pph_10 = buffer.data(pph + 10);
    const auto *pph_14 = buffer.data(pph + 14);
    const auto *pph_15 = buffer.data(pph + 15);
    const auto *pph_17 = buffer.data(pph + 17);
    const auto *pph_18 = buffer.data(pph + 18);
    const auto *pph_19 = buffer.data(pph + 19);
    const auto *pph_20 = buffer.data(pph + 20);
    const auto *pph_21 = buffer.data(pph + 21);
    const auto *pph_23 = buffer.data(pph + 23);
    const auto *pph_24 = buffer.data(pph + 24);
    const auto *pph_26 = buffer.data(pph + 26);
    const auto *pph_27 = buffer.data(pph + 27);
    const auto *pph_30 = buffer.data(pph + 30);
    const auto *pph_31 = buffer.data(pph + 31);
    const auto *pph_35 = buffer.data(pph + 35);
    const auto *pph_36 = buffer.data(pph + 36);
    const auto *pph_38 = buffer.data(pph + 38);
    const auto *pph_39 = buffer.data(pph + 39);
    const auto *pph_41 = buffer.data(pph + 41);
    const auto *pph_42 = buffer.data(pph + 42);
    const auto *pph_44 = buffer.data(pph + 44);
    const auto *pph_45 = buffer.data(pph + 45);
    const auto *pph_47 = buffer.data(pph + 47);
    const auto *pph_48 = buffer.data(pph + 48);
    const auto *pph_50 = buffer.data(pph + 50);
    const auto *pph_51 = buffer.data(pph + 51);
    const auto *pph_52 = buffer.data(pph + 52);
    const auto *pph_56 = buffer.data(pph + 56);
    const auto *pph_57 = buffer.data(pph + 57);
    const auto *pph_59 = buffer.data(pph + 59);
    const auto *pph_60 = buffer.data(pph + 60);
    const auto *pph_62 = buffer.data(pph + 62);
    const auto *pph_63 = buffer.data(pph + 63);
    const auto *pph_64 = buffer.data(pph + 64);
    const auto *pph_66 = buffer.data(pph + 66);
    const auto *pph_69 = buffer.data(pph + 69);
    const auto *pph_78 = buffer.data(pph + 78);
    const auto *pph_79 = buffer.data(pph + 79);
    const auto *pph_80 = buffer.data(pph + 80);
    const auto *pph_81 = buffer.data(pph + 81);
    const auto *pph_82 = buffer.data(pph + 82);
    const auto *pph_83 = buffer.data(pph + 83);
    const auto *pph_84 = buffer.data(pph + 84);
    const auto *pph_85 = buffer.data(pph + 85);
    const auto *pph_87 = buffer.data(pph + 87);
    const auto *pph_89 = buffer.data(pph + 89);
    const auto *pph_90 = buffer.data(pph + 90);
    const auto *pph_92 = buffer.data(pph + 92);
    const auto *pph_93 = buffer.data(pph + 93);
    const auto *pph_94 = buffer.data(pph + 94);
    const auto *pph_96 = buffer.data(pph + 96);
    const auto *pph_97 = buffer.data(pph + 97);
    const auto *pph_98 = buffer.data(pph + 98);
    const auto *pph_99 = buffer.data(pph + 99);
    const auto *pph_100 = buffer.data(pph + 100);
    const auto *pph_101 = buffer.data(pph + 101);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sph_0, psh_0, ppg0_0, \
                         ppg1_0, pph_0, pph_1, pph_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sph_0[k]
                 + f_0 * psh_0[k]
                 + f_1 * ppg0_0[k]
                 - f_2 * ppg1_0[k]
                 + f_3 * pc_x[k] * pph_0[k];

        t_1[k] = f_3 * pc_y[k] * pph_0[k];

        t_2[k] = f_3 * pc_z[k] * pph_0[k];

        t_3[k] = f_4 * ppg0_0[k]
                 - f_5 * ppg1_0[k]
                 + f_3 * pc_y[k] * pph_1[k];

        t_4[k] = f_3 * pc_y[k] * pph_2[k];

        t_5[k] = f_4 * ppg0_0[k]
                 - f_5 * ppg1_0[k]
                 + f_3 * pc_z[k] * pph_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, ppg0_1, ppg0_2, ppg0_3, ppg1_1, \
                         ppg1_2, ppg1_3, pph_3, pph_5, pph_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ppg0_1[k]
                 - f_7 * ppg1_1[k]
                 + f_3 * pc_y[k] * pph_3[k];

        t_7[k] = f_3 * pc_z[k] * pph_3[k];

        t_8[k] = f_3 * pc_y[k] * pph_5[k];

        t_9[k] = f_6 * ppg0_2[k]
                 - f_7 * ppg1_2[k]
                 + f_3 * pc_z[k] * pph_5[k];

        t_10[k] = f_8 * ppg0_3[k]
                  - f_9 * ppg1_3[k]
                  + f_3 * pc_y[k] * pph_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, sph_15, psh_15, \
                         ppg0_5, ppg1_5, pph_6, pph_8, pph_9, pph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * pph_6[k];

        t_12[k] = f_4 * ppg0_5[k]
                  - f_5 * ppg1_5[k]
                  + f_3 * pc_y[k] * pph_8[k];

        t_13[k] = f_3 * pc_y[k] * pph_9[k];

        t_14[k] = f_8 * ppg0_5[k]
                  - f_9 * ppg1_5[k]
                  + f_3 * pc_z[k] * pph_9[k];

        t_15[k] = f_0 * sph_15[k]
                  + f_0 * psh_15[k]
                  + f_3 * pc_x[k] * pph_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, sph_17, sph_18, psh_17, \
                         psh_18, pph_10, pph_14, pph_17, pph_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * pph_10[k];

        t_17[k] = f_0 * sph_17[k]
                  + f_0 * psh_17[k]
                  + f_3 * pc_x[k] * pph_17[k];

        t_18[k] = f_0 * sph_18[k]
                  + f_0 * psh_18[k]
                  + f_3 * pc_x[k] * pph_18[k];

        t_19[k] = f_3 * pc_y[k] * pph_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sph_20, psh_20, ppg0_10, \
                         ppg0_12, ppg1_10, ppg1_12, pph_15, pph_17, \
                         pph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * sph_20[k]
                  + f_0 * psh_20[k]
                  + f_3 * pc_x[k] * pph_20[k];

        t_21[k] = f_1 * ppg0_10[k]
                  - f_2 * ppg1_10[k]
                  + f_3 * pc_y[k] * pph_15[k];

        t_22[k] = f_3 * pc_z[k] * pph_15[k];

        t_23[k] = f_8 * ppg0_12[k]
                  - f_9 * ppg1_12[k]
                  + f_3 * pc_y[k] * pph_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, ppg0_13, ppg0_14, ppg1_13, \
                         ppg1_14, pph_18, pph_19, pph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * ppg0_13[k]
                  - f_7 * ppg1_13[k]
                  + f_3 * pc_y[k] * pph_18[k];

        t_25[k] = f_4 * ppg0_14[k]
                  - f_5 * ppg1_14[k]
                  + f_3 * pc_y[k] * pph_19[k];

        t_26[k] = f_3 * pc_y[k] * pph_20[k];

        t_27[k] = f_1 * ppg0_14[k]
                  - f_2 * ppg1_14[k]
                  + f_3 * pc_z[k] * pph_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_y, pc_x, pc_y, pc_z, spi0_31, \
                         sph_24, spi1_31, psi0_0, psh_0, psi1_0, \
                         pph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * psi0_0[k]
                  - f_10 * pc_y[k] * psi1_0[k];

        t_29[k] = f_0 * psh_0[k]
                  + f_3 * pc_y[k] * pph_21[k];

        t_30[k] = f_3 * pc_z[k] * pph_21[k];

        t_31[k] = pa_x[k] * spi0_31[k]
                  + f_11 * sph_24[k]
                  - f_10 * pc_x[k] * spi1_31[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pb_y, pc_x, pc_y, spi0_34, sph_27, spi1_34, \
                         psi0_5, psh_2, psi1_5, pph_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * psh_2[k]
                  + f_3 * pc_y[k] * pph_23[k];

        t_33[k] = pb_y[k] * psi0_5[k]
                  - f_10 * pc_y[k] * psi1_5[k];

        t_34[k] = pa_x[k] * spi0_34[k]
                  + f_12 * sph_27[k]
                  - f_10 * pc_x[k] * spi1_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_y, pc_z, psi0_9, psh_5, psi1_9, pph_24, \
                         pph_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * pc_z[k] * pph_24[k];

        t_36[k] = f_0 * psh_5[k]
                  + f_3 * pc_y[k] * pph_26[k];

        t_37[k] = pb_y[k] * psi0_9[k]
                  - f_10 * pc_y[k] * psi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pc_x, pc_z, spi0_38, spi0_40, sph_31, sph_33, \
                         spi1_38, spi1_40, pph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_x[k] * spi0_38[k]
                  + f_13 * sph_31[k]
                  - f_10 * pc_x[k] * spi1_38[k];

        t_39[k] = f_3 * pc_z[k] * pph_27[k];

        t_40[k] = pa_x[k] * spi0_40[k]
                  + f_13 * sph_33[k]
                  - f_10 * pc_x[k] * spi1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pb_y, pc_x, pc_y, pc_z, sph_36, psi0_14, \
                         psh_9, psi1_14, pph_30, pph_31, pph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_0 * psh_9[k]
                  + f_3 * pc_y[k] * pph_30[k];

        t_42[k] = pb_y[k] * psi0_14[k]
                  - f_10 * pc_y[k] * psi1_14[k];

        t_43[k] = f_0 * sph_36[k]
                  + f_3 * pc_x[k] * pph_36[k];

        t_44[k] = f_3 * pc_z[k] * pph_31[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, sph_38, sph_39, sph_41, psh_14, \
                         pph_35, pph_38, pph_39, pph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * sph_38[k]
                  + f_3 * pc_x[k] * pph_38[k];

        t_46[k] = f_0 * sph_39[k]
                  + f_3 * pc_x[k] * pph_39[k];

        t_47[k] = f_0 * psh_14[k]
                  + f_3 * pc_y[k] * pph_35[k];

        t_48[k] = f_0 * sph_41[k]
                  + f_3 * pc_x[k] * pph_41[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pc_x, pc_z, spi0_49, spi0_51, spi0_52, \
                         spi1_49, spi1_51, spi1_52, pph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_x[k] * spi0_49[k]
                  - f_10 * pc_x[k] * spi1_49[k];

        t_50[k] = f_3 * pc_z[k] * pph_36[k];

        t_51[k] = pa_x[k] * spi0_51[k]
                  - f_10 * pc_x[k] * spi1_51[k];

        t_52[k] = pa_x[k] * spi0_52[k]
                  - f_10 * pc_x[k] * spi1_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pc_x, pc_y, spi0_53, spi0_55, spi1_53, \
                         spi1_55, psh_20, pph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_x[k] * spi0_53[k]
                  - f_10 * pc_x[k] * spi1_53[k];

        t_54[k] = f_0 * psh_20[k]
                  + f_3 * pc_y[k] * pph_41[k];

        t_55[k] = pa_x[k] * spi0_55[k]
                  - f_10 * pc_x[k] * spi1_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, psi0_0, psi0_3, \
                         psh_0, psi1_0, psi1_3, pph_42, pph_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * psi0_0[k]
                  - f_10 * pc_z[k] * psi1_0[k];

        t_57[k] = f_3 * pc_y[k] * pph_42[k];

        t_58[k] = f_0 * psh_0[k]
                  + f_3 * pc_z[k] * pph_42[k];

        t_59[k] = pb_z[k] * psi0_3[k]
                  - f_10 * pc_z[k] * psi1_3[k];

        t_60[k] = f_3 * pc_y[k] * pph_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_x, pb_z, pc_x, pc_z, spi0_61, sph_47, spi1_61, \
                         psi0_6, psh_3, psi1_6, pph_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_x[k] * spi0_61[k]
                  + f_11 * sph_47[k]
                  - f_10 * pc_x[k] * spi1_61[k];

        t_62[k] = pb_z[k] * psi0_6[k]
                  - f_10 * pc_z[k] * psi1_6[k];

        t_63[k] = f_0 * psh_3[k]
                  + f_3 * pc_z[k] * pph_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_x, pb_z, pc_x, pc_y, pc_z, spi0_65, sph_51, \
                         spi1_65, psi0_10, psi1_10, pph_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * pc_y[k] * pph_47[k];

        t_65[k] = pa_x[k] * spi0_65[k]
                  + f_12 * sph_51[k]
                  - f_10 * pc_x[k] * spi1_65[k];

        t_66[k] = pb_z[k] * psi0_10[k]
                  - f_10 * pc_z[k] * psi1_10[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_y, pc_z, psh_6, ppg0_35, ppg1_35, pph_48, \
                         pph_50, pph_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * psh_6[k]
                  + f_3 * pc_z[k] * pph_48[k];

        t_68[k] = f_4 * ppg0_35[k]
                  - f_5 * ppg1_35[k]
                  + f_3 * pc_y[k] * pph_50[k];

        t_69[k] = f_3 * pc_y[k] * pph_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pc_x, pc_z, spi0_70, sph_56, sph_57, \
                         sph_59, spi1_70, psh_10, pph_52, pph_57, \
                         pph_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * spi0_70[k]
                  + f_13 * sph_56[k]
                  - f_10 * pc_x[k] * spi1_70[k];

        t_71[k] = f_0 * sph_57[k]
                  + f_3 * pc_x[k] * pph_57[k];

        t_72[k] = f_0 * psh_10[k]
                  + f_3 * pc_z[k] * pph_52[k];

        t_73[k] = f_0 * sph_59[k]
                  + f_3 * pc_x[k] * pph_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_y, spi0_77, sph_60, sph_62, \
                         spi1_77, pph_56, pph_60, pph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * sph_60[k]
                  + f_3 * pc_x[k] * pph_60[k];

        t_75[k] = f_3 * pc_y[k] * pph_56[k];

        t_76[k] = f_0 * sph_62[k]
                  + f_3 * pc_x[k] * pph_62[k];

        t_77[k] = pa_x[k] * spi0_77[k]
                  - f_10 * pc_x[k] * spi1_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pc_x, pc_z, spi0_79, spi0_80, spi0_81, \
                         spi1_79, spi1_80, spi1_81, psh_15, pph_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * psh_15[k]
                  + f_3 * pc_z[k] * pph_57[k];

        t_79[k] = pa_x[k] * spi0_79[k]
                  - f_10 * pc_x[k] * spi1_79[k];

        t_80[k] = pa_x[k] * spi0_80[k]
                  - f_10 * pc_x[k] * spi1_80[k];

        t_81[k] = pa_x[k] * spi0_81[k]
                  - f_10 * pc_x[k] * spi1_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, pc_x, pc_y, spi0_0, spi0_1, \
                         spi0_83, sph_0, spi1_0, spi1_1, spi1_83, \
                         pph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_y[k] * pph_62[k];

        t_83[k] = pa_x[k] * spi0_83[k]
                  - f_10 * pc_x[k] * spi1_83[k];

        t_84[k] = pa_y[k] * spi0_0[k]
                  - f_10 * pc_y[k] * spi1_0[k];

        t_85[k] = pa_y[k] * spi0_1[k]
                  + f_0 * sph_0[k]
                  - f_10 * pc_y[k] * spi1_1[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pc_y, pc_z, spi0_3, spi0_5, sph_1, \
                         spi1_3, spi1_5, pph_63, pph_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * pc_z[k] * pph_63[k];

        t_87[k] = pa_y[k] * spi0_3[k]
                  + f_13 * sph_1[k]
                  - f_10 * pc_y[k] * spi1_3[k];

        t_88[k] = f_3 * pc_z[k] * pph_64[k];

        t_89[k] = pa_y[k] * spi0_5[k]
                  - f_10 * pc_y[k] * spi1_5[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_y, pc_y, pc_z, spi0_6, spi0_8, spi0_9, \
                         sph_3, sph_5, spi1_6, spi1_8, spi1_9, pph_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_y[k] * spi0_6[k]
                  + f_12 * sph_3[k]
                  - f_10 * pc_y[k] * spi1_6[k];

        t_91[k] = f_3 * pc_z[k] * pph_66[k];

        t_92[k] = pa_y[k] * spi0_8[k]
                  + f_0 * sph_5[k]
                  - f_10 * pc_y[k] * spi1_8[k];

        t_93[k] = pa_y[k] * spi0_9[k]
                  - f_10 * pc_y[k] * spi1_9[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pc_y, pc_z, spi0_10, spi0_12, sph_6, sph_8, \
                         spi1_10, spi1_12, pph_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pa_y[k] * spi0_10[k]
                  + f_11 * sph_6[k]
                  - f_10 * pc_y[k] * spi1_10[k];

        t_95[k] = f_3 * pc_z[k] * pph_69[k];

        t_96[k] = pa_y[k] * spi0_12[k]
                  + f_13 * sph_8[k]
                  - f_10 * pc_y[k] * spi1_12[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pc_x, pc_y, spi0_13, spi0_14, sph_9, \
                         spi1_13, spi1_14, psh_36, psh_37, pph_78, \
                         pph_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = pa_y[k] * spi0_13[k]
                  + f_0 * sph_9[k]
                  - f_10 * pc_y[k] * spi1_13[k];

        t_98[k] = pa_y[k] * spi0_14[k]
                  - f_10 * pc_y[k] * spi1_14[k];

        t_99[k] = f_0 * psh_36[k]
                  + f_3 * pc_x[k] * pph_78[k];

        t_100[k] = f_0 * psh_37[k]
                   + f_3 * pc_x[k] * pph_79[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, psh_38, psh_39, psh_40, psh_41, \
                         pph_80, pph_81, pph_82, pph_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * psh_38[k]
                   + f_3 * pc_x[k] * pph_80[k];

        t_102[k] = f_0 * psh_39[k]
                   + f_3 * pc_x[k] * pph_81[k];

        t_103[k] = f_0 * psh_40[k]
                   + f_3 * pc_x[k] * pph_82[k];

        t_104[k] = f_0 * psh_41[k]
                   + f_3 * pc_x[k] * pph_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, pc_x, pc_z, psi0_49, psi0_51, \
                         psi0_52, psi1_49, psi1_51, psi1_52, pph_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_x[k] * psi0_49[k]
                   - f_10 * pc_x[k] * psi1_49[k];

        t_106[k] = f_3 * pc_z[k] * pph_78[k];

        t_107[k] = pb_x[k] * psi0_51[k]
                   - f_10 * pc_x[k] * psi1_51[k];

        t_108[k] = pb_x[k] * psi0_52[k]
                   - f_10 * pc_x[k] * psi1_52[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, pa_y, pb_x, pc_x, pc_y, spi0_27, sph_20, \
                         spi1_27, psi0_53, psi1_53, pph_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_x[k] * psi0_53[k]
                   - f_10 * pc_x[k] * psi1_53[k];

        t_110[k] = f_0 * sph_20[k]
                   + f_3 * pc_y[k] * pph_83[k];

        t_111[k] = pa_y[k] * spi0_27[k]
                   - f_10 * pc_y[k] * spi1_27[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, pc_x, pc_z, ppg0_60, ppg0_61, \
                         ppg0_63, ppg1_60, ppg1_61, ppg1_63, pph_84, pph_85, \
                         pph_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * ppg0_60[k]
                   - f_2 * ppg1_60[k]
                   + f_3 * pc_x[k] * pph_84[k];

        t_113[k] = f_14 * ppg0_61[k]
                   - f_15 * ppg1_61[k]
                   + f_3 * pc_x[k] * pph_85[k];

        t_114[k] = f_3 * pc_z[k] * pph_84[k];

        t_115[k] = f_8 * ppg0_63[k]
                   - f_9 * ppg1_63[k]
                   + f_3 * pc_x[k] * pph_87[k];

        t_116[k] = f_3 * pc_z[k] * pph_85[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pc_x, pc_z, ppg0_65, ppg0_66, ppg0_68, \
                         ppg1_65, ppg1_66, ppg1_68, pph_87, pph_89, pph_90, \
                         pph_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * ppg0_65[k]
                   - f_9 * ppg1_65[k]
                   + f_3 * pc_x[k] * pph_89[k];

        t_118[k] = f_6 * ppg0_66[k]
                   - f_7 * ppg1_66[k]
                   + f_3 * pc_x[k] * pph_90[k];

        t_119[k] = f_3 * pc_z[k] * pph_87[k];

        t_120[k] = f_6 * ppg0_68[k]
                   - f_7 * ppg1_68[k]
                   + f_3 * pc_x[k] * pph_92[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_z, ppg0_69, ppg0_70, ppg0_72, \
                         ppg1_69, ppg1_70, ppg1_72, pph_90, pph_93, pph_94, \
                         pph_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_6 * ppg0_69[k]
                   - f_7 * ppg1_69[k]
                   + f_3 * pc_x[k] * pph_93[k];

        t_122[k] = f_4 * ppg0_70[k]
                   - f_5 * ppg1_70[k]
                   + f_3 * pc_x[k] * pph_94[k];

        t_123[k] = f_3 * pc_z[k] * pph_90[k];

        t_124[k] = f_4 * ppg0_72[k]
                   - f_5 * ppg1_72[k]
                   + f_3 * pc_x[k] * pph_96[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pc_x, ppg0_73, ppg0_74, ppg1_73, \
                         ppg1_74, pph_97, pph_98, pph_99, pph_100, \
                         pph_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_4 * ppg0_73[k]
                   - f_5 * ppg1_73[k]
                   + f_3 * pc_x[k] * pph_97[k];

        t_126[k] = f_4 * ppg0_74[k]
                   - f_5 * ppg1_74[k]
                   + f_3 * pc_x[k] * pph_98[k];

        t_127[k] = f_3 * pc_x[k] * pph_99[k];

        t_128[k] = f_3 * pc_x[k] * pph_100[k];

        t_129[k] = f_3 * pc_x[k] * pph_101[k];
    }
}

static auto
compute_prim_ppi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t spi0, const size_t sph,
                                                          const size_t spi1, const size_t psi0,
                                                          const size_t psh, const size_t psi1,
                                                          const size_t ppg0, const size_t ppg1,
                                                          const size_t pph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
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
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 1.0 / q;
    const auto f_14 = 2.0 / gamma;
    const auto f_15 = 2.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spi0_0 = buffer.data(spi0 + 0);
    const auto *spi0_2 = buffer.data(spi0 + 2);
    const auto *spi0_3 = buffer.data(spi0 + 3);
    const auto *spi0_5 = buffer.data(spi0 + 5);
    const auto *spi0_6 = buffer.data(spi0 + 6);
    const auto *spi0_7 = buffer.data(spi0 + 7);
    const auto *spi0_9 = buffer.data(spi0 + 9);
    const auto *spi0_10 = buffer.data(spi0 + 10);
    const auto *spi0_11 = buffer.data(spi0 + 11);
    const auto *spi0_12 = buffer.data(spi0 + 12);
    const auto *spi0_14 = buffer.data(spi0 + 14);
    const auto *spi0_21 = buffer.data(spi0 + 21);
    const auto *spi0_28 = buffer.data(spi0 + 28);
    const auto *spi0_31 = buffer.data(spi0 + 31);
    const auto *spi0_34 = buffer.data(spi0 + 34);
    const auto *spi0_38 = buffer.data(spi0 + 38);
    const auto *spi0_49 = buffer.data(spi0 + 49);
    const auto *spi0_50 = buffer.data(spi0 + 50);
    const auto *spi0_51 = buffer.data(spi0 + 51);
    const auto *spi0_52 = buffer.data(spi0 + 52);
    const auto *spi0_53 = buffer.data(spi0 + 53);
    const auto *spi0_56 = buffer.data(spi0 + 56);
    const auto *spi0_61 = buffer.data(spi0 + 61);
    const auto *spi0_65 = buffer.data(spi0 + 65);
    const auto *spi0_70 = buffer.data(spi0 + 70);
    const auto *spi0_79 = buffer.data(spi0 + 79);
    const auto *spi0_80 = buffer.data(spi0 + 80);
    const auto *spi0_81 = buffer.data(spi0 + 81);
    const auto *spi0_83 = buffer.data(spi0 + 83);

    const auto *sph_0 = buffer.data(sph + 0);
    const auto *sph_2 = buffer.data(sph + 2);
    const auto *sph_3 = buffer.data(sph + 3);
    const auto *sph_5 = buffer.data(sph + 5);
    const auto *sph_6 = buffer.data(sph + 6);
    const auto *sph_7 = buffer.data(sph + 7);
    const auto *sph_9 = buffer.data(sph + 9);
    const auto *sph_36 = buffer.data(sph + 36);
    const auto *sph_37 = buffer.data(sph + 37);
    const auto *sph_38 = buffer.data(sph + 38);
    const auto *sph_39 = buffer.data(sph + 39);
    const auto *sph_41 = buffer.data(sph + 41);
    const auto *sph_59 = buffer.data(sph + 59);
    const auto *sph_60 = buffer.data(sph + 60);
    const auto *sph_61 = buffer.data(sph + 61);
    const auto *sph_62 = buffer.data(sph + 62);

    const auto *spi1_0 = buffer.data(spi1 + 0);
    const auto *spi1_2 = buffer.data(spi1 + 2);
    const auto *spi1_3 = buffer.data(spi1 + 3);
    const auto *spi1_5 = buffer.data(spi1 + 5);
    const auto *spi1_6 = buffer.data(spi1 + 6);
    const auto *spi1_7 = buffer.data(spi1 + 7);
    const auto *spi1_9 = buffer.data(spi1 + 9);
    const auto *spi1_10 = buffer.data(spi1 + 10);
    const auto *spi1_11 = buffer.data(spi1 + 11);
    const auto *spi1_12 = buffer.data(spi1 + 12);
    const auto *spi1_14 = buffer.data(spi1 + 14);
    const auto *spi1_21 = buffer.data(spi1 + 21);
    const auto *spi1_28 = buffer.data(spi1 + 28);
    const auto *spi1_31 = buffer.data(spi1 + 31);
    const auto *spi1_34 = buffer.data(spi1 + 34);
    const auto *spi1_38 = buffer.data(spi1 + 38);
    const auto *spi1_49 = buffer.data(spi1 + 49);
    const auto *spi1_50 = buffer.data(spi1 + 50);
    const auto *spi1_51 = buffer.data(spi1 + 51);
    const auto *spi1_52 = buffer.data(spi1 + 52);
    const auto *spi1_53 = buffer.data(spi1 + 53);
    const auto *spi1_56 = buffer.data(spi1 + 56);
    const auto *spi1_61 = buffer.data(spi1 + 61);
    const auto *spi1_65 = buffer.data(spi1 + 65);
    const auto *spi1_70 = buffer.data(spi1 + 70);
    const auto *spi1_79 = buffer.data(spi1 + 79);
    const auto *spi1_80 = buffer.data(spi1 + 80);
    const auto *spi1_81 = buffer.data(spi1 + 81);
    const auto *spi1_83 = buffer.data(spi1 + 83);

    const auto *psi0_29 = buffer.data(psi0 + 29);
    const auto *psi0_31 = buffer.data(psi0 + 31);
    const auto *psi0_34 = buffer.data(psi0 + 34);
    const auto *psi0_38 = buffer.data(psi0 + 38);
    const auto *psi0_49 = buffer.data(psi0 + 49);
    const auto *psi0_58 = buffer.data(psi0 + 58);
    const auto *psi0_61 = buffer.data(psi0 + 61);
    const auto *psi0_65 = buffer.data(psi0 + 65);
    const auto *psi0_70 = buffer.data(psi0 + 70);
    const auto *psi0_78 = buffer.data(psi0 + 78);
    const auto *psi0_79 = buffer.data(psi0 + 79);
    const auto *psi0_80 = buffer.data(psi0 + 80);
    const auto *psi0_81 = buffer.data(psi0 + 81);
    const auto *psi0_83 = buffer.data(psi0 + 83);

    const auto *psh_21 = buffer.data(psh + 21);
    const auto *psh_22 = buffer.data(psh + 22);
    const auto *psh_24 = buffer.data(psh + 24);
    const auto *psh_27 = buffer.data(psh + 27);
    const auto *psh_36 = buffer.data(psh + 36);
    const auto *psh_41 = buffer.data(psh + 41);
    const auto *psh_42 = buffer.data(psh + 42);
    const auto *psh_44 = buffer.data(psh + 44);
    const auto *psh_47 = buffer.data(psh + 47);
    const auto *psh_51 = buffer.data(psh + 51);
    const auto *psh_57 = buffer.data(psh + 57);
    const auto *psh_58 = buffer.data(psh + 58);
    const auto *psh_59 = buffer.data(psh + 59);
    const auto *psh_60 = buffer.data(psh + 60);
    const auto *psh_61 = buffer.data(psh + 61);
    const auto *psh_62 = buffer.data(psh + 62);

    const auto *psi1_29 = buffer.data(psi1 + 29);
    const auto *psi1_31 = buffer.data(psi1 + 31);
    const auto *psi1_34 = buffer.data(psi1 + 34);
    const auto *psi1_38 = buffer.data(psi1 + 38);
    const auto *psi1_49 = buffer.data(psi1 + 49);
    const auto *psi1_58 = buffer.data(psi1 + 58);
    const auto *psi1_61 = buffer.data(psi1 + 61);
    const auto *psi1_65 = buffer.data(psi1 + 65);
    const auto *psi1_70 = buffer.data(psi1 + 70);
    const auto *psi1_78 = buffer.data(psi1 + 78);
    const auto *psi1_79 = buffer.data(psi1 + 79);
    const auto *psi1_80 = buffer.data(psi1 + 80);
    const auto *psi1_81 = buffer.data(psi1 + 81);
    const auto *psi1_83 = buffer.data(psi1 + 83);

    const auto *ppg0_70 = buffer.data(ppg0 + 70);
    const auto *ppg0_71 = buffer.data(ppg0 + 71);
    const auto *ppg0_72 = buffer.data(ppg0 + 72);
    const auto *ppg0_74 = buffer.data(ppg0 + 74);
    const auto *ppg0_83 = buffer.data(ppg0 + 83);
    const auto *ppg0_87 = buffer.data(ppg0 + 87);
    const auto *ppg0_88 = buffer.data(ppg0 + 88);
    const auto *ppg0_112 = buffer.data(ppg0 + 112);
    const auto *ppg0_116 = buffer.data(ppg0 + 116);
    const auto *ppg0_117 = buffer.data(ppg0 + 117);
    const auto *ppg0_120 = buffer.data(ppg0 + 120);
    const auto *ppg0_122 = buffer.data(ppg0 + 122);
    const auto *ppg0_123 = buffer.data(ppg0 + 123);
    const auto *ppg0_125 = buffer.data(ppg0 + 125);
    const auto *ppg0_126 = buffer.data(ppg0 + 126);
    const auto *ppg0_127 = buffer.data(ppg0 + 127);
    const auto *ppg0_129 = buffer.data(ppg0 + 129);
    const auto *ppg0_130 = buffer.data(ppg0 + 130);
    const auto *ppg0_131 = buffer.data(ppg0 + 131);
    const auto *ppg0_132 = buffer.data(ppg0 + 132);
    const auto *ppg0_133 = buffer.data(ppg0 + 133);
    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppg1_70 = buffer.data(ppg1 + 70);
    const auto *ppg1_71 = buffer.data(ppg1 + 71);
    const auto *ppg1_72 = buffer.data(ppg1 + 72);
    const auto *ppg1_74 = buffer.data(ppg1 + 74);
    const auto *ppg1_83 = buffer.data(ppg1 + 83);
    const auto *ppg1_87 = buffer.data(ppg1 + 87);
    const auto *ppg1_88 = buffer.data(ppg1 + 88);
    const auto *ppg1_112 = buffer.data(ppg1 + 112);
    const auto *ppg1_116 = buffer.data(ppg1 + 116);
    const auto *ppg1_117 = buffer.data(ppg1 + 117);
    const auto *ppg1_120 = buffer.data(ppg1 + 120);
    const auto *ppg1_122 = buffer.data(ppg1 + 122);
    const auto *ppg1_123 = buffer.data(ppg1 + 123);
    const auto *ppg1_125 = buffer.data(ppg1 + 125);
    const auto *ppg1_126 = buffer.data(ppg1 + 126);
    const auto *ppg1_127 = buffer.data(ppg1 + 127);
    const auto *ppg1_129 = buffer.data(ppg1 + 129);
    const auto *ppg1_130 = buffer.data(ppg1 + 130);
    const auto *ppg1_131 = buffer.data(ppg1 + 131);
    const auto *ppg1_132 = buffer.data(ppg1 + 132);
    const auto *ppg1_133 = buffer.data(ppg1 + 133);
    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *pph_99 = buffer.data(pph + 99);
    const auto *pph_100 = buffer.data(pph + 100);
    const auto *pph_101 = buffer.data(pph + 101);
    const auto *pph_102 = buffer.data(pph + 102);
    const auto *pph_103 = buffer.data(pph + 103);
    const auto *pph_104 = buffer.data(pph + 104);
    const auto *pph_105 = buffer.data(pph + 105);
    const auto *pph_106 = buffer.data(pph + 106);
    const auto *pph_108 = buffer.data(pph + 108);
    const auto *pph_111 = buffer.data(pph + 111);
    const auto *pph_113 = buffer.data(pph + 113);
    const auto *pph_117 = buffer.data(pph + 117);
    const auto *pph_118 = buffer.data(pph + 118);
    const auto *pph_120 = buffer.data(pph + 120);
    const auto *pph_121 = buffer.data(pph + 121);
    const auto *pph_122 = buffer.data(pph + 122);
    const auto *pph_123 = buffer.data(pph + 123);
    const auto *pph_124 = buffer.data(pph + 124);
    const auto *pph_125 = buffer.data(pph + 125);
    const auto *pph_126 = buffer.data(pph + 126);
    const auto *pph_128 = buffer.data(pph + 128);
    const auto *pph_131 = buffer.data(pph + 131);
    const auto *pph_135 = buffer.data(pph + 135);
    const auto *pph_141 = buffer.data(pph + 141);
    const auto *pph_142 = buffer.data(pph + 142);
    const auto *pph_143 = buffer.data(pph + 143);
    const auto *pph_144 = buffer.data(pph + 144);
    const auto *pph_145 = buffer.data(pph + 145);
    const auto *pph_146 = buffer.data(pph + 146);
    const auto *pph_147 = buffer.data(pph + 147);
    const auto *pph_149 = buffer.data(pph + 149);
    const auto *pph_152 = buffer.data(pph + 152);
    const auto *pph_154 = buffer.data(pph + 154);
    const auto *pph_156 = buffer.data(pph + 156);
    const auto *pph_158 = buffer.data(pph + 158);
    const auto *pph_159 = buffer.data(pph + 159);
    const auto *pph_162 = buffer.data(pph + 162);
    const auto *pph_163 = buffer.data(pph + 163);
    const auto *pph_164 = buffer.data(pph + 164);
    const auto *pph_165 = buffer.data(pph + 165);
    const auto *pph_166 = buffer.data(pph + 166);
    const auto *pph_167 = buffer.data(pph + 167);
    const auto *pph_168 = buffer.data(pph + 168);
    const auto *pph_170 = buffer.data(pph + 170);
    const auto *pph_171 = buffer.data(pph + 171);
    const auto *pph_173 = buffer.data(pph + 173);
    const auto *pph_174 = buffer.data(pph + 174);
    const auto *pph_175 = buffer.data(pph + 175);
    const auto *pph_177 = buffer.data(pph + 177);
    const auto *pph_178 = buffer.data(pph + 178);
    const auto *pph_179 = buffer.data(pph + 179);
    const auto *pph_180 = buffer.data(pph + 180);
    const auto *pph_182 = buffer.data(pph + 182);
    const auto *pph_183 = buffer.data(pph + 183);
    const auto *pph_184 = buffer.data(pph + 184);
    const auto *pph_185 = buffer.data(pph + 185);
    const auto *pph_186 = buffer.data(pph + 186);
    const auto *pph_187 = buffer.data(pph + 187);
    const auto *pph_188 = buffer.data(pph + 188);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pc_x, pc_y, pc_z, sph_36, psh_36, \
                         ppg0_70, ppg1_70, pph_99, pph_102, pph_103, \
                         pph_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_x[k] * pph_102[k];

        t_131[k] = f_3 * pc_x[k] * pph_103[k];

        t_132[k] = f_3 * pc_x[k] * pph_104[k];

        t_133[k] = f_0 * sph_36[k]
                   + f_0 * psh_36[k]
                   + f_1 * ppg0_70[k]
                   - f_2 * ppg1_70[k]
                   + f_3 * pc_y[k] * pph_99[k];

        t_134[k] = f_3 * pc_z[k] * pph_99[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pc_z, ppg0_70, ppg0_71, ppg0_72, ppg1_70, \
                         ppg1_71, ppg1_72, pph_100, pph_101, pph_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_4 * ppg0_70[k]
                   - f_5 * ppg1_70[k]
                   + f_3 * pc_z[k] * pph_100[k];

        t_136[k] = f_6 * ppg0_71[k]
                   - f_7 * ppg1_71[k]
                   + f_3 * pc_z[k] * pph_101[k];

        t_137[k] = f_8 * ppg0_72[k]
                   - f_9 * ppg1_72[k]
                   + f_3 * pc_z[k] * pph_102[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pa_y, pc_y, pc_z, spi0_56, sph_41, spi1_56, \
                         psh_41, ppg0_74, ppg1_74, pph_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * sph_41[k]
                   + f_0 * psh_41[k]
                   + f_3 * pc_y[k] * pph_104[k];

        t_139[k] = f_1 * ppg0_74[k]
                   - f_2 * ppg1_74[k]
                   + f_3 * pc_z[k] * pph_104[k];

        t_140[k] = pa_y[k] * spi0_56[k]
                   - f_10 * pc_y[k] * spi1_56[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_z, pc_z, psi0_29, psi0_31, psh_21, \
                         psh_22, psi1_29, psi1_31, pph_105, pph_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_z[k] * psi0_29[k]
                   - f_10 * pc_z[k] * psi1_29[k];

        t_142[k] = f_0 * psh_21[k]
                   + f_3 * pc_z[k] * pph_105[k];

        t_143[k] = pb_z[k] * psi0_31[k]
                   - f_10 * pc_z[k] * psi1_31[k];

        t_144[k] = f_0 * psh_22[k]
                   + f_3 * pc_z[k] * pph_106[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pa_y, pb_z, pc_y, pc_z, spi0_61, spi1_61, \
                         psi0_34, psh_24, psi1_34, pph_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_y[k] * spi0_61[k]
                   - f_10 * pc_y[k] * spi1_61[k];

        t_146[k] = pb_z[k] * psi0_34[k]
                   - f_10 * pc_z[k] * psi1_34[k];

        t_147[k] = f_0 * psh_24[k]
                   + f_3 * pc_z[k] * pph_108[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_y, pb_z, pc_x, pc_y, pc_z, spi0_65, spi1_65, \
                         psi0_38, psi1_38, ppg0_83, ppg1_83, pph_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_6 * ppg0_83[k]
                   - f_7 * ppg1_83[k]
                   + f_3 * pc_x[k] * pph_113[k];

        t_149[k] = pa_y[k] * spi0_65[k]
                   - f_10 * pc_y[k] * spi1_65[k];

        t_150[k] = pb_z[k] * psi0_38[k]
                   - f_10 * pc_z[k] * psi1_38[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pc_x, pc_z, psh_27, ppg0_87, ppg0_88, ppg1_87, \
                         ppg1_88, pph_111, pph_117, pph_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_0 * psh_27[k]
                   + f_3 * pc_z[k] * pph_111[k];

        t_152[k] = f_4 * ppg0_87[k]
                   - f_5 * ppg1_87[k]
                   + f_3 * pc_x[k] * pph_117[k];

        t_153[k] = f_4 * ppg0_88[k]
                   - f_5 * ppg1_88[k]
                   + f_3 * pc_x[k] * pph_118[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, t_159, pa_y, pc_x, pc_y, spi0_70, \
                         spi1_70, pph_120, pph_121, pph_122, pph_123, \
                         pph_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_y[k] * spi0_70[k]
                   - f_10 * pc_y[k] * spi1_70[k];

        t_155[k] = f_3 * pc_x[k] * pph_120[k];

        t_156[k] = f_3 * pc_x[k] * pph_121[k];

        t_157[k] = f_3 * pc_x[k] * pph_122[k];

        t_158[k] = f_3 * pc_x[k] * pph_123[k];

        t_159[k] = f_3 * pc_x[k] * pph_124[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pb_z, pc_x, pc_z, psi0_49, psh_36, psi1_49, \
                         pph_120, pph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_3 * pc_x[k] * pph_125[k];

        t_161[k] = pb_z[k] * psi0_49[k]
                   - f_10 * pc_z[k] * psi1_49[k];

        t_162[k] = f_0 * psh_36[k]
                   + f_3 * pc_z[k] * pph_120[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_y, pc_y, spi0_79, spi0_80, spi0_81, sph_59, \
                         sph_60, sph_61, spi1_79, spi1_80, spi1_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_y[k] * spi0_79[k]
                   + f_11 * sph_59[k]
                   - f_10 * pc_y[k] * spi1_79[k];

        t_164[k] = pa_y[k] * spi0_80[k]
                   + f_12 * sph_60[k]
                   - f_10 * pc_y[k] * spi1_80[k];

        t_165[k] = pa_y[k] * spi0_81[k]
                   + f_13 * sph_61[k]
                   - f_10 * pc_y[k] * spi1_81[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_y, pa_z, pc_y, pc_z, spi0_0, spi0_83, \
                         sph_62, spi1_0, spi1_83, pph_125, pph_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_0 * sph_62[k]
                   + f_3 * pc_y[k] * pph_125[k];

        t_167[k] = pa_y[k] * spi0_83[k]
                   - f_10 * pc_y[k] * spi1_83[k];

        t_168[k] = pa_z[k] * spi0_0[k]
                   - f_10 * pc_z[k] * spi1_0[k];

        t_169[k] = f_3 * pc_y[k] * pph_126[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_z, pc_y, pc_z, spi0_2, spi0_3, spi0_5, \
                         sph_0, sph_2, spi1_2, spi1_3, spi1_5, \
                         pph_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pa_z[k] * spi0_2[k]
                   + f_0 * sph_0[k]
                   - f_10 * pc_z[k] * spi1_2[k];

        t_171[k] = pa_z[k] * spi0_3[k]
                   - f_10 * pc_z[k] * spi1_3[k];

        t_172[k] = f_3 * pc_y[k] * pph_128[k];

        t_173[k] = pa_z[k] * spi0_5[k]
                   + f_13 * sph_2[k]
                   - f_10 * pc_z[k] * spi1_5[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_z, pc_y, pc_z, spi0_6, spi0_7, spi0_9, \
                         sph_3, sph_5, spi1_6, spi1_7, spi1_9, \
                         pph_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_z[k] * spi0_6[k]
                   - f_10 * pc_z[k] * spi1_6[k];

        t_175[k] = pa_z[k] * spi0_7[k]
                   + f_0 * sph_3[k]
                   - f_10 * pc_z[k] * spi1_7[k];

        t_176[k] = f_3 * pc_y[k] * pph_131[k];

        t_177[k] = pa_z[k] * spi0_9[k]
                   + f_12 * sph_5[k]
                   - f_10 * pc_z[k] * spi1_9[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_z, pc_y, pc_z, spi0_10, spi0_11, \
                         spi0_12, sph_6, sph_7, spi1_10, spi1_11, spi1_12, \
                         pph_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = pa_z[k] * spi0_10[k]
                   - f_10 * pc_z[k] * spi1_10[k];

        t_179[k] = pa_z[k] * spi0_11[k]
                   + f_0 * sph_6[k]
                   - f_10 * pc_z[k] * spi1_11[k];

        t_180[k] = pa_z[k] * spi0_12[k]
                   + f_13 * sph_7[k]
                   - f_10 * pc_z[k] * spi1_12[k];

        t_181[k] = f_3 * pc_y[k] * pph_135[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pa_z, pc_x, pc_z, spi0_14, sph_9, \
                         spi1_14, psh_57, psh_58, psh_59, pph_141, pph_142, \
                         pph_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = pa_z[k] * spi0_14[k]
                   + f_11 * sph_9[k]
                   - f_10 * pc_z[k] * spi1_14[k];

        t_183[k] = f_0 * psh_57[k]
                   + f_3 * pc_x[k] * pph_141[k];

        t_184[k] = f_0 * psh_58[k]
                   + f_3 * pc_x[k] * pph_142[k];

        t_185[k] = f_0 * psh_59[k]
                   + f_3 * pc_x[k] * pph_143[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pa_z, pc_x, pc_z, spi0_21, spi1_21, \
                         psh_60, psh_61, psh_62, pph_144, pph_145, \
                         pph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_0 * psh_60[k]
                   + f_3 * pc_x[k] * pph_144[k];

        t_187[k] = f_0 * psh_61[k]
                   + f_3 * pc_x[k] * pph_145[k];

        t_188[k] = f_0 * psh_62[k]
                   + f_3 * pc_x[k] * pph_146[k];

        t_189[k] = pa_z[k] * spi0_21[k]
                   - f_10 * pc_z[k] * spi1_21[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pb_x, pc_x, psi0_78, psi0_79, psi0_80, \
                         psi0_81, psi1_78, psi1_79, psi1_80, psi1_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pb_x[k] * psi0_78[k]
                   - f_10 * pc_x[k] * psi1_78[k];

        t_191[k] = pb_x[k] * psi0_79[k]
                   - f_10 * pc_x[k] * psi1_79[k];

        t_192[k] = pb_x[k] * psi0_80[k]
                   - f_10 * pc_x[k] * psi1_80[k];

        t_193[k] = pb_x[k] * psi0_81[k]
                   - f_10 * pc_x[k] * psi1_81[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pa_z, pb_x, pc_x, pc_y, pc_z, spi0_28, \
                         spi1_28, psi0_83, psh_42, psi1_83, pph_146, \
                         pph_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_3 * pc_y[k] * pph_146[k];

        t_195[k] = pb_x[k] * psi0_83[k]
                   - f_10 * pc_x[k] * psi1_83[k];

        t_196[k] = pa_z[k] * spi0_28[k]
                   - f_10 * pc_z[k] * spi1_28[k];

        t_197[k] = f_0 * psh_42[k]
                   + f_3 * pc_y[k] * pph_147[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_z, pb_y, pc_y, pc_z, spi0_31, spi1_31, \
                         psi0_58, psi0_61, psh_44, psi1_58, psi1_61, \
                         pph_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * psi0_58[k]
                   - f_10 * pc_y[k] * psi1_58[k];

        t_199[k] = pa_z[k] * spi0_31[k]
                   - f_10 * pc_z[k] * spi1_31[k];

        t_200[k] = f_0 * psh_44[k]
                   + f_3 * pc_y[k] * pph_149[k];

        t_201[k] = pb_y[k] * psi0_61[k]
                   - f_10 * pc_y[k] * psi1_61[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pa_z, pc_x, pc_y, pc_z, spi0_34, spi1_34, \
                         psh_47, ppg0_112, ppg1_112, pph_152, pph_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = pa_z[k] * spi0_34[k]
                   - f_10 * pc_z[k] * spi1_34[k];

        t_203[k] = f_6 * ppg0_112[k]
                   - f_7 * ppg1_112[k]
                   + f_3 * pc_x[k] * pph_154[k];

        t_204[k] = f_0 * psh_47[k]
                   + f_3 * pc_y[k] * pph_152[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pa_z, pb_y, pc_x, pc_y, pc_z, spi0_38, spi1_38, \
                         psi0_65, psi1_65, ppg0_116, ppg1_116, \
                         pph_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = pb_y[k] * psi0_65[k]
                   - f_10 * pc_y[k] * psi1_65[k];

        t_206[k] = pa_z[k] * spi0_38[k]
                   - f_10 * pc_z[k] * spi1_38[k];

        t_207[k] = f_4 * ppg0_116[k]
                   - f_5 * ppg1_116[k]
                   + f_3 * pc_x[k] * pph_158[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, psi0_70, psh_51, \
                         psi1_70, ppg0_117, ppg1_117, pph_156, pph_159, \
                         pph_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_4 * ppg0_117[k]
                   - f_5 * ppg1_117[k]
                   + f_3 * pc_x[k] * pph_159[k];

        t_209[k] = f_0 * psh_51[k]
                   + f_3 * pc_y[k] * pph_156[k];

        t_210[k] = pb_y[k] * psi0_70[k]
                   - f_10 * pc_y[k] * psi1_70[k];

        t_211[k] = f_3 * pc_x[k] * pph_162[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, spi0_49, \
                         spi1_49, pph_163, pph_164, pph_165, pph_166, \
                         pph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_3 * pc_x[k] * pph_163[k];

        t_213[k] = f_3 * pc_x[k] * pph_164[k];

        t_214[k] = f_3 * pc_x[k] * pph_165[k];

        t_215[k] = f_3 * pc_x[k] * pph_166[k];

        t_216[k] = f_3 * pc_x[k] * pph_167[k];

        t_217[k] = pa_z[k] * spi0_49[k]
                   - f_10 * pc_z[k] * spi1_49[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pa_z, pc_z, spi0_50, spi0_51, spi0_52, sph_36, \
                         sph_37, sph_38, spi1_50, spi1_51, spi1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = pa_z[k] * spi0_50[k]
                   + f_0 * sph_36[k]
                   - f_10 * pc_z[k] * spi1_50[k];

        t_219[k] = pa_z[k] * spi0_51[k]
                   + f_13 * sph_37[k]
                   - f_10 * pc_z[k] * spi1_51[k];

        t_220[k] = pa_z[k] * spi0_52[k]
                   + f_12 * sph_38[k]
                   - f_10 * pc_z[k] * spi1_52[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pa_z, pb_y, pc_y, pc_z, spi0_53, sph_39, \
                         spi1_53, psi0_83, psh_62, psi1_83, pph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pa_z[k] * spi0_53[k]
                   + f_11 * sph_39[k]
                   - f_10 * pc_z[k] * spi1_53[k];

        t_222[k] = f_0 * psh_62[k]
                   + f_3 * pc_y[k] * pph_167[k];

        t_223[k] = pb_y[k] * psi0_83[k]
                   - f_10 * pc_y[k] * psi1_83[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pc_x, pc_y, ppg0_120, ppg0_122, \
                         ppg0_123, ppg1_120, ppg1_122, ppg1_123, pph_168, pph_170, \
                         pph_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_1 * ppg0_120[k]
                   - f_2 * ppg1_120[k]
                   + f_3 * pc_x[k] * pph_168[k];

        t_225[k] = f_3 * pc_y[k] * pph_168[k];

        t_226[k] = f_14 * ppg0_122[k]
                   - f_15 * ppg1_122[k]
                   + f_3 * pc_x[k] * pph_170[k];

        t_227[k] = f_8 * ppg0_123[k]
                   - f_9 * ppg1_123[k]
                   + f_3 * pc_x[k] * pph_171[k];

        t_228[k] = f_3 * pc_y[k] * pph_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pc_x, pc_y, ppg0_125, ppg0_126, ppg0_127, \
                         ppg1_125, ppg1_126, ppg1_127, pph_173, pph_174, \
                         pph_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_8 * ppg0_125[k]
                   - f_9 * ppg1_125[k]
                   + f_3 * pc_x[k] * pph_173[k];

        t_230[k] = f_6 * ppg0_126[k]
                   - f_7 * ppg1_126[k]
                   + f_3 * pc_x[k] * pph_174[k];

        t_231[k] = f_6 * ppg0_127[k]
                   - f_7 * ppg1_127[k]
                   + f_3 * pc_x[k] * pph_175[k];

        t_232[k] = f_3 * pc_y[k] * pph_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, ppg0_129, ppg0_130, ppg0_131, ppg1_129, \
                         ppg1_130, ppg1_131, pph_177, pph_178, \
                         pph_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_6 * ppg0_129[k]
                   - f_7 * ppg1_129[k]
                   + f_3 * pc_x[k] * pph_177[k];

        t_234[k] = f_4 * ppg0_130[k]
                   - f_5 * ppg1_130[k]
                   + f_3 * pc_x[k] * pph_178[k];

        t_235[k] = f_4 * ppg0_131[k]
                   - f_5 * ppg1_131[k]
                   + f_3 * pc_x[k] * pph_179[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, t_240, pc_x, pc_y, ppg0_132, ppg0_134, \
                         ppg1_132, ppg1_134, pph_177, pph_180, pph_182, pph_183, \
                         pph_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_4 * ppg0_132[k]
                   - f_5 * ppg1_132[k]
                   + f_3 * pc_x[k] * pph_180[k];

        t_237[k] = f_3 * pc_y[k] * pph_177[k];

        t_238[k] = f_4 * ppg0_134[k]
                   - f_5 * ppg1_134[k]
                   + f_3 * pc_x[k] * pph_182[k];

        t_239[k] = f_3 * pc_x[k] * pph_183[k];

        t_240[k] = f_3 * pc_x[k] * pph_184[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pc_x, pc_y, ppg0_130, ppg1_130, \
                         pph_183, pph_185, pph_186, pph_187, pph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_3 * pc_x[k] * pph_185[k];

        t_242[k] = f_3 * pc_x[k] * pph_186[k];

        t_243[k] = f_3 * pc_x[k] * pph_187[k];

        t_244[k] = f_3 * pc_x[k] * pph_188[k];

        t_245[k] = f_1 * ppg0_130[k]
                   - f_2 * ppg1_130[k]
                   + f_3 * pc_y[k] * pph_183[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pc_y, ppg0_131, ppg0_132, ppg0_133, ppg1_131, \
                         ppg1_132, ppg1_133, pph_184, pph_185, \
                         pph_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_14 * ppg0_131[k]
                   - f_15 * ppg1_131[k]
                   + f_3 * pc_y[k] * pph_184[k];

        t_247[k] = f_8 * ppg0_132[k]
                   - f_9 * ppg1_132[k]
                   + f_3 * pc_y[k] * pph_185[k];

        t_248[k] = f_6 * ppg0_133[k]
                   - f_7 * ppg1_133[k]
                   + f_3 * pc_y[k] * pph_186[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, pc_z, sph_62, psh_62, ppg0_134, ppg1_134, \
                         pph_187, pph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_4 * ppg0_134[k]
                   - f_5 * ppg1_134[k]
                   + f_3 * pc_y[k] * pph_187[k];

        t_250[k] = f_3 * pc_y[k] * pph_188[k];

        t_251[k] = f_0 * sph_62[k]
                   + f_0 * psh_62[k]
                   + f_1 * ppg0_134[k]
                   - f_2 * ppg1_134[k]
                   + f_3 * pc_z[k] * pph_188[k];
    }
}

auto
compute_prim_ppi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t spi0,
                                                   const size_t sph, const size_t spi1,
                                                   const size_t psi0, const size_t psh,
                                                   const size_t psi1, const size_t ppg0,
                                                   const size_t ppg1, const size_t pph,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ppi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, spi0,
                                                              sph, spi1, psi0, psh, psi1, ppg0,
                                                              ppg1, pph, ncols, gamma, p, q);

    compute_prim_ppi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, spi0,
                                                              sph, spi1, psi0, psh, psi1, ppg0,
                                                              ppg1, pph, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
