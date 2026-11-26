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


#include "SimdThreeCenterElectronRepulsionVrrRecPDI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pdi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdi0, const size_t sdh,
                                                          const size_t sdi1, const size_t ppi0,
                                                          const size_t pph, const size_t ppi1,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdi0_84 = buffer.data(sdi0 + 84);
    const auto *sdi0_87 = buffer.data(sdi0 + 87);
    const auto *sdi0_90 = buffer.data(sdi0 + 90);
    const auto *sdi0_94 = buffer.data(sdi0 + 94);
    const auto *sdi0_96 = buffer.data(sdi0 + 96);
    const auto *sdi0_105 = buffer.data(sdi0 + 105);
    const auto *sdi0_107 = buffer.data(sdi0 + 107);
    const auto *sdi0_108 = buffer.data(sdi0 + 108);
    const auto *sdi0_109 = buffer.data(sdi0 + 109);
    const auto *sdi0_111 = buffer.data(sdi0 + 111);
    const auto *sdi0_124 = buffer.data(sdi0 + 124);

    const auto *sdh_0 = buffer.data(sdh + 0);
    const auto *sdh_15 = buffer.data(sdh + 15);
    const auto *sdh_17 = buffer.data(sdh + 17);
    const auto *sdh_18 = buffer.data(sdh + 18);
    const auto *sdh_20 = buffer.data(sdh + 20);
    const auto *sdh_36 = buffer.data(sdh + 36);
    const auto *sdh_38 = buffer.data(sdh + 38);
    const auto *sdh_39 = buffer.data(sdh + 39);
    const auto *sdh_59 = buffer.data(sdh + 59);
    const auto *sdh_60 = buffer.data(sdh + 60);
    const auto *sdh_62 = buffer.data(sdh + 62);
    const auto *sdh_63 = buffer.data(sdh + 63);
    const auto *sdh_66 = buffer.data(sdh + 66);
    const auto *sdh_69 = buffer.data(sdh + 69);
    const auto *sdh_73 = buffer.data(sdh + 73);
    const auto *sdh_75 = buffer.data(sdh + 75);
    const auto *sdh_78 = buffer.data(sdh + 78);
    const auto *sdh_80 = buffer.data(sdh + 80);
    const auto *sdh_81 = buffer.data(sdh + 81);
    const auto *sdh_83 = buffer.data(sdh + 83);
    const auto *sdh_96 = buffer.data(sdh + 96);

    const auto *sdi1_84 = buffer.data(sdi1 + 84);
    const auto *sdi1_87 = buffer.data(sdi1 + 87);
    const auto *sdi1_90 = buffer.data(sdi1 + 90);
    const auto *sdi1_94 = buffer.data(sdi1 + 94);
    const auto *sdi1_96 = buffer.data(sdi1 + 96);
    const auto *sdi1_105 = buffer.data(sdi1 + 105);
    const auto *sdi1_107 = buffer.data(sdi1 + 107);
    const auto *sdi1_108 = buffer.data(sdi1 + 108);
    const auto *sdi1_109 = buffer.data(sdi1 + 109);
    const auto *sdi1_111 = buffer.data(sdi1 + 111);
    const auto *sdi1_124 = buffer.data(sdi1 + 124);

    const auto *ppi0_0 = buffer.data(ppi0 + 0);
    const auto *ppi0_3 = buffer.data(ppi0 + 3);
    const auto *ppi0_5 = buffer.data(ppi0 + 5);
    const auto *ppi0_6 = buffer.data(ppi0 + 6);
    const auto *ppi0_9 = buffer.data(ppi0 + 9);
    const auto *ppi0_10 = buffer.data(ppi0 + 10);
    const auto *ppi0_14 = buffer.data(ppi0 + 14);
    const auto *ppi0_15 = buffer.data(ppi0 + 15);
    const auto *ppi0_20 = buffer.data(ppi0 + 20);
    const auto *ppi0_31 = buffer.data(ppi0 + 31);
    const auto *ppi0_34 = buffer.data(ppi0 + 34);
    const auto *ppi0_38 = buffer.data(ppi0 + 38);
    const auto *ppi0_56 = buffer.data(ppi0 + 56);
    const auto *ppi0_61 = buffer.data(ppi0 + 61);
    const auto *ppi0_65 = buffer.data(ppi0 + 65);

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
    const auto *pph_35 = buffer.data(pph + 35);
    const auto *pph_36 = buffer.data(pph + 36);
    const auto *pph_38 = buffer.data(pph + 38);
    const auto *pph_39 = buffer.data(pph + 39);
    const auto *pph_41 = buffer.data(pph + 41);
    const auto *pph_42 = buffer.data(pph + 42);
    const auto *pph_44 = buffer.data(pph + 44);
    const auto *pph_47 = buffer.data(pph + 47);
    const auto *pph_51 = buffer.data(pph + 51);
    const auto *pph_59 = buffer.data(pph + 59);
    const auto *pph_60 = buffer.data(pph + 60);
    const auto *pph_62 = buffer.data(pph + 62);

    const auto *ppi1_0 = buffer.data(ppi1 + 0);
    const auto *ppi1_3 = buffer.data(ppi1 + 3);
    const auto *ppi1_5 = buffer.data(ppi1 + 5);
    const auto *ppi1_6 = buffer.data(ppi1 + 6);
    const auto *ppi1_9 = buffer.data(ppi1 + 9);
    const auto *ppi1_10 = buffer.data(ppi1 + 10);
    const auto *ppi1_14 = buffer.data(ppi1 + 14);
    const auto *ppi1_15 = buffer.data(ppi1 + 15);
    const auto *ppi1_20 = buffer.data(ppi1 + 20);
    const auto *ppi1_31 = buffer.data(ppi1 + 31);
    const auto *ppi1_34 = buffer.data(ppi1 + 34);
    const auto *ppi1_38 = buffer.data(ppi1 + 38);
    const auto *ppi1_56 = buffer.data(ppi1 + 56);
    const auto *ppi1_61 = buffer.data(ppi1 + 61);
    const auto *ppi1_65 = buffer.data(ppi1 + 65);

    const auto *pdg0_0 = buffer.data(pdg0 + 0);
    const auto *pdg0_1 = buffer.data(pdg0 + 1);
    const auto *pdg0_2 = buffer.data(pdg0 + 2);
    const auto *pdg0_3 = buffer.data(pdg0 + 3);
    const auto *pdg0_5 = buffer.data(pdg0 + 5);
    const auto *pdg0_10 = buffer.data(pdg0 + 10);
    const auto *pdg0_12 = buffer.data(pdg0 + 12);
    const auto *pdg0_13 = buffer.data(pdg0 + 13);
    const auto *pdg0_14 = buffer.data(pdg0 + 14);
    const auto *pdg0_15 = buffer.data(pdg0 + 15);
    const auto *pdg0_16 = buffer.data(pdg0 + 16);
    const auto *pdg0_18 = buffer.data(pdg0 + 18);
    const auto *pdg0_20 = buffer.data(pdg0 + 20);
    const auto *pdg0_25 = buffer.data(pdg0 + 25);
    const auto *pdg0_27 = buffer.data(pdg0 + 27);
    const auto *pdg0_28 = buffer.data(pdg0 + 28);
    const auto *pdg0_29 = buffer.data(pdg0 + 29);
    const auto *pdg0_30 = buffer.data(pdg0 + 30);
    const auto *pdg0_32 = buffer.data(pdg0 + 32);
    const auto *pdg0_35 = buffer.data(pdg0 + 35);
    const auto *pdg0_40 = buffer.data(pdg0 + 40);
    const auto *pdg0_42 = buffer.data(pdg0 + 42);
    const auto *pdg0_43 = buffer.data(pdg0 + 43);
    const auto *pdg0_44 = buffer.data(pdg0 + 44);
    const auto *pdg0_45 = buffer.data(pdg0 + 45);
    const auto *pdg0_47 = buffer.data(pdg0 + 47);
    const auto *pdg0_50 = buffer.data(pdg0 + 50);

    const auto *pdg1_0 = buffer.data(pdg1 + 0);
    const auto *pdg1_1 = buffer.data(pdg1 + 1);
    const auto *pdg1_2 = buffer.data(pdg1 + 2);
    const auto *pdg1_3 = buffer.data(pdg1 + 3);
    const auto *pdg1_5 = buffer.data(pdg1 + 5);
    const auto *pdg1_10 = buffer.data(pdg1 + 10);
    const auto *pdg1_12 = buffer.data(pdg1 + 12);
    const auto *pdg1_13 = buffer.data(pdg1 + 13);
    const auto *pdg1_14 = buffer.data(pdg1 + 14);
    const auto *pdg1_15 = buffer.data(pdg1 + 15);
    const auto *pdg1_16 = buffer.data(pdg1 + 16);
    const auto *pdg1_18 = buffer.data(pdg1 + 18);
    const auto *pdg1_20 = buffer.data(pdg1 + 20);
    const auto *pdg1_25 = buffer.data(pdg1 + 25);
    const auto *pdg1_27 = buffer.data(pdg1 + 27);
    const auto *pdg1_28 = buffer.data(pdg1 + 28);
    const auto *pdg1_29 = buffer.data(pdg1 + 29);
    const auto *pdg1_30 = buffer.data(pdg1 + 30);
    const auto *pdg1_32 = buffer.data(pdg1 + 32);
    const auto *pdg1_35 = buffer.data(pdg1 + 35);
    const auto *pdg1_40 = buffer.data(pdg1 + 40);
    const auto *pdg1_42 = buffer.data(pdg1 + 42);
    const auto *pdg1_43 = buffer.data(pdg1 + 43);
    const auto *pdg1_44 = buffer.data(pdg1 + 44);
    const auto *pdg1_45 = buffer.data(pdg1 + 45);
    const auto *pdg1_47 = buffer.data(pdg1 + 47);
    const auto *pdg1_50 = buffer.data(pdg1 + 50);

    const auto *pdh_0 = buffer.data(pdh + 0);
    const auto *pdh_1 = buffer.data(pdh + 1);
    const auto *pdh_2 = buffer.data(pdh + 2);
    const auto *pdh_3 = buffer.data(pdh + 3);
    const auto *pdh_5 = buffer.data(pdh + 5);
    const auto *pdh_6 = buffer.data(pdh + 6);
    const auto *pdh_8 = buffer.data(pdh + 8);
    const auto *pdh_9 = buffer.data(pdh + 9);
    const auto *pdh_10 = buffer.data(pdh + 10);
    const auto *pdh_14 = buffer.data(pdh + 14);
    const auto *pdh_15 = buffer.data(pdh + 15);
    const auto *pdh_17 = buffer.data(pdh + 17);
    const auto *pdh_18 = buffer.data(pdh + 18);
    const auto *pdh_19 = buffer.data(pdh + 19);
    const auto *pdh_20 = buffer.data(pdh + 20);
    const auto *pdh_21 = buffer.data(pdh + 21);
    const auto *pdh_22 = buffer.data(pdh + 22);
    const auto *pdh_23 = buffer.data(pdh + 23);
    const auto *pdh_24 = buffer.data(pdh + 24);
    const auto *pdh_26 = buffer.data(pdh + 26);
    const auto *pdh_27 = buffer.data(pdh + 27);
    const auto *pdh_29 = buffer.data(pdh + 29);
    const auto *pdh_30 = buffer.data(pdh + 30);
    const auto *pdh_31 = buffer.data(pdh + 31);
    const auto *pdh_35 = buffer.data(pdh + 35);
    const auto *pdh_36 = buffer.data(pdh + 36);
    const auto *pdh_38 = buffer.data(pdh + 38);
    const auto *pdh_39 = buffer.data(pdh + 39);
    const auto *pdh_40 = buffer.data(pdh + 40);
    const auto *pdh_41 = buffer.data(pdh + 41);
    const auto *pdh_42 = buffer.data(pdh + 42);
    const auto *pdh_44 = buffer.data(pdh + 44);
    const auto *pdh_45 = buffer.data(pdh + 45);
    const auto *pdh_47 = buffer.data(pdh + 47);
    const auto *pdh_48 = buffer.data(pdh + 48);
    const auto *pdh_50 = buffer.data(pdh + 50);
    const auto *pdh_51 = buffer.data(pdh + 51);
    const auto *pdh_52 = buffer.data(pdh + 52);
    const auto *pdh_56 = buffer.data(pdh + 56);
    const auto *pdh_57 = buffer.data(pdh + 57);
    const auto *pdh_59 = buffer.data(pdh + 59);
    const auto *pdh_60 = buffer.data(pdh + 60);
    const auto *pdh_61 = buffer.data(pdh + 61);
    const auto *pdh_62 = buffer.data(pdh + 62);
    const auto *pdh_63 = buffer.data(pdh + 63);
    const auto *pdh_65 = buffer.data(pdh + 65);
    const auto *pdh_66 = buffer.data(pdh + 66);
    const auto *pdh_68 = buffer.data(pdh + 68);
    const auto *pdh_69 = buffer.data(pdh + 69);
    const auto *pdh_72 = buffer.data(pdh + 72);
    const auto *pdh_73 = buffer.data(pdh + 73);
    const auto *pdh_77 = buffer.data(pdh + 77);
    const auto *pdh_78 = buffer.data(pdh + 78);
    const auto *pdh_80 = buffer.data(pdh + 80);
    const auto *pdh_81 = buffer.data(pdh + 81);
    const auto *pdh_83 = buffer.data(pdh + 83);
    const auto *pdh_84 = buffer.data(pdh + 84);
    const auto *pdh_86 = buffer.data(pdh + 86);
    const auto *pdh_87 = buffer.data(pdh + 87);
    const auto *pdh_89 = buffer.data(pdh + 89);
    const auto *pdh_90 = buffer.data(pdh + 90);
    const auto *pdh_93 = buffer.data(pdh + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sdh_0, pph_0, pdg0_0, \
                         pdg1_0, pdh_0, pdh_1, pdh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdh_0[k]
                 + f_1 * pph_0[k]
                 + f_2 * pdg0_0[k]
                 - f_3 * pdg1_0[k]
                 + f_4 * pc_x[k] * pdh_0[k];

        t_1[k] = f_4 * pc_y[k] * pdh_0[k];

        t_2[k] = f_4 * pc_z[k] * pdh_0[k];

        t_3[k] = f_5 * pdg0_0[k]
                 - f_6 * pdg1_0[k]
                 + f_4 * pc_y[k] * pdh_1[k];

        t_4[k] = f_4 * pc_y[k] * pdh_2[k];

        t_5[k] = f_5 * pdg0_0[k]
                 - f_6 * pdg1_0[k]
                 + f_4 * pc_z[k] * pdh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, pdg0_1, pdg0_2, pdg0_3, pdg1_1, \
                         pdg1_2, pdg1_3, pdh_3, pdh_5, pdh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * pdg0_1[k]
                 - f_8 * pdg1_1[k]
                 + f_4 * pc_y[k] * pdh_3[k];

        t_7[k] = f_4 * pc_z[k] * pdh_3[k];

        t_8[k] = f_4 * pc_y[k] * pdh_5[k];

        t_9[k] = f_7 * pdg0_2[k]
                 - f_8 * pdg1_2[k]
                 + f_4 * pc_z[k] * pdh_5[k];

        t_10[k] = f_9 * pdg0_3[k]
                  - f_10 * pdg1_3[k]
                  + f_4 * pc_y[k] * pdh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, sdh_15, pph_15, \
                         pdg0_5, pdg1_5, pdh_6, pdh_8, pdh_9, pdh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * pc_z[k] * pdh_6[k];

        t_12[k] = f_5 * pdg0_5[k]
                  - f_6 * pdg1_5[k]
                  + f_4 * pc_y[k] * pdh_8[k];

        t_13[k] = f_4 * pc_y[k] * pdh_9[k];

        t_14[k] = f_9 * pdg0_5[k]
                  - f_10 * pdg1_5[k]
                  + f_4 * pc_z[k] * pdh_9[k];

        t_15[k] = f_0 * sdh_15[k]
                  + f_1 * pph_15[k]
                  + f_4 * pc_x[k] * pdh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, sdh_17, sdh_18, pph_17, \
                         pph_18, pdh_10, pdh_14, pdh_17, pdh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * pc_z[k] * pdh_10[k];

        t_17[k] = f_0 * sdh_17[k]
                  + f_1 * pph_17[k]
                  + f_4 * pc_x[k] * pdh_17[k];

        t_18[k] = f_0 * sdh_18[k]
                  + f_1 * pph_18[k]
                  + f_4 * pc_x[k] * pdh_18[k];

        t_19[k] = f_4 * pc_y[k] * pdh_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sdh_20, pph_20, pdg0_10, \
                         pdg0_12, pdg1_10, pdg1_12, pdh_15, pdh_17, \
                         pdh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * sdh_20[k]
                  + f_1 * pph_20[k]
                  + f_4 * pc_x[k] * pdh_20[k];

        t_21[k] = f_2 * pdg0_10[k]
                  - f_3 * pdg1_10[k]
                  + f_4 * pc_y[k] * pdh_15[k];

        t_22[k] = f_4 * pc_z[k] * pdh_15[k];

        t_23[k] = f_9 * pdg0_12[k]
                  - f_10 * pdg1_12[k]
                  + f_4 * pc_y[k] * pdh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, pdg0_13, pdg0_14, pdg1_13, \
                         pdg1_14, pdh_18, pdh_19, pdh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * pdg0_13[k]
                  - f_8 * pdg1_13[k]
                  + f_4 * pc_y[k] * pdh_18[k];

        t_25[k] = f_5 * pdg0_14[k]
                  - f_6 * pdg1_14[k]
                  + f_4 * pc_y[k] * pdh_19[k];

        t_26[k] = f_4 * pc_y[k] * pdh_20[k];

        t_27[k] = f_2 * pdg0_14[k]
                  - f_3 * pdg1_14[k]
                  + f_4 * pc_z[k] * pdh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_y, pc_z, ppi0_0, pph_0, pph_1, \
                         ppi1_0, pdg0_15, pdg1_15, pdh_21, pdh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * ppi0_0[k]
                  - f_11 * pc_y[k] * ppi1_0[k];

        t_29[k] = f_0 * pph_0[k]
                  + f_4 * pc_y[k] * pdh_21[k];

        t_30[k] = f_4 * pc_z[k] * pdh_21[k];

        t_31[k] = f_0 * pph_1[k]
                  + f_5 * pdg0_15[k]
                  - f_6 * pdg1_15[k]
                  + f_4 * pc_y[k] * pdh_22[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pc_y, pc_z, ppi0_5, pph_2, pph_3, \
                         ppi1_5, pdg0_16, pdg1_16, pdh_23, pdh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pph_2[k]
                  + f_4 * pc_y[k] * pdh_23[k];

        t_33[k] = pb_y[k] * ppi0_5[k]
                  - f_11 * pc_y[k] * ppi1_5[k];

        t_34[k] = f_0 * pph_3[k]
                  + f_7 * pdg0_16[k]
                  - f_8 * pdg1_16[k]
                  + f_4 * pc_y[k] * pdh_24[k];

        t_35[k] = f_4 * pc_z[k] * pdh_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, ppi0_9, pph_5, pph_6, \
                         ppi1_9, pdg0_18, pdg1_18, pdh_26, pdh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * pph_5[k]
                  + f_4 * pc_y[k] * pdh_26[k];

        t_37[k] = pb_y[k] * ppi0_9[k]
                  - f_11 * pc_y[k] * ppi1_9[k];

        t_38[k] = f_0 * pph_6[k]
                  + f_9 * pdg0_18[k]
                  - f_10 * pdg1_18[k]
                  + f_4 * pc_y[k] * pdh_27[k];

        t_39[k] = f_4 * pc_z[k] * pdh_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pc_y, ppi0_14, pph_8, pph_9, ppi1_14, \
                         pdg0_20, pdg1_20, pdh_29, pdh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * pph_8[k]
                  + f_5 * pdg0_20[k]
                  - f_6 * pdg1_20[k]
                  + f_4 * pc_y[k] * pdh_29[k];

        t_41[k] = f_0 * pph_9[k]
                  + f_4 * pc_y[k] * pdh_30[k];

        t_42[k] = pb_y[k] * ppi0_14[k]
                  - f_11 * pc_y[k] * ppi1_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pc_x, pc_z, sdh_36, sdh_38, sdh_39, pph_36, \
                         pph_38, pph_39, pdh_31, pdh_36, pdh_38, \
                         pdh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * sdh_36[k]
                  + f_0 * pph_36[k]
                  + f_4 * pc_x[k] * pdh_36[k];

        t_44[k] = f_4 * pc_z[k] * pdh_31[k];

        t_45[k] = f_0 * sdh_38[k]
                  + f_0 * pph_38[k]
                  + f_4 * pc_x[k] * pdh_38[k];

        t_46[k] = f_0 * sdh_39[k]
                  + f_0 * pph_39[k]
                  + f_4 * pc_x[k] * pdh_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pc_y, pc_z, ppi0_20, pph_14, pph_15, \
                         ppi1_20, pdg0_25, pdg1_25, pdh_35, pdh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * pph_14[k]
                  + f_4 * pc_y[k] * pdh_35[k];

        t_48[k] = pb_y[k] * ppi0_20[k]
                  - f_11 * pc_y[k] * ppi1_20[k];

        t_49[k] = f_0 * pph_15[k]
                  + f_2 * pdg0_25[k]
                  - f_3 * pdg1_25[k]
                  + f_4 * pc_y[k] * pdh_36[k];

        t_50[k] = f_4 * pc_z[k] * pdh_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pc_y, pph_17, pph_18, pph_19, pdg0_27, pdg0_28, \
                         pdg0_29, pdg1_27, pdg1_28, pdg1_29, pdh_38, pdh_39, \
                         pdh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * pph_17[k]
                  + f_9 * pdg0_27[k]
                  - f_10 * pdg1_27[k]
                  + f_4 * pc_y[k] * pdh_38[k];

        t_52[k] = f_0 * pph_18[k]
                  + f_7 * pdg0_28[k]
                  - f_8 * pdg1_28[k]
                  + f_4 * pc_y[k] * pdh_39[k];

        t_53[k] = f_0 * pph_19[k]
                  + f_5 * pdg0_29[k]
                  - f_6 * pdg1_29[k]
                  + f_4 * pc_y[k] * pdh_40[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pb_z, pc_y, pc_z, ppi0_0, pph_0, \
                         pph_20, ppi1_0, pdg0_29, pdg1_29, pdh_41, \
                         pdh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * pph_20[k]
                  + f_4 * pc_y[k] * pdh_41[k];

        t_55[k] = f_2 * pdg0_29[k]
                  - f_3 * pdg1_29[k]
                  + f_4 * pc_z[k] * pdh_41[k];

        t_56[k] = pb_z[k] * ppi0_0[k]
                  - f_11 * pc_z[k] * ppi1_0[k];

        t_57[k] = f_4 * pc_y[k] * pdh_42[k];

        t_58[k] = f_0 * pph_0[k]
                  + f_4 * pc_z[k] * pdh_42[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pb_z, pc_y, pc_z, ppi0_3, ppi0_6, pph_2, \
                         ppi1_3, ppi1_6, pdg0_30, pdg1_30, pdh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pb_z[k] * ppi0_3[k]
                  - f_11 * pc_z[k] * ppi1_3[k];

        t_60[k] = f_4 * pc_y[k] * pdh_44[k];

        t_61[k] = f_0 * pph_2[k]
                  + f_5 * pdg0_30[k]
                  - f_6 * pdg1_30[k]
                  + f_4 * pc_z[k] * pdh_44[k];

        t_62[k] = pb_z[k] * ppi0_6[k]
                  - f_11 * pc_z[k] * ppi1_6[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pb_z, pc_y, pc_z, ppi0_10, pph_3, pph_5, \
                         ppi1_10, pdg0_32, pdg1_32, pdh_45, pdh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * pph_3[k]
                  + f_4 * pc_z[k] * pdh_45[k];

        t_64[k] = f_4 * pc_y[k] * pdh_47[k];

        t_65[k] = f_0 * pph_5[k]
                  + f_7 * pdg0_32[k]
                  - f_8 * pdg1_32[k]
                  + f_4 * pc_z[k] * pdh_47[k];

        t_66[k] = pb_z[k] * ppi0_10[k]
                  - f_11 * pc_z[k] * ppi1_10[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pc_y, pc_z, pph_6, pph_9, pdg0_35, pdg1_35, \
                         pdh_48, pdh_50, pdh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_0 * pph_6[k]
                  + f_4 * pc_z[k] * pdh_48[k];

        t_68[k] = f_5 * pdg0_35[k]
                  - f_6 * pdg1_35[k]
                  + f_4 * pc_y[k] * pdh_50[k];

        t_69[k] = f_4 * pc_y[k] * pdh_51[k];

        t_70[k] = f_0 * pph_9[k]
                  + f_9 * pdg0_35[k]
                  - f_10 * pdg1_35[k]
                  + f_4 * pc_z[k] * pdh_51[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pb_z, pc_x, pc_z, sdh_59, ppi0_15, pph_10, pph_59, \
                         ppi1_15, pdh_52, pdh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_z[k] * ppi0_15[k]
                  - f_11 * pc_z[k] * ppi1_15[k];

        t_72[k] = f_0 * pph_10[k]
                  + f_4 * pc_z[k] * pdh_52[k];

        t_73[k] = f_0 * sdh_59[k]
                  + f_0 * pph_59[k]
                  + f_4 * pc_x[k] * pdh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, sdh_60, sdh_62, pph_60, pph_62, \
                         pdg0_40, pdg1_40, pdh_56, pdh_57, pdh_60, \
                         pdh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * sdh_60[k]
                  + f_0 * pph_60[k]
                  + f_4 * pc_x[k] * pdh_60[k];

        t_75[k] = f_4 * pc_y[k] * pdh_56[k];

        t_76[k] = f_0 * sdh_62[k]
                  + f_0 * pph_62[k]
                  + f_4 * pc_x[k] * pdh_62[k];

        t_77[k] = f_2 * pdg0_40[k]
                  - f_3 * pdg1_40[k]
                  + f_4 * pc_y[k] * pdh_57[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, pph_15, pdg0_42, pdg0_43, pdg1_42, \
                         pdg1_43, pdh_57, pdh_59, pdh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * pph_15[k]
                  + f_4 * pc_z[k] * pdh_57[k];

        t_79[k] = f_9 * pdg0_42[k]
                  - f_10 * pdg1_42[k]
                  + f_4 * pc_y[k] * pdh_59[k];

        t_80[k] = f_7 * pdg0_43[k]
                  - f_8 * pdg1_43[k]
                  + f_4 * pc_y[k] * pdh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_x, pc_x, pc_y, pc_z, sdi0_84, sdh_63, \
                         sdi1_84, pph_20, pdg0_44, pdg1_44, pdh_61, \
                         pdh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * pdg0_44[k]
                  - f_6 * pdg1_44[k]
                  + f_4 * pc_y[k] * pdh_61[k];

        t_82[k] = f_4 * pc_y[k] * pdh_62[k];

        t_83[k] = f_0 * pph_20[k]
                  + f_2 * pdg0_44[k]
                  - f_3 * pdg1_44[k]
                  + f_4 * pc_z[k] * pdh_62[k];

        t_84[k] = pa_x[k] * sdi0_84[k]
                  + f_12 * sdh_63[k]
                  - f_11 * pc_x[k] * sdi1_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_x, pc_x, pc_y, pc_z, sdi0_87, sdh_66, \
                         sdi1_87, pph_21, pph_23, pdh_63, pdh_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * pph_21[k]
                  + f_4 * pc_y[k] * pdh_63[k];

        t_86[k] = f_4 * pc_z[k] * pdh_63[k];

        t_87[k] = pa_x[k] * sdi0_87[k]
                  + f_13 * sdh_66[k]
                  - f_11 * pc_x[k] * sdi1_87[k];

        t_88[k] = f_1 * pph_23[k]
                  + f_4 * pc_y[k] * pdh_65[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pc_x, pc_z, sdi0_90, sdh_69, sdi1_90, \
                         pdg0_45, pdg1_45, pdh_65, pdh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_5 * pdg0_45[k]
                  - f_6 * pdg1_45[k]
                  + f_4 * pc_z[k] * pdh_65[k];

        t_90[k] = pa_x[k] * sdi0_90[k]
                  + f_14 * sdh_69[k]
                  - f_11 * pc_x[k] * sdi1_90[k];

        t_91[k] = f_4 * pc_z[k] * pdh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pc_x, pc_y, pc_z, sdi0_94, sdh_73, \
                         sdi1_94, pph_26, pdg0_47, pdg1_47, pdh_68, \
                         pdh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_1 * pph_26[k]
                  + f_4 * pc_y[k] * pdh_68[k];

        t_93[k] = f_7 * pdg0_47[k]
                  - f_8 * pdg1_47[k]
                  + f_4 * pc_z[k] * pdh_68[k];

        t_94[k] = pa_x[k] * sdi0_94[k]
                  + f_1 * sdh_73[k]
                  - f_11 * pc_x[k] * sdi1_94[k];

        t_95[k] = f_4 * pc_z[k] * pdh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_x, pc_x, pc_y, pc_z, sdi0_96, sdh_75, sdi1_96, \
                         pph_30, pdg0_50, pdg1_50, pdh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_x[k] * sdi0_96[k]
                  + f_1 * sdh_75[k]
                  - f_11 * pc_x[k] * sdi1_96[k];

        t_97[k] = f_1 * pph_30[k]
                  + f_4 * pc_y[k] * pdh_72[k];

        t_98[k] = f_9 * pdg0_50[k]
                  - f_10 * pdg1_50[k]
                  + f_4 * pc_z[k] * pdh_72[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_z, sdh_78, sdh_80, sdh_81, \
                         pdh_73, pdh_78, pdh_80, pdh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_0 * sdh_78[k]
                  + f_4 * pc_x[k] * pdh_78[k];

        t_100[k] = f_4 * pc_z[k] * pdh_73[k];

        t_101[k] = f_0 * sdh_80[k]
                   + f_4 * pc_x[k] * pdh_80[k];

        t_102[k] = f_0 * sdh_81[k]
                   + f_4 * pc_x[k] * pdh_81[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_x, pc_x, pc_y, pc_z, sdi0_105, sdh_83, \
                         sdi1_105, pph_35, pdh_77, pdh_78, pdh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * pph_35[k]
                   + f_4 * pc_y[k] * pdh_77[k];

        t_104[k] = f_0 * sdh_83[k]
                   + f_4 * pc_x[k] * pdh_83[k];

        t_105[k] = pa_x[k] * sdi0_105[k]
                   - f_11 * pc_x[k] * sdi1_105[k];

        t_106[k] = f_4 * pc_z[k] * pdh_78[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pc_x, pc_y, sdi0_107, sdi0_108, \
                         sdi0_109, sdi1_107, sdi1_108, sdi1_109, pph_41, \
                         pdh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * sdi0_107[k]
                   - f_11 * pc_x[k] * sdi1_107[k];

        t_108[k] = pa_x[k] * sdi0_108[k]
                   - f_11 * pc_x[k] * sdi1_108[k];

        t_109[k] = pa_x[k] * sdi0_109[k]
                   - f_11 * pc_x[k] * sdi1_109[k];

        t_110[k] = f_1 * pph_41[k]
                   + f_4 * pc_y[k] * pdh_83[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pb_y, pc_x, pc_y, pc_z, sdi0_111, \
                         sdi1_111, ppi0_56, pph_21, pph_42, ppi1_56, \
                         pdh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_x[k] * sdi0_111[k]
                   - f_11 * pc_x[k] * sdi1_111[k];

        t_112[k] = pb_y[k] * ppi0_56[k]
                   - f_11 * pc_y[k] * ppi1_56[k];

        t_113[k] = f_0 * pph_42[k]
                   + f_4 * pc_y[k] * pdh_84[k];

        t_114[k] = f_0 * pph_21[k]
                   + f_4 * pc_z[k] * pdh_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_y, pb_z, pc_y, pc_z, ppi0_31, ppi0_34, \
                         ppi0_61, pph_44, ppi1_31, ppi1_34, ppi1_61, \
                         pdh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * ppi0_31[k]
                   - f_11 * pc_z[k] * ppi1_31[k];

        t_116[k] = f_0 * pph_44[k]
                   + f_4 * pc_y[k] * pdh_86[k];

        t_117[k] = pb_y[k] * ppi0_61[k]
                   - f_11 * pc_y[k] * ppi1_61[k];

        t_118[k] = pb_z[k] * ppi0_34[k]
                   - f_11 * pc_z[k] * ppi1_34[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_y, pb_z, pc_y, pc_z, ppi0_38, ppi0_65, \
                         pph_24, pph_47, ppi1_38, ppi1_65, pdh_87, \
                         pdh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_0 * pph_24[k]
                   + f_4 * pc_z[k] * pdh_87[k];

        t_120[k] = f_0 * pph_47[k]
                   + f_4 * pc_y[k] * pdh_89[k];

        t_121[k] = pb_y[k] * ppi0_65[k]
                   - f_11 * pc_y[k] * ppi1_65[k];

        t_122[k] = pb_z[k] * ppi0_38[k]
                   - f_11 * pc_z[k] * ppi1_38[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pa_x, pc_x, pc_y, pc_z, sdi0_124, sdh_96, \
                         sdi1_124, pph_27, pph_51, pdh_90, pdh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * pph_27[k]
                   + f_4 * pc_z[k] * pdh_90[k];

        t_124[k] = pa_x[k] * sdi0_124[k]
                   + f_1 * sdh_96[k]
                   - f_11 * pc_x[k] * sdi1_124[k];

        t_125[k] = f_0 * pph_51[k]
                   + f_4 * pc_y[k] * pdh_93[k];
    }
}

static auto
compute_prim_pdi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdi0, const size_t sdh,
                                                          const size_t sdi1, const size_t ppi0,
                                                          const size_t pph, const size_t ppi1,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdi0_0 = buffer.data(sdi0 + 0);
    const auto *sdi0_1 = buffer.data(sdi0 + 1);
    const auto *sdi0_3 = buffer.data(sdi0 + 3);
    const auto *sdi0_5 = buffer.data(sdi0 + 5);
    const auto *sdi0_6 = buffer.data(sdi0 + 6);
    const auto *sdi0_8 = buffer.data(sdi0 + 8);
    const auto *sdi0_9 = buffer.data(sdi0 + 9);
    const auto *sdi0_10 = buffer.data(sdi0 + 10);
    const auto *sdi0_12 = buffer.data(sdi0 + 12);
    const auto *sdi0_13 = buffer.data(sdi0 + 13);
    const auto *sdi0_14 = buffer.data(sdi0 + 14);
    const auto *sdi0_21 = buffer.data(sdi0 + 21);
    const auto *sdi0_27 = buffer.data(sdi0 + 27);
    const auto *sdi0_56 = buffer.data(sdi0 + 56);
    const auto *sdi0_61 = buffer.data(sdi0 + 61);
    const auto *sdi0_64 = buffer.data(sdi0 + 64);
    const auto *sdi0_65 = buffer.data(sdi0 + 65);
    const auto *sdi0_68 = buffer.data(sdi0 + 68);
    const auto *sdi0_69 = buffer.data(sdi0 + 69);
    const auto *sdi0_70 = buffer.data(sdi0 + 70);
    const auto *sdi0_133 = buffer.data(sdi0 + 133);
    const auto *sdi0_135 = buffer.data(sdi0 + 135);
    const auto *sdi0_136 = buffer.data(sdi0 + 136);
    const auto *sdi0_137 = buffer.data(sdi0 + 137);
    const auto *sdi0_139 = buffer.data(sdi0 + 139);
    const auto *sdi0_140 = buffer.data(sdi0 + 140);
    const auto *sdi0_145 = buffer.data(sdi0 + 145);
    const auto *sdi0_149 = buffer.data(sdi0 + 149);
    const auto *sdi0_154 = buffer.data(sdi0 + 154);
    const auto *sdi0_161 = buffer.data(sdi0 + 161);
    const auto *sdi0_163 = buffer.data(sdi0 + 163);
    const auto *sdi0_164 = buffer.data(sdi0 + 164);
    const auto *sdi0_165 = buffer.data(sdi0 + 165);
    const auto *sdi0_167 = buffer.data(sdi0 + 167);

    const auto *sdh_0 = buffer.data(sdh + 0);
    const auto *sdh_1 = buffer.data(sdh + 1);
    const auto *sdh_3 = buffer.data(sdh + 3);
    const auto *sdh_5 = buffer.data(sdh + 5);
    const auto *sdh_6 = buffer.data(sdh + 6);
    const auto *sdh_8 = buffer.data(sdh + 8);
    const auto *sdh_9 = buffer.data(sdh + 9);
    const auto *sdh_15 = buffer.data(sdh + 15);
    const auto *sdh_20 = buffer.data(sdh + 20);
    const auto *sdh_47 = buffer.data(sdh + 47);
    const auto *sdh_50 = buffer.data(sdh + 50);
    const auto *sdh_51 = buffer.data(sdh + 51);
    const auto *sdh_99 = buffer.data(sdh + 99);
    const auto *sdh_101 = buffer.data(sdh + 101);
    const auto *sdh_102 = buffer.data(sdh + 102);
    const auto *sdh_104 = buffer.data(sdh + 104);
    const auto *sdh_105 = buffer.data(sdh + 105);
    const auto *sdh_110 = buffer.data(sdh + 110);
    const auto *sdh_114 = buffer.data(sdh + 114);
    const auto *sdh_119 = buffer.data(sdh + 119);
    const auto *sdh_120 = buffer.data(sdh + 120);
    const auto *sdh_122 = buffer.data(sdh + 122);
    const auto *sdh_123 = buffer.data(sdh + 123);
    const auto *sdh_125 = buffer.data(sdh + 125);

    const auto *sdi1_0 = buffer.data(sdi1 + 0);
    const auto *sdi1_1 = buffer.data(sdi1 + 1);
    const auto *sdi1_3 = buffer.data(sdi1 + 3);
    const auto *sdi1_5 = buffer.data(sdi1 + 5);
    const auto *sdi1_6 = buffer.data(sdi1 + 6);
    const auto *sdi1_8 = buffer.data(sdi1 + 8);
    const auto *sdi1_9 = buffer.data(sdi1 + 9);
    const auto *sdi1_10 = buffer.data(sdi1 + 10);
    const auto *sdi1_12 = buffer.data(sdi1 + 12);
    const auto *sdi1_13 = buffer.data(sdi1 + 13);
    const auto *sdi1_14 = buffer.data(sdi1 + 14);
    const auto *sdi1_21 = buffer.data(sdi1 + 21);
    const auto *sdi1_27 = buffer.data(sdi1 + 27);
    const auto *sdi1_56 = buffer.data(sdi1 + 56);
    const auto *sdi1_61 = buffer.data(sdi1 + 61);
    const auto *sdi1_64 = buffer.data(sdi1 + 64);
    const auto *sdi1_65 = buffer.data(sdi1 + 65);
    const auto *sdi1_68 = buffer.data(sdi1 + 68);
    const auto *sdi1_69 = buffer.data(sdi1 + 69);
    const auto *sdi1_70 = buffer.data(sdi1 + 70);
    const auto *sdi1_133 = buffer.data(sdi1 + 133);
    const auto *sdi1_135 = buffer.data(sdi1 + 135);
    const auto *sdi1_136 = buffer.data(sdi1 + 136);
    const auto *sdi1_137 = buffer.data(sdi1 + 137);
    const auto *sdi1_139 = buffer.data(sdi1 + 139);
    const auto *sdi1_140 = buffer.data(sdi1 + 140);
    const auto *sdi1_145 = buffer.data(sdi1 + 145);
    const auto *sdi1_149 = buffer.data(sdi1 + 149);
    const auto *sdi1_154 = buffer.data(sdi1 + 154);
    const auto *sdi1_161 = buffer.data(sdi1 + 161);
    const auto *sdi1_163 = buffer.data(sdi1 + 163);
    const auto *sdi1_164 = buffer.data(sdi1 + 164);
    const auto *sdi1_165 = buffer.data(sdi1 + 165);
    const auto *sdi1_167 = buffer.data(sdi1 + 167);

    const auto *ppi0_70 = buffer.data(ppi0 + 70);
    const auto *ppi0_85 = buffer.data(ppi0 + 85);
    const auto *ppi0_87 = buffer.data(ppi0 + 87);
    const auto *ppi0_90 = buffer.data(ppi0 + 90);
    const auto *ppi0_94 = buffer.data(ppi0 + 94);
    const auto *ppi0_113 = buffer.data(ppi0 + 113);
    const auto *ppi0_115 = buffer.data(ppi0 + 115);
    const auto *ppi0_118 = buffer.data(ppi0 + 118);
    const auto *ppi0_120 = buffer.data(ppi0 + 120);
    const auto *ppi0_122 = buffer.data(ppi0 + 122);
    const auto *ppi0_124 = buffer.data(ppi0 + 124);
    const auto *ppi0_125 = buffer.data(ppi0 + 125);
    const auto *ppi0_133 = buffer.data(ppi0 + 133);
    const auto *ppi0_135 = buffer.data(ppi0 + 135);
    const auto *ppi0_136 = buffer.data(ppi0 + 136);
    const auto *ppi0_137 = buffer.data(ppi0 + 137);
    const auto *ppi0_138 = buffer.data(ppi0 + 138);
    const auto *ppi0_139 = buffer.data(ppi0 + 139);
    const auto *ppi0_161 = buffer.data(ppi0 + 161);
    const auto *ppi0_163 = buffer.data(ppi0 + 163);
    const auto *ppi0_164 = buffer.data(ppi0 + 164);
    const auto *ppi0_165 = buffer.data(ppi0 + 165);

    const auto *pph_31 = buffer.data(pph + 31);
    const auto *pph_36 = buffer.data(pph + 36);
    const auto *pph_42 = buffer.data(pph + 42);
    const auto *pph_45 = buffer.data(pph + 45);
    const auto *pph_48 = buffer.data(pph + 48);
    const auto *pph_52 = buffer.data(pph + 52);
    const auto *pph_56 = buffer.data(pph + 56);
    const auto *pph_57 = buffer.data(pph + 57);
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
    const auto *pph_102 = buffer.data(pph + 102);
    const auto *pph_103 = buffer.data(pph + 103);
    const auto *pph_104 = buffer.data(pph + 104);
    const auto *pph_120 = buffer.data(pph + 120);
    const auto *pph_121 = buffer.data(pph + 121);
    const auto *pph_122 = buffer.data(pph + 122);
    const auto *pph_123 = buffer.data(pph + 123);
    const auto *pph_124 = buffer.data(pph + 124);
    const auto *pph_125 = buffer.data(pph + 125);

    const auto *ppi1_70 = buffer.data(ppi1 + 70);
    const auto *ppi1_85 = buffer.data(ppi1 + 85);
    const auto *ppi1_87 = buffer.data(ppi1 + 87);
    const auto *ppi1_90 = buffer.data(ppi1 + 90);
    const auto *ppi1_94 = buffer.data(ppi1 + 94);
    const auto *ppi1_113 = buffer.data(ppi1 + 113);
    const auto *ppi1_115 = buffer.data(ppi1 + 115);
    const auto *ppi1_118 = buffer.data(ppi1 + 118);
    const auto *ppi1_120 = buffer.data(ppi1 + 120);
    const auto *ppi1_122 = buffer.data(ppi1 + 122);
    const auto *ppi1_124 = buffer.data(ppi1 + 124);
    const auto *ppi1_125 = buffer.data(ppi1 + 125);
    const auto *ppi1_133 = buffer.data(ppi1 + 133);
    const auto *ppi1_135 = buffer.data(ppi1 + 135);
    const auto *ppi1_136 = buffer.data(ppi1 + 136);
    const auto *ppi1_137 = buffer.data(ppi1 + 137);
    const auto *ppi1_138 = buffer.data(ppi1 + 138);
    const auto *ppi1_139 = buffer.data(ppi1 + 139);
    const auto *ppi1_161 = buffer.data(ppi1 + 161);
    const auto *ppi1_163 = buffer.data(ppi1 + 163);
    const auto *ppi1_164 = buffer.data(ppi1 + 164);
    const auto *ppi1_165 = buffer.data(ppi1 + 165);

    const auto *pdg0_75 = buffer.data(pdg0 + 75);
    const auto *pdg0_76 = buffer.data(pdg0 + 76);
    const auto *pdg0_78 = buffer.data(pdg0 + 78);
    const auto *pdg0_80 = buffer.data(pdg0 + 80);
    const auto *pdg0_100 = buffer.data(pdg0 + 100);
    const auto *pdg0_101 = buffer.data(pdg0 + 101);
    const auto *pdg0_102 = buffer.data(pdg0 + 102);
    const auto *pdg0_105 = buffer.data(pdg0 + 105);
    const auto *pdg0_110 = buffer.data(pdg0 + 110);
    const auto *pdg0_114 = buffer.data(pdg0 + 114);
    const auto *pdg0_119 = buffer.data(pdg0 + 119);

    const auto *pdg1_75 = buffer.data(pdg1 + 75);
    const auto *pdg1_76 = buffer.data(pdg1 + 76);
    const auto *pdg1_78 = buffer.data(pdg1 + 78);
    const auto *pdg1_80 = buffer.data(pdg1 + 80);
    const auto *pdg1_100 = buffer.data(pdg1 + 100);
    const auto *pdg1_101 = buffer.data(pdg1 + 101);
    const auto *pdg1_102 = buffer.data(pdg1 + 102);
    const auto *pdg1_105 = buffer.data(pdg1 + 105);
    const auto *pdg1_110 = buffer.data(pdg1 + 110);
    const auto *pdg1_114 = buffer.data(pdg1 + 114);
    const auto *pdg1_119 = buffer.data(pdg1 + 119);

    const auto *pdh_94 = buffer.data(pdh + 94);
    const auto *pdh_98 = buffer.data(pdh + 98);
    const auto *pdh_99 = buffer.data(pdh + 99);
    const auto *pdh_101 = buffer.data(pdh + 101);
    const auto *pdh_102 = buffer.data(pdh + 102);
    const auto *pdh_104 = buffer.data(pdh + 104);
    const auto *pdh_105 = buffer.data(pdh + 105);
    const auto *pdh_106 = buffer.data(pdh + 106);
    const auto *pdh_107 = buffer.data(pdh + 107);
    const auto *pdh_108 = buffer.data(pdh + 108);
    const auto *pdh_110 = buffer.data(pdh + 110);
    const auto *pdh_111 = buffer.data(pdh + 111);
    const auto *pdh_113 = buffer.data(pdh + 113);
    const auto *pdh_114 = buffer.data(pdh + 114);
    const auto *pdh_115 = buffer.data(pdh + 115);
    const auto *pdh_119 = buffer.data(pdh + 119);
    const auto *pdh_120 = buffer.data(pdh + 120);
    const auto *pdh_122 = buffer.data(pdh + 122);
    const auto *pdh_123 = buffer.data(pdh + 123);
    const auto *pdh_125 = buffer.data(pdh + 125);
    const auto *pdh_126 = buffer.data(pdh + 126);
    const auto *pdh_127 = buffer.data(pdh + 127);
    const auto *pdh_129 = buffer.data(pdh + 129);
    const auto *pdh_132 = buffer.data(pdh + 132);
    const auto *pdh_141 = buffer.data(pdh + 141);
    const auto *pdh_142 = buffer.data(pdh + 142);
    const auto *pdh_143 = buffer.data(pdh + 143);
    const auto *pdh_144 = buffer.data(pdh + 144);
    const auto *pdh_145 = buffer.data(pdh + 145);
    const auto *pdh_146 = buffer.data(pdh + 146);
    const auto *pdh_147 = buffer.data(pdh + 147);
    const auto *pdh_148 = buffer.data(pdh + 148);
    const auto *pdh_150 = buffer.data(pdh + 150);
    const auto *pdh_152 = buffer.data(pdh + 152);
    const auto *pdh_153 = buffer.data(pdh + 153);
    const auto *pdh_156 = buffer.data(pdh + 156);
    const auto *pdh_161 = buffer.data(pdh + 161);
    const auto *pdh_162 = buffer.data(pdh + 162);
    const auto *pdh_163 = buffer.data(pdh + 163);
    const auto *pdh_164 = buffer.data(pdh + 164);
    const auto *pdh_165 = buffer.data(pdh + 165);
    const auto *pdh_166 = buffer.data(pdh + 166);
    const auto *pdh_167 = buffer.data(pdh + 167);
    const auto *pdh_168 = buffer.data(pdh + 168);
    const auto *pdh_169 = buffer.data(pdh + 169);
    const auto *pdh_171 = buffer.data(pdh + 171);
    const auto *pdh_174 = buffer.data(pdh + 174);
    const auto *pdh_183 = buffer.data(pdh + 183);
    const auto *pdh_184 = buffer.data(pdh + 184);
    const auto *pdh_185 = buffer.data(pdh + 185);
    const auto *pdh_186 = buffer.data(pdh + 186);
    const auto *pdh_187 = buffer.data(pdh + 187);
    const auto *pdh_188 = buffer.data(pdh + 188);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pb_y, pc_x, pc_y, pc_z, sdh_99, sdh_101, \
                         ppi0_70, pph_31, ppi1_70, pdh_94, pdh_99, \
                         pdh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pb_y[k] * ppi0_70[k]
                   - f_11 * pc_y[k] * ppi1_70[k];

        t_127[k] = f_0 * sdh_99[k]
                   + f_4 * pc_x[k] * pdh_99[k];

        t_128[k] = f_0 * pph_31[k]
                   + f_4 * pc_z[k] * pdh_94[k];

        t_129[k] = f_0 * sdh_101[k]
                   + f_4 * pc_x[k] * pdh_101[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_x, pc_x, pc_y, sdi0_133, sdh_102, \
                         sdh_104, sdi1_133, pph_56, pdh_98, pdh_102, \
                         pdh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_0 * sdh_102[k]
                   + f_4 * pc_x[k] * pdh_102[k];

        t_131[k] = f_0 * pph_56[k]
                   + f_4 * pc_y[k] * pdh_98[k];

        t_132[k] = f_0 * sdh_104[k]
                   + f_4 * pc_x[k] * pdh_104[k];

        t_133[k] = pa_x[k] * sdi0_133[k]
                   - f_11 * pc_x[k] * sdi1_133[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pc_x, pc_z, sdi0_135, sdi0_136, \
                         sdi0_137, sdi1_135, sdi1_136, sdi1_137, pph_36, \
                         pdh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_0 * pph_36[k]
                   + f_4 * pc_z[k] * pdh_99[k];

        t_135[k] = pa_x[k] * sdi0_135[k]
                   - f_11 * pc_x[k] * sdi1_135[k];

        t_136[k] = pa_x[k] * sdi0_136[k]
                   - f_11 * pc_x[k] * sdi1_136[k];

        t_137[k] = pa_x[k] * sdi0_137[k]
                   - f_11 * pc_x[k] * sdi1_137[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_x, pc_x, pc_y, sdi0_139, sdi0_140, \
                         sdh_105, sdi1_139, sdi1_140, pph_62, pdh_104, \
                         pdh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * pph_62[k]
                   + f_4 * pc_y[k] * pdh_104[k];

        t_139[k] = pa_x[k] * sdi0_139[k]
                   - f_11 * pc_x[k] * sdi1_139[k];

        t_140[k] = pa_x[k] * sdi0_140[k]
                   + f_12 * sdh_105[k]
                   - f_11 * pc_x[k] * sdi1_140[k];

        t_141[k] = f_4 * pc_y[k] * pdh_105[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pc_y, pc_z, pph_42, pdg0_75, pdg1_75, pdh_105, \
                         pdh_106, pdh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_1 * pph_42[k]
                   + f_4 * pc_z[k] * pdh_105[k];

        t_143[k] = f_5 * pdg0_75[k]
                   - f_6 * pdg1_75[k]
                   + f_4 * pc_y[k] * pdh_106[k];

        t_144[k] = f_4 * pc_y[k] * pdh_107[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_x, pc_x, pc_y, pc_z, sdi0_145, \
                         sdh_110, sdi1_145, pph_45, pdg0_76, pdg1_76, pdh_108, \
                         pdh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_x[k] * sdi0_145[k]
                   + f_13 * sdh_110[k]
                   - f_11 * pc_x[k] * sdi1_145[k];

        t_146[k] = f_7 * pdg0_76[k]
                   - f_8 * pdg1_76[k]
                   + f_4 * pc_y[k] * pdh_108[k];

        t_147[k] = f_1 * pph_45[k]
                   + f_4 * pc_z[k] * pdh_108[k];

        t_148[k] = f_4 * pc_y[k] * pdh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_x, pc_x, pc_y, pc_z, sdi0_149, sdh_114, \
                         sdi1_149, pph_48, pdg0_78, pdg1_78, pdh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_x[k] * sdi0_149[k]
                   + f_14 * sdh_114[k]
                   - f_11 * pc_x[k] * sdi1_149[k];

        t_150[k] = f_9 * pdg0_78[k]
                   - f_10 * pdg1_78[k]
                   + f_4 * pc_y[k] * pdh_111[k];

        t_151[k] = f_1 * pph_48[k]
                   + f_4 * pc_z[k] * pdh_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_x, pc_x, pc_y, sdi0_154, sdh_119, \
                         sdh_120, sdi1_154, pdg0_80, pdg1_80, pdh_113, pdh_114, \
                         pdh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * pdg0_80[k]
                   - f_6 * pdg1_80[k]
                   + f_4 * pc_y[k] * pdh_113[k];

        t_153[k] = f_4 * pc_y[k] * pdh_114[k];

        t_154[k] = pa_x[k] * sdi0_154[k]
                   + f_1 * sdh_119[k]
                   - f_11 * pc_x[k] * sdi1_154[k];

        t_155[k] = f_0 * sdh_120[k]
                   + f_4 * pc_x[k] * pdh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, pc_y, pc_z, sdh_122, sdh_123, \
                         pph_52, pdh_115, pdh_119, pdh_122, pdh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * pph_52[k]
                   + f_4 * pc_z[k] * pdh_115[k];

        t_157[k] = f_0 * sdh_122[k]
                   + f_4 * pc_x[k] * pdh_122[k];

        t_158[k] = f_0 * sdh_123[k]
                   + f_4 * pc_x[k] * pdh_123[k];

        t_159[k] = f_4 * pc_y[k] * pdh_119[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pc_x, pc_z, sdi0_161, sdi0_163, \
                         sdh_125, sdi1_161, sdi1_163, pph_57, pdh_120, \
                         pdh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_0 * sdh_125[k]
                   + f_4 * pc_x[k] * pdh_125[k];

        t_161[k] = pa_x[k] * sdi0_161[k]
                   - f_11 * pc_x[k] * sdi1_161[k];

        t_162[k] = f_1 * pph_57[k]
                   + f_4 * pc_z[k] * pdh_120[k];

        t_163[k] = pa_x[k] * sdi0_163[k]
                   - f_11 * pc_x[k] * sdi1_163[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_x, pc_x, pc_y, sdi0_164, sdi0_165, \
                         sdi0_167, sdi1_164, sdi1_165, sdi1_167, \
                         pdh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_x[k] * sdi0_164[k]
                   - f_11 * pc_x[k] * sdi1_164[k];

        t_165[k] = pa_x[k] * sdi0_165[k]
                   - f_11 * pc_x[k] * sdi1_165[k];

        t_166[k] = f_4 * pc_y[k] * pdh_125[k];

        t_167[k] = pa_x[k] * sdi0_167[k]
                   - f_11 * pc_x[k] * sdi1_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pc_y, pc_z, sdi0_0, sdi0_1, sdi0_3, \
                         sdh_0, sdh_1, sdi1_0, sdi1_1, sdi1_3, \
                         pdh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * sdi0_0[k]
                   - f_11 * pc_y[k] * sdi1_0[k];

        t_169[k] = pa_y[k] * sdi0_1[k]
                   + f_0 * sdh_0[k]
                   - f_11 * pc_y[k] * sdi1_1[k];

        t_170[k] = f_4 * pc_z[k] * pdh_126[k];

        t_171[k] = pa_y[k] * sdi0_3[k]
                   + f_1 * sdh_1[k]
                   - f_11 * pc_y[k] * sdi1_3[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pc_y, pc_z, sdi0_5, sdi0_6, sdh_3, \
                         sdi1_5, sdi1_6, pdh_127, pdh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_4 * pc_z[k] * pdh_127[k];

        t_173[k] = pa_y[k] * sdi0_5[k]
                   - f_11 * pc_y[k] * sdi1_5[k];

        t_174[k] = pa_y[k] * sdi0_6[k]
                   + f_14 * sdh_3[k]
                   - f_11 * pc_y[k] * sdi1_6[k];

        t_175[k] = f_4 * pc_z[k] * pdh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pc_y, pc_z, sdi0_8, sdi0_9, \
                         sdi0_10, sdh_5, sdh_6, sdi1_8, sdi1_9, sdi1_10, \
                         pdh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_y[k] * sdi0_8[k]
                   + f_0 * sdh_5[k]
                   - f_11 * pc_y[k] * sdi1_8[k];

        t_177[k] = pa_y[k] * sdi0_9[k]
                   - f_11 * pc_y[k] * sdi1_9[k];

        t_178[k] = pa_y[k] * sdi0_10[k]
                   + f_13 * sdh_6[k]
                   - f_11 * pc_y[k] * sdi1_10[k];

        t_179[k] = f_4 * pc_z[k] * pdh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_y, pc_y, sdi0_12, sdi0_13, sdi0_14, sdh_8, \
                         sdh_9, sdi1_12, sdi1_13, sdi1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * sdi0_12[k]
                   + f_1 * sdh_8[k]
                   - f_11 * pc_y[k] * sdi1_12[k];

        t_181[k] = pa_y[k] * sdi0_13[k]
                   + f_0 * sdh_9[k]
                   - f_11 * pc_y[k] * sdi1_13[k];

        t_182[k] = pa_y[k] * sdi0_14[k]
                   - f_11 * pc_y[k] * sdi1_14[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pc_x, pph_78, pph_79, pph_80, \
                         pph_81, pph_82, pdh_141, pdh_142, pdh_143, pdh_144, \
                         pdh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_1 * pph_78[k]
                   + f_4 * pc_x[k] * pdh_141[k];

        t_184[k] = f_1 * pph_79[k]
                   + f_4 * pc_x[k] * pdh_142[k];

        t_185[k] = f_1 * pph_80[k]
                   + f_4 * pc_x[k] * pdh_143[k];

        t_186[k] = f_1 * pph_81[k]
                   + f_4 * pc_x[k] * pdh_144[k];

        t_187[k] = f_1 * pph_82[k]
                   + f_4 * pc_x[k] * pdh_145[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_y, pc_x, pc_y, pc_z, sdi0_21, sdh_15, \
                         sdi1_21, pph_83, pdh_141, pdh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_1 * pph_83[k]
                   + f_4 * pc_x[k] * pdh_146[k];

        t_189[k] = pa_y[k] * sdi0_21[k]
                   + f_12 * sdh_15[k]
                   - f_11 * pc_y[k] * sdi1_21[k];

        t_190[k] = f_4 * pc_z[k] * pdh_141[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_z, pdg0_100, pdg0_101, pdg0_102, pdg1_100, \
                         pdg1_101, pdg1_102, pdh_142, pdh_143, \
                         pdh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_5 * pdg0_100[k]
                   - f_6 * pdg1_100[k]
                   + f_4 * pc_z[k] * pdh_142[k];

        t_192[k] = f_7 * pdg0_101[k]
                   - f_8 * pdg1_101[k]
                   + f_4 * pc_z[k] * pdh_143[k];

        t_193[k] = f_9 * pdg0_102[k]
                   - f_10 * pdg1_102[k]
                   + f_4 * pc_z[k] * pdh_144[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_y, pc_x, pc_y, sdi0_27, sdh_20, sdi1_27, \
                         pph_84, pdg0_105, pdg1_105, pdh_146, pdh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_0 * sdh_20[k]
                   + f_4 * pc_y[k] * pdh_146[k];

        t_195[k] = pa_y[k] * sdi0_27[k]
                   - f_11 * pc_y[k] * sdi1_27[k];

        t_196[k] = f_0 * pph_84[k]
                   + f_2 * pdg0_105[k]
                   - f_3 * pdg1_105[k]
                   + f_4 * pc_x[k] * pdh_147[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_x, pc_x, pc_z, ppi0_113, ppi0_115, \
                         pph_85, pph_87, ppi1_113, ppi1_115, pdh_147, \
                         pdh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pb_x[k] * ppi0_113[k]
                   + f_15 * pph_85[k]
                   - f_11 * pc_x[k] * ppi1_113[k];

        t_198[k] = f_4 * pc_z[k] * pdh_147[k];

        t_199[k] = pb_x[k] * ppi0_115[k]
                   + f_13 * pph_87[k]
                   - f_11 * pc_x[k] * ppi1_115[k];

        t_200[k] = f_4 * pc_z[k] * pdh_148[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_x, pc_x, pc_z, ppi0_118, pph_89, pph_90, \
                         ppi1_118, pdg0_110, pdg1_110, pdh_150, \
                         pdh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_0 * pph_89[k]
                   + f_9 * pdg0_110[k]
                   - f_10 * pdg1_110[k]
                   + f_4 * pc_x[k] * pdh_152[k];

        t_202[k] = pb_x[k] * ppi0_118[k]
                   + f_14 * pph_90[k]
                   - f_11 * pc_x[k] * ppi1_118[k];

        t_203[k] = f_4 * pc_z[k] * pdh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pb_x, pc_x, ppi0_120, ppi0_122, pph_92, pph_93, \
                         pph_94, ppi1_120, ppi1_122, pdg0_114, pdg1_114, \
                         pdh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_x[k] * ppi0_120[k]
                   + f_14 * pph_92[k]
                   - f_11 * pc_x[k] * ppi1_120[k];

        t_205[k] = f_0 * pph_93[k]
                   + f_7 * pdg0_114[k]
                   - f_8 * pdg1_114[k]
                   + f_4 * pc_x[k] * pdh_156[k];

        t_206[k] = pb_x[k] * ppi0_122[k]
                   + f_1 * pph_94[k]
                   - f_11 * pc_x[k] * ppi1_122[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pb_x, pc_x, pc_z, ppi0_124, ppi0_125, pph_96, \
                         pph_97, ppi1_124, ppi1_125, pdh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_4 * pc_z[k] * pdh_153[k];

        t_208[k] = pb_x[k] * ppi0_124[k]
                   + f_1 * pph_96[k]
                   - f_11 * pc_x[k] * ppi1_124[k];

        t_209[k] = pb_x[k] * ppi0_125[k]
                   + f_1 * pph_97[k]
                   - f_11 * pc_x[k] * ppi1_125[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pph_98, pph_99, pph_100, pph_101, \
                         pdg0_119, pdg1_119, pdh_161, pdh_162, pdh_163, \
                         pdh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * pph_98[k]
                   + f_5 * pdg0_119[k]
                   - f_6 * pdg1_119[k]
                   + f_4 * pc_x[k] * pdh_161[k];

        t_211[k] = f_0 * pph_99[k]
                   + f_4 * pc_x[k] * pdh_162[k];

        t_212[k] = f_0 * pph_100[k]
                   + f_4 * pc_x[k] * pdh_163[k];

        t_213[k] = f_0 * pph_101[k]
                   + f_4 * pc_x[k] * pdh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pb_x, pc_x, ppi0_133, pph_102, pph_103, \
                         pph_104, ppi1_133, pdh_165, pdh_166, pdh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_0 * pph_102[k]
                   + f_4 * pc_x[k] * pdh_165[k];

        t_215[k] = f_0 * pph_103[k]
                   + f_4 * pc_x[k] * pdh_166[k];

        t_216[k] = f_0 * pph_104[k]
                   + f_4 * pc_x[k] * pdh_167[k];

        t_217[k] = pb_x[k] * ppi0_133[k]
                   - f_11 * pc_x[k] * ppi1_133[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pb_x, pc_x, pc_z, ppi0_135, ppi0_136, \
                         ppi0_137, ppi1_135, ppi1_136, ppi1_137, \
                         pdh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_4 * pc_z[k] * pdh_162[k];

        t_219[k] = pb_x[k] * ppi0_135[k]
                   - f_11 * pc_x[k] * ppi1_135[k];

        t_220[k] = pb_x[k] * ppi0_136[k]
                   - f_11 * pc_x[k] * ppi1_136[k];

        t_221[k] = pb_x[k] * ppi0_137[k]
                   - f_11 * pc_x[k] * ppi1_137[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_y, pb_x, pc_x, pc_y, sdi0_56, sdi1_56, \
                         ppi0_138, ppi0_139, ppi1_138, ppi1_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = pb_x[k] * ppi0_138[k]
                   - f_11 * pc_x[k] * ppi1_138[k];

        t_223[k] = pb_x[k] * ppi0_139[k]
                   - f_11 * pc_x[k] * ppi1_139[k];

        t_224[k] = pa_y[k] * sdi0_56[k]
                   - f_11 * pc_y[k] * sdi1_56[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_z, pc_z, ppi0_85, ppi0_87, pph_63, \
                         pph_64, ppi1_85, ppi1_87, pdh_168, pdh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pb_z[k] * ppi0_85[k]
                   - f_11 * pc_z[k] * ppi1_85[k];

        t_226[k] = f_0 * pph_63[k]
                   + f_4 * pc_z[k] * pdh_168[k];

        t_227[k] = pb_z[k] * ppi0_87[k]
                   - f_11 * pc_z[k] * ppi1_87[k];

        t_228[k] = f_0 * pph_64[k]
                   + f_4 * pc_z[k] * pdh_169[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pa_y, pb_z, pc_y, pc_z, sdi0_61, sdi1_61, \
                         ppi0_90, pph_66, ppi1_90, pdh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * sdi0_61[k]
                   - f_11 * pc_y[k] * sdi1_61[k];

        t_230[k] = pb_z[k] * ppi0_90[k]
                   - f_11 * pc_z[k] * ppi1_90[k];

        t_231[k] = f_0 * pph_66[k]
                   + f_4 * pc_z[k] * pdh_171[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pa_y, pb_z, pc_y, pc_z, sdi0_64, sdi0_65, \
                         sdh_47, sdi1_64, sdi1_65, ppi0_94, ppi1_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pa_y[k] * sdi0_64[k]
                   + f_0 * sdh_47[k]
                   - f_11 * pc_y[k] * sdi1_64[k];

        t_233[k] = pa_y[k] * sdi0_65[k]
                   - f_11 * pc_y[k] * sdi1_65[k];

        t_234[k] = pb_z[k] * ppi0_94[k]
                   - f_11 * pc_z[k] * ppi1_94[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_y, pc_y, pc_z, sdi0_68, sdi0_69, sdh_50, \
                         sdh_51, sdi1_68, sdi1_69, pph_69, pdh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_0 * pph_69[k]
                   + f_4 * pc_z[k] * pdh_174[k];

        t_236[k] = pa_y[k] * sdi0_68[k]
                   + f_1 * sdh_50[k]
                   - f_11 * pc_y[k] * sdi1_68[k];

        t_237[k] = pa_y[k] * sdi0_69[k]
                   + f_0 * sdh_51[k]
                   - f_11 * pc_y[k] * sdi1_69[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_y, pc_x, pc_y, sdi0_70, sdi1_70, \
                         pph_120, pph_121, pph_122, pdh_183, pdh_184, \
                         pdh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_y[k] * sdi0_70[k]
                   - f_11 * pc_y[k] * sdi1_70[k];

        t_239[k] = f_0 * pph_120[k]
                   + f_4 * pc_x[k] * pdh_183[k];

        t_240[k] = f_0 * pph_121[k]
                   + f_4 * pc_x[k] * pdh_184[k];

        t_241[k] = f_0 * pph_122[k]
                   + f_4 * pc_x[k] * pdh_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pb_x, pc_x, ppi0_161, pph_123, pph_124, \
                         pph_125, ppi1_161, pdh_186, pdh_187, pdh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * pph_123[k]
                   + f_4 * pc_x[k] * pdh_186[k];

        t_243[k] = f_0 * pph_124[k]
                   + f_4 * pc_x[k] * pdh_187[k];

        t_244[k] = f_0 * pph_125[k]
                   + f_4 * pc_x[k] * pdh_188[k];

        t_245[k] = pb_x[k] * ppi0_161[k]
                   - f_11 * pc_x[k] * ppi1_161[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pb_x, pc_x, pc_z, ppi0_163, ppi0_164, \
                         ppi0_165, pph_78, ppi1_163, ppi1_164, ppi1_165, \
                         pdh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_0 * pph_78[k]
                   + f_4 * pc_z[k] * pdh_183[k];

        t_247[k] = pb_x[k] * ppi0_163[k]
                   - f_11 * pc_x[k] * ppi1_163[k];

        t_248[k] = pb_x[k] * ppi0_164[k]
                   - f_11 * pc_x[k] * ppi1_164[k];

        t_249[k] = pb_x[k] * ppi0_165[k]
                   - f_11 * pc_x[k] * ppi1_165[k];
    }
}

static auto
compute_prim_pdi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdi0, const size_t sdh,
                                                          const size_t sdi1, const size_t ppi0,
                                                          const size_t pph, const size_t ppi1,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdi0_0 = buffer.data(sdi0 + 0);
    const auto *sdi0_2 = buffer.data(sdi0 + 2);
    const auto *sdi0_3 = buffer.data(sdi0 + 3);
    const auto *sdi0_5 = buffer.data(sdi0 + 5);
    const auto *sdi0_6 = buffer.data(sdi0 + 6);
    const auto *sdi0_7 = buffer.data(sdi0 + 7);
    const auto *sdi0_9 = buffer.data(sdi0 + 9);
    const auto *sdi0_10 = buffer.data(sdi0 + 10);
    const auto *sdi0_11 = buffer.data(sdi0 + 11);
    const auto *sdi0_12 = buffer.data(sdi0 + 12);
    const auto *sdi0_14 = buffer.data(sdi0 + 14);
    const auto *sdi0_21 = buffer.data(sdi0 + 21);
    const auto *sdi0_27 = buffer.data(sdi0 + 27);
    const auto *sdi0_28 = buffer.data(sdi0 + 28);
    const auto *sdi0_31 = buffer.data(sdi0 + 31);
    const auto *sdi0_34 = buffer.data(sdi0 + 34);
    const auto *sdi0_35 = buffer.data(sdi0 + 35);
    const auto *sdi0_83 = buffer.data(sdi0 + 83);
    const auto *sdi0_140 = buffer.data(sdi0 + 140);
    const auto *sdi0_145 = buffer.data(sdi0 + 145);
    const auto *sdi0_149 = buffer.data(sdi0 + 149);
    const auto *sdi0_154 = buffer.data(sdi0 + 154);
    const auto *sdi0_161 = buffer.data(sdi0 + 161);
    const auto *sdi0_163 = buffer.data(sdi0 + 163);
    const auto *sdi0_164 = buffer.data(sdi0 + 164);
    const auto *sdi0_165 = buffer.data(sdi0 + 165);
    const auto *sdi0_167 = buffer.data(sdi0 + 167);

    const auto *sdh_0 = buffer.data(sdh + 0);
    const auto *sdh_2 = buffer.data(sdh + 2);
    const auto *sdh_3 = buffer.data(sdh + 3);
    const auto *sdh_5 = buffer.data(sdh + 5);
    const auto *sdh_6 = buffer.data(sdh + 6);
    const auto *sdh_7 = buffer.data(sdh + 7);
    const auto *sdh_9 = buffer.data(sdh + 9);
    const auto *sdh_20 = buffer.data(sdh + 20);
    const auto *sdh_24 = buffer.data(sdh + 24);
    const auto *sdh_62 = buffer.data(sdh + 62);
    const auto *sdh_78 = buffer.data(sdh + 78);
    const auto *sdh_83 = buffer.data(sdh + 83);
    const auto *sdh_104 = buffer.data(sdh + 104);
    const auto *sdh_120 = buffer.data(sdh + 120);
    const auto *sdh_122 = buffer.data(sdh + 122);
    const auto *sdh_123 = buffer.data(sdh + 123);
    const auto *sdh_124 = buffer.data(sdh + 124);
    const auto *sdh_125 = buffer.data(sdh + 125);

    const auto *sdi1_0 = buffer.data(sdi1 + 0);
    const auto *sdi1_2 = buffer.data(sdi1 + 2);
    const auto *sdi1_3 = buffer.data(sdi1 + 3);
    const auto *sdi1_5 = buffer.data(sdi1 + 5);
    const auto *sdi1_6 = buffer.data(sdi1 + 6);
    const auto *sdi1_7 = buffer.data(sdi1 + 7);
    const auto *sdi1_9 = buffer.data(sdi1 + 9);
    const auto *sdi1_10 = buffer.data(sdi1 + 10);
    const auto *sdi1_11 = buffer.data(sdi1 + 11);
    const auto *sdi1_12 = buffer.data(sdi1 + 12);
    const auto *sdi1_14 = buffer.data(sdi1 + 14);
    const auto *sdi1_21 = buffer.data(sdi1 + 21);
    const auto *sdi1_27 = buffer.data(sdi1 + 27);
    const auto *sdi1_28 = buffer.data(sdi1 + 28);
    const auto *sdi1_31 = buffer.data(sdi1 + 31);
    const auto *sdi1_34 = buffer.data(sdi1 + 34);
    const auto *sdi1_35 = buffer.data(sdi1 + 35);
    const auto *sdi1_83 = buffer.data(sdi1 + 83);
    const auto *sdi1_140 = buffer.data(sdi1 + 140);
    const auto *sdi1_145 = buffer.data(sdi1 + 145);
    const auto *sdi1_149 = buffer.data(sdi1 + 149);
    const auto *sdi1_154 = buffer.data(sdi1 + 154);
    const auto *sdi1_161 = buffer.data(sdi1 + 161);
    const auto *sdi1_163 = buffer.data(sdi1 + 163);
    const auto *sdi1_164 = buffer.data(sdi1 + 164);
    const auto *sdi1_165 = buffer.data(sdi1 + 165);
    const auto *sdi1_167 = buffer.data(sdi1 + 167);

    const auto *ppi0_113 = buffer.data(ppi0 + 113);
    const auto *ppi0_115 = buffer.data(ppi0 + 115);
    const auto *ppi0_118 = buffer.data(ppi0 + 118);
    const auto *ppi0_122 = buffer.data(ppi0 + 122);
    const auto *ppi0_133 = buffer.data(ppi0 + 133);
    const auto *ppi0_170 = buffer.data(ppi0 + 170);
    const auto *ppi0_173 = buffer.data(ppi0 + 173);

    const auto *pph_84 = buffer.data(pph + 84);
    const auto *pph_85 = buffer.data(pph + 85);
    const auto *pph_87 = buffer.data(pph + 87);
    const auto *pph_90 = buffer.data(pph + 90);
    const auto *pph_99 = buffer.data(pph + 99);
    const auto *pph_100 = buffer.data(pph + 100);
    const auto *pph_101 = buffer.data(pph + 101);
    const auto *pph_102 = buffer.data(pph + 102);
    const auto *pph_104 = buffer.data(pph + 104);
    const auto *pph_105 = buffer.data(pph + 105);
    const auto *pph_106 = buffer.data(pph + 106);
    const auto *pph_108 = buffer.data(pph + 108);
    const auto *pph_111 = buffer.data(pph + 111);
    const auto *pph_120 = buffer.data(pph + 120);
    const auto *pph_125 = buffer.data(pph + 125);
    const auto *pph_126 = buffer.data(pph + 126);
    const auto *pph_128 = buffer.data(pph + 128);
    const auto *pph_141 = buffer.data(pph + 141);
    const auto *pph_142 = buffer.data(pph + 142);
    const auto *pph_143 = buffer.data(pph + 143);
    const auto *pph_144 = buffer.data(pph + 144);
    const auto *pph_145 = buffer.data(pph + 145);
    const auto *pph_146 = buffer.data(pph + 146);

    const auto *ppi1_113 = buffer.data(ppi1 + 113);
    const auto *ppi1_115 = buffer.data(ppi1 + 115);
    const auto *ppi1_118 = buffer.data(ppi1 + 118);
    const auto *ppi1_122 = buffer.data(ppi1 + 122);
    const auto *ppi1_133 = buffer.data(ppi1 + 133);
    const auto *ppi1_170 = buffer.data(ppi1 + 170);
    const auto *ppi1_173 = buffer.data(ppi1 + 173);

    const auto *pdg0_135 = buffer.data(pdg0 + 135);
    const auto *pdg0_136 = buffer.data(pdg0 + 136);
    const auto *pdg0_138 = buffer.data(pdg0 + 138);
    const auto *pdg0_140 = buffer.data(pdg0 + 140);
    const auto *pdg0_141 = buffer.data(pdg0 + 141);
    const auto *pdg0_143 = buffer.data(pdg0 + 143);
    const auto *pdg0_144 = buffer.data(pdg0 + 144);
    const auto *pdg0_145 = buffer.data(pdg0 + 145);
    const auto *pdg0_146 = buffer.data(pdg0 + 146);
    const auto *pdg0_147 = buffer.data(pdg0 + 147);
    const auto *pdg0_148 = buffer.data(pdg0 + 148);
    const auto *pdg0_149 = buffer.data(pdg0 + 149);
    const auto *pdg0_150 = buffer.data(pdg0 + 150);
    const auto *pdg0_155 = buffer.data(pdg0 + 155);
    const auto *pdg0_158 = buffer.data(pdg0 + 158);
    const auto *pdg0_159 = buffer.data(pdg0 + 159);
    const auto *pdg0_160 = buffer.data(pdg0 + 160);
    const auto *pdg0_161 = buffer.data(pdg0 + 161);
    const auto *pdg0_162 = buffer.data(pdg0 + 162);
    const auto *pdg0_163 = buffer.data(pdg0 + 163);
    const auto *pdg0_164 = buffer.data(pdg0 + 164);
    const auto *pdg0_166 = buffer.data(pdg0 + 166);
    const auto *pdg0_168 = buffer.data(pdg0 + 168);
    const auto *pdg0_171 = buffer.data(pdg0 + 171);
    const auto *pdg0_173 = buffer.data(pdg0 + 173);
    const auto *pdg0_175 = buffer.data(pdg0 + 175);
    const auto *pdg0_177 = buffer.data(pdg0 + 177);
    const auto *pdg0_178 = buffer.data(pdg0 + 178);
    const auto *pdg0_191 = buffer.data(pdg0 + 191);
    const auto *pdg0_192 = buffer.data(pdg0 + 192);
    const auto *pdg0_193 = buffer.data(pdg0 + 193);
    const auto *pdg0_194 = buffer.data(pdg0 + 194);

    const auto *pdg1_135 = buffer.data(pdg1 + 135);
    const auto *pdg1_136 = buffer.data(pdg1 + 136);
    const auto *pdg1_138 = buffer.data(pdg1 + 138);
    const auto *pdg1_140 = buffer.data(pdg1 + 140);
    const auto *pdg1_141 = buffer.data(pdg1 + 141);
    const auto *pdg1_143 = buffer.data(pdg1 + 143);
    const auto *pdg1_144 = buffer.data(pdg1 + 144);
    const auto *pdg1_145 = buffer.data(pdg1 + 145);
    const auto *pdg1_146 = buffer.data(pdg1 + 146);
    const auto *pdg1_147 = buffer.data(pdg1 + 147);
    const auto *pdg1_148 = buffer.data(pdg1 + 148);
    const auto *pdg1_149 = buffer.data(pdg1 + 149);
    const auto *pdg1_150 = buffer.data(pdg1 + 150);
    const auto *pdg1_155 = buffer.data(pdg1 + 155);
    const auto *pdg1_158 = buffer.data(pdg1 + 158);
    const auto *pdg1_159 = buffer.data(pdg1 + 159);
    const auto *pdg1_160 = buffer.data(pdg1 + 160);
    const auto *pdg1_161 = buffer.data(pdg1 + 161);
    const auto *pdg1_162 = buffer.data(pdg1 + 162);
    const auto *pdg1_163 = buffer.data(pdg1 + 163);
    const auto *pdg1_164 = buffer.data(pdg1 + 164);
    const auto *pdg1_166 = buffer.data(pdg1 + 166);
    const auto *pdg1_168 = buffer.data(pdg1 + 168);
    const auto *pdg1_171 = buffer.data(pdg1 + 171);
    const auto *pdg1_173 = buffer.data(pdg1 + 173);
    const auto *pdg1_175 = buffer.data(pdg1 + 175);
    const auto *pdg1_177 = buffer.data(pdg1 + 177);
    const auto *pdg1_178 = buffer.data(pdg1 + 178);
    const auto *pdg1_191 = buffer.data(pdg1 + 191);
    const auto *pdg1_192 = buffer.data(pdg1 + 192);
    const auto *pdg1_193 = buffer.data(pdg1 + 193);
    const auto *pdg1_194 = buffer.data(pdg1 + 194);

    const auto *pdh_188 = buffer.data(pdh + 188);
    const auto *pdh_189 = buffer.data(pdh + 189);
    const auto *pdh_190 = buffer.data(pdh + 190);
    const auto *pdh_192 = buffer.data(pdh + 192);
    const auto *pdh_194 = buffer.data(pdh + 194);
    const auto *pdh_195 = buffer.data(pdh + 195);
    const auto *pdh_197 = buffer.data(pdh + 197);
    const auto *pdh_198 = buffer.data(pdh + 198);
    const auto *pdh_199 = buffer.data(pdh + 199);
    const auto *pdh_201 = buffer.data(pdh + 201);
    const auto *pdh_202 = buffer.data(pdh + 202);
    const auto *pdh_203 = buffer.data(pdh + 203);
    const auto *pdh_204 = buffer.data(pdh + 204);
    const auto *pdh_205 = buffer.data(pdh + 205);
    const auto *pdh_206 = buffer.data(pdh + 206);
    const auto *pdh_207 = buffer.data(pdh + 207);
    const auto *pdh_208 = buffer.data(pdh + 208);
    const auto *pdh_209 = buffer.data(pdh + 209);
    const auto *pdh_210 = buffer.data(pdh + 210);
    const auto *pdh_211 = buffer.data(pdh + 211);
    const auto *pdh_213 = buffer.data(pdh + 213);
    const auto *pdh_215 = buffer.data(pdh + 215);
    const auto *pdh_216 = buffer.data(pdh + 216);
    const auto *pdh_218 = buffer.data(pdh + 218);
    const auto *pdh_219 = buffer.data(pdh + 219);
    const auto *pdh_222 = buffer.data(pdh + 222);
    const auto *pdh_223 = buffer.data(pdh + 223);
    const auto *pdh_224 = buffer.data(pdh + 224);
    const auto *pdh_225 = buffer.data(pdh + 225);
    const auto *pdh_226 = buffer.data(pdh + 226);
    const auto *pdh_227 = buffer.data(pdh + 227);
    const auto *pdh_228 = buffer.data(pdh + 228);
    const auto *pdh_229 = buffer.data(pdh + 229);
    const auto *pdh_230 = buffer.data(pdh + 230);
    const auto *pdh_231 = buffer.data(pdh + 231);
    const auto *pdh_232 = buffer.data(pdh + 232);
    const auto *pdh_234 = buffer.data(pdh + 234);
    const auto *pdh_237 = buffer.data(pdh + 237);
    const auto *pdh_239 = buffer.data(pdh + 239);
    const auto *pdh_241 = buffer.data(pdh + 241);
    const auto *pdh_243 = buffer.data(pdh + 243);
    const auto *pdh_244 = buffer.data(pdh + 244);
    const auto *pdh_246 = buffer.data(pdh + 246);
    const auto *pdh_247 = buffer.data(pdh + 247);
    const auto *pdh_248 = buffer.data(pdh + 248);
    const auto *pdh_249 = buffer.data(pdh + 249);
    const auto *pdh_250 = buffer.data(pdh + 250);
    const auto *pdh_251 = buffer.data(pdh + 251);
    const auto *pdh_252 = buffer.data(pdh + 252);
    const auto *pdh_254 = buffer.data(pdh + 254);
    const auto *pdh_257 = buffer.data(pdh + 257);
    const auto *pdh_261 = buffer.data(pdh + 261);
    const auto *pdh_267 = buffer.data(pdh + 267);
    const auto *pdh_268 = buffer.data(pdh + 268);
    const auto *pdh_269 = buffer.data(pdh + 269);
    const auto *pdh_270 = buffer.data(pdh + 270);
    const auto *pdh_271 = buffer.data(pdh + 271);
    const auto *pdh_272 = buffer.data(pdh + 272);
    const auto *pdh_273 = buffer.data(pdh + 273);
    const auto *pdh_275 = buffer.data(pdh + 275);

#pragma omp simd aligned(t_250, t_251, t_252, pa_y, pc_x, pc_y, sdi0_83, sdh_62, sdi1_83, \
                         pdg0_135, pdg1_135, pdh_188, pdh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * sdh_62[k]
                   + f_4 * pc_y[k] * pdh_188[k];

        t_251[k] = pa_y[k] * sdi0_83[k]
                   - f_11 * pc_y[k] * sdi1_83[k];

        t_252[k] = f_2 * pdg0_135[k]
                   - f_3 * pdg1_135[k]
                   + f_4 * pc_x[k] * pdh_189[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_z, pdg0_136, pdg0_138, pdg1_136, \
                         pdg1_138, pdh_189, pdh_190, pdh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_16 * pdg0_136[k]
                   - f_17 * pdg1_136[k]
                   + f_4 * pc_x[k] * pdh_190[k];

        t_254[k] = f_4 * pc_z[k] * pdh_189[k];

        t_255[k] = f_9 * pdg0_138[k]
                   - f_10 * pdg1_138[k]
                   + f_4 * pc_x[k] * pdh_192[k];

        t_256[k] = f_4 * pc_z[k] * pdh_190[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pc_x, pc_z, pdg0_140, pdg0_141, pdg0_143, \
                         pdg1_140, pdg1_141, pdg1_143, pdh_192, pdh_194, pdh_195, \
                         pdh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_9 * pdg0_140[k]
                   - f_10 * pdg1_140[k]
                   + f_4 * pc_x[k] * pdh_194[k];

        t_258[k] = f_7 * pdg0_141[k]
                   - f_8 * pdg1_141[k]
                   + f_4 * pc_x[k] * pdh_195[k];

        t_259[k] = f_4 * pc_z[k] * pdh_192[k];

        t_260[k] = f_7 * pdg0_143[k]
                   - f_8 * pdg1_143[k]
                   + f_4 * pc_x[k] * pdh_197[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pc_x, pc_z, pdg0_144, pdg0_145, pdg0_147, \
                         pdg1_144, pdg1_145, pdg1_147, pdh_195, pdh_198, pdh_199, \
                         pdh_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_7 * pdg0_144[k]
                   - f_8 * pdg1_144[k]
                   + f_4 * pc_x[k] * pdh_198[k];

        t_262[k] = f_5 * pdg0_145[k]
                   - f_6 * pdg1_145[k]
                   + f_4 * pc_x[k] * pdh_199[k];

        t_263[k] = f_4 * pc_z[k] * pdh_195[k];

        t_264[k] = f_5 * pdg0_147[k]
                   - f_6 * pdg1_147[k]
                   + f_4 * pc_x[k] * pdh_201[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pc_x, pdg0_148, pdg0_149, \
                         pdg1_148, pdg1_149, pdh_202, pdh_203, pdh_204, pdh_205, \
                         pdh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_5 * pdg0_148[k]
                   - f_6 * pdg1_148[k]
                   + f_4 * pc_x[k] * pdh_202[k];

        t_266[k] = f_5 * pdg0_149[k]
                   - f_6 * pdg1_149[k]
                   + f_4 * pc_x[k] * pdh_203[k];

        t_267[k] = f_4 * pc_x[k] * pdh_204[k];

        t_268[k] = f_4 * pc_x[k] * pdh_205[k];

        t_269[k] = f_4 * pc_x[k] * pdh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, sdh_78, pph_99, \
                         pdg0_145, pdg1_145, pdh_204, pdh_207, pdh_208, \
                         pdh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_4 * pc_x[k] * pdh_207[k];

        t_271[k] = f_4 * pc_x[k] * pdh_208[k];

        t_272[k] = f_4 * pc_x[k] * pdh_209[k];

        t_273[k] = f_0 * sdh_78[k]
                   + f_1 * pph_99[k]
                   + f_2 * pdg0_145[k]
                   - f_3 * pdg1_145[k]
                   + f_4 * pc_y[k] * pdh_204[k];

        t_274[k] = f_4 * pc_z[k] * pdh_204[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_z, pdg0_145, pdg0_146, pdg0_147, pdg1_145, \
                         pdg1_146, pdg1_147, pdh_205, pdh_206, \
                         pdh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_5 * pdg0_145[k]
                   - f_6 * pdg1_145[k]
                   + f_4 * pc_z[k] * pdh_205[k];

        t_276[k] = f_7 * pdg0_146[k]
                   - f_8 * pdg1_146[k]
                   + f_4 * pc_z[k] * pdh_206[k];

        t_277[k] = f_9 * pdg0_147[k]
                   - f_10 * pdg1_147[k]
                   + f_4 * pc_z[k] * pdh_207[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, pc_z, sdh_83, pph_104, pdg0_149, \
                         pdg0_150, pdg1_149, pdg1_150, pdh_209, \
                         pdh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * sdh_83[k]
                   + f_1 * pph_104[k]
                   + f_4 * pc_y[k] * pdh_209[k];

        t_279[k] = f_2 * pdg0_149[k]
                   - f_3 * pdg1_149[k]
                   + f_4 * pc_z[k] * pdh_209[k];

        t_280[k] = f_2 * pdg0_150[k]
                   - f_3 * pdg1_150[k]
                   + f_4 * pc_x[k] * pdh_210[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pb_z, pc_z, ppi0_113, ppi0_115, pph_84, \
                         pph_85, ppi1_113, ppi1_115, pdh_210, pdh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = pb_z[k] * ppi0_113[k]
                   - f_11 * pc_z[k] * ppi1_113[k];

        t_282[k] = f_0 * pph_84[k]
                   + f_4 * pc_z[k] * pdh_210[k];

        t_283[k] = pb_z[k] * ppi0_115[k]
                   - f_11 * pc_z[k] * ppi1_115[k];

        t_284[k] = f_0 * pph_85[k]
                   + f_4 * pc_z[k] * pdh_211[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pb_z, pc_x, pc_z, ppi0_118, pph_87, ppi1_118, \
                         pdg0_155, pdg1_155, pdh_213, pdh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_9 * pdg0_155[k]
                   - f_10 * pdg1_155[k]
                   + f_4 * pc_x[k] * pdh_215[k];

        t_286[k] = pb_z[k] * ppi0_118[k]
                   - f_11 * pc_z[k] * ppi1_118[k];

        t_287[k] = f_0 * pph_87[k]
                   + f_4 * pc_z[k] * pdh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pb_z, pc_x, pc_z, ppi0_122, ppi1_122, pdg0_158, \
                         pdg0_159, pdg1_158, pdg1_159, pdh_218, \
                         pdh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_7 * pdg0_158[k]
                   - f_8 * pdg1_158[k]
                   + f_4 * pc_x[k] * pdh_218[k];

        t_289[k] = f_7 * pdg0_159[k]
                   - f_8 * pdg1_159[k]
                   + f_4 * pc_x[k] * pdh_219[k];

        t_290[k] = pb_z[k] * ppi0_122[k]
                   - f_11 * pc_z[k] * ppi1_122[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pc_x, pc_z, pph_90, pdg0_162, pdg0_163, \
                         pdg1_162, pdg1_163, pdh_216, pdh_222, \
                         pdh_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_0 * pph_90[k]
                   + f_4 * pc_z[k] * pdh_216[k];

        t_292[k] = f_5 * pdg0_162[k]
                   - f_6 * pdg1_162[k]
                   + f_4 * pc_x[k] * pdh_222[k];

        t_293[k] = f_5 * pdg0_163[k]
                   - f_6 * pdg1_163[k]
                   + f_4 * pc_x[k] * pdh_223[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, t_299, pc_x, pdg0_164, pdg1_164, \
                         pdh_224, pdh_225, pdh_226, pdh_227, pdh_228, \
                         pdh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_5 * pdg0_164[k]
                   - f_6 * pdg1_164[k]
                   + f_4 * pc_x[k] * pdh_224[k];

        t_295[k] = f_4 * pc_x[k] * pdh_225[k];

        t_296[k] = f_4 * pc_x[k] * pdh_226[k];

        t_297[k] = f_4 * pc_x[k] * pdh_227[k];

        t_298[k] = f_4 * pc_x[k] * pdh_228[k];

        t_299[k] = f_4 * pc_x[k] * pdh_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pb_z, pc_x, pc_z, ppi0_133, pph_99, \
                         pph_100, ppi1_133, pdg0_160, pdg1_160, pdh_225, pdh_226, \
                         pdh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_4 * pc_x[k] * pdh_230[k];

        t_301[k] = pb_z[k] * ppi0_133[k]
                   - f_11 * pc_z[k] * ppi1_133[k];

        t_302[k] = f_0 * pph_99[k]
                   + f_4 * pc_z[k] * pdh_225[k];

        t_303[k] = f_0 * pph_100[k]
                   + f_5 * pdg0_160[k]
                   - f_6 * pdg1_160[k]
                   + f_4 * pc_z[k] * pdh_226[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, pc_y, pc_z, sdh_104, pph_101, pph_102, pph_125, \
                         pdg0_161, pdg0_162, pdg1_161, pdg1_162, pdh_227, pdh_228, \
                         pdh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_0 * pph_101[k]
                   + f_7 * pdg0_161[k]
                   - f_8 * pdg1_161[k]
                   + f_4 * pc_z[k] * pdh_227[k];

        t_305[k] = f_0 * pph_102[k]
                   + f_9 * pdg0_162[k]
                   - f_10 * pdg1_162[k]
                   + f_4 * pc_z[k] * pdh_228[k];

        t_306[k] = f_0 * sdh_104[k]
                   + f_0 * pph_125[k]
                   + f_4 * pc_y[k] * pdh_230[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, pa_y, pc_x, pc_y, pc_z, sdi0_140, sdi1_140, \
                         pph_104, pdg0_164, pdg0_166, pdg1_164, pdg1_166, pdh_230, \
                         pdh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * pph_104[k]
                   + f_2 * pdg0_164[k]
                   - f_3 * pdg1_164[k]
                   + f_4 * pc_z[k] * pdh_230[k];

        t_308[k] = pa_y[k] * sdi0_140[k]
                   - f_11 * pc_y[k] * sdi1_140[k];

        t_309[k] = f_16 * pdg0_166[k]
                   - f_17 * pdg1_166[k]
                   + f_4 * pc_x[k] * pdh_232[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pc_x, pc_z, pph_105, pph_106, pdg0_168, \
                         pdg1_168, pdh_231, pdh_232, pdh_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_1 * pph_105[k]
                   + f_4 * pc_z[k] * pdh_231[k];

        t_311[k] = f_9 * pdg0_168[k]
                   - f_10 * pdg1_168[k]
                   + f_4 * pc_x[k] * pdh_234[k];

        t_312[k] = f_1 * pph_106[k]
                   + f_4 * pc_z[k] * pdh_232[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_y, pc_x, pc_y, pc_z, sdi0_145, sdi1_145, \
                         pph_108, pdg0_171, pdg1_171, pdh_234, \
                         pdh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pa_y[k] * sdi0_145[k]
                   - f_11 * pc_y[k] * sdi1_145[k];

        t_314[k] = f_7 * pdg0_171[k]
                   - f_8 * pdg1_171[k]
                   + f_4 * pc_x[k] * pdh_237[k];

        t_315[k] = f_1 * pph_108[k]
                   + f_4 * pc_z[k] * pdh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_y, pc_x, pc_y, sdi0_149, sdi1_149, pdg0_173, \
                         pdg0_175, pdg1_173, pdg1_175, pdh_239, \
                         pdh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_7 * pdg0_173[k]
                   - f_8 * pdg1_173[k]
                   + f_4 * pc_x[k] * pdh_239[k];

        t_317[k] = pa_y[k] * sdi0_149[k]
                   - f_11 * pc_y[k] * sdi1_149[k];

        t_318[k] = f_5 * pdg0_175[k]
                   - f_6 * pdg1_175[k]
                   + f_4 * pc_x[k] * pdh_241[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pc_x, pc_z, pph_111, pdg0_177, pdg0_178, \
                         pdg1_177, pdg1_178, pdh_237, pdh_243, \
                         pdh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_1 * pph_111[k]
                   + f_4 * pc_z[k] * pdh_237[k];

        t_320[k] = f_5 * pdg0_177[k]
                   - f_6 * pdg1_177[k]
                   + f_4 * pc_x[k] * pdh_243[k];

        t_321[k] = f_5 * pdg0_178[k]
                   - f_6 * pdg1_178[k]
                   + f_4 * pc_x[k] * pdh_244[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, t_327, pa_y, pc_x, pc_y, sdi0_154, \
                         sdi1_154, pdh_246, pdh_247, pdh_248, pdh_249, \
                         pdh_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = pa_y[k] * sdi0_154[k]
                   - f_11 * pc_y[k] * sdi1_154[k];

        t_323[k] = f_4 * pc_x[k] * pdh_246[k];

        t_324[k] = f_4 * pc_x[k] * pdh_247[k];

        t_325[k] = f_4 * pc_x[k] * pdh_248[k];

        t_326[k] = f_4 * pc_x[k] * pdh_249[k];

        t_327[k] = f_4 * pc_x[k] * pdh_250[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_y, pc_x, pc_y, pc_z, sdi0_161, sdh_120, \
                         sdi1_161, pph_120, pdh_246, pdh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_4 * pc_x[k] * pdh_251[k];

        t_329[k] = pa_y[k] * sdi0_161[k]
                   + f_12 * sdh_120[k]
                   - f_11 * pc_y[k] * sdi1_161[k];

        t_330[k] = f_1 * pph_120[k]
                   + f_4 * pc_z[k] * pdh_246[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_y, pc_y, sdi0_163, sdi0_164, sdi0_165, \
                         sdh_122, sdh_123, sdh_124, sdi1_163, sdi1_164, \
                         sdi1_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_y[k] * sdi0_163[k]
                   + f_13 * sdh_122[k]
                   - f_11 * pc_y[k] * sdi1_163[k];

        t_332[k] = pa_y[k] * sdi0_164[k]
                   + f_14 * sdh_123[k]
                   - f_11 * pc_y[k] * sdi1_164[k];

        t_333[k] = pa_y[k] * sdi0_165[k]
                   + f_1 * sdh_124[k]
                   - f_11 * pc_y[k] * sdi1_165[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_y, pa_z, pc_y, pc_z, sdi0_0, sdi0_167, \
                         sdh_125, sdi1_0, sdi1_167, pdh_251, pdh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * sdh_125[k]
                   + f_4 * pc_y[k] * pdh_251[k];

        t_335[k] = pa_y[k] * sdi0_167[k]
                   - f_11 * pc_y[k] * sdi1_167[k];

        t_336[k] = pa_z[k] * sdi0_0[k]
                   - f_11 * pc_z[k] * sdi1_0[k];

        t_337[k] = f_4 * pc_y[k] * pdh_252[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_z, pc_y, pc_z, sdi0_2, sdi0_3, sdi0_5, \
                         sdh_0, sdh_2, sdi1_2, sdi1_3, sdi1_5, \
                         pdh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pa_z[k] * sdi0_2[k]
                   + f_0 * sdh_0[k]
                   - f_11 * pc_z[k] * sdi1_2[k];

        t_339[k] = pa_z[k] * sdi0_3[k]
                   - f_11 * pc_z[k] * sdi1_3[k];

        t_340[k] = f_4 * pc_y[k] * pdh_254[k];

        t_341[k] = pa_z[k] * sdi0_5[k]
                   + f_1 * sdh_2[k]
                   - f_11 * pc_z[k] * sdi1_5[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pa_z, pc_y, pc_z, sdi0_6, sdi0_7, sdi0_9, \
                         sdh_3, sdh_5, sdi1_6, sdi1_7, sdi1_9, \
                         pdh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * sdi0_6[k]
                   - f_11 * pc_z[k] * sdi1_6[k];

        t_343[k] = pa_z[k] * sdi0_7[k]
                   + f_0 * sdh_3[k]
                   - f_11 * pc_z[k] * sdi1_7[k];

        t_344[k] = f_4 * pc_y[k] * pdh_257[k];

        t_345[k] = pa_z[k] * sdi0_9[k]
                   + f_14 * sdh_5[k]
                   - f_11 * pc_z[k] * sdi1_9[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_z, pc_y, pc_z, sdi0_10, sdi0_11, \
                         sdi0_12, sdh_6, sdh_7, sdi1_10, sdi1_11, sdi1_12, \
                         pdh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pa_z[k] * sdi0_10[k]
                   - f_11 * pc_z[k] * sdi1_10[k];

        t_347[k] = pa_z[k] * sdi0_11[k]
                   + f_0 * sdh_6[k]
                   - f_11 * pc_z[k] * sdi1_11[k];

        t_348[k] = pa_z[k] * sdi0_12[k]
                   + f_1 * sdh_7[k]
                   - f_11 * pc_z[k] * sdi1_12[k];

        t_349[k] = f_4 * pc_y[k] * pdh_261[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_z, pc_x, pc_z, sdi0_14, sdh_9, \
                         sdi1_14, pph_141, pph_142, pph_143, pdh_267, pdh_268, \
                         pdh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pa_z[k] * sdi0_14[k]
                   + f_13 * sdh_9[k]
                   - f_11 * pc_z[k] * sdi1_14[k];

        t_351[k] = f_1 * pph_141[k]
                   + f_4 * pc_x[k] * pdh_267[k];

        t_352[k] = f_1 * pph_142[k]
                   + f_4 * pc_x[k] * pdh_268[k];

        t_353[k] = f_1 * pph_143[k]
                   + f_4 * pc_x[k] * pdh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_z, pc_x, pc_z, sdi0_21, sdi1_21, \
                         pph_144, pph_145, pph_146, pdh_270, pdh_271, \
                         pdh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_1 * pph_144[k]
                   + f_4 * pc_x[k] * pdh_270[k];

        t_355[k] = f_1 * pph_145[k]
                   + f_4 * pc_x[k] * pdh_271[k];

        t_356[k] = f_1 * pph_146[k]
                   + f_4 * pc_x[k] * pdh_272[k];

        t_357[k] = pa_z[k] * sdi0_21[k]
                   - f_11 * pc_z[k] * sdi1_21[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pc_y, pdg0_191, pdg0_192, pdg0_193, pdg1_191, \
                         pdg1_192, pdg1_193, pdh_268, pdh_269, \
                         pdh_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_16 * pdg0_191[k]
                   - f_17 * pdg1_191[k]
                   + f_4 * pc_y[k] * pdh_268[k];

        t_359[k] = f_9 * pdg0_192[k]
                   - f_10 * pdg1_192[k]
                   + f_4 * pc_y[k] * pdh_269[k];

        t_360[k] = f_7 * pdg0_193[k]
                   - f_8 * pdg1_193[k]
                   + f_4 * pc_y[k] * pdh_270[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pa_z, pc_y, pc_z, sdi0_27, sdi0_28, \
                         sdh_20, sdi1_27, sdi1_28, pdg0_194, pdg1_194, pdh_271, \
                         pdh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_5 * pdg0_194[k]
                   - f_6 * pdg1_194[k]
                   + f_4 * pc_y[k] * pdh_271[k];

        t_362[k] = f_4 * pc_y[k] * pdh_272[k];

        t_363[k] = pa_z[k] * sdi0_27[k]
                   + f_12 * sdh_20[k]
                   - f_11 * pc_z[k] * sdi1_27[k];

        t_364[k] = pa_z[k] * sdi0_28[k]
                   - f_11 * pc_z[k] * sdi1_28[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pa_z, pb_y, pc_y, pc_z, sdi0_31, sdi1_31, \
                         ppi0_170, pph_126, pph_128, ppi1_170, pdh_273, \
                         pdh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_0 * pph_126[k]
                   + f_4 * pc_y[k] * pdh_273[k];

        t_366[k] = pb_y[k] * ppi0_170[k]
                   - f_11 * pc_y[k] * ppi1_170[k];

        t_367[k] = pa_z[k] * sdi0_31[k]
                   - f_11 * pc_z[k] * sdi1_31[k];

        t_368[k] = f_0 * pph_128[k]
                   + f_4 * pc_y[k] * pdh_275[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_z, pb_y, pc_y, pc_z, sdi0_34, sdi0_35, \
                         sdh_24, sdi1_34, sdi1_35, ppi0_173, ppi1_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_y[k] * ppi0_173[k]
                   - f_11 * pc_y[k] * ppi1_173[k];

        t_370[k] = pa_z[k] * sdi0_34[k]
                   - f_11 * pc_z[k] * sdi1_34[k];

        t_371[k] = pa_z[k] * sdi0_35[k]
                   + f_0 * sdh_24[k]
                   - f_11 * pc_z[k] * sdi1_35[k];
    }
}

static auto
compute_prim_pdi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdi0, const size_t sdh,
                                                          const size_t sdi1, const size_t ppi0,
                                                          const size_t pph, const size_t ppi1,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_11 = gamma / q;
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdi0_38 = buffer.data(sdi0 + 38);
    const auto *sdi0_39 = buffer.data(sdi0 + 39);
    const auto *sdi0_40 = buffer.data(sdi0 + 40);
    const auto *sdi0_49 = buffer.data(sdi0 + 49);
    const auto *sdi0_84 = buffer.data(sdi0 + 84);
    const auto *sdi0_87 = buffer.data(sdi0 + 87);
    const auto *sdi0_90 = buffer.data(sdi0 + 90);
    const auto *sdi0_94 = buffer.data(sdi0 + 94);
    const auto *sdi0_105 = buffer.data(sdi0 + 105);
    const auto *sdi0_106 = buffer.data(sdi0 + 106);
    const auto *sdi0_107 = buffer.data(sdi0 + 107);
    const auto *sdi0_108 = buffer.data(sdi0 + 108);
    const auto *sdi0_109 = buffer.data(sdi0 + 109);
    const auto *sdi0_111 = buffer.data(sdi0 + 111);

    const auto *sdh_27 = buffer.data(sdh + 27);
    const auto *sdh_28 = buffer.data(sdh + 28);
    const auto *sdh_78 = buffer.data(sdh + 78);
    const auto *sdh_79 = buffer.data(sdh + 79);
    const auto *sdh_80 = buffer.data(sdh + 80);
    const auto *sdh_81 = buffer.data(sdh + 81);
    const auto *sdh_83 = buffer.data(sdh + 83);

    const auto *sdi1_38 = buffer.data(sdi1 + 38);
    const auto *sdi1_39 = buffer.data(sdi1 + 39);
    const auto *sdi1_40 = buffer.data(sdi1 + 40);
    const auto *sdi1_49 = buffer.data(sdi1 + 49);
    const auto *sdi1_84 = buffer.data(sdi1 + 84);
    const auto *sdi1_87 = buffer.data(sdi1 + 87);
    const auto *sdi1_90 = buffer.data(sdi1 + 90);
    const auto *sdi1_94 = buffer.data(sdi1 + 94);
    const auto *sdi1_105 = buffer.data(sdi1 + 105);
    const auto *sdi1_106 = buffer.data(sdi1 + 106);
    const auto *sdi1_107 = buffer.data(sdi1 + 107);
    const auto *sdi1_108 = buffer.data(sdi1 + 108);
    const auto *sdi1_109 = buffer.data(sdi1 + 109);
    const auto *sdi1_111 = buffer.data(sdi1 + 111);

    const auto *ppi0_177 = buffer.data(ppi0 + 177);
    const auto *ppi0_182 = buffer.data(ppi0 + 182);
    const auto *ppi0_218 = buffer.data(ppi0 + 218);
    const auto *ppi0_219 = buffer.data(ppi0 + 219);
    const auto *ppi0_220 = buffer.data(ppi0 + 220);
    const auto *ppi0_221 = buffer.data(ppi0 + 221);
    const auto *ppi0_223 = buffer.data(ppi0 + 223);
    const auto *ppi0_224 = buffer.data(ppi0 + 224);
    const auto *ppi0_226 = buffer.data(ppi0 + 226);
    const auto *ppi0_229 = buffer.data(ppi0 + 229);
    const auto *ppi0_231 = buffer.data(ppi0 + 231);
    const auto *ppi0_233 = buffer.data(ppi0 + 233);
    const auto *ppi0_235 = buffer.data(ppi0 + 235);
    const auto *ppi0_236 = buffer.data(ppi0 + 236);
    const auto *ppi0_238 = buffer.data(ppi0 + 238);
    const auto *ppi0_245 = buffer.data(ppi0 + 245);
    const auto *ppi0_246 = buffer.data(ppi0 + 246);
    const auto *ppi0_247 = buffer.data(ppi0 + 247);
    const auto *ppi0_248 = buffer.data(ppi0 + 248);
    const auto *ppi0_249 = buffer.data(ppi0 + 249);
    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *pph_131 = buffer.data(pph + 131);
    const auto *pph_135 = buffer.data(pph + 135);
    const auto *pph_146 = buffer.data(pph + 146);
    const auto *pph_147 = buffer.data(pph + 147);
    const auto *pph_149 = buffer.data(pph + 149);
    const auto *pph_152 = buffer.data(pph + 152);
    const auto *pph_156 = buffer.data(pph + 156);
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

    const auto *ppi1_177 = buffer.data(ppi1 + 177);
    const auto *ppi1_182 = buffer.data(ppi1 + 182);
    const auto *ppi1_218 = buffer.data(ppi1 + 218);
    const auto *ppi1_219 = buffer.data(ppi1 + 219);
    const auto *ppi1_220 = buffer.data(ppi1 + 220);
    const auto *ppi1_221 = buffer.data(ppi1 + 221);
    const auto *ppi1_223 = buffer.data(ppi1 + 223);
    const auto *ppi1_224 = buffer.data(ppi1 + 224);
    const auto *ppi1_226 = buffer.data(ppi1 + 226);
    const auto *ppi1_229 = buffer.data(ppi1 + 229);
    const auto *ppi1_231 = buffer.data(ppi1 + 231);
    const auto *ppi1_233 = buffer.data(ppi1 + 233);
    const auto *ppi1_235 = buffer.data(ppi1 + 235);
    const auto *ppi1_236 = buffer.data(ppi1 + 236);
    const auto *ppi1_238 = buffer.data(ppi1 + 238);
    const auto *ppi1_245 = buffer.data(ppi1 + 245);
    const auto *ppi1_246 = buffer.data(ppi1 + 246);
    const auto *ppi1_247 = buffer.data(ppi1 + 247);
    const auto *ppi1_248 = buffer.data(ppi1 + 248);
    const auto *ppi1_249 = buffer.data(ppi1 + 249);
    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *pdg0_210 = buffer.data(pdg0 + 210);
    const auto *pdg0_213 = buffer.data(pdg0 + 213);
    const auto *pdg0_216 = buffer.data(pdg0 + 216);
    const auto *pdg0_220 = buffer.data(pdg0 + 220);
    const auto *pdg0_227 = buffer.data(pdg0 + 227);
    const auto *pdg0_230 = buffer.data(pdg0 + 230);
    const auto *pdg0_232 = buffer.data(pdg0 + 232);
    const auto *pdg0_234 = buffer.data(pdg0 + 234);
    const auto *pdg0_236 = buffer.data(pdg0 + 236);
    const auto *pdg0_237 = buffer.data(pdg0 + 237);
    const auto *pdg0_239 = buffer.data(pdg0 + 239);
    const auto *pdg0_243 = buffer.data(pdg0 + 243);
    const auto *pdg0_246 = buffer.data(pdg0 + 246);
    const auto *pdg0_247 = buffer.data(pdg0 + 247);
    const auto *pdg0_250 = buffer.data(pdg0 + 250);
    const auto *pdg0_251 = buffer.data(pdg0 + 251);
    const auto *pdg0_252 = buffer.data(pdg0 + 252);
    const auto *pdg0_255 = buffer.data(pdg0 + 255);
    const auto *pdg0_257 = buffer.data(pdg0 + 257);
    const auto *pdg0_258 = buffer.data(pdg0 + 258);
    const auto *pdg0_260 = buffer.data(pdg0 + 260);
    const auto *pdg0_261 = buffer.data(pdg0 + 261);
    const auto *pdg0_262 = buffer.data(pdg0 + 262);
    const auto *pdg0_264 = buffer.data(pdg0 + 264);
    const auto *pdg0_265 = buffer.data(pdg0 + 265);
    const auto *pdg0_266 = buffer.data(pdg0 + 266);
    const auto *pdg0_267 = buffer.data(pdg0 + 267);
    const auto *pdg0_269 = buffer.data(pdg0 + 269);

    const auto *pdg1_210 = buffer.data(pdg1 + 210);
    const auto *pdg1_213 = buffer.data(pdg1 + 213);
    const auto *pdg1_216 = buffer.data(pdg1 + 216);
    const auto *pdg1_220 = buffer.data(pdg1 + 220);
    const auto *pdg1_227 = buffer.data(pdg1 + 227);
    const auto *pdg1_230 = buffer.data(pdg1 + 230);
    const auto *pdg1_232 = buffer.data(pdg1 + 232);
    const auto *pdg1_234 = buffer.data(pdg1 + 234);
    const auto *pdg1_236 = buffer.data(pdg1 + 236);
    const auto *pdg1_237 = buffer.data(pdg1 + 237);
    const auto *pdg1_239 = buffer.data(pdg1 + 239);
    const auto *pdg1_243 = buffer.data(pdg1 + 243);
    const auto *pdg1_246 = buffer.data(pdg1 + 246);
    const auto *pdg1_247 = buffer.data(pdg1 + 247);
    const auto *pdg1_250 = buffer.data(pdg1 + 250);
    const auto *pdg1_251 = buffer.data(pdg1 + 251);
    const auto *pdg1_252 = buffer.data(pdg1 + 252);
    const auto *pdg1_255 = buffer.data(pdg1 + 255);
    const auto *pdg1_257 = buffer.data(pdg1 + 257);
    const auto *pdg1_258 = buffer.data(pdg1 + 258);
    const auto *pdg1_260 = buffer.data(pdg1 + 260);
    const auto *pdg1_261 = buffer.data(pdg1 + 261);
    const auto *pdg1_262 = buffer.data(pdg1 + 262);
    const auto *pdg1_264 = buffer.data(pdg1 + 264);
    const auto *pdg1_265 = buffer.data(pdg1 + 265);
    const auto *pdg1_266 = buffer.data(pdg1 + 266);
    const auto *pdg1_267 = buffer.data(pdg1 + 267);
    const auto *pdg1_269 = buffer.data(pdg1 + 269);

    const auto *pdh_278 = buffer.data(pdh + 278);
    const auto *pdh_282 = buffer.data(pdh + 282);
    const auto *pdh_288 = buffer.data(pdh + 288);
    const auto *pdh_289 = buffer.data(pdh + 289);
    const auto *pdh_290 = buffer.data(pdh + 290);
    const auto *pdh_291 = buffer.data(pdh + 291);
    const auto *pdh_292 = buffer.data(pdh + 292);
    const auto *pdh_293 = buffer.data(pdh + 293);
    const auto *pdh_294 = buffer.data(pdh + 294);
    const auto *pdh_296 = buffer.data(pdh + 296);
    const auto *pdh_297 = buffer.data(pdh + 297);
    const auto *pdh_299 = buffer.data(pdh + 299);
    const auto *pdh_300 = buffer.data(pdh + 300);
    const auto *pdh_303 = buffer.data(pdh + 303);
    const auto *pdh_304 = buffer.data(pdh + 304);
    const auto *pdh_309 = buffer.data(pdh + 309);
    const auto *pdh_310 = buffer.data(pdh + 310);
    const auto *pdh_311 = buffer.data(pdh + 311);
    const auto *pdh_312 = buffer.data(pdh + 312);
    const auto *pdh_313 = buffer.data(pdh + 313);
    const auto *pdh_314 = buffer.data(pdh + 314);
    const auto *pdh_315 = buffer.data(pdh + 315);
    const auto *pdh_317 = buffer.data(pdh + 317);
    const auto *pdh_320 = buffer.data(pdh + 320);
    const auto *pdh_322 = buffer.data(pdh + 322);
    const auto *pdh_324 = buffer.data(pdh + 324);
    const auto *pdh_326 = buffer.data(pdh + 326);
    const auto *pdh_327 = buffer.data(pdh + 327);
    const auto *pdh_329 = buffer.data(pdh + 329);
    const auto *pdh_330 = buffer.data(pdh + 330);
    const auto *pdh_331 = buffer.data(pdh + 331);
    const auto *pdh_332 = buffer.data(pdh + 332);
    const auto *pdh_333 = buffer.data(pdh + 333);
    const auto *pdh_334 = buffer.data(pdh + 334);
    const auto *pdh_335 = buffer.data(pdh + 335);
    const auto *pdh_336 = buffer.data(pdh + 336);
    const auto *pdh_338 = buffer.data(pdh + 338);
    const auto *pdh_339 = buffer.data(pdh + 339);
    const auto *pdh_341 = buffer.data(pdh + 341);
    const auto *pdh_342 = buffer.data(pdh + 342);
    const auto *pdh_343 = buffer.data(pdh + 343);
    const auto *pdh_345 = buffer.data(pdh + 345);
    const auto *pdh_346 = buffer.data(pdh + 346);
    const auto *pdh_347 = buffer.data(pdh + 347);
    const auto *pdh_348 = buffer.data(pdh + 348);
    const auto *pdh_351 = buffer.data(pdh + 351);
    const auto *pdh_352 = buffer.data(pdh + 352);
    const auto *pdh_353 = buffer.data(pdh + 353);
    const auto *pdh_354 = buffer.data(pdh + 354);
    const auto *pdh_355 = buffer.data(pdh + 355);
    const auto *pdh_356 = buffer.data(pdh + 356);
    const auto *pdh_357 = buffer.data(pdh + 357);
    const auto *pdh_359 = buffer.data(pdh + 359);
    const auto *pdh_360 = buffer.data(pdh + 360);
    const auto *pdh_362 = buffer.data(pdh + 362);
    const auto *pdh_363 = buffer.data(pdh + 363);
    const auto *pdh_364 = buffer.data(pdh + 364);
    const auto *pdh_366 = buffer.data(pdh + 366);
    const auto *pdh_367 = buffer.data(pdh + 367);
    const auto *pdh_368 = buffer.data(pdh + 368);
    const auto *pdh_369 = buffer.data(pdh + 369);
    const auto *pdh_371 = buffer.data(pdh + 371);
    const auto *pdh_372 = buffer.data(pdh + 372);
    const auto *pdh_373 = buffer.data(pdh + 373);
    const auto *pdh_374 = buffer.data(pdh + 374);
    const auto *pdh_375 = buffer.data(pdh + 375);
    const auto *pdh_376 = buffer.data(pdh + 376);
    const auto *pdh_377 = buffer.data(pdh + 377);

#pragma omp simd aligned(t_372, t_373, t_374, pa_z, pb_y, pc_y, pc_z, sdi0_38, sdi1_38, \
                         ppi0_177, pph_131, ppi1_177, pdh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_0 * pph_131[k]
                   + f_4 * pc_y[k] * pdh_278[k];

        t_373[k] = pb_y[k] * ppi0_177[k]
                   - f_11 * pc_y[k] * ppi1_177[k];

        t_374[k] = pa_z[k] * sdi0_38[k]
                   - f_11 * pc_z[k] * sdi1_38[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_z, pc_y, pc_z, sdi0_39, sdi0_40, sdh_27, \
                         sdh_28, sdi1_39, sdi1_40, pph_135, pdh_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pa_z[k] * sdi0_39[k]
                   + f_0 * sdh_27[k]
                   - f_11 * pc_z[k] * sdi1_39[k];

        t_376[k] = pa_z[k] * sdi0_40[k]
                   + f_1 * sdh_28[k]
                   - f_11 * pc_z[k] * sdi1_40[k];

        t_377[k] = f_0 * pph_135[k]
                   + f_4 * pc_y[k] * pdh_282[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_x, pc_y, ppi0_182, pph_162, \
                         pph_163, pph_164, ppi1_182, pdh_288, pdh_289, \
                         pdh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * ppi0_182[k]
                   - f_11 * pc_y[k] * ppi1_182[k];

        t_379[k] = f_0 * pph_162[k]
                   + f_4 * pc_x[k] * pdh_288[k];

        t_380[k] = f_0 * pph_163[k]
                   + f_4 * pc_x[k] * pdh_289[k];

        t_381[k] = f_0 * pph_164[k]
                   + f_4 * pc_x[k] * pdh_290[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_z, pc_x, pc_z, sdi0_49, sdi1_49, \
                         pph_165, pph_166, pph_167, pdh_291, pdh_292, \
                         pdh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_0 * pph_165[k]
                   + f_4 * pc_x[k] * pdh_291[k];

        t_383[k] = f_0 * pph_166[k]
                   + f_4 * pc_x[k] * pdh_292[k];

        t_384[k] = f_0 * pph_167[k]
                   + f_4 * pc_x[k] * pdh_293[k];

        t_385[k] = pa_z[k] * sdi0_49[k]
                   - f_11 * pc_z[k] * sdi1_49[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pb_x, pc_x, ppi0_218, ppi0_219, ppi0_220, \
                         ppi0_221, ppi1_218, ppi1_219, ppi1_220, \
                         ppi1_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = pb_x[k] * ppi0_218[k]
                   - f_11 * pc_x[k] * ppi1_218[k];

        t_387[k] = pb_x[k] * ppi0_219[k]
                   - f_11 * pc_x[k] * ppi1_219[k];

        t_388[k] = pb_x[k] * ppi0_220[k]
                   - f_11 * pc_x[k] * ppi1_220[k];

        t_389[k] = pb_x[k] * ppi0_221[k]
                   - f_11 * pc_x[k] * ppi1_221[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pb_x, pc_x, pc_y, ppi0_223, pph_146, \
                         pph_168, ppi1_223, pdg0_210, pdg1_210, pdh_293, \
                         pdh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_0 * pph_146[k]
                   + f_4 * pc_y[k] * pdh_293[k];

        t_391[k] = pb_x[k] * ppi0_223[k]
                   - f_11 * pc_x[k] * ppi1_223[k];

        t_392[k] = f_0 * pph_168[k]
                   + f_2 * pdg0_210[k]
                   - f_3 * pdg1_210[k]
                   + f_4 * pc_x[k] * pdh_294[k];

        t_393[k] = f_4 * pc_y[k] * pdh_294[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pb_x, pc_x, pc_y, ppi0_226, pph_170, pph_171, \
                         ppi1_226, pdg0_213, pdg1_213, pdh_296, \
                         pdh_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = pb_x[k] * ppi0_226[k]
                   + f_15 * pph_170[k]
                   - f_11 * pc_x[k] * ppi1_226[k];

        t_395[k] = f_0 * pph_171[k]
                   + f_9 * pdg0_213[k]
                   - f_10 * pdg1_213[k]
                   + f_4 * pc_x[k] * pdh_297[k];

        t_396[k] = f_4 * pc_y[k] * pdh_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pb_x, pc_x, ppi0_229, ppi0_231, pph_173, \
                         pph_174, pph_175, ppi1_229, ppi1_231, pdg0_216, pdg1_216, \
                         pdh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pb_x[k] * ppi0_229[k]
                   + f_13 * pph_173[k]
                   - f_11 * pc_x[k] * ppi1_229[k];

        t_398[k] = f_0 * pph_174[k]
                   + f_7 * pdg0_216[k]
                   - f_8 * pdg1_216[k]
                   + f_4 * pc_x[k] * pdh_300[k];

        t_399[k] = pb_x[k] * ppi0_231[k]
                   + f_14 * pph_175[k]
                   - f_11 * pc_x[k] * ppi1_231[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pb_x, pc_x, pc_y, ppi0_233, pph_177, pph_178, \
                         ppi1_233, pdg0_220, pdg1_220, pdh_299, \
                         pdh_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_4 * pc_y[k] * pdh_299[k];

        t_401[k] = pb_x[k] * ppi0_233[k]
                   + f_14 * pph_177[k]
                   - f_11 * pc_x[k] * ppi1_233[k];

        t_402[k] = f_0 * pph_178[k]
                   + f_5 * pdg0_220[k]
                   - f_6 * pdg1_220[k]
                   + f_4 * pc_x[k] * pdh_304[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pb_x, pc_x, pc_y, ppi0_235, ppi0_236, pph_179, \
                         pph_180, ppi1_235, ppi1_236, pdh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = pb_x[k] * ppi0_235[k]
                   + f_1 * pph_179[k]
                   - f_11 * pc_x[k] * ppi1_235[k];

        t_404[k] = pb_x[k] * ppi0_236[k]
                   + f_1 * pph_180[k]
                   - f_11 * pc_x[k] * ppi1_236[k];

        t_405[k] = f_4 * pc_y[k] * pdh_303[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pb_x, pc_x, ppi0_238, pph_182, pph_183, \
                         pph_184, pph_185, ppi1_238, pdh_309, pdh_310, \
                         pdh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pb_x[k] * ppi0_238[k]
                   + f_1 * pph_182[k]
                   - f_11 * pc_x[k] * ppi1_238[k];

        t_407[k] = f_0 * pph_183[k]
                   + f_4 * pc_x[k] * pdh_309[k];

        t_408[k] = f_0 * pph_184[k]
                   + f_4 * pc_x[k] * pdh_310[k];

        t_409[k] = f_0 * pph_185[k]
                   + f_4 * pc_x[k] * pdh_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, pc_x, ppi0_245, pph_186, pph_187, \
                         pph_188, ppi1_245, pdh_312, pdh_313, pdh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_0 * pph_186[k]
                   + f_4 * pc_x[k] * pdh_312[k];

        t_411[k] = f_0 * pph_187[k]
                   + f_4 * pc_x[k] * pdh_313[k];

        t_412[k] = f_0 * pph_188[k]
                   + f_4 * pc_x[k] * pdh_314[k];

        t_413[k] = pb_x[k] * ppi0_245[k]
                   - f_11 * pc_x[k] * ppi1_245[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pc_x, ppi0_246, ppi0_247, ppi0_248, \
                         ppi0_249, ppi1_246, ppi1_247, ppi1_248, \
                         ppi1_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pb_x[k] * ppi0_246[k]
                   - f_11 * pc_x[k] * ppi1_246[k];

        t_415[k] = pb_x[k] * ppi0_247[k]
                   - f_11 * pc_x[k] * ppi1_247[k];

        t_416[k] = pb_x[k] * ppi0_248[k]
                   - f_11 * pc_x[k] * ppi1_248[k];

        t_417[k] = pb_x[k] * ppi0_249[k]
                   - f_11 * pc_x[k] * ppi1_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pa_z, pb_x, pc_x, pc_y, pc_z, sdi0_84, \
                         sdi1_84, ppi0_251, pph_147, ppi1_251, pdh_314, \
                         pdh_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_4 * pc_y[k] * pdh_314[k];

        t_419[k] = pb_x[k] * ppi0_251[k]
                   - f_11 * pc_x[k] * ppi1_251[k];

        t_420[k] = pa_z[k] * sdi0_84[k]
                   - f_11 * pc_z[k] * sdi1_84[k];

        t_421[k] = f_1 * pph_147[k]
                   + f_4 * pc_y[k] * pdh_315[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, pa_z, pc_x, pc_y, pc_z, sdi0_87, sdi1_87, \
                         pph_149, pdg0_227, pdg1_227, pdh_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_16 * pdg0_227[k]
                   - f_17 * pdg1_227[k]
                   + f_4 * pc_x[k] * pdh_317[k];

        t_423[k] = pa_z[k] * sdi0_87[k]
                   - f_11 * pc_z[k] * sdi1_87[k];

        t_424[k] = f_1 * pph_149[k]
                   + f_4 * pc_y[k] * pdh_317[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, pa_z, pc_x, pc_z, sdi0_90, sdi1_90, pdg0_230, \
                         pdg0_232, pdg1_230, pdg1_232, pdh_320, \
                         pdh_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_9 * pdg0_230[k]
                   - f_10 * pdg1_230[k]
                   + f_4 * pc_x[k] * pdh_320[k];

        t_426[k] = pa_z[k] * sdi0_90[k]
                   - f_11 * pc_z[k] * sdi1_90[k];

        t_427[k] = f_7 * pdg0_232[k]
                   - f_8 * pdg1_232[k]
                   + f_4 * pc_x[k] * pdh_322[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, pa_z, pc_x, pc_y, pc_z, sdi0_94, sdi1_94, \
                         pph_152, pdg0_234, pdg1_234, pdh_320, \
                         pdh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_1 * pph_152[k]
                   + f_4 * pc_y[k] * pdh_320[k];

        t_429[k] = f_7 * pdg0_234[k]
                   - f_8 * pdg1_234[k]
                   + f_4 * pc_x[k] * pdh_324[k];

        t_430[k] = pa_z[k] * sdi0_94[k]
                   - f_11 * pc_z[k] * sdi1_94[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, pc_x, pc_y, pph_156, pdg0_236, pdg0_237, \
                         pdg1_236, pdg1_237, pdh_324, pdh_326, \
                         pdh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_5 * pdg0_236[k]
                   - f_6 * pdg1_236[k]
                   + f_4 * pc_x[k] * pdh_326[k];

        t_432[k] = f_5 * pdg0_237[k]
                   - f_6 * pdg1_237[k]
                   + f_4 * pc_x[k] * pdh_327[k];

        t_433[k] = f_1 * pph_156[k]
                   + f_4 * pc_y[k] * pdh_324[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, t_439, pc_x, pdg0_239, pdg1_239, \
                         pdh_329, pdh_330, pdh_331, pdh_332, pdh_333, \
                         pdh_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_5 * pdg0_239[k]
                   - f_6 * pdg1_239[k]
                   + f_4 * pc_x[k] * pdh_329[k];

        t_435[k] = f_4 * pc_x[k] * pdh_330[k];

        t_436[k] = f_4 * pc_x[k] * pdh_331[k];

        t_437[k] = f_4 * pc_x[k] * pdh_332[k];

        t_438[k] = f_4 * pc_x[k] * pdh_333[k];

        t_439[k] = f_4 * pc_x[k] * pdh_334[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pa_z, pc_x, pc_z, sdi0_105, sdi0_106, \
                         sdi0_107, sdh_78, sdh_79, sdi1_105, sdi1_106, sdi1_107, \
                         pdh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_4 * pc_x[k] * pdh_335[k];

        t_441[k] = pa_z[k] * sdi0_105[k]
                   - f_11 * pc_z[k] * sdi1_105[k];

        t_442[k] = pa_z[k] * sdi0_106[k]
                   + f_0 * sdh_78[k]
                   - f_11 * pc_z[k] * sdi1_106[k];

        t_443[k] = pa_z[k] * sdi0_107[k]
                   + f_1 * sdh_79[k]
                   - f_11 * pc_z[k] * sdi1_107[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pa_z, pc_y, pc_z, sdi0_108, sdi0_109, sdh_80, \
                         sdh_81, sdi1_108, sdi1_109, pph_167, pdh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = pa_z[k] * sdi0_108[k]
                   + f_14 * sdh_80[k]
                   - f_11 * pc_z[k] * sdi1_108[k];

        t_445[k] = pa_z[k] * sdi0_109[k]
                   + f_13 * sdh_81[k]
                   - f_11 * pc_z[k] * sdi1_109[k];

        t_446[k] = f_1 * pph_167[k]
                   + f_4 * pc_y[k] * pdh_335[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pa_z, pb_y, pc_y, pc_z, sdi0_111, sdh_83, \
                         sdi1_111, ppi0_224, pph_168, ppi1_224, \
                         pdh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = pa_z[k] * sdi0_111[k]
                   + f_12 * sdh_83[k]
                   - f_11 * pc_z[k] * sdi1_111[k];

        t_448[k] = pb_y[k] * ppi0_224[k]
                   - f_11 * pc_y[k] * ppi1_224[k];

        t_449[k] = f_0 * pph_168[k]
                   + f_4 * pc_y[k] * pdh_336[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pb_y, pc_x, pc_y, ppi0_226, ppi0_229, \
                         pph_170, ppi1_226, ppi1_229, pdg0_243, pdg1_243, pdh_338, \
                         pdh_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = pb_y[k] * ppi0_226[k]
                   - f_11 * pc_y[k] * ppi1_226[k];

        t_451[k] = f_9 * pdg0_243[k]
                   - f_10 * pdg1_243[k]
                   + f_4 * pc_x[k] * pdh_339[k];

        t_452[k] = f_0 * pph_170[k]
                   + f_4 * pc_y[k] * pdh_338[k];

        t_453[k] = pb_y[k] * ppi0_229[k]
                   - f_11 * pc_y[k] * ppi1_229[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, pc_x, pc_y, pph_173, pdg0_246, pdg0_247, \
                         pdg1_246, pdg1_247, pdh_341, pdh_342, \
                         pdh_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_7 * pdg0_246[k]
                   - f_8 * pdg1_246[k]
                   + f_4 * pc_x[k] * pdh_342[k];

        t_455[k] = f_7 * pdg0_247[k]
                   - f_8 * pdg1_247[k]
                   + f_4 * pc_x[k] * pdh_343[k];

        t_456[k] = f_0 * pph_173[k]
                   + f_4 * pc_y[k] * pdh_341[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pb_y, pc_x, pc_y, ppi0_233, ppi1_233, pdg0_250, \
                         pdg0_251, pdg1_250, pdg1_251, pdh_346, \
                         pdh_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = pb_y[k] * ppi0_233[k]
                   - f_11 * pc_y[k] * ppi1_233[k];

        t_458[k] = f_5 * pdg0_250[k]
                   - f_6 * pdg1_250[k]
                   + f_4 * pc_x[k] * pdh_346[k];

        t_459[k] = f_5 * pdg0_251[k]
                   - f_6 * pdg1_251[k]
                   + f_4 * pc_x[k] * pdh_347[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_y, pc_x, pc_y, ppi0_238, pph_177, \
                         ppi1_238, pdg0_252, pdg1_252, pdh_345, pdh_348, \
                         pdh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_5 * pdg0_252[k]
                   - f_6 * pdg1_252[k]
                   + f_4 * pc_x[k] * pdh_348[k];

        t_461[k] = f_0 * pph_177[k]
                   + f_4 * pc_y[k] * pdh_345[k];

        t_462[k] = pb_y[k] * ppi0_238[k]
                   - f_11 * pc_y[k] * ppi1_238[k];

        t_463[k] = f_4 * pc_x[k] * pdh_351[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, pc_x, pdh_352, pdh_353, pdh_354, \
                         pdh_355, pdh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_4 * pc_x[k] * pdh_352[k];

        t_465[k] = f_4 * pc_x[k] * pdh_353[k];

        t_466[k] = f_4 * pc_x[k] * pdh_354[k];

        t_467[k] = f_4 * pc_x[k] * pdh_355[k];

        t_468[k] = f_4 * pc_x[k] * pdh_356[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, pb_y, pc_y, ppi0_246, ppi0_247, pph_183, \
                         pph_184, pph_185, ppi1_246, ppi1_247, pdg0_250, pdg1_250, \
                         pdh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_0 * pph_183[k]
                   + f_2 * pdg0_250[k]
                   - f_3 * pdg1_250[k]
                   + f_4 * pc_y[k] * pdh_351[k];

        t_470[k] = pb_y[k] * ppi0_246[k]
                   + f_15 * pph_184[k]
                   - f_11 * pc_y[k] * ppi1_246[k];

        t_471[k] = pb_y[k] * ppi0_247[k]
                   + f_13 * pph_185[k]
                   - f_11 * pc_y[k] * ppi1_247[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, pb_y, pc_y, ppi0_248, ppi0_249, ppi0_251, \
                         pph_186, pph_187, pph_188, ppi1_248, ppi1_249, ppi1_251, \
                         pdh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = pb_y[k] * ppi0_248[k]
                   + f_14 * pph_186[k]
                   - f_11 * pc_y[k] * ppi1_248[k];

        t_473[k] = pb_y[k] * ppi0_249[k]
                   + f_1 * pph_187[k]
                   - f_11 * pc_y[k] * ppi1_249[k];

        t_474[k] = f_0 * pph_188[k]
                   + f_4 * pc_y[k] * pdh_356[k];

        t_475[k] = pb_y[k] * ppi0_251[k]
                   - f_11 * pc_y[k] * ppi1_251[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, t_480, pc_x, pc_y, pdg0_255, pdg0_257, \
                         pdg0_258, pdg1_255, pdg1_257, pdg1_258, pdh_357, pdh_359, \
                         pdh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_2 * pdg0_255[k]
                   - f_3 * pdg1_255[k]
                   + f_4 * pc_x[k] * pdh_357[k];

        t_477[k] = f_4 * pc_y[k] * pdh_357[k];

        t_478[k] = f_16 * pdg0_257[k]
                   - f_17 * pdg1_257[k]
                   + f_4 * pc_x[k] * pdh_359[k];

        t_479[k] = f_9 * pdg0_258[k]
                   - f_10 * pdg1_258[k]
                   + f_4 * pc_x[k] * pdh_360[k];

        t_480[k] = f_4 * pc_y[k] * pdh_359[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pdg0_260, pdg0_261, pdg0_262, \
                         pdg1_260, pdg1_261, pdg1_262, pdh_362, pdh_363, \
                         pdh_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_9 * pdg0_260[k]
                   - f_10 * pdg1_260[k]
                   + f_4 * pc_x[k] * pdh_362[k];

        t_482[k] = f_7 * pdg0_261[k]
                   - f_8 * pdg1_261[k]
                   + f_4 * pc_x[k] * pdh_363[k];

        t_483[k] = f_7 * pdg0_262[k]
                   - f_8 * pdg1_262[k]
                   + f_4 * pc_x[k] * pdh_364[k];

        t_484[k] = f_4 * pc_y[k] * pdh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pdg0_264, pdg0_265, pdg0_266, pdg1_264, \
                         pdg1_265, pdg1_266, pdh_366, pdh_367, \
                         pdh_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_7 * pdg0_264[k]
                   - f_8 * pdg1_264[k]
                   + f_4 * pc_x[k] * pdh_366[k];

        t_486[k] = f_5 * pdg0_265[k]
                   - f_6 * pdg1_265[k]
                   + f_4 * pc_x[k] * pdh_367[k];

        t_487[k] = f_5 * pdg0_266[k]
                   - f_6 * pdg1_266[k]
                   + f_4 * pc_x[k] * pdh_368[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, pc_x, pc_y, pdg0_267, pdg0_269, \
                         pdg1_267, pdg1_269, pdh_366, pdh_369, pdh_371, pdh_372, \
                         pdh_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * pdg0_267[k]
                   - f_6 * pdg1_267[k]
                   + f_4 * pc_x[k] * pdh_369[k];

        t_489[k] = f_4 * pc_y[k] * pdh_366[k];

        t_490[k] = f_5 * pdg0_269[k]
                   - f_6 * pdg1_269[k]
                   + f_4 * pc_x[k] * pdh_371[k];

        t_491[k] = f_4 * pc_x[k] * pdh_372[k];

        t_492[k] = f_4 * pc_x[k] * pdh_373[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, pc_x, pc_y, pdg0_265, pdg1_265, \
                         pdh_372, pdh_374, pdh_375, pdh_376, pdh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_4 * pc_x[k] * pdh_374[k];

        t_494[k] = f_4 * pc_x[k] * pdh_375[k];

        t_495[k] = f_4 * pc_x[k] * pdh_376[k];

        t_496[k] = f_4 * pc_x[k] * pdh_377[k];

        t_497[k] = f_2 * pdg0_265[k]
                   - f_3 * pdg1_265[k]
                   + f_4 * pc_y[k] * pdh_372[k];
    }
}

static auto
compute_prim_pdi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sdh, const size_t pph,
                                                          const size_t pdg0, const size_t pdg1,
                                                          const size_t pdh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.5 / gamma;
    const auto f_3 = 2.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 1.5 / gamma;
    const auto f_10 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh_125 = buffer.data(sdh + 125);

    const auto *pph_188 = buffer.data(pph + 188);

    const auto *pdg0_266 = buffer.data(pdg0 + 266);
    const auto *pdg0_267 = buffer.data(pdg0 + 267);
    const auto *pdg0_268 = buffer.data(pdg0 + 268);
    const auto *pdg0_269 = buffer.data(pdg0 + 269);

    const auto *pdg1_266 = buffer.data(pdg1 + 266);
    const auto *pdg1_267 = buffer.data(pdg1 + 267);
    const auto *pdg1_268 = buffer.data(pdg1 + 268);
    const auto *pdg1_269 = buffer.data(pdg1 + 269);

    const auto *pdh_373 = buffer.data(pdh + 373);
    const auto *pdh_374 = buffer.data(pdh + 374);
    const auto *pdh_375 = buffer.data(pdh + 375);
    const auto *pdh_376 = buffer.data(pdh + 376);
    const auto *pdh_377 = buffer.data(pdh + 377);

#pragma omp simd aligned(t_498, t_499, t_500, pc_y, pdg0_266, pdg0_267, pdg0_268, pdg1_266, \
                         pdg1_267, pdg1_268, pdh_373, pdh_374, \
                         pdh_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_16 * pdg0_266[k]
                   - f_17 * pdg1_266[k]
                   + f_4 * pc_y[k] * pdh_373[k];

        t_499[k] = f_9 * pdg0_267[k]
                   - f_10 * pdg1_267[k]
                   + f_4 * pc_y[k] * pdh_374[k];

        t_500[k] = f_7 * pdg0_268[k]
                   - f_8 * pdg1_268[k]
                   + f_4 * pc_y[k] * pdh_375[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pc_y, pc_z, sdh_125, pph_188, pdg0_269, \
                         pdg1_269, pdh_376, pdh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_5 * pdg0_269[k]
                   - f_6 * pdg1_269[k]
                   + f_4 * pc_y[k] * pdh_376[k];

        t_502[k] = f_4 * pc_y[k] * pdh_377[k];

        t_503[k] = f_0 * sdh_125[k]
                   + f_1 * pph_188[k]
                   + f_2 * pdg0_269[k]
                   - f_3 * pdg1_269[k]
                   + f_4 * pc_z[k] * pdh_377[k];
    }
}

auto
compute_prim_pdi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sdi0,
                                                   const size_t sdh, const size_t sdi1,
                                                   const size_t ppi0, const size_t pph,
                                                   const size_t ppi1, const size_t pdg0,
                                                   const size_t pdg1, const size_t pdh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pdi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sdi0,
                                                              sdh, sdi1, ppi0, pph, ppi1, pdg0,
                                                              pdg1, pdh, ncols, gamma, p, q);

    compute_prim_pdi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sdi0,
                                                              sdh, sdi1, ppi0, pph, ppi1, pdg0,
                                                              pdg1, pdh, ncols, gamma, p, q);

    compute_prim_pdi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sdi0,
                                                              sdh, sdi1, ppi0, pph, ppi1, pdg0,
                                                              pdg1, pdh, ncols, gamma, p, q);

    compute_prim_pdi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sdi0,
                                                              sdh, sdi1, ppi0, pph, ppi1, pdg0,
                                                              pdg1, pdh, ncols, gamma, p, q);

    compute_prim_pdi_three_center_electron_repulsion_0_piece4(buffer, target, pc, sdh, pph,
                                                              pdg0, pdg1, pdh, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
