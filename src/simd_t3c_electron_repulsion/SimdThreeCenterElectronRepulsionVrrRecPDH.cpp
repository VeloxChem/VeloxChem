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


#include "SimdThreeCenterElectronRepulsionVrrRecPDH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pdh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdh0, const size_t sdg,
                                                          const size_t sdh1, const size_t pph0,
                                                          const size_t ppg, const size_t pph1,
                                                          const size_t pdf0, const size_t pdf1,
                                                          const size_t pdg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 1.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh0_0 = buffer.data(sdh0 + 0);
    const auto *sdh0_1 = buffer.data(sdh0 + 1);
    const auto *sdh0_63 = buffer.data(sdh0 + 63);
    const auto *sdh0_66 = buffer.data(sdh0 + 66);
    const auto *sdh0_69 = buffer.data(sdh0 + 69);
    const auto *sdh0_78 = buffer.data(sdh0 + 78);
    const auto *sdh0_80 = buffer.data(sdh0 + 80);
    const auto *sdh0_81 = buffer.data(sdh0 + 81);
    const auto *sdh0_83 = buffer.data(sdh0 + 83);
    const auto *sdh0_99 = buffer.data(sdh0 + 99);
    const auto *sdh0_101 = buffer.data(sdh0 + 101);
    const auto *sdh0_102 = buffer.data(sdh0 + 102);
    const auto *sdh0_104 = buffer.data(sdh0 + 104);
    const auto *sdh0_105 = buffer.data(sdh0 + 105);
    const auto *sdh0_110 = buffer.data(sdh0 + 110);
    const auto *sdh0_114 = buffer.data(sdh0 + 114);
    const auto *sdh0_120 = buffer.data(sdh0 + 120);
    const auto *sdh0_122 = buffer.data(sdh0 + 122);
    const auto *sdh0_123 = buffer.data(sdh0 + 123);
    const auto *sdh0_125 = buffer.data(sdh0 + 125);

    const auto *sdg_0 = buffer.data(sdg + 0);
    const auto *sdg_10 = buffer.data(sdg + 10);
    const auto *sdg_12 = buffer.data(sdg + 12);
    const auto *sdg_14 = buffer.data(sdg + 14);
    const auto *sdg_25 = buffer.data(sdg + 25);
    const auto *sdg_27 = buffer.data(sdg + 27);
    const auto *sdg_42 = buffer.data(sdg + 42);
    const auto *sdg_44 = buffer.data(sdg + 44);
    const auto *sdg_45 = buffer.data(sdg + 45);
    const auto *sdg_48 = buffer.data(sdg + 48);
    const auto *sdg_51 = buffer.data(sdg + 51);
    const auto *sdg_55 = buffer.data(sdg + 55);
    const auto *sdg_57 = buffer.data(sdg + 57);
    const auto *sdg_59 = buffer.data(sdg + 59);
    const auto *sdg_70 = buffer.data(sdg + 70);
    const auto *sdg_72 = buffer.data(sdg + 72);
    const auto *sdg_74 = buffer.data(sdg + 74);
    const auto *sdg_75 = buffer.data(sdg + 75);
    const auto *sdg_80 = buffer.data(sdg + 80);
    const auto *sdg_84 = buffer.data(sdg + 84);
    const auto *sdg_85 = buffer.data(sdg + 85);
    const auto *sdg_87 = buffer.data(sdg + 87);
    const auto *sdg_89 = buffer.data(sdg + 89);

    const auto *sdh1_0 = buffer.data(sdh1 + 0);
    const auto *sdh1_1 = buffer.data(sdh1 + 1);
    const auto *sdh1_63 = buffer.data(sdh1 + 63);
    const auto *sdh1_66 = buffer.data(sdh1 + 66);
    const auto *sdh1_69 = buffer.data(sdh1 + 69);
    const auto *sdh1_78 = buffer.data(sdh1 + 78);
    const auto *sdh1_80 = buffer.data(sdh1 + 80);
    const auto *sdh1_81 = buffer.data(sdh1 + 81);
    const auto *sdh1_83 = buffer.data(sdh1 + 83);
    const auto *sdh1_99 = buffer.data(sdh1 + 99);
    const auto *sdh1_101 = buffer.data(sdh1 + 101);
    const auto *sdh1_102 = buffer.data(sdh1 + 102);
    const auto *sdh1_104 = buffer.data(sdh1 + 104);
    const auto *sdh1_105 = buffer.data(sdh1 + 105);
    const auto *sdh1_110 = buffer.data(sdh1 + 110);
    const auto *sdh1_114 = buffer.data(sdh1 + 114);
    const auto *sdh1_120 = buffer.data(sdh1 + 120);
    const auto *sdh1_122 = buffer.data(sdh1 + 122);
    const auto *sdh1_123 = buffer.data(sdh1 + 123);
    const auto *sdh1_125 = buffer.data(sdh1 + 125);

    const auto *pph0_0 = buffer.data(pph0 + 0);
    const auto *pph0_3 = buffer.data(pph0 + 3);
    const auto *pph0_5 = buffer.data(pph0 + 5);
    const auto *pph0_6 = buffer.data(pph0 + 6);
    const auto *pph0_9 = buffer.data(pph0 + 9);
    const auto *pph0_10 = buffer.data(pph0 + 10);
    const auto *pph0_14 = buffer.data(pph0 + 14);
    const auto *pph0_24 = buffer.data(pph0 + 24);
    const auto *pph0_27 = buffer.data(pph0 + 27);
    const auto *pph0_42 = buffer.data(pph0 + 42);
    const auto *pph0_47 = buffer.data(pph0 + 47);
    const auto *pph0_51 = buffer.data(pph0 + 51);

    const auto *ppg_0 = buffer.data(ppg + 0);
    const auto *ppg_1 = buffer.data(ppg + 1);
    const auto *ppg_2 = buffer.data(ppg + 2);
    const auto *ppg_3 = buffer.data(ppg + 3);
    const auto *ppg_5 = buffer.data(ppg + 5);
    const auto *ppg_6 = buffer.data(ppg + 6);
    const auto *ppg_9 = buffer.data(ppg + 9);
    const auto *ppg_10 = buffer.data(ppg + 10);
    const auto *ppg_12 = buffer.data(ppg + 12);
    const auto *ppg_13 = buffer.data(ppg + 13);
    const auto *ppg_14 = buffer.data(ppg + 14);
    const auto *ppg_15 = buffer.data(ppg + 15);
    const auto *ppg_17 = buffer.data(ppg + 17);
    const auto *ppg_18 = buffer.data(ppg + 18);
    const auto *ppg_20 = buffer.data(ppg + 20);
    const auto *ppg_21 = buffer.data(ppg + 21);
    const auto *ppg_24 = buffer.data(ppg + 24);
    const auto *ppg_25 = buffer.data(ppg + 25);
    const auto *ppg_27 = buffer.data(ppg + 27);
    const auto *ppg_29 = buffer.data(ppg + 29);
    const auto *ppg_30 = buffer.data(ppg + 30);
    const auto *ppg_32 = buffer.data(ppg + 32);
    const auto *ppg_33 = buffer.data(ppg + 33);
    const auto *ppg_35 = buffer.data(ppg + 35);
    const auto *ppg_36 = buffer.data(ppg + 36);
    const auto *ppg_39 = buffer.data(ppg + 39);
    const auto *ppg_40 = buffer.data(ppg + 40);
    const auto *ppg_42 = buffer.data(ppg + 42);
    const auto *ppg_44 = buffer.data(ppg + 44);

    const auto *pph1_0 = buffer.data(pph1 + 0);
    const auto *pph1_3 = buffer.data(pph1 + 3);
    const auto *pph1_5 = buffer.data(pph1 + 5);
    const auto *pph1_6 = buffer.data(pph1 + 6);
    const auto *pph1_9 = buffer.data(pph1 + 9);
    const auto *pph1_10 = buffer.data(pph1 + 10);
    const auto *pph1_14 = buffer.data(pph1 + 14);
    const auto *pph1_24 = buffer.data(pph1 + 24);
    const auto *pph1_27 = buffer.data(pph1 + 27);
    const auto *pph1_42 = buffer.data(pph1 + 42);
    const auto *pph1_47 = buffer.data(pph1 + 47);
    const auto *pph1_51 = buffer.data(pph1 + 51);

    const auto *pdf0_0 = buffer.data(pdf0 + 0);
    const auto *pdf0_1 = buffer.data(pdf0 + 1);
    const auto *pdf0_2 = buffer.data(pdf0 + 2);
    const auto *pdf0_6 = buffer.data(pdf0 + 6);
    const auto *pdf0_8 = buffer.data(pdf0 + 8);
    const auto *pdf0_9 = buffer.data(pdf0 + 9);
    const auto *pdf0_10 = buffer.data(pdf0 + 10);
    const auto *pdf0_11 = buffer.data(pdf0 + 11);
    const auto *pdf0_16 = buffer.data(pdf0 + 16);
    const auto *pdf0_18 = buffer.data(pdf0 + 18);
    const auto *pdf0_19 = buffer.data(pdf0 + 19);
    const auto *pdf0_20 = buffer.data(pdf0 + 20);
    const auto *pdf0_22 = buffer.data(pdf0 + 22);
    const auto *pdf0_26 = buffer.data(pdf0 + 26);
    const auto *pdf0_28 = buffer.data(pdf0 + 28);
    const auto *pdf0_29 = buffer.data(pdf0 + 29);
    const auto *pdf0_30 = buffer.data(pdf0 + 30);
    const auto *pdf0_32 = buffer.data(pdf0 + 32);
    const auto *pdf0_50 = buffer.data(pdf0 + 50);
    const auto *pdf0_51 = buffer.data(pdf0 + 51);

    const auto *pdf1_0 = buffer.data(pdf1 + 0);
    const auto *pdf1_1 = buffer.data(pdf1 + 1);
    const auto *pdf1_2 = buffer.data(pdf1 + 2);
    const auto *pdf1_6 = buffer.data(pdf1 + 6);
    const auto *pdf1_8 = buffer.data(pdf1 + 8);
    const auto *pdf1_9 = buffer.data(pdf1 + 9);
    const auto *pdf1_10 = buffer.data(pdf1 + 10);
    const auto *pdf1_11 = buffer.data(pdf1 + 11);
    const auto *pdf1_16 = buffer.data(pdf1 + 16);
    const auto *pdf1_18 = buffer.data(pdf1 + 18);
    const auto *pdf1_19 = buffer.data(pdf1 + 19);
    const auto *pdf1_20 = buffer.data(pdf1 + 20);
    const auto *pdf1_22 = buffer.data(pdf1 + 22);
    const auto *pdf1_26 = buffer.data(pdf1 + 26);
    const auto *pdf1_28 = buffer.data(pdf1 + 28);
    const auto *pdf1_29 = buffer.data(pdf1 + 29);
    const auto *pdf1_30 = buffer.data(pdf1 + 30);
    const auto *pdf1_32 = buffer.data(pdf1 + 32);
    const auto *pdf1_50 = buffer.data(pdf1 + 50);
    const auto *pdf1_51 = buffer.data(pdf1 + 51);

    const auto *pdg_0 = buffer.data(pdg + 0);
    const auto *pdg_1 = buffer.data(pdg + 1);
    const auto *pdg_2 = buffer.data(pdg + 2);
    const auto *pdg_3 = buffer.data(pdg + 3);
    const auto *pdg_5 = buffer.data(pdg + 5);
    const auto *pdg_6 = buffer.data(pdg + 6);
    const auto *pdg_9 = buffer.data(pdg + 9);
    const auto *pdg_10 = buffer.data(pdg + 10);
    const auto *pdg_12 = buffer.data(pdg + 12);
    const auto *pdg_13 = buffer.data(pdg + 13);
    const auto *pdg_14 = buffer.data(pdg + 14);
    const auto *pdg_15 = buffer.data(pdg + 15);
    const auto *pdg_16 = buffer.data(pdg + 16);
    const auto *pdg_17 = buffer.data(pdg + 17);
    const auto *pdg_18 = buffer.data(pdg + 18);
    const auto *pdg_20 = buffer.data(pdg + 20);
    const auto *pdg_21 = buffer.data(pdg + 21);
    const auto *pdg_24 = buffer.data(pdg + 24);
    const auto *pdg_25 = buffer.data(pdg + 25);
    const auto *pdg_27 = buffer.data(pdg + 27);
    const auto *pdg_28 = buffer.data(pdg + 28);
    const auto *pdg_29 = buffer.data(pdg + 29);
    const auto *pdg_30 = buffer.data(pdg + 30);
    const auto *pdg_32 = buffer.data(pdg + 32);
    const auto *pdg_33 = buffer.data(pdg + 33);
    const auto *pdg_35 = buffer.data(pdg + 35);
    const auto *pdg_36 = buffer.data(pdg + 36);
    const auto *pdg_39 = buffer.data(pdg + 39);
    const auto *pdg_40 = buffer.data(pdg + 40);
    const auto *pdg_42 = buffer.data(pdg + 42);
    const auto *pdg_43 = buffer.data(pdg + 43);
    const auto *pdg_44 = buffer.data(pdg + 44);
    const auto *pdg_45 = buffer.data(pdg + 45);
    const auto *pdg_47 = buffer.data(pdg + 47);
    const auto *pdg_48 = buffer.data(pdg + 48);
    const auto *pdg_50 = buffer.data(pdg + 50);
    const auto *pdg_51 = buffer.data(pdg + 51);
    const auto *pdg_54 = buffer.data(pdg + 54);
    const auto *pdg_55 = buffer.data(pdg + 55);
    const auto *pdg_57 = buffer.data(pdg + 57);
    const auto *pdg_59 = buffer.data(pdg + 59);
    const auto *pdg_60 = buffer.data(pdg + 60);
    const auto *pdg_62 = buffer.data(pdg + 62);
    const auto *pdg_63 = buffer.data(pdg + 63);
    const auto *pdg_65 = buffer.data(pdg + 65);
    const auto *pdg_66 = buffer.data(pdg + 66);
    const auto *pdg_69 = buffer.data(pdg + 69);
    const auto *pdg_70 = buffer.data(pdg + 70);
    const auto *pdg_72 = buffer.data(pdg + 72);
    const auto *pdg_74 = buffer.data(pdg + 74);
    const auto *pdg_75 = buffer.data(pdg + 75);
    const auto *pdg_76 = buffer.data(pdg + 76);
    const auto *pdg_77 = buffer.data(pdg + 77);
    const auto *pdg_78 = buffer.data(pdg + 78);
    const auto *pdg_80 = buffer.data(pdg + 80);
    const auto *pdg_81 = buffer.data(pdg + 81);
    const auto *pdg_84 = buffer.data(pdg + 84);
    const auto *pdg_85 = buffer.data(pdg + 85);
    const auto *pdg_87 = buffer.data(pdg + 87);
    const auto *pdg_89 = buffer.data(pdg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sdg_0, ppg_0, pdf0_0, \
                         pdf1_0, pdg_0, pdg_1, pdg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdg_0[k]
                 + f_1 * ppg_0[k]
                 + f_2 * pdf0_0[k]
                 - f_3 * pdf1_0[k]
                 + f_4 * pc_x[k] * pdg_0[k];

        t_1[k] = f_4 * pc_y[k] * pdg_0[k];

        t_2[k] = f_4 * pc_z[k] * pdg_0[k];

        t_3[k] = f_5 * pdf0_0[k]
                 - f_6 * pdf1_0[k]
                 + f_4 * pc_y[k] * pdg_1[k];

        t_4[k] = f_4 * pc_y[k] * pdg_2[k];

        t_5[k] = f_5 * pdf0_0[k]
                 - f_6 * pdf1_0[k]
                 + f_4 * pc_z[k] * pdg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, pdf0_1, pdf0_2, pdf1_1, pdf1_2, \
                         pdg_3, pdg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * pdf0_1[k]
                 - f_8 * pdf1_1[k]
                 + f_4 * pc_y[k] * pdg_3[k];

        t_7[k] = f_4 * pc_z[k] * pdg_3[k];

        t_8[k] = f_4 * pc_y[k] * pdg_5[k];

        t_9[k] = f_7 * pdf0_2[k]
                 - f_8 * pdf1_2[k]
                 + f_4 * pc_z[k] * pdg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, sdg_10, sdg_12, ppg_10, \
                         ppg_12, pdg_6, pdg_9, pdg_10, pdg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sdg_10[k]
                  + f_1 * ppg_10[k]
                  + f_4 * pc_x[k] * pdg_10[k];

        t_11[k] = f_4 * pc_z[k] * pdg_6[k];

        t_12[k] = f_0 * sdg_12[k]
                  + f_1 * ppg_12[k]
                  + f_4 * pc_x[k] * pdg_12[k];

        t_13[k] = f_4 * pc_y[k] * pdg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sdg_14, ppg_14, pdf0_6, \
                         pdf0_8, pdf1_6, pdf1_8, pdg_10, pdg_12, \
                         pdg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * sdg_14[k]
                  + f_1 * ppg_14[k]
                  + f_4 * pc_x[k] * pdg_14[k];

        t_15[k] = f_2 * pdf0_6[k]
                  - f_3 * pdf1_6[k]
                  + f_4 * pc_y[k] * pdg_10[k];

        t_16[k] = f_4 * pc_z[k] * pdg_10[k];

        t_17[k] = f_7 * pdf0_8[k]
                  - f_8 * pdf1_8[k]
                  + f_4 * pc_y[k] * pdg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, pph0_0, ppg_0, \
                         pph1_0, pdf0_9, pdf1_9, pdg_13, pdg_14, \
                         pdg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * pdf0_9[k]
                  - f_6 * pdf1_9[k]
                  + f_4 * pc_y[k] * pdg_13[k];

        t_19[k] = f_4 * pc_y[k] * pdg_14[k];

        t_20[k] = f_2 * pdf0_9[k]
                  - f_3 * pdf1_9[k]
                  + f_4 * pc_z[k] * pdg_14[k];

        t_21[k] = pb_y[k] * pph0_0[k]
                  - f_9 * pc_y[k] * pph1_0[k];

        t_22[k] = f_0 * ppg_0[k]
                  + f_4 * pc_y[k] * pdg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, pph0_5, ppg_1, ppg_2, \
                         pph1_5, pdf0_10, pdf1_10, pdg_15, pdg_16, \
                         pdg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * pc_z[k] * pdg_15[k];

        t_24[k] = f_0 * ppg_1[k]
                  + f_5 * pdf0_10[k]
                  - f_6 * pdf1_10[k]
                  + f_4 * pc_y[k] * pdg_16[k];

        t_25[k] = f_0 * ppg_2[k]
                  + f_4 * pc_y[k] * pdg_17[k];

        t_26[k] = pb_y[k] * pph0_5[k]
                  - f_9 * pc_y[k] * pph1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_y, pc_z, pph0_9, ppg_3, ppg_5, \
                         pph1_9, pdf0_11, pdf1_11, pdg_18, pdg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * ppg_3[k]
                  + f_7 * pdf0_11[k]
                  - f_8 * pdf1_11[k]
                  + f_4 * pc_y[k] * pdg_18[k];

        t_28[k] = f_4 * pc_z[k] * pdg_18[k];

        t_29[k] = f_0 * ppg_5[k]
                  + f_4 * pc_y[k] * pdg_20[k];

        t_30[k] = pb_y[k] * pph0_9[k]
                  - f_9 * pc_y[k] * pph1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, sdg_25, sdg_27, ppg_9, \
                         ppg_25, ppg_27, pdg_21, pdg_24, pdg_25, \
                         pdg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * sdg_25[k]
                  + f_0 * ppg_25[k]
                  + f_4 * pc_x[k] * pdg_25[k];

        t_32[k] = f_4 * pc_z[k] * pdg_21[k];

        t_33[k] = f_0 * sdg_27[k]
                  + f_0 * ppg_27[k]
                  + f_4 * pc_x[k] * pdg_27[k];

        t_34[k] = f_0 * ppg_9[k]
                  + f_4 * pc_y[k] * pdg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_y, pc_z, pph0_14, ppg_10, pph1_14, \
                         pdf0_16, pdf1_16, pdg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_y[k] * pph0_14[k]
                  - f_9 * pc_y[k] * pph1_14[k];

        t_36[k] = f_0 * ppg_10[k]
                  + f_2 * pdf0_16[k]
                  - f_3 * pdf1_16[k]
                  + f_4 * pc_y[k] * pdg_25[k];

        t_37[k] = f_4 * pc_z[k] * pdg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_y, pc_z, ppg_12, ppg_13, ppg_14, pdf0_18, \
                         pdf0_19, pdf1_18, pdf1_19, pdg_27, pdg_28, \
                         pdg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * ppg_12[k]
                  + f_7 * pdf0_18[k]
                  - f_8 * pdf1_18[k]
                  + f_4 * pc_y[k] * pdg_27[k];

        t_39[k] = f_0 * ppg_13[k]
                  + f_5 * pdf0_19[k]
                  - f_6 * pdf1_19[k]
                  + f_4 * pc_y[k] * pdg_28[k];

        t_40[k] = f_0 * ppg_14[k]
                  + f_4 * pc_y[k] * pdg_29[k];

        t_41[k] = f_2 * pdf0_19[k]
                  - f_3 * pdf1_19[k]
                  + f_4 * pc_z[k] * pdg_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, pph0_0, pph0_3, \
                         ppg_0, pph1_0, pph1_3, pdg_30, pdg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * pph0_0[k]
                  - f_9 * pc_z[k] * pph1_0[k];

        t_43[k] = f_4 * pc_y[k] * pdg_30[k];

        t_44[k] = f_0 * ppg_0[k]
                  + f_4 * pc_z[k] * pdg_30[k];

        t_45[k] = pb_z[k] * pph0_3[k]
                  - f_9 * pc_z[k] * pph1_3[k];

        t_46[k] = f_4 * pc_y[k] * pdg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, pph0_6, ppg_2, ppg_3, \
                         pph1_6, pdf0_20, pdf1_20, pdg_32, pdg_33, \
                         pdg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * ppg_2[k]
                  + f_5 * pdf0_20[k]
                  - f_6 * pdf1_20[k]
                  + f_4 * pc_z[k] * pdg_32[k];

        t_48[k] = pb_z[k] * pph0_6[k]
                  - f_9 * pc_z[k] * pph1_6[k];

        t_49[k] = f_0 * ppg_3[k]
                  + f_4 * pc_z[k] * pdg_33[k];

        t_50[k] = f_4 * pc_y[k] * pdg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_z, pc_z, pph0_10, ppg_5, ppg_6, pph1_10, \
                         pdf0_22, pdf1_22, pdg_35, pdg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * ppg_5[k]
                  + f_7 * pdf0_22[k]
                  - f_8 * pdf1_22[k]
                  + f_4 * pc_z[k] * pdg_35[k];

        t_52[k] = pb_z[k] * pph0_10[k]
                  - f_9 * pc_z[k] * pph1_10[k];

        t_53[k] = f_0 * ppg_6[k]
                  + f_4 * pc_z[k] * pdg_36[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pc_x, pc_y, sdg_42, sdg_44, ppg_42, ppg_44, \
                         pdf0_26, pdf1_26, pdg_39, pdg_40, pdg_42, \
                         pdg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * sdg_42[k]
                  + f_0 * ppg_42[k]
                  + f_4 * pc_x[k] * pdg_42[k];

        t_55[k] = f_4 * pc_y[k] * pdg_39[k];

        t_56[k] = f_0 * sdg_44[k]
                  + f_0 * ppg_44[k]
                  + f_4 * pc_x[k] * pdg_44[k];

        t_57[k] = f_2 * pdf0_26[k]
                  - f_3 * pdf1_26[k]
                  + f_4 * pc_y[k] * pdg_40[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pc_y, pc_z, ppg_10, pdf0_28, pdf0_29, \
                         pdf1_28, pdf1_29, pdg_40, pdg_42, pdg_43, \
                         pdg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * ppg_10[k]
                  + f_4 * pc_z[k] * pdg_40[k];

        t_59[k] = f_7 * pdf0_28[k]
                  - f_8 * pdf1_28[k]
                  + f_4 * pc_y[k] * pdg_42[k];

        t_60[k] = f_5 * pdf0_29[k]
                  - f_6 * pdf1_29[k]
                  + f_4 * pc_y[k] * pdg_43[k];

        t_61[k] = f_4 * pc_y[k] * pdg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_x, pc_x, pc_y, pc_z, sdh0_63, sdg_45, sdh1_63, \
                         ppg_14, ppg_15, pdf0_29, pdf1_29, pdg_44, \
                         pdg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_0 * ppg_14[k]
                  + f_2 * pdf0_29[k]
                  - f_3 * pdf1_29[k]
                  + f_4 * pc_z[k] * pdg_44[k];

        t_63[k] = pa_x[k] * sdh0_63[k]
                  + f_10 * sdg_45[k]
                  - f_9 * pc_x[k] * sdh1_63[k];

        t_64[k] = f_1 * ppg_15[k]
                  + f_4 * pc_y[k] * pdg_45[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pc_x, pc_y, pc_z, sdh0_66, sdg_48, \
                         sdh1_66, ppg_17, pdf0_30, pdf1_30, pdg_45, \
                         pdg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_4 * pc_z[k] * pdg_45[k];

        t_66[k] = pa_x[k] * sdh0_66[k]
                  + f_11 * sdg_48[k]
                  - f_9 * pc_x[k] * sdh1_66[k];

        t_67[k] = f_1 * ppg_17[k]
                  + f_4 * pc_y[k] * pdg_47[k];

        t_68[k] = f_5 * pdf0_30[k]
                  - f_6 * pdf1_30[k]
                  + f_4 * pc_z[k] * pdg_47[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pc_x, pc_y, pc_z, sdh0_69, sdg_51, \
                         sdh1_69, ppg_20, pdf0_32, pdf1_32, pdg_48, \
                         pdg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_x[k] * sdh0_69[k]
                  + f_1 * sdg_51[k]
                  - f_9 * pc_x[k] * sdh1_69[k];

        t_70[k] = f_4 * pc_z[k] * pdg_48[k];

        t_71[k] = f_1 * ppg_20[k]
                  + f_4 * pc_y[k] * pdg_50[k];

        t_72[k] = f_7 * pdf0_32[k]
                  - f_8 * pdf1_32[k]
                  + f_4 * pc_z[k] * pdg_50[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, pc_y, pc_z, sdg_55, sdg_57, ppg_24, \
                         pdg_51, pdg_54, pdg_55, pdg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * sdg_55[k]
                  + f_4 * pc_x[k] * pdg_55[k];

        t_74[k] = f_4 * pc_z[k] * pdg_51[k];

        t_75[k] = f_0 * sdg_57[k]
                  + f_4 * pc_x[k] * pdg_57[k];

        t_76[k] = f_1 * ppg_24[k]
                  + f_4 * pc_y[k] * pdg_54[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pc_x, pc_z, sdh0_78, sdh0_80, sdg_59, \
                         sdh1_78, sdh1_80, pdg_55, pdg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * sdg_59[k]
                  + f_4 * pc_x[k] * pdg_59[k];

        t_78[k] = pa_x[k] * sdh0_78[k]
                  - f_9 * pc_x[k] * sdh1_78[k];

        t_79[k] = f_4 * pc_z[k] * pdg_55[k];

        t_80[k] = pa_x[k] * sdh0_80[k]
                  - f_9 * pc_x[k] * sdh1_80[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_x, pb_y, pc_x, pc_y, sdh0_81, sdh0_83, \
                         sdh1_81, sdh1_83, pph0_42, ppg_29, pph1_42, \
                         pdg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = pa_x[k] * sdh0_81[k]
                  - f_9 * pc_x[k] * sdh1_81[k];

        t_82[k] = f_1 * ppg_29[k]
                  + f_4 * pc_y[k] * pdg_59[k];

        t_83[k] = pa_x[k] * sdh0_83[k]
                  - f_9 * pc_x[k] * sdh1_83[k];

        t_84[k] = pb_y[k] * pph0_42[k]
                  - f_9 * pc_y[k] * pph1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, pph0_24, ppg_15, ppg_30, \
                         ppg_32, pph1_24, pdg_60, pdg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_0 * ppg_30[k]
                  + f_4 * pc_y[k] * pdg_60[k];

        t_86[k] = f_0 * ppg_15[k]
                  + f_4 * pc_z[k] * pdg_60[k];

        t_87[k] = pb_z[k] * pph0_24[k]
                  - f_9 * pc_z[k] * pph1_24[k];

        t_88[k] = f_0 * ppg_32[k]
                  + f_4 * pc_y[k] * pdg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, pph0_27, pph0_47, \
                         ppg_18, ppg_35, pph1_27, pph1_47, pdg_63, \
                         pdg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * pph0_47[k]
                  - f_9 * pc_y[k] * pph1_47[k];

        t_90[k] = pb_z[k] * pph0_27[k]
                  - f_9 * pc_z[k] * pph1_27[k];

        t_91[k] = f_0 * ppg_18[k]
                  + f_4 * pc_z[k] * pdg_63[k];

        t_92[k] = f_0 * ppg_35[k]
                  + f_4 * pc_y[k] * pdg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, pc_z, sdg_70, sdg_72, \
                         pph0_51, ppg_21, pph1_51, pdg_66, pdg_70, \
                         pdg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * pph0_51[k]
                  - f_9 * pc_y[k] * pph1_51[k];

        t_94[k] = f_0 * sdg_70[k]
                  + f_4 * pc_x[k] * pdg_70[k];

        t_95[k] = f_0 * ppg_21[k]
                  + f_4 * pc_z[k] * pdg_66[k];

        t_96[k] = f_0 * sdg_72[k]
                  + f_4 * pc_x[k] * pdg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pc_x, pc_y, pc_z, sdh0_99, sdg_74, \
                         sdh1_99, ppg_25, ppg_39, pdg_69, pdg_70, \
                         pdg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_0 * ppg_39[k]
                  + f_4 * pc_y[k] * pdg_69[k];

        t_98[k] = f_0 * sdg_74[k]
                  + f_4 * pc_x[k] * pdg_74[k];

        t_99[k] = pa_x[k] * sdh0_99[k]
                  - f_9 * pc_x[k] * sdh1_99[k];

        t_100[k] = f_0 * ppg_25[k]
                   + f_4 * pc_z[k] * pdg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pc_x, pc_y, sdh0_101, sdh0_102, \
                         sdh0_104, sdh1_101, sdh1_102, sdh1_104, ppg_44, \
                         pdg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pa_x[k] * sdh0_101[k]
                   - f_9 * pc_x[k] * sdh1_101[k];

        t_102[k] = pa_x[k] * sdh0_102[k]
                   - f_9 * pc_x[k] * sdh1_102[k];

        t_103[k] = f_0 * ppg_44[k]
                   + f_4 * pc_y[k] * pdg_74[k];

        t_104[k] = pa_x[k] * sdh0_104[k]
                   - f_9 * pc_x[k] * sdh1_104[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pc_x, pc_y, pc_z, sdh0_105, sdg_75, \
                         sdh1_105, ppg_30, pdf0_50, pdf1_50, pdg_75, \
                         pdg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_x[k] * sdh0_105[k]
                   + f_10 * sdg_75[k]
                   - f_9 * pc_x[k] * sdh1_105[k];

        t_106[k] = f_4 * pc_y[k] * pdg_75[k];

        t_107[k] = f_1 * ppg_30[k]
                   + f_4 * pc_z[k] * pdg_75[k];

        t_108[k] = f_5 * pdf0_50[k]
                   - f_6 * pdf1_50[k]
                   + f_4 * pc_y[k] * pdg_76[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pc_x, pc_y, pc_z, sdh0_110, sdg_80, \
                         sdh1_110, ppg_33, pdf0_51, pdf1_51, pdg_77, \
                         pdg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_4 * pc_y[k] * pdg_77[k];

        t_110[k] = pa_x[k] * sdh0_110[k]
                   + f_11 * sdg_80[k]
                   - f_9 * pc_x[k] * sdh1_110[k];

        t_111[k] = f_7 * pdf0_51[k]
                   - f_8 * pdf1_51[k]
                   + f_4 * pc_y[k] * pdg_78[k];

        t_112[k] = f_1 * ppg_33[k]
                   + f_4 * pc_z[k] * pdg_78[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pc_x, pc_y, pc_z, sdh0_114, sdg_84, \
                         sdg_85, sdh1_114, ppg_36, pdg_80, pdg_81, \
                         pdg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_4 * pc_y[k] * pdg_80[k];

        t_114[k] = pa_x[k] * sdh0_114[k]
                   + f_1 * sdg_84[k]
                   - f_9 * pc_x[k] * sdh1_114[k];

        t_115[k] = f_0 * sdg_85[k]
                   + f_4 * pc_x[k] * pdg_85[k];

        t_116[k] = f_1 * ppg_36[k]
                   + f_4 * pc_z[k] * pdg_81[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pc_x, pc_y, sdh0_120, sdg_87, \
                         sdg_89, sdh1_120, pdg_84, pdg_87, pdg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * sdg_87[k]
                   + f_4 * pc_x[k] * pdg_87[k];

        t_118[k] = f_4 * pc_y[k] * pdg_84[k];

        t_119[k] = f_0 * sdg_89[k]
                   + f_4 * pc_x[k] * pdg_89[k];

        t_120[k] = pa_x[k] * sdh0_120[k]
                   - f_9 * pc_x[k] * sdh1_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pc_x, pc_y, pc_z, sdh0_122, \
                         sdh0_123, sdh1_122, sdh1_123, ppg_40, pdg_85, \
                         pdg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * ppg_40[k]
                   + f_4 * pc_z[k] * pdg_85[k];

        t_122[k] = pa_x[k] * sdh0_122[k]
                   - f_9 * pc_x[k] * sdh1_122[k];

        t_123[k] = pa_x[k] * sdh0_123[k]
                   - f_9 * pc_x[k] * sdh1_123[k];

        t_124[k] = f_4 * pc_y[k] * pdg_89[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pa_y, pc_x, pc_y, sdh0_0, sdh0_1, \
                         sdh0_125, sdg_0, sdh1_0, sdh1_1, sdh1_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_x[k] * sdh0_125[k]
                   - f_9 * pc_x[k] * sdh1_125[k];

        t_126[k] = pa_y[k] * sdh0_0[k]
                   - f_9 * pc_y[k] * sdh1_0[k];

        t_127[k] = pa_y[k] * sdh0_1[k]
                   + f_0 * sdg_0[k]
                   - f_9 * pc_y[k] * sdh1_1[k];
    }
}

static auto
compute_prim_pdh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdh0, const size_t sdg,
                                                          const size_t sdh1, const size_t pph0,
                                                          const size_t ppg, const size_t pph1,
                                                          const size_t pdf0, const size_t pdf1,
                                                          const size_t pdg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);

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
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh0_0 = buffer.data(sdh0 + 0);
    const auto *sdh0_2 = buffer.data(sdh0 + 2);
    const auto *sdh0_3 = buffer.data(sdh0 + 3);
    const auto *sdh0_5 = buffer.data(sdh0 + 5);
    const auto *sdh0_6 = buffer.data(sdh0 + 6);
    const auto *sdh0_8 = buffer.data(sdh0 + 8);
    const auto *sdh0_9 = buffer.data(sdh0 + 9);
    const auto *sdh0_15 = buffer.data(sdh0 + 15);
    const auto *sdh0_20 = buffer.data(sdh0 + 20);
    const auto *sdh0_42 = buffer.data(sdh0 + 42);
    const auto *sdh0_47 = buffer.data(sdh0 + 47);
    const auto *sdh0_50 = buffer.data(sdh0 + 50);
    const auto *sdh0_51 = buffer.data(sdh0 + 51);
    const auto *sdh0_62 = buffer.data(sdh0 + 62);
    const auto *sdh0_105 = buffer.data(sdh0 + 105);
    const auto *sdh0_110 = buffer.data(sdh0 + 110);
    const auto *sdh0_114 = buffer.data(sdh0 + 114);
    const auto *sdh0_120 = buffer.data(sdh0 + 120);
    const auto *sdh0_122 = buffer.data(sdh0 + 122);
    const auto *sdh0_123 = buffer.data(sdh0 + 123);
    const auto *sdh0_125 = buffer.data(sdh0 + 125);

    const auto *sdg_0 = buffer.data(sdg + 0);
    const auto *sdg_1 = buffer.data(sdg + 1);
    const auto *sdg_3 = buffer.data(sdg + 3);
    const auto *sdg_5 = buffer.data(sdg + 5);
    const auto *sdg_10 = buffer.data(sdg + 10);
    const auto *sdg_14 = buffer.data(sdg + 14);
    const auto *sdg_35 = buffer.data(sdg + 35);
    const auto *sdg_44 = buffer.data(sdg + 44);
    const auto *sdg_55 = buffer.data(sdg + 55);
    const auto *sdg_59 = buffer.data(sdg + 59);
    const auto *sdg_74 = buffer.data(sdg + 74);
    const auto *sdg_85 = buffer.data(sdg + 85);
    const auto *sdg_87 = buffer.data(sdg + 87);
    const auto *sdg_88 = buffer.data(sdg + 88);
    const auto *sdg_89 = buffer.data(sdg + 89);

    const auto *sdh1_0 = buffer.data(sdh1 + 0);
    const auto *sdh1_2 = buffer.data(sdh1 + 2);
    const auto *sdh1_3 = buffer.data(sdh1 + 3);
    const auto *sdh1_5 = buffer.data(sdh1 + 5);
    const auto *sdh1_6 = buffer.data(sdh1 + 6);
    const auto *sdh1_8 = buffer.data(sdh1 + 8);
    const auto *sdh1_9 = buffer.data(sdh1 + 9);
    const auto *sdh1_15 = buffer.data(sdh1 + 15);
    const auto *sdh1_20 = buffer.data(sdh1 + 20);
    const auto *sdh1_42 = buffer.data(sdh1 + 42);
    const auto *sdh1_47 = buffer.data(sdh1 + 47);
    const auto *sdh1_50 = buffer.data(sdh1 + 50);
    const auto *sdh1_51 = buffer.data(sdh1 + 51);
    const auto *sdh1_62 = buffer.data(sdh1 + 62);
    const auto *sdh1_105 = buffer.data(sdh1 + 105);
    const auto *sdh1_110 = buffer.data(sdh1 + 110);
    const auto *sdh1_114 = buffer.data(sdh1 + 114);
    const auto *sdh1_120 = buffer.data(sdh1 + 120);
    const auto *sdh1_122 = buffer.data(sdh1 + 122);
    const auto *sdh1_123 = buffer.data(sdh1 + 123);
    const auto *sdh1_125 = buffer.data(sdh1 + 125);

    const auto *pph0_64 = buffer.data(pph0 + 64);
    const auto *pph0_66 = buffer.data(pph0 + 66);
    const auto *pph0_69 = buffer.data(pph0 + 69);
    const auto *pph0_85 = buffer.data(pph0 + 85);
    const auto *pph0_87 = buffer.data(pph0 + 87);
    const auto *pph0_90 = buffer.data(pph0 + 90);
    const auto *pph0_92 = buffer.data(pph0 + 92);
    const auto *pph0_99 = buffer.data(pph0 + 99);
    const auto *pph0_101 = buffer.data(pph0 + 101);
    const auto *pph0_102 = buffer.data(pph0 + 102);
    const auto *pph0_103 = buffer.data(pph0 + 103);
    const auto *pph0_104 = buffer.data(pph0 + 104);
    const auto *pph0_120 = buffer.data(pph0 + 120);
    const auto *pph0_122 = buffer.data(pph0 + 122);
    const auto *pph0_123 = buffer.data(pph0 + 123);

    const auto *ppg_45 = buffer.data(ppg + 45);
    const auto *ppg_46 = buffer.data(ppg + 46);
    const auto *ppg_48 = buffer.data(ppg + 48);
    const auto *ppg_55 = buffer.data(ppg + 55);
    const auto *ppg_56 = buffer.data(ppg + 56);
    const auto *ppg_57 = buffer.data(ppg + 57);
    const auto *ppg_58 = buffer.data(ppg + 58);
    const auto *ppg_59 = buffer.data(ppg + 59);
    const auto *ppg_60 = buffer.data(ppg + 60);
    const auto *ppg_61 = buffer.data(ppg + 61);
    const auto *ppg_63 = buffer.data(ppg + 63);
    const auto *ppg_65 = buffer.data(ppg + 65);
    const auto *ppg_66 = buffer.data(ppg + 66);
    const auto *ppg_68 = buffer.data(ppg + 68);
    const auto *ppg_69 = buffer.data(ppg + 69);
    const auto *ppg_70 = buffer.data(ppg + 70);
    const auto *ppg_71 = buffer.data(ppg + 71);
    const auto *ppg_72 = buffer.data(ppg + 72);
    const auto *ppg_73 = buffer.data(ppg + 73);
    const auto *ppg_74 = buffer.data(ppg + 74);
    const auto *ppg_75 = buffer.data(ppg + 75);
    const auto *ppg_76 = buffer.data(ppg + 76);
    const auto *ppg_78 = buffer.data(ppg + 78);
    const auto *ppg_85 = buffer.data(ppg + 85);
    const auto *ppg_86 = buffer.data(ppg + 86);
    const auto *ppg_87 = buffer.data(ppg + 87);
    const auto *ppg_88 = buffer.data(ppg + 88);
    const auto *ppg_89 = buffer.data(ppg + 89);

    const auto *pph1_64 = buffer.data(pph1 + 64);
    const auto *pph1_66 = buffer.data(pph1 + 66);
    const auto *pph1_69 = buffer.data(pph1 + 69);
    const auto *pph1_85 = buffer.data(pph1 + 85);
    const auto *pph1_87 = buffer.data(pph1 + 87);
    const auto *pph1_90 = buffer.data(pph1 + 90);
    const auto *pph1_92 = buffer.data(pph1 + 92);
    const auto *pph1_99 = buffer.data(pph1 + 99);
    const auto *pph1_101 = buffer.data(pph1 + 101);
    const auto *pph1_102 = buffer.data(pph1 + 102);
    const auto *pph1_103 = buffer.data(pph1 + 103);
    const auto *pph1_104 = buffer.data(pph1 + 104);
    const auto *pph1_120 = buffer.data(pph1 + 120);
    const auto *pph1_122 = buffer.data(pph1 + 122);
    const auto *pph1_123 = buffer.data(pph1 + 123);

    const auto *pdf0_66 = buffer.data(pdf0 + 66);
    const auto *pdf0_67 = buffer.data(pdf0 + 67);
    const auto *pdf0_70 = buffer.data(pdf0 + 70);
    const auto *pdf0_75 = buffer.data(pdf0 + 75);
    const auto *pdf0_79 = buffer.data(pdf0 + 79);
    const auto *pdf0_90 = buffer.data(pdf0 + 90);
    const auto *pdf0_91 = buffer.data(pdf0 + 91);
    const auto *pdf0_93 = buffer.data(pdf0 + 93);
    const auto *pdf0_95 = buffer.data(pdf0 + 95);
    const auto *pdf0_96 = buffer.data(pdf0 + 96);
    const auto *pdf0_97 = buffer.data(pdf0 + 97);
    const auto *pdf0_98 = buffer.data(pdf0 + 98);
    const auto *pdf0_99 = buffer.data(pdf0 + 99);
    const auto *pdf0_100 = buffer.data(pdf0 + 100);
    const auto *pdf0_105 = buffer.data(pdf0 + 105);
    const auto *pdf0_106 = buffer.data(pdf0 + 106);
    const auto *pdf0_107 = buffer.data(pdf0 + 107);
    const auto *pdf0_108 = buffer.data(pdf0 + 108);
    const auto *pdf0_109 = buffer.data(pdf0 + 109);
    const auto *pdf0_111 = buffer.data(pdf0 + 111);
    const auto *pdf0_113 = buffer.data(pdf0 + 113);
    const auto *pdf0_116 = buffer.data(pdf0 + 116);
    const auto *pdf0_118 = buffer.data(pdf0 + 118);

    const auto *pdf1_66 = buffer.data(pdf1 + 66);
    const auto *pdf1_67 = buffer.data(pdf1 + 67);
    const auto *pdf1_70 = buffer.data(pdf1 + 70);
    const auto *pdf1_75 = buffer.data(pdf1 + 75);
    const auto *pdf1_79 = buffer.data(pdf1 + 79);
    const auto *pdf1_90 = buffer.data(pdf1 + 90);
    const auto *pdf1_91 = buffer.data(pdf1 + 91);
    const auto *pdf1_93 = buffer.data(pdf1 + 93);
    const auto *pdf1_95 = buffer.data(pdf1 + 95);
    const auto *pdf1_96 = buffer.data(pdf1 + 96);
    const auto *pdf1_97 = buffer.data(pdf1 + 97);
    const auto *pdf1_98 = buffer.data(pdf1 + 98);
    const auto *pdf1_99 = buffer.data(pdf1 + 99);
    const auto *pdf1_100 = buffer.data(pdf1 + 100);
    const auto *pdf1_105 = buffer.data(pdf1 + 105);
    const auto *pdf1_106 = buffer.data(pdf1 + 106);
    const auto *pdf1_107 = buffer.data(pdf1 + 107);
    const auto *pdf1_108 = buffer.data(pdf1 + 108);
    const auto *pdf1_109 = buffer.data(pdf1 + 109);
    const auto *pdf1_111 = buffer.data(pdf1 + 111);
    const auto *pdf1_113 = buffer.data(pdf1 + 113);
    const auto *pdf1_116 = buffer.data(pdf1 + 116);
    const auto *pdf1_118 = buffer.data(pdf1 + 118);

    const auto *pdg_90 = buffer.data(pdg + 90);
    const auto *pdg_91 = buffer.data(pdg + 91);
    const auto *pdg_93 = buffer.data(pdg + 93);
    const auto *pdg_100 = buffer.data(pdg + 100);
    const auto *pdg_101 = buffer.data(pdg + 101);
    const auto *pdg_102 = buffer.data(pdg + 102);
    const auto *pdg_103 = buffer.data(pdg + 103);
    const auto *pdg_104 = buffer.data(pdg + 104);
    const auto *pdg_105 = buffer.data(pdg + 105);
    const auto *pdg_106 = buffer.data(pdg + 106);
    const auto *pdg_108 = buffer.data(pdg + 108);
    const auto *pdg_110 = buffer.data(pdg + 110);
    const auto *pdg_114 = buffer.data(pdg + 114);
    const auto *pdg_115 = buffer.data(pdg + 115);
    const auto *pdg_116 = buffer.data(pdg + 116);
    const auto *pdg_117 = buffer.data(pdg + 117);
    const auto *pdg_118 = buffer.data(pdg + 118);
    const auto *pdg_119 = buffer.data(pdg + 119);
    const auto *pdg_120 = buffer.data(pdg + 120);
    const auto *pdg_121 = buffer.data(pdg + 121);
    const auto *pdg_123 = buffer.data(pdg + 123);
    const auto *pdg_130 = buffer.data(pdg + 130);
    const auto *pdg_131 = buffer.data(pdg + 131);
    const auto *pdg_132 = buffer.data(pdg + 132);
    const auto *pdg_133 = buffer.data(pdg + 133);
    const auto *pdg_134 = buffer.data(pdg + 134);
    const auto *pdg_135 = buffer.data(pdg + 135);
    const auto *pdg_136 = buffer.data(pdg + 136);
    const auto *pdg_138 = buffer.data(pdg + 138);
    const auto *pdg_140 = buffer.data(pdg + 140);
    const auto *pdg_141 = buffer.data(pdg + 141);
    const auto *pdg_143 = buffer.data(pdg + 143);
    const auto *pdg_144 = buffer.data(pdg + 144);
    const auto *pdg_145 = buffer.data(pdg + 145);
    const auto *pdg_146 = buffer.data(pdg + 146);
    const auto *pdg_147 = buffer.data(pdg + 147);
    const auto *pdg_148 = buffer.data(pdg + 148);
    const auto *pdg_149 = buffer.data(pdg + 149);
    const auto *pdg_150 = buffer.data(pdg + 150);
    const auto *pdg_151 = buffer.data(pdg + 151);
    const auto *pdg_153 = buffer.data(pdg + 153);
    const auto *pdg_155 = buffer.data(pdg + 155);
    const auto *pdg_158 = buffer.data(pdg + 158);
    const auto *pdg_159 = buffer.data(pdg + 159);
    const auto *pdg_160 = buffer.data(pdg + 160);
    const auto *pdg_161 = buffer.data(pdg + 161);
    const auto *pdg_162 = buffer.data(pdg + 162);
    const auto *pdg_163 = buffer.data(pdg + 163);
    const auto *pdg_164 = buffer.data(pdg + 164);
    const auto *pdg_165 = buffer.data(pdg + 165);
    const auto *pdg_166 = buffer.data(pdg + 166);
    const auto *pdg_168 = buffer.data(pdg + 168);
    const auto *pdg_171 = buffer.data(pdg + 171);
    const auto *pdg_173 = buffer.data(pdg + 173);
    const auto *pdg_175 = buffer.data(pdg + 175);
    const auto *pdg_176 = buffer.data(pdg + 176);
    const auto *pdg_177 = buffer.data(pdg + 177);
    const auto *pdg_178 = buffer.data(pdg + 178);
    const auto *pdg_179 = buffer.data(pdg + 179);
    const auto *pdg_180 = buffer.data(pdg + 180);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_y, pc_y, pc_z, sdh0_3, sdh0_5, sdg_1, \
                         sdh1_3, sdh1_5, pdg_90, pdg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_4 * pc_z[k] * pdg_90[k];

        t_129[k] = pa_y[k] * sdh0_3[k]
                   + f_1 * sdg_1[k]
                   - f_9 * pc_y[k] * sdh1_3[k];

        t_130[k] = f_4 * pc_z[k] * pdg_91[k];

        t_131[k] = pa_y[k] * sdh0_5[k]
                   - f_9 * pc_y[k] * sdh1_5[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_y, pc_y, pc_z, sdh0_6, sdh0_8, sdh0_9, \
                         sdg_3, sdg_5, sdh1_6, sdh1_8, sdh1_9, pdg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * sdh0_6[k]
                   + f_11 * sdg_3[k]
                   - f_9 * pc_y[k] * sdh1_6[k];

        t_133[k] = f_4 * pc_z[k] * pdg_93[k];

        t_134[k] = pa_y[k] * sdh0_8[k]
                   + f_0 * sdg_5[k]
                   - f_9 * pc_y[k] * sdh1_8[k];

        t_135[k] = pa_y[k] * sdh0_9[k]
                   - f_9 * pc_y[k] * sdh1_9[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pc_x, ppg_55, ppg_56, ppg_57, \
                         ppg_58, ppg_59, pdg_100, pdg_101, pdg_102, pdg_103, \
                         pdg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_1 * ppg_55[k]
                   + f_4 * pc_x[k] * pdg_100[k];

        t_137[k] = f_1 * ppg_56[k]
                   + f_4 * pc_x[k] * pdg_101[k];

        t_138[k] = f_1 * ppg_57[k]
                   + f_4 * pc_x[k] * pdg_102[k];

        t_139[k] = f_1 * ppg_58[k]
                   + f_4 * pc_x[k] * pdg_103[k];

        t_140[k] = f_1 * ppg_59[k]
                   + f_4 * pc_x[k] * pdg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pa_y, pc_y, pc_z, sdh0_15, sdg_10, sdh1_15, \
                         pdf0_66, pdf1_66, pdg_100, pdg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pa_y[k] * sdh0_15[k]
                   + f_10 * sdg_10[k]
                   - f_9 * pc_y[k] * sdh1_15[k];

        t_142[k] = f_4 * pc_z[k] * pdg_100[k];

        t_143[k] = f_5 * pdf0_66[k]
                   - f_6 * pdf1_66[k]
                   + f_4 * pc_z[k] * pdg_101[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pa_y, pc_y, pc_z, sdh0_20, sdg_14, sdh1_20, \
                         pdf0_67, pdf1_67, pdg_102, pdg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_7 * pdf0_67[k]
                   - f_8 * pdf1_67[k]
                   + f_4 * pc_z[k] * pdg_102[k];

        t_145[k] = f_0 * sdg_14[k]
                   + f_4 * pc_y[k] * pdg_104[k];

        t_146[k] = pa_y[k] * sdh0_20[k]
                   - f_9 * pc_y[k] * sdh1_20[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_x, pc_x, pc_z, pph0_85, ppg_60, ppg_61, \
                         pph1_85, pdf0_70, pdf1_70, pdg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_0 * ppg_60[k]
                   + f_2 * pdf0_70[k]
                   - f_3 * pdf1_70[k]
                   + f_4 * pc_x[k] * pdg_105[k];

        t_148[k] = pb_x[k] * pph0_85[k]
                   + f_12 * ppg_61[k]
                   - f_9 * pc_x[k] * pph1_85[k];

        t_149[k] = f_4 * pc_z[k] * pdg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_x, pc_x, pc_z, pph0_87, ppg_63, ppg_65, \
                         pph1_87, pdf0_75, pdf1_75, pdg_106, pdg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_x[k] * pph0_87[k]
                   + f_11 * ppg_63[k]
                   - f_9 * pc_x[k] * pph1_87[k];

        t_151[k] = f_4 * pc_z[k] * pdg_106[k];

        t_152[k] = f_0 * ppg_65[k]
                   + f_7 * pdf0_75[k]
                   - f_8 * pdf1_75[k]
                   + f_4 * pc_x[k] * pdg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_x, pc_x, pc_z, pph0_90, pph0_92, ppg_66, \
                         ppg_68, pph1_90, pph1_92, pdg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pb_x[k] * pph0_90[k]
                   + f_1 * ppg_66[k]
                   - f_9 * pc_x[k] * pph1_90[k];

        t_154[k] = f_4 * pc_z[k] * pdg_108[k];

        t_155[k] = pb_x[k] * pph0_92[k]
                   + f_1 * ppg_68[k]
                   - f_9 * pc_x[k] * pph1_92[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, ppg_69, ppg_70, ppg_71, ppg_72, \
                         pdf0_79, pdf1_79, pdg_114, pdg_115, pdg_116, \
                         pdg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * ppg_69[k]
                   + f_5 * pdf0_79[k]
                   - f_6 * pdf1_79[k]
                   + f_4 * pc_x[k] * pdg_114[k];

        t_157[k] = f_0 * ppg_70[k]
                   + f_4 * pc_x[k] * pdg_115[k];

        t_158[k] = f_0 * ppg_71[k]
                   + f_4 * pc_x[k] * pdg_116[k];

        t_159[k] = f_0 * ppg_72[k]
                   + f_4 * pc_x[k] * pdg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pb_x, pc_x, pc_z, pph0_99, ppg_73, \
                         ppg_74, pph1_99, pdg_115, pdg_118, pdg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_0 * ppg_73[k]
                   + f_4 * pc_x[k] * pdg_118[k];

        t_161[k] = f_0 * ppg_74[k]
                   + f_4 * pc_x[k] * pdg_119[k];

        t_162[k] = pb_x[k] * pph0_99[k]
                   - f_9 * pc_x[k] * pph1_99[k];

        t_163[k] = f_4 * pc_z[k] * pdg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pb_x, pc_x, pph0_101, pph0_102, pph0_103, \
                         pph0_104, pph1_101, pph1_102, pph1_103, \
                         pph1_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pb_x[k] * pph0_101[k]
                   - f_9 * pc_x[k] * pph1_101[k];

        t_165[k] = pb_x[k] * pph0_102[k]
                   - f_9 * pc_x[k] * pph1_102[k];

        t_166[k] = pb_x[k] * pph0_103[k]
                   - f_9 * pc_x[k] * pph1_103[k];

        t_167[k] = pb_x[k] * pph0_104[k]
                   - f_9 * pc_x[k] * pph1_104[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pb_z, pc_y, pc_z, sdh0_42, sdh1_42, \
                         pph0_64, pph0_66, ppg_45, pph1_64, pph1_66, \
                         pdg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * sdh0_42[k]
                   - f_9 * pc_y[k] * sdh1_42[k];

        t_169[k] = pb_z[k] * pph0_64[k]
                   - f_9 * pc_z[k] * pph1_64[k];

        t_170[k] = f_0 * ppg_45[k]
                   + f_4 * pc_z[k] * pdg_120[k];

        t_171[k] = pb_z[k] * pph0_66[k]
                   - f_9 * pc_z[k] * pph1_66[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_z, pc_y, pc_z, sdh0_47, sdh1_47, \
                         pph0_69, ppg_46, ppg_48, pph1_69, pdg_121, \
                         pdg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_0 * ppg_46[k]
                   + f_4 * pc_z[k] * pdg_121[k];

        t_173[k] = pa_y[k] * sdh0_47[k]
                   - f_9 * pc_y[k] * sdh1_47[k];

        t_174[k] = pb_z[k] * pph0_69[k]
                   - f_9 * pc_z[k] * pph1_69[k];

        t_175[k] = f_0 * ppg_48[k]
                   + f_4 * pc_z[k] * pdg_123[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pc_x, pc_y, sdh0_50, sdh0_51, \
                         sdg_35, sdh1_50, sdh1_51, ppg_85, ppg_86, pdg_130, \
                         pdg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_y[k] * sdh0_50[k]
                   + f_0 * sdg_35[k]
                   - f_9 * pc_y[k] * sdh1_50[k];

        t_177[k] = pa_y[k] * sdh0_51[k]
                   - f_9 * pc_y[k] * sdh1_51[k];

        t_178[k] = f_0 * ppg_85[k]
                   + f_4 * pc_x[k] * pdg_130[k];

        t_179[k] = f_0 * ppg_86[k]
                   + f_4 * pc_x[k] * pdg_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pb_x, pc_x, pph0_120, ppg_87, ppg_88, \
                         ppg_89, pph1_120, pdg_132, pdg_133, pdg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_0 * ppg_87[k]
                   + f_4 * pc_x[k] * pdg_132[k];

        t_181[k] = f_0 * ppg_88[k]
                   + f_4 * pc_x[k] * pdg_133[k];

        t_182[k] = f_0 * ppg_89[k]
                   + f_4 * pc_x[k] * pdg_134[k];

        t_183[k] = pb_x[k] * pph0_120[k]
                   - f_9 * pc_x[k] * pph1_120[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pc_x, pc_y, pc_z, sdg_44, pph0_122, \
                         pph0_123, ppg_55, pph1_122, pph1_123, pdg_130, \
                         pdg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_0 * ppg_55[k]
                   + f_4 * pc_z[k] * pdg_130[k];

        t_185[k] = pb_x[k] * pph0_122[k]
                   - f_9 * pc_x[k] * pph1_122[k];

        t_186[k] = pb_x[k] * pph0_123[k]
                   - f_9 * pc_x[k] * pph1_123[k];

        t_187[k] = f_0 * sdg_44[k]
                   + f_4 * pc_y[k] * pdg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pc_x, pc_y, pc_z, sdh0_62, sdh1_62, \
                         pdf0_90, pdf0_91, pdf1_90, pdf1_91, pdg_135, \
                         pdg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * sdh0_62[k]
                   - f_9 * pc_y[k] * sdh1_62[k];

        t_189[k] = f_2 * pdf0_90[k]
                   - f_3 * pdf1_90[k]
                   + f_4 * pc_x[k] * pdg_135[k];

        t_190[k] = f_13 * pdf0_91[k]
                   - f_14 * pdf1_91[k]
                   + f_4 * pc_x[k] * pdg_136[k];

        t_191[k] = f_4 * pc_z[k] * pdg_135[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_x, pc_z, pdf0_93, pdf0_95, pdf0_96, \
                         pdf1_93, pdf1_95, pdf1_96, pdg_136, pdg_138, pdg_140, \
                         pdg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_7 * pdf0_93[k]
                   - f_8 * pdf1_93[k]
                   + f_4 * pc_x[k] * pdg_138[k];

        t_193[k] = f_4 * pc_z[k] * pdg_136[k];

        t_194[k] = f_7 * pdf0_95[k]
                   - f_8 * pdf1_95[k]
                   + f_4 * pc_x[k] * pdg_140[k];

        t_195[k] = f_5 * pdf0_96[k]
                   - f_6 * pdf1_96[k]
                   + f_4 * pc_x[k] * pdg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pc_x, pc_z, pdf0_98, pdf0_99, \
                         pdf1_98, pdf1_99, pdg_138, pdg_143, pdg_144, pdg_145, \
                         pdg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_4 * pc_z[k] * pdg_138[k];

        t_197[k] = f_5 * pdf0_98[k]
                   - f_6 * pdf1_98[k]
                   + f_4 * pc_x[k] * pdg_143[k];

        t_198[k] = f_5 * pdf0_99[k]
                   - f_6 * pdf1_99[k]
                   + f_4 * pc_x[k] * pdg_144[k];

        t_199[k] = f_4 * pc_x[k] * pdg_145[k];

        t_200[k] = f_4 * pc_x[k] * pdg_146[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pc_x, pc_y, pc_z, sdg_55, ppg_70, \
                         pdf0_96, pdf1_96, pdg_145, pdg_147, pdg_148, \
                         pdg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_4 * pc_x[k] * pdg_147[k];

        t_202[k] = f_4 * pc_x[k] * pdg_148[k];

        t_203[k] = f_4 * pc_x[k] * pdg_149[k];

        t_204[k] = f_0 * sdg_55[k]
                   + f_1 * ppg_70[k]
                   + f_2 * pdf0_96[k]
                   - f_3 * pdf1_96[k]
                   + f_4 * pc_y[k] * pdg_145[k];

        t_205[k] = f_4 * pc_z[k] * pdg_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pc_y, pc_z, sdg_59, ppg_74, pdf0_96, pdf0_97, \
                         pdf1_96, pdf1_97, pdg_146, pdg_147, pdg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_5 * pdf0_96[k]
                   - f_6 * pdf1_96[k]
                   + f_4 * pc_z[k] * pdg_146[k];

        t_207[k] = f_7 * pdf0_97[k]
                   - f_8 * pdf1_97[k]
                   + f_4 * pc_z[k] * pdg_147[k];

        t_208[k] = f_0 * sdg_59[k]
                   + f_1 * ppg_74[k]
                   + f_4 * pc_y[k] * pdg_149[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_z, pc_x, pc_z, pph0_85, ppg_60, \
                         pph1_85, pdf0_99, pdf0_100, pdf1_99, pdf1_100, pdg_149, \
                         pdg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_2 * pdf0_99[k]
                   - f_3 * pdf1_99[k]
                   + f_4 * pc_z[k] * pdg_149[k];

        t_210[k] = f_2 * pdf0_100[k]
                   - f_3 * pdf1_100[k]
                   + f_4 * pc_x[k] * pdg_150[k];

        t_211[k] = pb_z[k] * pph0_85[k]
                   - f_9 * pc_z[k] * pph1_85[k];

        t_212[k] = f_0 * ppg_60[k]
                   + f_4 * pc_z[k] * pdg_150[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, pph0_87, pph0_90, \
                         ppg_61, pph1_87, pph1_90, pdf0_105, pdf1_105, pdg_151, \
                         pdg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = pb_z[k] * pph0_87[k]
                   - f_9 * pc_z[k] * pph1_87[k];

        t_214[k] = f_0 * ppg_61[k]
                   + f_4 * pc_z[k] * pdg_151[k];

        t_215[k] = f_7 * pdf0_105[k]
                   - f_8 * pdf1_105[k]
                   + f_4 * pc_x[k] * pdg_155[k];

        t_216[k] = pb_z[k] * pph0_90[k]
                   - f_9 * pc_z[k] * pph1_90[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pc_x, pc_z, ppg_63, pdf0_108, pdf0_109, \
                         pdf1_108, pdf1_109, pdg_153, pdg_158, pdg_159, \
                         pdg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_0 * ppg_63[k]
                   + f_4 * pc_z[k] * pdg_153[k];

        t_218[k] = f_5 * pdf0_108[k]
                   - f_6 * pdf1_108[k]
                   + f_4 * pc_x[k] * pdg_158[k];

        t_219[k] = f_5 * pdf0_109[k]
                   - f_6 * pdf1_109[k]
                   + f_4 * pc_x[k] * pdg_159[k];

        t_220[k] = f_4 * pc_x[k] * pdg_160[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, pb_z, pc_x, pc_z, pph0_99, \
                         pph1_99, pdg_161, pdg_162, pdg_163, pdg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_4 * pc_x[k] * pdg_161[k];

        t_222[k] = f_4 * pc_x[k] * pdg_162[k];

        t_223[k] = f_4 * pc_x[k] * pdg_163[k];

        t_224[k] = f_4 * pc_x[k] * pdg_164[k];

        t_225[k] = pb_z[k] * pph0_99[k]
                   - f_9 * pc_z[k] * pph1_99[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_z, ppg_70, ppg_71, ppg_72, pdf0_106, \
                         pdf0_107, pdf1_106, pdf1_107, pdg_160, pdg_161, \
                         pdg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_0 * ppg_70[k]
                   + f_4 * pc_z[k] * pdg_160[k];

        t_227[k] = f_0 * ppg_71[k]
                   + f_5 * pdf0_106[k]
                   - f_6 * pdf1_106[k]
                   + f_4 * pc_z[k] * pdg_161[k];

        t_228[k] = f_0 * ppg_72[k]
                   + f_7 * pdf0_107[k]
                   - f_8 * pdf1_107[k]
                   + f_4 * pc_z[k] * pdg_162[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pa_y, pc_y, pc_z, sdh0_105, sdg_74, sdh1_105, \
                         ppg_74, ppg_89, pdf0_109, pdf1_109, pdg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * sdg_74[k]
                   + f_0 * ppg_89[k]
                   + f_4 * pc_y[k] * pdg_164[k];

        t_230[k] = f_0 * ppg_74[k]
                   + f_2 * pdf0_109[k]
                   - f_3 * pdf1_109[k]
                   + f_4 * pc_z[k] * pdg_164[k];

        t_231[k] = pa_y[k] * sdh0_105[k]
                   - f_9 * pc_y[k] * sdh1_105[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_z, ppg_75, ppg_76, pdf0_111, \
                         pdf0_113, pdf1_111, pdf1_113, pdg_165, pdg_166, \
                         pdg_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_13 * pdf0_111[k]
                   - f_14 * pdf1_111[k]
                   + f_4 * pc_x[k] * pdg_166[k];

        t_233[k] = f_1 * ppg_75[k]
                   + f_4 * pc_z[k] * pdg_165[k];

        t_234[k] = f_7 * pdf0_113[k]
                   - f_8 * pdf1_113[k]
                   + f_4 * pc_x[k] * pdg_168[k];

        t_235[k] = f_1 * ppg_76[k]
                   + f_4 * pc_z[k] * pdg_166[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_y, pc_x, pc_y, pc_z, sdh0_110, sdh1_110, \
                         ppg_78, pdf0_116, pdf1_116, pdg_168, pdg_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * sdh0_110[k]
                   - f_9 * pc_y[k] * sdh1_110[k];

        t_237[k] = f_5 * pdf0_116[k]
                   - f_6 * pdf1_116[k]
                   + f_4 * pc_x[k] * pdg_171[k];

        t_238[k] = f_1 * ppg_78[k]
                   + f_4 * pc_z[k] * pdg_168[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_y, pc_x, pc_y, sdh0_114, \
                         sdh1_114, pdf0_118, pdf1_118, pdg_173, pdg_175, pdg_176, \
                         pdg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_5 * pdf0_118[k]
                   - f_6 * pdf1_118[k]
                   + f_4 * pc_x[k] * pdg_173[k];

        t_240[k] = pa_y[k] * sdh0_114[k]
                   - f_9 * pc_y[k] * sdh1_114[k];

        t_241[k] = f_4 * pc_x[k] * pdg_175[k];

        t_242[k] = f_4 * pc_x[k] * pdg_176[k];

        t_243[k] = f_4 * pc_x[k] * pdg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_y, pc_x, pc_y, pc_z, sdh0_120, sdg_85, \
                         sdh1_120, ppg_85, pdg_175, pdg_178, pdg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_4 * pc_x[k] * pdg_178[k];

        t_245[k] = f_4 * pc_x[k] * pdg_179[k];

        t_246[k] = pa_y[k] * sdh0_120[k]
                   + f_10 * sdg_85[k]
                   - f_9 * pc_y[k] * sdh1_120[k];

        t_247[k] = f_1 * ppg_85[k]
                   + f_4 * pc_z[k] * pdg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_y, pc_y, sdh0_122, sdh0_123, sdh0_125, \
                         sdg_87, sdg_88, sdg_89, sdh1_122, sdh1_123, sdh1_125, \
                         pdg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pa_y[k] * sdh0_122[k]
                   + f_11 * sdg_87[k]
                   - f_9 * pc_y[k] * sdh1_122[k];

        t_249[k] = pa_y[k] * sdh0_123[k]
                   + f_1 * sdg_88[k]
                   - f_9 * pc_y[k] * sdh1_123[k];

        t_250[k] = f_0 * sdg_89[k]
                   + f_4 * pc_y[k] * pdg_179[k];

        t_251[k] = pa_y[k] * sdh0_125[k]
                   - f_9 * pc_y[k] * sdh1_125[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_z, pc_y, pc_z, sdh0_0, sdh0_2, sdh0_3, \
                         sdg_0, sdh1_0, sdh1_2, sdh1_3, pdg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pa_z[k] * sdh0_0[k]
                   - f_9 * pc_z[k] * sdh1_0[k];

        t_253[k] = f_4 * pc_y[k] * pdg_180[k];

        t_254[k] = pa_z[k] * sdh0_2[k]
                   + f_0 * sdg_0[k]
                   - f_9 * pc_z[k] * sdh1_2[k];

        t_255[k] = pa_z[k] * sdh0_3[k]
                   - f_9 * pc_z[k] * sdh1_3[k];
    }
}

static auto
compute_prim_pdh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdh0, const size_t sdg,
                                                          const size_t sdh1, const size_t pph0,
                                                          const size_t ppg, const size_t pph1,
                                                          const size_t pdf0, const size_t pdf1,
                                                          const size_t pdg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);

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
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdh0_5 = buffer.data(sdh0 + 5);
    const auto *sdh0_6 = buffer.data(sdh0 + 6);
    const auto *sdh0_7 = buffer.data(sdh0 + 7);
    const auto *sdh0_9 = buffer.data(sdh0 + 9);
    const auto *sdh0_15 = buffer.data(sdh0 + 15);
    const auto *sdh0_20 = buffer.data(sdh0 + 20);
    const auto *sdh0_21 = buffer.data(sdh0 + 21);
    const auto *sdh0_24 = buffer.data(sdh0 + 24);
    const auto *sdh0_27 = buffer.data(sdh0 + 27);
    const auto *sdh0_28 = buffer.data(sdh0 + 28);
    const auto *sdh0_36 = buffer.data(sdh0 + 36);
    const auto *sdh0_63 = buffer.data(sdh0 + 63);
    const auto *sdh0_66 = buffer.data(sdh0 + 66);
    const auto *sdh0_69 = buffer.data(sdh0 + 69);
    const auto *sdh0_78 = buffer.data(sdh0 + 78);
    const auto *sdh0_79 = buffer.data(sdh0 + 79);
    const auto *sdh0_80 = buffer.data(sdh0 + 80);
    const auto *sdh0_81 = buffer.data(sdh0 + 81);
    const auto *sdh0_83 = buffer.data(sdh0 + 83);

    const auto *sdg_2 = buffer.data(sdg + 2);
    const auto *sdg_3 = buffer.data(sdg + 3);
    const auto *sdg_5 = buffer.data(sdg + 5);
    const auto *sdg_14 = buffer.data(sdg + 14);
    const auto *sdg_18 = buffer.data(sdg + 18);
    const auto *sdg_55 = buffer.data(sdg + 55);
    const auto *sdg_56 = buffer.data(sdg + 56);
    const auto *sdg_57 = buffer.data(sdg + 57);
    const auto *sdg_59 = buffer.data(sdg + 59);
    const auto *sdg_89 = buffer.data(sdg + 89);

    const auto *sdh1_5 = buffer.data(sdh1 + 5);
    const auto *sdh1_6 = buffer.data(sdh1 + 6);
    const auto *sdh1_7 = buffer.data(sdh1 + 7);
    const auto *sdh1_9 = buffer.data(sdh1 + 9);
    const auto *sdh1_15 = buffer.data(sdh1 + 15);
    const auto *sdh1_20 = buffer.data(sdh1 + 20);
    const auto *sdh1_21 = buffer.data(sdh1 + 21);
    const auto *sdh1_24 = buffer.data(sdh1 + 24);
    const auto *sdh1_27 = buffer.data(sdh1 + 27);
    const auto *sdh1_28 = buffer.data(sdh1 + 28);
    const auto *sdh1_36 = buffer.data(sdh1 + 36);
    const auto *sdh1_63 = buffer.data(sdh1 + 63);
    const auto *sdh1_66 = buffer.data(sdh1 + 66);
    const auto *sdh1_69 = buffer.data(sdh1 + 69);
    const auto *sdh1_78 = buffer.data(sdh1 + 78);
    const auto *sdh1_79 = buffer.data(sdh1 + 79);
    const auto *sdh1_80 = buffer.data(sdh1 + 80);
    const auto *sdh1_81 = buffer.data(sdh1 + 81);
    const auto *sdh1_83 = buffer.data(sdh1 + 83);

    const auto *pph0_128 = buffer.data(pph0 + 128);
    const auto *pph0_131 = buffer.data(pph0 + 131);
    const auto *pph0_135 = buffer.data(pph0 + 135);
    const auto *pph0_163 = buffer.data(pph0 + 163);
    const auto *pph0_164 = buffer.data(pph0 + 164);
    const auto *pph0_165 = buffer.data(pph0 + 165);
    const auto *pph0_167 = buffer.data(pph0 + 167);
    const auto *pph0_168 = buffer.data(pph0 + 168);
    const auto *pph0_170 = buffer.data(pph0 + 170);
    const auto *pph0_173 = buffer.data(pph0 + 173);
    const auto *pph0_175 = buffer.data(pph0 + 175);
    const auto *pph0_177 = buffer.data(pph0 + 177);
    const auto *pph0_183 = buffer.data(pph0 + 183);
    const auto *pph0_184 = buffer.data(pph0 + 184);
    const auto *pph0_185 = buffer.data(pph0 + 185);
    const auto *pph0_186 = buffer.data(pph0 + 186);
    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *ppg_90 = buffer.data(ppg + 90);
    const auto *ppg_92 = buffer.data(ppg + 92);
    const auto *ppg_95 = buffer.data(ppg + 95);
    const auto *ppg_100 = buffer.data(ppg + 100);
    const auto *ppg_101 = buffer.data(ppg + 101);
    const auto *ppg_102 = buffer.data(ppg + 102);
    const auto *ppg_103 = buffer.data(ppg + 103);
    const auto *ppg_104 = buffer.data(ppg + 104);
    const auto *ppg_105 = buffer.data(ppg + 105);
    const auto *ppg_107 = buffer.data(ppg + 107);
    const auto *ppg_110 = buffer.data(ppg + 110);
    const auto *ppg_115 = buffer.data(ppg + 115);
    const auto *ppg_116 = buffer.data(ppg + 116);
    const auto *ppg_117 = buffer.data(ppg + 117);
    const auto *ppg_118 = buffer.data(ppg + 118);
    const auto *ppg_119 = buffer.data(ppg + 119);
    const auto *ppg_120 = buffer.data(ppg + 120);
    const auto *ppg_122 = buffer.data(ppg + 122);
    const auto *ppg_123 = buffer.data(ppg + 123);
    const auto *ppg_125 = buffer.data(ppg + 125);
    const auto *ppg_126 = buffer.data(ppg + 126);
    const auto *ppg_127 = buffer.data(ppg + 127);
    const auto *ppg_129 = buffer.data(ppg + 129);
    const auto *ppg_130 = buffer.data(ppg + 130);
    const auto *ppg_131 = buffer.data(ppg + 131);
    const auto *ppg_132 = buffer.data(ppg + 132);
    const auto *ppg_133 = buffer.data(ppg + 133);
    const auto *ppg_134 = buffer.data(ppg + 134);

    const auto *pph1_128 = buffer.data(pph1 + 128);
    const auto *pph1_131 = buffer.data(pph1 + 131);
    const auto *pph1_135 = buffer.data(pph1 + 135);
    const auto *pph1_163 = buffer.data(pph1 + 163);
    const auto *pph1_164 = buffer.data(pph1 + 164);
    const auto *pph1_165 = buffer.data(pph1 + 165);
    const auto *pph1_167 = buffer.data(pph1 + 167);
    const auto *pph1_168 = buffer.data(pph1 + 168);
    const auto *pph1_170 = buffer.data(pph1 + 170);
    const auto *pph1_173 = buffer.data(pph1 + 173);
    const auto *pph1_175 = buffer.data(pph1 + 175);
    const auto *pph1_177 = buffer.data(pph1 + 177);
    const auto *pph1_183 = buffer.data(pph1 + 183);
    const auto *pph1_184 = buffer.data(pph1 + 184);
    const auto *pph1_185 = buffer.data(pph1 + 185);
    const auto *pph1_186 = buffer.data(pph1 + 186);
    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *pdf0_127 = buffer.data(pdf0 + 127);
    const auto *pdf0_128 = buffer.data(pdf0 + 128);
    const auto *pdf0_129 = buffer.data(pdf0 + 129);
    const auto *pdf0_140 = buffer.data(pdf0 + 140);
    const auto *pdf0_143 = buffer.data(pdf0 + 143);
    const auto *pdf0_146 = buffer.data(pdf0 + 146);
    const auto *pdf0_152 = buffer.data(pdf0 + 152);
    const auto *pdf0_155 = buffer.data(pdf0 + 155);
    const auto *pdf0_157 = buffer.data(pdf0 + 157);
    const auto *pdf0_159 = buffer.data(pdf0 + 159);
    const auto *pdf0_163 = buffer.data(pdf0 + 163);
    const auto *pdf0_166 = buffer.data(pdf0 + 166);
    const auto *pdf0_167 = buffer.data(pdf0 + 167);
    const auto *pdf0_170 = buffer.data(pdf0 + 170);
    const auto *pdf0_172 = buffer.data(pdf0 + 172);
    const auto *pdf0_173 = buffer.data(pdf0 + 173);
    const auto *pdf0_175 = buffer.data(pdf0 + 175);
    const auto *pdf0_176 = buffer.data(pdf0 + 176);
    const auto *pdf0_177 = buffer.data(pdf0 + 177);
    const auto *pdf0_178 = buffer.data(pdf0 + 178);
    const auto *pdf0_179 = buffer.data(pdf0 + 179);

    const auto *pdf1_127 = buffer.data(pdf1 + 127);
    const auto *pdf1_128 = buffer.data(pdf1 + 128);
    const auto *pdf1_129 = buffer.data(pdf1 + 129);
    const auto *pdf1_140 = buffer.data(pdf1 + 140);
    const auto *pdf1_143 = buffer.data(pdf1 + 143);
    const auto *pdf1_146 = buffer.data(pdf1 + 146);
    const auto *pdf1_152 = buffer.data(pdf1 + 152);
    const auto *pdf1_155 = buffer.data(pdf1 + 155);
    const auto *pdf1_157 = buffer.data(pdf1 + 157);
    const auto *pdf1_159 = buffer.data(pdf1 + 159);
    const auto *pdf1_163 = buffer.data(pdf1 + 163);
    const auto *pdf1_166 = buffer.data(pdf1 + 166);
    const auto *pdf1_167 = buffer.data(pdf1 + 167);
    const auto *pdf1_170 = buffer.data(pdf1 + 170);
    const auto *pdf1_172 = buffer.data(pdf1 + 172);
    const auto *pdf1_173 = buffer.data(pdf1 + 173);
    const auto *pdf1_175 = buffer.data(pdf1 + 175);
    const auto *pdf1_176 = buffer.data(pdf1 + 176);
    const auto *pdf1_177 = buffer.data(pdf1 + 177);
    const auto *pdf1_178 = buffer.data(pdf1 + 178);
    const auto *pdf1_179 = buffer.data(pdf1 + 179);

    const auto *pdg_182 = buffer.data(pdg + 182);
    const auto *pdg_185 = buffer.data(pdg + 185);
    const auto *pdg_190 = buffer.data(pdg + 190);
    const auto *pdg_191 = buffer.data(pdg + 191);
    const auto *pdg_192 = buffer.data(pdg + 192);
    const auto *pdg_193 = buffer.data(pdg + 193);
    const auto *pdg_194 = buffer.data(pdg + 194);
    const auto *pdg_195 = buffer.data(pdg + 195);
    const auto *pdg_197 = buffer.data(pdg + 197);
    const auto *pdg_200 = buffer.data(pdg + 200);
    const auto *pdg_205 = buffer.data(pdg + 205);
    const auto *pdg_206 = buffer.data(pdg + 206);
    const auto *pdg_207 = buffer.data(pdg + 207);
    const auto *pdg_208 = buffer.data(pdg + 208);
    const auto *pdg_209 = buffer.data(pdg + 209);
    const auto *pdg_210 = buffer.data(pdg + 210);
    const auto *pdg_212 = buffer.data(pdg + 212);
    const auto *pdg_213 = buffer.data(pdg + 213);
    const auto *pdg_215 = buffer.data(pdg + 215);
    const auto *pdg_216 = buffer.data(pdg + 216);
    const auto *pdg_220 = buffer.data(pdg + 220);
    const auto *pdg_221 = buffer.data(pdg + 221);
    const auto *pdg_222 = buffer.data(pdg + 222);
    const auto *pdg_223 = buffer.data(pdg + 223);
    const auto *pdg_224 = buffer.data(pdg + 224);
    const auto *pdg_225 = buffer.data(pdg + 225);
    const auto *pdg_227 = buffer.data(pdg + 227);
    const auto *pdg_230 = buffer.data(pdg + 230);
    const auto *pdg_232 = buffer.data(pdg + 232);
    const auto *pdg_234 = buffer.data(pdg + 234);
    const auto *pdg_235 = buffer.data(pdg + 235);
    const auto *pdg_236 = buffer.data(pdg + 236);
    const auto *pdg_237 = buffer.data(pdg + 237);
    const auto *pdg_238 = buffer.data(pdg + 238);
    const auto *pdg_239 = buffer.data(pdg + 239);
    const auto *pdg_240 = buffer.data(pdg + 240);
    const auto *pdg_242 = buffer.data(pdg + 242);
    const auto *pdg_243 = buffer.data(pdg + 243);
    const auto *pdg_245 = buffer.data(pdg + 245);
    const auto *pdg_246 = buffer.data(pdg + 246);
    const auto *pdg_247 = buffer.data(pdg + 247);
    const auto *pdg_250 = buffer.data(pdg + 250);
    const auto *pdg_251 = buffer.data(pdg + 251);
    const auto *pdg_252 = buffer.data(pdg + 252);
    const auto *pdg_253 = buffer.data(pdg + 253);
    const auto *pdg_254 = buffer.data(pdg + 254);
    const auto *pdg_255 = buffer.data(pdg + 255);
    const auto *pdg_257 = buffer.data(pdg + 257);
    const auto *pdg_258 = buffer.data(pdg + 258);
    const auto *pdg_260 = buffer.data(pdg + 260);
    const auto *pdg_261 = buffer.data(pdg + 261);
    const auto *pdg_262 = buffer.data(pdg + 262);
    const auto *pdg_264 = buffer.data(pdg + 264);
    const auto *pdg_265 = buffer.data(pdg + 265);
    const auto *pdg_266 = buffer.data(pdg + 266);
    const auto *pdg_267 = buffer.data(pdg + 267);
    const auto *pdg_268 = buffer.data(pdg + 268);
    const auto *pdg_269 = buffer.data(pdg + 269);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_z, pc_y, pc_z, sdh0_5, sdh0_6, sdh0_7, \
                         sdg_2, sdg_3, sdh1_5, sdh1_6, sdh1_7, \
                         pdg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_4 * pc_y[k] * pdg_182[k];

        t_257[k] = pa_z[k] * sdh0_5[k]
                   + f_1 * sdg_2[k]
                   - f_9 * pc_z[k] * sdh1_5[k];

        t_258[k] = pa_z[k] * sdh0_6[k]
                   - f_9 * pc_z[k] * sdh1_6[k];

        t_259[k] = pa_z[k] * sdh0_7[k]
                   + f_0 * sdg_3[k]
                   - f_9 * pc_z[k] * sdh1_7[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_z, pc_x, pc_y, pc_z, sdh0_9, sdg_5, \
                         sdh1_9, ppg_100, ppg_101, pdg_185, pdg_190, \
                         pdg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_4 * pc_y[k] * pdg_185[k];

        t_261[k] = pa_z[k] * sdh0_9[k]
                   + f_11 * sdg_5[k]
                   - f_9 * pc_z[k] * sdh1_9[k];

        t_262[k] = f_1 * ppg_100[k]
                   + f_4 * pc_x[k] * pdg_190[k];

        t_263[k] = f_1 * ppg_101[k]
                   + f_4 * pc_x[k] * pdg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_z, pc_x, pc_z, sdh0_15, sdh1_15, \
                         ppg_102, ppg_103, ppg_104, pdg_192, pdg_193, \
                         pdg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_1 * ppg_102[k]
                   + f_4 * pc_x[k] * pdg_192[k];

        t_265[k] = f_1 * ppg_103[k]
                   + f_4 * pc_x[k] * pdg_193[k];

        t_266[k] = f_1 * ppg_104[k]
                   + f_4 * pc_x[k] * pdg_194[k];

        t_267[k] = pa_z[k] * sdh0_15[k]
                   - f_9 * pc_z[k] * sdh1_15[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pc_y, pdf0_127, pdf0_128, pdf0_129, \
                         pdf1_127, pdf1_128, pdf1_129, pdg_191, pdg_192, pdg_193, \
                         pdg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_13 * pdf0_127[k]
                   - f_14 * pdf1_127[k]
                   + f_4 * pc_y[k] * pdg_191[k];

        t_269[k] = f_7 * pdf0_128[k]
                   - f_8 * pdf1_128[k]
                   + f_4 * pc_y[k] * pdg_192[k];

        t_270[k] = f_5 * pdf0_129[k]
                   - f_6 * pdf1_129[k]
                   + f_4 * pc_y[k] * pdg_193[k];

        t_271[k] = f_4 * pc_y[k] * pdg_194[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pa_z, pc_y, pc_z, sdh0_20, sdh0_21, sdg_14, \
                         sdh1_20, sdh1_21, ppg_90, pdg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = pa_z[k] * sdh0_20[k]
                   + f_10 * sdg_14[k]
                   - f_9 * pc_z[k] * sdh1_20[k];

        t_273[k] = pa_z[k] * sdh0_21[k]
                   - f_9 * pc_z[k] * sdh1_21[k];

        t_274[k] = f_0 * ppg_90[k]
                   + f_4 * pc_y[k] * pdg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_z, pb_y, pc_y, pc_z, sdh0_24, sdh1_24, \
                         pph0_128, pph0_131, ppg_92, pph1_128, pph1_131, \
                         pdg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = pb_y[k] * pph0_128[k]
                   - f_9 * pc_y[k] * pph1_128[k];

        t_276[k] = pa_z[k] * sdh0_24[k]
                   - f_9 * pc_z[k] * sdh1_24[k];

        t_277[k] = f_0 * ppg_92[k]
                   + f_4 * pc_y[k] * pdg_197[k];

        t_278[k] = pb_y[k] * pph0_131[k]
                   - f_9 * pc_y[k] * pph1_131[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pa_z, pc_y, pc_z, sdh0_27, sdh0_28, sdg_18, \
                         sdh1_27, sdh1_28, ppg_95, pdg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_z[k] * sdh0_27[k]
                   - f_9 * pc_z[k] * sdh1_27[k];

        t_280[k] = pa_z[k] * sdh0_28[k]
                   + f_0 * sdg_18[k]
                   - f_9 * pc_z[k] * sdh1_28[k];

        t_281[k] = f_0 * ppg_95[k]
                   + f_4 * pc_y[k] * pdg_200[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_x, pc_y, pph0_135, ppg_115, \
                         ppg_116, ppg_117, pph1_135, pdg_205, pdg_206, \
                         pdg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pb_y[k] * pph0_135[k]
                   - f_9 * pc_y[k] * pph1_135[k];

        t_283[k] = f_0 * ppg_115[k]
                   + f_4 * pc_x[k] * pdg_205[k];

        t_284[k] = f_0 * ppg_116[k]
                   + f_4 * pc_x[k] * pdg_206[k];

        t_285[k] = f_0 * ppg_117[k]
                   + f_4 * pc_x[k] * pdg_207[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_z, pb_x, pc_x, pc_z, sdh0_36, sdh1_36, \
                         pph0_163, ppg_118, ppg_119, pph1_163, pdg_208, \
                         pdg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * ppg_118[k]
                   + f_4 * pc_x[k] * pdg_208[k];

        t_287[k] = f_0 * ppg_119[k]
                   + f_4 * pc_x[k] * pdg_209[k];

        t_288[k] = pa_z[k] * sdh0_36[k]
                   - f_9 * pc_z[k] * sdh1_36[k];

        t_289[k] = pb_x[k] * pph0_163[k]
                   - f_9 * pc_x[k] * pph1_163[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_x, pc_x, pc_y, pph0_164, pph0_165, \
                         pph0_167, ppg_104, pph1_164, pph1_165, pph1_167, \
                         pdg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * pph0_164[k]
                   - f_9 * pc_x[k] * pph1_164[k];

        t_291[k] = pb_x[k] * pph0_165[k]
                   - f_9 * pc_x[k] * pph1_165[k];

        t_292[k] = f_0 * ppg_104[k]
                   + f_4 * pc_y[k] * pdg_209[k];

        t_293[k] = pb_x[k] * pph0_167[k]
                   - f_9 * pc_x[k] * pph1_167[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pb_x, pc_x, pc_y, pph0_170, ppg_120, ppg_122, \
                         pph1_170, pdf0_140, pdf1_140, pdg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * ppg_120[k]
                   + f_2 * pdf0_140[k]
                   - f_3 * pdf1_140[k]
                   + f_4 * pc_x[k] * pdg_210[k];

        t_295[k] = f_4 * pc_y[k] * pdg_210[k];

        t_296[k] = pb_x[k] * pph0_170[k]
                   + f_12 * ppg_122[k]
                   - f_9 * pc_x[k] * pph1_170[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pc_x, pc_y, pph0_173, ppg_123, ppg_125, \
                         pph1_173, pdf0_143, pdf1_143, pdg_212, \
                         pdg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_0 * ppg_123[k]
                   + f_7 * pdf0_143[k]
                   - f_8 * pdf1_143[k]
                   + f_4 * pc_x[k] * pdg_213[k];

        t_298[k] = f_4 * pc_y[k] * pdg_212[k];

        t_299[k] = pb_x[k] * pph0_173[k]
                   + f_11 * ppg_125[k]
                   - f_9 * pc_x[k] * pph1_173[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pb_x, pc_x, pc_y, pph0_175, ppg_126, ppg_127, \
                         pph1_175, pdf0_146, pdf1_146, pdg_215, \
                         pdg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_0 * ppg_126[k]
                   + f_5 * pdf0_146[k]
                   - f_6 * pdf1_146[k]
                   + f_4 * pc_x[k] * pdg_216[k];

        t_301[k] = pb_x[k] * pph0_175[k]
                   + f_1 * ppg_127[k]
                   - f_9 * pc_x[k] * pph1_175[k];

        t_302[k] = f_4 * pc_y[k] * pdg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pb_x, pc_x, pph0_177, ppg_129, ppg_130, \
                         ppg_131, ppg_132, pph1_177, pdg_220, pdg_221, \
                         pdg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pb_x[k] * pph0_177[k]
                   + f_1 * ppg_129[k]
                   - f_9 * pc_x[k] * pph1_177[k];

        t_304[k] = f_0 * ppg_130[k]
                   + f_4 * pc_x[k] * pdg_220[k];

        t_305[k] = f_0 * ppg_131[k]
                   + f_4 * pc_x[k] * pdg_221[k];

        t_306[k] = f_0 * ppg_132[k]
                   + f_4 * pc_x[k] * pdg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pb_x, pc_x, pph0_183, pph0_184, ppg_133, \
                         ppg_134, pph1_183, pph1_184, pdg_223, \
                         pdg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_0 * ppg_133[k]
                   + f_4 * pc_x[k] * pdg_223[k];

        t_308[k] = f_0 * ppg_134[k]
                   + f_4 * pc_x[k] * pdg_224[k];

        t_309[k] = pb_x[k] * pph0_183[k]
                   - f_9 * pc_x[k] * pph1_183[k];

        t_310[k] = pb_x[k] * pph0_184[k]
                   - f_9 * pc_x[k] * pph1_184[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_x, pc_x, pc_y, pph0_185, pph0_186, \
                         pph0_188, pph1_185, pph1_186, pph1_188, \
                         pdg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = pb_x[k] * pph0_185[k]
                   - f_9 * pc_x[k] * pph1_185[k];

        t_312[k] = pb_x[k] * pph0_186[k]
                   - f_9 * pc_x[k] * pph1_186[k];

        t_313[k] = f_4 * pc_y[k] * pdg_224[k];

        t_314[k] = pb_x[k] * pph0_188[k]
                   - f_9 * pc_x[k] * pph1_188[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_z, pc_x, pc_y, pc_z, sdh0_63, sdh1_63, \
                         ppg_105, pdf0_152, pdf1_152, pdg_225, \
                         pdg_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * sdh0_63[k]
                   - f_9 * pc_z[k] * sdh1_63[k];

        t_316[k] = f_1 * ppg_105[k]
                   + f_4 * pc_y[k] * pdg_225[k];

        t_317[k] = f_13 * pdf0_152[k]
                   - f_14 * pdf1_152[k]
                   + f_4 * pc_x[k] * pdg_227[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_z, pc_x, pc_y, pc_z, sdh0_66, sdh1_66, \
                         ppg_107, pdf0_155, pdf1_155, pdg_227, \
                         pdg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * sdh0_66[k]
                   - f_9 * pc_z[k] * sdh1_66[k];

        t_319[k] = f_1 * ppg_107[k]
                   + f_4 * pc_y[k] * pdg_227[k];

        t_320[k] = f_7 * pdf0_155[k]
                   - f_8 * pdf1_155[k]
                   + f_4 * pc_x[k] * pdg_230[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_z, pc_x, pc_y, pc_z, sdh0_69, sdh1_69, \
                         ppg_110, pdf0_157, pdf1_157, pdg_230, \
                         pdg_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = pa_z[k] * sdh0_69[k]
                   - f_9 * pc_z[k] * sdh1_69[k];

        t_322[k] = f_5 * pdf0_157[k]
                   - f_6 * pdf1_157[k]
                   + f_4 * pc_x[k] * pdg_232[k];

        t_323[k] = f_1 * ppg_110[k]
                   + f_4 * pc_y[k] * pdg_230[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pc_x, pdf0_159, pdf1_159, \
                         pdg_234, pdg_235, pdg_236, pdg_237, pdg_238, \
                         pdg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_5 * pdf0_159[k]
                   - f_6 * pdf1_159[k]
                   + f_4 * pc_x[k] * pdg_234[k];

        t_325[k] = f_4 * pc_x[k] * pdg_235[k];

        t_326[k] = f_4 * pc_x[k] * pdg_236[k];

        t_327[k] = f_4 * pc_x[k] * pdg_237[k];

        t_328[k] = f_4 * pc_x[k] * pdg_238[k];

        t_329[k] = f_4 * pc_x[k] * pdg_239[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pa_z, pc_z, sdh0_78, sdh0_79, sdh0_80, sdg_55, \
                         sdg_56, sdh1_78, sdh1_79, sdh1_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = pa_z[k] * sdh0_78[k]
                   - f_9 * pc_z[k] * sdh1_78[k];

        t_331[k] = pa_z[k] * sdh0_79[k]
                   + f_0 * sdg_55[k]
                   - f_9 * pc_z[k] * sdh1_79[k];

        t_332[k] = pa_z[k] * sdh0_80[k]
                   + f_1 * sdg_56[k]
                   - f_9 * pc_z[k] * sdh1_80[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_z, pc_y, pc_z, sdh0_81, sdh0_83, sdg_57, \
                         sdg_59, sdh1_81, sdh1_83, ppg_119, pdg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * sdh0_81[k]
                   + f_11 * sdg_57[k]
                   - f_9 * pc_z[k] * sdh1_81[k];

        t_334[k] = f_1 * ppg_119[k]
                   + f_4 * pc_y[k] * pdg_239[k];

        t_335[k] = pa_z[k] * sdh0_83[k]
                   + f_10 * sdg_59[k]
                   - f_9 * pc_z[k] * sdh1_83[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pb_y, pc_x, pc_y, pph0_168, pph0_170, \
                         ppg_120, pph1_168, pph1_170, pdf0_163, pdf1_163, pdg_240, \
                         pdg_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pb_y[k] * pph0_168[k]
                   - f_9 * pc_y[k] * pph1_168[k];

        t_337[k] = f_0 * ppg_120[k]
                   + f_4 * pc_y[k] * pdg_240[k];

        t_338[k] = pb_y[k] * pph0_170[k]
                   - f_9 * pc_y[k] * pph1_170[k];

        t_339[k] = f_7 * pdf0_163[k]
                   - f_8 * pdf1_163[k]
                   + f_4 * pc_x[k] * pdg_243[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pb_y, pc_x, pc_y, pph0_173, ppg_122, pph1_173, \
                         pdf0_166, pdf1_166, pdg_242, pdg_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_0 * ppg_122[k]
                   + f_4 * pc_y[k] * pdg_242[k];

        t_341[k] = pb_y[k] * pph0_173[k]
                   - f_9 * pc_y[k] * pph1_173[k];

        t_342[k] = f_5 * pdf0_166[k]
                   - f_6 * pdf1_166[k]
                   + f_4 * pc_x[k] * pdg_246[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pb_y, pc_x, pc_y, pph0_177, ppg_125, \
                         pph1_177, pdf0_167, pdf1_167, pdg_245, pdg_247, \
                         pdg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_5 * pdf0_167[k]
                   - f_6 * pdf1_167[k]
                   + f_4 * pc_x[k] * pdg_247[k];

        t_344[k] = f_0 * ppg_125[k]
                   + f_4 * pc_y[k] * pdg_245[k];

        t_345[k] = pb_y[k] * pph0_177[k]
                   - f_9 * pc_y[k] * pph1_177[k];

        t_346[k] = f_4 * pc_x[k] * pdg_250[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pc_x, pc_y, ppg_130, pdf0_166, \
                         pdf1_166, pdg_250, pdg_251, pdg_252, pdg_253, \
                         pdg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_4 * pc_x[k] * pdg_251[k];

        t_348[k] = f_4 * pc_x[k] * pdg_252[k];

        t_349[k] = f_4 * pc_x[k] * pdg_253[k];

        t_350[k] = f_4 * pc_x[k] * pdg_254[k];

        t_351[k] = f_0 * ppg_130[k]
                   + f_2 * pdf0_166[k]
                   - f_3 * pdf1_166[k]
                   + f_4 * pc_y[k] * pdg_250[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pb_y, pc_y, pph0_184, pph0_185, pph0_186, \
                         ppg_131, ppg_132, ppg_133, pph1_184, pph1_185, \
                         pph1_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = pb_y[k] * pph0_184[k]
                   + f_12 * ppg_131[k]
                   - f_9 * pc_y[k] * pph1_184[k];

        t_353[k] = pb_y[k] * pph0_185[k]
                   + f_11 * ppg_132[k]
                   - f_9 * pc_y[k] * pph1_185[k];

        t_354[k] = pb_y[k] * pph0_186[k]
                   + f_1 * ppg_133[k]
                   - f_9 * pc_y[k] * pph1_186[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pb_y, pc_x, pc_y, pph0_188, ppg_134, \
                         pph1_188, pdf0_170, pdf1_170, pdg_254, \
                         pdg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_0 * ppg_134[k]
                   + f_4 * pc_y[k] * pdg_254[k];

        t_356[k] = pb_y[k] * pph0_188[k]
                   - f_9 * pc_y[k] * pph1_188[k];

        t_357[k] = f_2 * pdf0_170[k]
                   - f_3 * pdf1_170[k]
                   + f_4 * pc_x[k] * pdg_255[k];

        t_358[k] = f_4 * pc_y[k] * pdg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pc_x, pc_y, pdf0_172, pdf0_173, pdf0_175, \
                         pdf1_172, pdf1_173, pdf1_175, pdg_257, pdg_258, \
                         pdg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_13 * pdf0_172[k]
                   - f_14 * pdf1_172[k]
                   + f_4 * pc_x[k] * pdg_257[k];

        t_360[k] = f_7 * pdf0_173[k]
                   - f_8 * pdf1_173[k]
                   + f_4 * pc_x[k] * pdg_258[k];

        t_361[k] = f_4 * pc_y[k] * pdg_257[k];

        t_362[k] = f_7 * pdf0_175[k]
                   - f_8 * pdf1_175[k]
                   + f_4 * pc_x[k] * pdg_260[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pc_x, pc_y, pdf0_176, pdf0_177, pdf0_179, \
                         pdf1_176, pdf1_177, pdf1_179, pdg_260, pdg_261, pdg_262, \
                         pdg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_5 * pdf0_176[k]
                   - f_6 * pdf1_176[k]
                   + f_4 * pc_x[k] * pdg_261[k];

        t_364[k] = f_5 * pdf0_177[k]
                   - f_6 * pdf1_177[k]
                   + f_4 * pc_x[k] * pdg_262[k];

        t_365[k] = f_4 * pc_y[k] * pdg_260[k];

        t_366[k] = f_5 * pdf0_179[k]
                   - f_6 * pdf1_179[k]
                   + f_4 * pc_x[k] * pdg_264[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, t_372, pc_x, pc_y, pdf0_176, \
                         pdf1_176, pdg_265, pdg_266, pdg_267, pdg_268, \
                         pdg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_4 * pc_x[k] * pdg_265[k];

        t_368[k] = f_4 * pc_x[k] * pdg_266[k];

        t_369[k] = f_4 * pc_x[k] * pdg_267[k];

        t_370[k] = f_4 * pc_x[k] * pdg_268[k];

        t_371[k] = f_4 * pc_x[k] * pdg_269[k];

        t_372[k] = f_2 * pdf0_176[k]
                   - f_3 * pdf1_176[k]
                   + f_4 * pc_y[k] * pdg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pc_y, pdf0_177, pdf0_178, pdf0_179, \
                         pdf1_177, pdf1_178, pdf1_179, pdg_266, pdg_267, pdg_268, \
                         pdg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_13 * pdf0_177[k]
                   - f_14 * pdf1_177[k]
                   + f_4 * pc_y[k] * pdg_266[k];

        t_374[k] = f_7 * pdf0_178[k]
                   - f_8 * pdf1_178[k]
                   + f_4 * pc_y[k] * pdg_267[k];

        t_375[k] = f_5 * pdf0_179[k]
                   - f_6 * pdf1_179[k]
                   + f_4 * pc_y[k] * pdg_268[k];

        t_376[k] = f_4 * pc_y[k] * pdg_269[k];
    }

#pragma omp simd aligned(t_377, pc_z, sdg_89, ppg_134, pdf0_179, pdf1_179, \
                         pdg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_0 * sdg_89[k]
                   + f_1 * ppg_134[k]
                   + f_2 * pdf0_179[k]
                   - f_3 * pdf1_179[k]
                   + f_4 * pc_z[k] * pdg_269[k];
    }
}

auto
compute_prim_pdh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sdh0,
                                                   const size_t sdg, const size_t sdh1,
                                                   const size_t pph0, const size_t ppg,
                                                   const size_t pph1, const size_t pdf0,
                                                   const size_t pdf1, const size_t pdg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pdh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sdh0,
                                                              sdg, sdh1, pph0, ppg, pph1, pdf0,
                                                              pdf1, pdg, ncols, gamma, p, q);

    compute_prim_pdh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sdh0,
                                                              sdg, sdh1, pph0, ppg, pph1, pdf0,
                                                              pdf1, pdg, ncols, gamma, p, q);

    compute_prim_pdh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sdh0,
                                                              sdg, sdh1, pph0, ppg, pph1, pdf0,
                                                              pdf1, pdg, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
