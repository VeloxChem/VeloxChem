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


#include "SimdThreeCenterElectronRepulsionVrrRecPDG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pdg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdg0, const size_t sdf,
                                                          const size_t sdg1, const size_t ppg0,
                                                          const size_t ppf, const size_t ppg1,
                                                          const size_t pdd0, const size_t pdd1,
                                                          const size_t pdf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 2.0 / q;
    const auto f_9 = 1.5 / q;

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

    const auto *sdg0_0 = buffer.data(sdg0 + 0);
    const auto *sdg0_1 = buffer.data(sdg0 + 1);
    const auto *sdg0_3 = buffer.data(sdg0 + 3);
    const auto *sdg0_5 = buffer.data(sdg0 + 5);
    const auto *sdg0_10 = buffer.data(sdg0 + 10);
    const auto *sdg0_14 = buffer.data(sdg0 + 14);
    const auto *sdg0_30 = buffer.data(sdg0 + 30);
    const auto *sdg0_35 = buffer.data(sdg0 + 35);
    const auto *sdg0_45 = buffer.data(sdg0 + 45);
    const auto *sdg0_48 = buffer.data(sdg0 + 48);
    const auto *sdg0_55 = buffer.data(sdg0 + 55);
    const auto *sdg0_57 = buffer.data(sdg0 + 57);
    const auto *sdg0_59 = buffer.data(sdg0 + 59);
    const auto *sdg0_70 = buffer.data(sdg0 + 70);
    const auto *sdg0_72 = buffer.data(sdg0 + 72);
    const auto *sdg0_74 = buffer.data(sdg0 + 74);
    const auto *sdg0_75 = buffer.data(sdg0 + 75);
    const auto *sdg0_80 = buffer.data(sdg0 + 80);
    const auto *sdg0_85 = buffer.data(sdg0 + 85);
    const auto *sdg0_87 = buffer.data(sdg0 + 87);
    const auto *sdg0_89 = buffer.data(sdg0 + 89);

    const auto *sdf_0 = buffer.data(sdf + 0);
    const auto *sdf_1 = buffer.data(sdf + 1);
    const auto *sdf_6 = buffer.data(sdf + 6);
    const auto *sdf_9 = buffer.data(sdf + 9);
    const auto *sdf_16 = buffer.data(sdf + 16);
    const auto *sdf_29 = buffer.data(sdf + 29);
    const auto *sdf_30 = buffer.data(sdf + 30);
    const auto *sdf_33 = buffer.data(sdf + 33);
    const auto *sdf_36 = buffer.data(sdf + 36);
    const auto *sdf_39 = buffer.data(sdf + 39);
    const auto *sdf_46 = buffer.data(sdf + 46);
    const auto *sdf_49 = buffer.data(sdf + 49);
    const auto *sdf_50 = buffer.data(sdf + 50);
    const auto *sdf_55 = buffer.data(sdf + 55);
    const auto *sdf_56 = buffer.data(sdf + 56);
    const auto *sdf_59 = buffer.data(sdf + 59);

    const auto *sdg1_0 = buffer.data(sdg1 + 0);
    const auto *sdg1_1 = buffer.data(sdg1 + 1);
    const auto *sdg1_3 = buffer.data(sdg1 + 3);
    const auto *sdg1_5 = buffer.data(sdg1 + 5);
    const auto *sdg1_10 = buffer.data(sdg1 + 10);
    const auto *sdg1_14 = buffer.data(sdg1 + 14);
    const auto *sdg1_30 = buffer.data(sdg1 + 30);
    const auto *sdg1_35 = buffer.data(sdg1 + 35);
    const auto *sdg1_45 = buffer.data(sdg1 + 45);
    const auto *sdg1_48 = buffer.data(sdg1 + 48);
    const auto *sdg1_55 = buffer.data(sdg1 + 55);
    const auto *sdg1_57 = buffer.data(sdg1 + 57);
    const auto *sdg1_59 = buffer.data(sdg1 + 59);
    const auto *sdg1_70 = buffer.data(sdg1 + 70);
    const auto *sdg1_72 = buffer.data(sdg1 + 72);
    const auto *sdg1_74 = buffer.data(sdg1 + 74);
    const auto *sdg1_75 = buffer.data(sdg1 + 75);
    const auto *sdg1_80 = buffer.data(sdg1 + 80);
    const auto *sdg1_85 = buffer.data(sdg1 + 85);
    const auto *sdg1_87 = buffer.data(sdg1 + 87);
    const auto *sdg1_89 = buffer.data(sdg1 + 89);

    const auto *ppg0_0 = buffer.data(ppg0 + 0);
    const auto *ppg0_3 = buffer.data(ppg0 + 3);
    const auto *ppg0_5 = buffer.data(ppg0 + 5);
    const auto *ppg0_6 = buffer.data(ppg0 + 6);
    const auto *ppg0_9 = buffer.data(ppg0 + 9);
    const auto *ppg0_18 = buffer.data(ppg0 + 18);
    const auto *ppg0_30 = buffer.data(ppg0 + 30);
    const auto *ppg0_35 = buffer.data(ppg0 + 35);
    const auto *ppg0_46 = buffer.data(ppg0 + 46);
    const auto *ppg0_48 = buffer.data(ppg0 + 48);
    const auto *ppg0_61 = buffer.data(ppg0 + 61);
    const auto *ppg0_63 = buffer.data(ppg0 + 63);
    const auto *ppg0_70 = buffer.data(ppg0 + 70);
    const auto *ppg0_72 = buffer.data(ppg0 + 72);
    const auto *ppg0_73 = buffer.data(ppg0 + 73);
    const auto *ppg0_74 = buffer.data(ppg0 + 74);

    const auto *ppf_0 = buffer.data(ppf + 0);
    const auto *ppf_1 = buffer.data(ppf + 1);
    const auto *ppf_2 = buffer.data(ppf + 2);
    const auto *ppf_3 = buffer.data(ppf + 3);
    const auto *ppf_5 = buffer.data(ppf + 5);
    const auto *ppf_6 = buffer.data(ppf + 6);
    const auto *ppf_8 = buffer.data(ppf + 8);
    const auto *ppf_9 = buffer.data(ppf + 9);
    const auto *ppf_10 = buffer.data(ppf + 10);
    const auto *ppf_12 = buffer.data(ppf + 12);
    const auto *ppf_13 = buffer.data(ppf + 13);
    const auto *ppf_15 = buffer.data(ppf + 15);
    const auto *ppf_16 = buffer.data(ppf + 16);
    const auto *ppf_19 = buffer.data(ppf + 19);
    const auto *ppf_20 = buffer.data(ppf + 20);
    const auto *ppf_22 = buffer.data(ppf + 22);
    const auto *ppf_23 = buffer.data(ppf + 23);
    const auto *ppf_25 = buffer.data(ppf + 25);
    const auto *ppf_26 = buffer.data(ppf + 26);
    const auto *ppf_29 = buffer.data(ppf + 29);
    const auto *ppf_30 = buffer.data(ppf + 30);
    const auto *ppf_31 = buffer.data(ppf + 31);
    const auto *ppf_36 = buffer.data(ppf + 36);
    const auto *ppf_37 = buffer.data(ppf + 37);
    const auto *ppf_38 = buffer.data(ppf + 38);
    const auto *ppf_39 = buffer.data(ppf + 39);
    const auto *ppf_40 = buffer.data(ppf + 40);
    const auto *ppf_41 = buffer.data(ppf + 41);
    const auto *ppf_43 = buffer.data(ppf + 43);
    const auto *ppf_45 = buffer.data(ppf + 45);
    const auto *ppf_46 = buffer.data(ppf + 46);
    const auto *ppf_47 = buffer.data(ppf + 47);
    const auto *ppf_48 = buffer.data(ppf + 48);
    const auto *ppf_49 = buffer.data(ppf + 49);
    const auto *ppf_56 = buffer.data(ppf + 56);
    const auto *ppf_57 = buffer.data(ppf + 57);
    const auto *ppf_58 = buffer.data(ppf + 58);
    const auto *ppf_59 = buffer.data(ppf + 59);

    const auto *ppg1_0 = buffer.data(ppg1 + 0);
    const auto *ppg1_3 = buffer.data(ppg1 + 3);
    const auto *ppg1_5 = buffer.data(ppg1 + 5);
    const auto *ppg1_6 = buffer.data(ppg1 + 6);
    const auto *ppg1_9 = buffer.data(ppg1 + 9);
    const auto *ppg1_18 = buffer.data(ppg1 + 18);
    const auto *ppg1_30 = buffer.data(ppg1 + 30);
    const auto *ppg1_35 = buffer.data(ppg1 + 35);
    const auto *ppg1_46 = buffer.data(ppg1 + 46);
    const auto *ppg1_48 = buffer.data(ppg1 + 48);
    const auto *ppg1_61 = buffer.data(ppg1 + 61);
    const auto *ppg1_63 = buffer.data(ppg1 + 63);
    const auto *ppg1_70 = buffer.data(ppg1 + 70);
    const auto *ppg1_72 = buffer.data(ppg1 + 72);
    const auto *ppg1_73 = buffer.data(ppg1 + 73);
    const auto *ppg1_74 = buffer.data(ppg1 + 74);

    const auto *pdd0_0 = buffer.data(pdd0 + 0);
    const auto *pdd0_3 = buffer.data(pdd0 + 3);
    const auto *pdd0_5 = buffer.data(pdd0 + 5);
    const auto *pdd0_6 = buffer.data(pdd0 + 6);
    const auto *pdd0_9 = buffer.data(pdd0 + 9);
    const auto *pdd0_11 = buffer.data(pdd0 + 11);
    const auto *pdd0_12 = buffer.data(pdd0 + 12);
    const auto *pdd0_15 = buffer.data(pdd0 + 15);
    const auto *pdd0_17 = buffer.data(pdd0 + 17);
    const auto *pdd0_18 = buffer.data(pdd0 + 18);
    const auto *pdd0_30 = buffer.data(pdd0 + 30);
    const auto *pdd0_39 = buffer.data(pdd0 + 39);
    const auto *pdd0_42 = buffer.data(pdd0 + 42);
    const auto *pdd0_47 = buffer.data(pdd0 + 47);

    const auto *pdd1_0 = buffer.data(pdd1 + 0);
    const auto *pdd1_3 = buffer.data(pdd1 + 3);
    const auto *pdd1_5 = buffer.data(pdd1 + 5);
    const auto *pdd1_6 = buffer.data(pdd1 + 6);
    const auto *pdd1_9 = buffer.data(pdd1 + 9);
    const auto *pdd1_11 = buffer.data(pdd1 + 11);
    const auto *pdd1_12 = buffer.data(pdd1 + 12);
    const auto *pdd1_15 = buffer.data(pdd1 + 15);
    const auto *pdd1_17 = buffer.data(pdd1 + 17);
    const auto *pdd1_18 = buffer.data(pdd1 + 18);
    const auto *pdd1_30 = buffer.data(pdd1 + 30);
    const auto *pdd1_39 = buffer.data(pdd1 + 39);
    const auto *pdd1_42 = buffer.data(pdd1 + 42);
    const auto *pdd1_47 = buffer.data(pdd1 + 47);

    const auto *pdf_0 = buffer.data(pdf + 0);
    const auto *pdf_1 = buffer.data(pdf + 1);
    const auto *pdf_2 = buffer.data(pdf + 2);
    const auto *pdf_3 = buffer.data(pdf + 3);
    const auto *pdf_5 = buffer.data(pdf + 5);
    const auto *pdf_6 = buffer.data(pdf + 6);
    const auto *pdf_8 = buffer.data(pdf + 8);
    const auto *pdf_9 = buffer.data(pdf + 9);
    const auto *pdf_10 = buffer.data(pdf + 10);
    const auto *pdf_11 = buffer.data(pdf + 11);
    const auto *pdf_12 = buffer.data(pdf + 12);
    const auto *pdf_13 = buffer.data(pdf + 13);
    const auto *pdf_15 = buffer.data(pdf + 15);
    const auto *pdf_16 = buffer.data(pdf + 16);
    const auto *pdf_18 = buffer.data(pdf + 18);
    const auto *pdf_19 = buffer.data(pdf + 19);
    const auto *pdf_20 = buffer.data(pdf + 20);
    const auto *pdf_22 = buffer.data(pdf + 22);
    const auto *pdf_23 = buffer.data(pdf + 23);
    const auto *pdf_25 = buffer.data(pdf + 25);
    const auto *pdf_26 = buffer.data(pdf + 26);
    const auto *pdf_28 = buffer.data(pdf + 28);
    const auto *pdf_29 = buffer.data(pdf + 29);
    const auto *pdf_30 = buffer.data(pdf + 30);
    const auto *pdf_32 = buffer.data(pdf + 32);
    const auto *pdf_33 = buffer.data(pdf + 33);
    const auto *pdf_35 = buffer.data(pdf + 35);
    const auto *pdf_36 = buffer.data(pdf + 36);
    const auto *pdf_39 = buffer.data(pdf + 39);
    const auto *pdf_40 = buffer.data(pdf + 40);
    const auto *pdf_42 = buffer.data(pdf + 42);
    const auto *pdf_43 = buffer.data(pdf + 43);
    const auto *pdf_45 = buffer.data(pdf + 45);
    const auto *pdf_46 = buffer.data(pdf + 46);
    const auto *pdf_49 = buffer.data(pdf + 49);
    const auto *pdf_50 = buffer.data(pdf + 50);
    const auto *pdf_51 = buffer.data(pdf + 51);
    const auto *pdf_52 = buffer.data(pdf + 52);
    const auto *pdf_53 = buffer.data(pdf + 53);
    const auto *pdf_55 = buffer.data(pdf + 55);
    const auto *pdf_56 = buffer.data(pdf + 56);
    const auto *pdf_59 = buffer.data(pdf + 59);
    const auto *pdf_60 = buffer.data(pdf + 60);
    const auto *pdf_61 = buffer.data(pdf + 61);
    const auto *pdf_66 = buffer.data(pdf + 66);
    const auto *pdf_67 = buffer.data(pdf + 67);
    const auto *pdf_68 = buffer.data(pdf + 68);
    const auto *pdf_69 = buffer.data(pdf + 69);
    const auto *pdf_70 = buffer.data(pdf + 70);
    const auto *pdf_71 = buffer.data(pdf + 71);
    const auto *pdf_75 = buffer.data(pdf + 75);
    const auto *pdf_76 = buffer.data(pdf + 76);
    const auto *pdf_77 = buffer.data(pdf + 77);
    const auto *pdf_78 = buffer.data(pdf + 78);
    const auto *pdf_79 = buffer.data(pdf + 79);
    const auto *pdf_80 = buffer.data(pdf + 80);
    const auto *pdf_81 = buffer.data(pdf + 81);
    const auto *pdf_86 = buffer.data(pdf + 86);
    const auto *pdf_87 = buffer.data(pdf + 87);
    const auto *pdf_88 = buffer.data(pdf + 88);
    const auto *pdf_89 = buffer.data(pdf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sdf_0, ppf_0, pdd0_0, \
                         pdd1_0, pdf_0, pdf_1, pdf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdf_0[k]
                 + f_1 * ppf_0[k]
                 + f_2 * pdd0_0[k]
                 - f_3 * pdd1_0[k]
                 + f_4 * pc_x[k] * pdf_0[k];

        t_1[k] = f_4 * pc_y[k] * pdf_0[k];

        t_2[k] = f_4 * pc_z[k] * pdf_0[k];

        t_3[k] = f_5 * pdd0_0[k]
                 - f_6 * pdd1_0[k]
                 + f_4 * pc_y[k] * pdf_1[k];

        t_4[k] = f_4 * pc_y[k] * pdf_2[k];

        t_5[k] = f_5 * pdd0_0[k]
                 - f_6 * pdd1_0[k]
                 + f_4 * pc_z[k] * pdf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sdf_6, sdf_9, ppf_6, ppf_9, \
                         pdf_3, pdf_5, pdf_6, pdf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * sdf_6[k]
                 + f_1 * ppf_6[k]
                 + f_4 * pc_x[k] * pdf_6[k];

        t_7[k] = f_4 * pc_z[k] * pdf_3[k];

        t_8[k] = f_4 * pc_y[k] * pdf_5[k];

        t_9[k] = f_0 * sdf_9[k]
                 + f_1 * ppf_9[k]
                 + f_4 * pc_x[k] * pdf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, pdd0_3, pdd0_5, pdd1_3, \
                         pdd1_5, pdf_6, pdf_8, pdf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * pdd0_3[k]
                  - f_3 * pdd1_3[k]
                  + f_4 * pc_y[k] * pdf_6[k];

        t_11[k] = f_4 * pc_z[k] * pdf_6[k];

        t_12[k] = f_5 * pdd0_5[k]
                  - f_6 * pdd1_5[k]
                  + f_4 * pc_y[k] * pdf_8[k];

        t_13[k] = f_4 * pc_y[k] * pdf_9[k];

        t_14[k] = f_2 * pdd0_5[k]
                  - f_3 * pdd1_5[k]
                  + f_4 * pc_z[k] * pdf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pc_y, pc_z, ppg0_0, ppf_0, ppf_1, \
                         ppg1_0, pdd0_6, pdd1_6, pdf_10, pdf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * ppg0_0[k]
                  - f_7 * pc_y[k] * ppg1_0[k];

        t_16[k] = f_0 * ppf_0[k]
                  + f_4 * pc_y[k] * pdf_10[k];

        t_17[k] = f_4 * pc_z[k] * pdf_10[k];

        t_18[k] = f_0 * ppf_1[k]
                  + f_5 * pdd0_6[k]
                  - f_6 * pdd1_6[k]
                  + f_4 * pc_y[k] * pdf_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, sdf_16, ppg0_5, \
                         ppf_2, ppf_16, ppg1_5, pdf_12, pdf_13, \
                         pdf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * ppf_2[k]
                  + f_4 * pc_y[k] * pdf_12[k];

        t_20[k] = pb_y[k] * ppg0_5[k]
                  - f_7 * pc_y[k] * ppg1_5[k];

        t_21[k] = f_0 * sdf_16[k]
                  + f_0 * ppf_16[k]
                  + f_4 * pc_x[k] * pdf_16[k];

        t_22[k] = f_4 * pc_z[k] * pdf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, ppg0_9, ppf_5, ppf_6, \
                         ppg1_9, pdd0_9, pdd1_9, pdf_15, pdf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * ppf_5[k]
                  + f_4 * pc_y[k] * pdf_15[k];

        t_24[k] = pb_y[k] * ppg0_9[k]
                  - f_7 * pc_y[k] * ppg1_9[k];

        t_25[k] = f_0 * ppf_6[k]
                  + f_2 * pdd0_9[k]
                  - f_3 * pdd1_9[k]
                  + f_4 * pc_y[k] * pdf_16[k];

        t_26[k] = f_4 * pc_z[k] * pdf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_z, pc_y, pc_z, ppg0_0, ppf_8, ppf_9, \
                         ppg1_0, pdd0_11, pdd1_11, pdf_18, pdf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * ppf_8[k]
                  + f_5 * pdd0_11[k]
                  - f_6 * pdd1_11[k]
                  + f_4 * pc_y[k] * pdf_18[k];

        t_28[k] = f_0 * ppf_9[k]
                  + f_4 * pc_y[k] * pdf_19[k];

        t_29[k] = f_2 * pdd0_11[k]
                  - f_3 * pdd1_11[k]
                  + f_4 * pc_z[k] * pdf_19[k];

        t_30[k] = pb_z[k] * ppg0_0[k]
                  - f_7 * pc_z[k] * ppg1_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_z, pc_y, pc_z, ppg0_3, ppf_0, ppf_2, \
                         ppg1_3, pdd0_12, pdd1_12, pdf_20, pdf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_4 * pc_y[k] * pdf_20[k];

        t_32[k] = f_0 * ppf_0[k]
                  + f_4 * pc_z[k] * pdf_20[k];

        t_33[k] = pb_z[k] * ppg0_3[k]
                  - f_7 * pc_z[k] * ppg1_3[k];

        t_34[k] = f_4 * pc_y[k] * pdf_22[k];

        t_35[k] = f_0 * ppf_2[k]
                  + f_5 * pdd0_12[k]
                  - f_6 * pdd1_12[k]
                  + f_4 * pc_z[k] * pdf_22[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_z, pc_x, pc_y, pc_z, sdf_29, ppg0_6, \
                         ppf_3, ppf_29, ppg1_6, pdf_23, pdf_25, \
                         pdf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_z[k] * ppg0_6[k]
                  - f_7 * pc_z[k] * ppg1_6[k];

        t_37[k] = f_0 * ppf_3[k]
                  + f_4 * pc_z[k] * pdf_23[k];

        t_38[k] = f_4 * pc_y[k] * pdf_25[k];

        t_39[k] = f_0 * sdf_29[k]
                  + f_0 * ppf_29[k]
                  + f_4 * pc_x[k] * pdf_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pc_y, pc_z, ppf_6, ppf_9, pdd0_15, \
                         pdd0_17, pdd1_15, pdd1_17, pdf_26, pdf_28, \
                         pdf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * pdd0_15[k]
                  - f_3 * pdd1_15[k]
                  + f_4 * pc_y[k] * pdf_26[k];

        t_41[k] = f_0 * ppf_6[k]
                  + f_4 * pc_z[k] * pdf_26[k];

        t_42[k] = f_5 * pdd0_17[k]
                  - f_6 * pdd1_17[k]
                  + f_4 * pc_y[k] * pdf_28[k];

        t_43[k] = f_4 * pc_y[k] * pdf_29[k];

        t_44[k] = f_0 * ppf_9[k]
                  + f_2 * pdd0_17[k]
                  - f_3 * pdd1_17[k]
                  + f_4 * pc_z[k] * pdf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pc_x, pc_y, pc_z, sdg0_45, sdg0_48, \
                         sdf_30, sdf_33, sdg1_45, sdg1_48, ppf_10, \
                         pdf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_x[k] * sdg0_45[k]
                  + f_8 * sdf_30[k]
                  - f_7 * pc_x[k] * sdg1_45[k];

        t_46[k] = f_1 * ppf_10[k]
                  + f_4 * pc_y[k] * pdf_30[k];

        t_47[k] = f_4 * pc_z[k] * pdf_30[k];

        t_48[k] = pa_x[k] * sdg0_48[k]
                  + f_1 * sdf_33[k]
                  - f_7 * pc_x[k] * sdg1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, pc_z, sdf_36, ppf_12, pdd0_18, \
                         pdd1_18, pdf_32, pdf_33, pdf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_1 * ppf_12[k]
                  + f_4 * pc_y[k] * pdf_32[k];

        t_50[k] = f_5 * pdd0_18[k]
                  - f_6 * pdd1_18[k]
                  + f_4 * pc_z[k] * pdf_32[k];

        t_51[k] = f_0 * sdf_36[k]
                  + f_4 * pc_x[k] * pdf_36[k];

        t_52[k] = f_4 * pc_z[k] * pdf_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pc_x, pc_y, pc_z, sdg0_55, sdf_39, \
                         sdg1_55, ppf_15, pdf_35, pdf_36, pdf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_1 * ppf_15[k]
                  + f_4 * pc_y[k] * pdf_35[k];

        t_54[k] = f_0 * sdf_39[k]
                  + f_4 * pc_x[k] * pdf_39[k];

        t_55[k] = pa_x[k] * sdg0_55[k]
                  - f_7 * pc_x[k] * sdg1_55[k];

        t_56[k] = f_4 * pc_z[k] * pdf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_y, pc_x, pc_y, sdg0_57, sdg0_59, \
                         sdg1_57, sdg1_59, ppg0_30, ppf_19, ppg1_30, \
                         pdf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_x[k] * sdg0_57[k]
                  - f_7 * pc_x[k] * sdg1_57[k];

        t_58[k] = f_1 * ppf_19[k]
                  + f_4 * pc_y[k] * pdf_39[k];

        t_59[k] = pa_x[k] * sdg0_59[k]
                  - f_7 * pc_x[k] * sdg1_59[k];

        t_60[k] = pb_y[k] * ppg0_30[k]
                  - f_7 * pc_y[k] * ppg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, ppg0_18, ppf_10, ppf_20, \
                         ppf_22, ppg1_18, pdf_40, pdf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * ppf_20[k]
                  + f_4 * pc_y[k] * pdf_40[k];

        t_62[k] = f_0 * ppf_10[k]
                  + f_4 * pc_z[k] * pdf_40[k];

        t_63[k] = pb_z[k] * ppg0_18[k]
                  - f_7 * pc_z[k] * ppg1_18[k];

        t_64[k] = f_0 * ppf_22[k]
                  + f_4 * pc_y[k] * pdf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, pc_z, sdf_46, ppg0_35, \
                         ppf_13, ppf_25, ppg1_35, pdf_43, pdf_45, \
                         pdf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * ppg0_35[k]
                  - f_7 * pc_y[k] * ppg1_35[k];

        t_66[k] = f_0 * sdf_46[k]
                  + f_4 * pc_x[k] * pdf_46[k];

        t_67[k] = f_0 * ppf_13[k]
                  + f_4 * pc_z[k] * pdf_43[k];

        t_68[k] = f_0 * ppf_25[k]
                  + f_4 * pc_y[k] * pdf_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pc_x, pc_z, sdg0_70, sdg0_72, sdf_49, \
                         sdg1_70, sdg1_72, ppf_16, pdf_46, pdf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_0 * sdf_49[k]
                  + f_4 * pc_x[k] * pdf_49[k];

        t_70[k] = pa_x[k] * sdg0_70[k]
                  - f_7 * pc_x[k] * sdg1_70[k];

        t_71[k] = f_0 * ppf_16[k]
                  + f_4 * pc_z[k] * pdf_46[k];

        t_72[k] = pa_x[k] * sdg0_72[k]
                  - f_7 * pc_x[k] * sdg1_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pc_x, pc_y, sdg0_74, sdg0_75, sdf_50, \
                         sdg1_74, sdg1_75, ppf_29, pdf_49, pdf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * ppf_29[k]
                  + f_4 * pc_y[k] * pdf_49[k];

        t_74[k] = pa_x[k] * sdg0_74[k]
                  - f_7 * pc_x[k] * sdg1_74[k];

        t_75[k] = pa_x[k] * sdg0_75[k]
                  + f_8 * sdf_50[k]
                  - f_7 * pc_x[k] * sdg1_75[k];

        t_76[k] = f_4 * pc_y[k] * pdf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, ppf_20, pdd0_30, pdd1_30, pdf_50, \
                         pdf_51, pdf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * ppf_20[k]
                  + f_4 * pc_z[k] * pdf_50[k];

        t_78[k] = f_5 * pdd0_30[k]
                  - f_6 * pdd1_30[k]
                  + f_4 * pc_y[k] * pdf_51[k];

        t_79[k] = f_4 * pc_y[k] * pdf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_x, pc_x, pc_y, pc_z, sdg0_80, sdf_55, \
                         sdf_56, sdg1_80, ppf_23, pdf_53, pdf_55, \
                         pdf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_x[k] * sdg0_80[k]
                  + f_1 * sdf_55[k]
                  - f_7 * pc_x[k] * sdg1_80[k];

        t_81[k] = f_0 * sdf_56[k]
                  + f_4 * pc_x[k] * pdf_56[k];

        t_82[k] = f_1 * ppf_23[k]
                  + f_4 * pc_z[k] * pdf_53[k];

        t_83[k] = f_4 * pc_y[k] * pdf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pc_x, pc_z, sdg0_85, sdg0_87, sdf_59, \
                         sdg1_85, sdg1_87, ppf_26, pdf_56, pdf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_0 * sdf_59[k]
                  + f_4 * pc_x[k] * pdf_59[k];

        t_85[k] = pa_x[k] * sdg0_85[k]
                  - f_7 * pc_x[k] * sdg1_85[k];

        t_86[k] = f_1 * ppf_26[k]
                  + f_4 * pc_z[k] * pdf_56[k];

        t_87[k] = pa_x[k] * sdg0_87[k]
                  - f_7 * pc_x[k] * sdg1_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pc_x, pc_y, sdg0_0, sdg0_1, \
                         sdg0_89, sdf_0, sdg1_0, sdg1_1, sdg1_89, \
                         pdf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_4 * pc_y[k] * pdf_59[k];

        t_89[k] = pa_x[k] * sdg0_89[k]
                  - f_7 * pc_x[k] * sdg1_89[k];

        t_90[k] = pa_y[k] * sdg0_0[k]
                  - f_7 * pc_y[k] * sdg1_0[k];

        t_91[k] = pa_y[k] * sdg0_1[k]
                  + f_0 * sdf_0[k]
                  - f_7 * pc_y[k] * sdg1_1[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pc_y, pc_z, sdg0_3, sdg0_5, sdf_1, \
                         sdg1_3, sdg1_5, pdf_60, pdf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * pc_z[k] * pdf_60[k];

        t_93[k] = pa_y[k] * sdg0_3[k]
                  + f_1 * sdf_1[k]
                  - f_7 * pc_y[k] * sdg1_3[k];

        t_94[k] = f_4 * pc_z[k] * pdf_61[k];

        t_95[k] = pa_y[k] * sdg0_5[k]
                  - f_7 * pc_y[k] * sdg1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, ppf_36, ppf_37, ppf_38, ppf_39, pdf_66, \
                         pdf_67, pdf_68, pdf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * ppf_36[k]
                  + f_4 * pc_x[k] * pdf_66[k];

        t_97[k] = f_1 * ppf_37[k]
                  + f_4 * pc_x[k] * pdf_67[k];

        t_98[k] = f_1 * ppf_38[k]
                  + f_4 * pc_x[k] * pdf_68[k];

        t_99[k] = f_1 * ppf_39[k]
                  + f_4 * pc_x[k] * pdf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pc_y, pc_z, sdg0_10, sdf_6, sdf_9, \
                         sdg1_10, pdd0_39, pdd1_39, pdf_66, pdf_67, \
                         pdf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_y[k] * sdg0_10[k]
                   + f_8 * sdf_6[k]
                   - f_7 * pc_y[k] * sdg1_10[k];

        t_101[k] = f_4 * pc_z[k] * pdf_66[k];

        t_102[k] = f_5 * pdd0_39[k]
                   - f_6 * pdd1_39[k]
                   + f_4 * pc_z[k] * pdf_67[k];

        t_103[k] = f_0 * sdf_9[k]
                   + f_4 * pc_y[k] * pdf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pb_x, pc_x, pc_y, sdg0_14, sdg1_14, \
                         ppg0_61, ppf_40, ppf_41, ppg1_61, pdd0_42, pdd1_42, \
                         pdf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * sdg0_14[k]
                   - f_7 * pc_y[k] * sdg1_14[k];

        t_105[k] = f_0 * ppf_40[k]
                   + f_2 * pdd0_42[k]
                   - f_3 * pdd1_42[k]
                   + f_4 * pc_x[k] * pdf_70[k];

        t_106[k] = pb_x[k] * ppg0_61[k]
                   + f_9 * ppf_41[k]
                   - f_7 * pc_x[k] * ppg1_61[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_x, pc_x, pc_z, ppg0_63, ppf_43, \
                         ppf_45, ppg1_63, pdd0_47, pdd1_47, pdf_70, pdf_71, \
                         pdf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_4 * pc_z[k] * pdf_70[k];

        t_108[k] = pb_x[k] * ppg0_63[k]
                   + f_1 * ppf_43[k]
                   - f_7 * pc_x[k] * ppg1_63[k];

        t_109[k] = f_4 * pc_z[k] * pdf_71[k];

        t_110[k] = f_0 * ppf_45[k]
                   + f_5 * pdd0_47[k]
                   - f_6 * pdd1_47[k]
                   + f_4 * pc_x[k] * pdf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, ppf_46, ppf_47, ppf_48, ppf_49, \
                         pdf_76, pdf_77, pdf_78, pdf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * ppf_46[k]
                   + f_4 * pc_x[k] * pdf_76[k];

        t_112[k] = f_0 * ppf_47[k]
                   + f_4 * pc_x[k] * pdf_77[k];

        t_113[k] = f_0 * ppf_48[k]
                   + f_4 * pc_x[k] * pdf_78[k];

        t_114[k] = f_0 * ppf_49[k]
                   + f_4 * pc_x[k] * pdf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_x, pc_x, pc_z, ppg0_70, ppg0_72, \
                         ppg0_73, ppg1_70, ppg1_72, ppg1_73, pdf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_x[k] * ppg0_70[k]
                   - f_7 * pc_x[k] * ppg1_70[k];

        t_116[k] = f_4 * pc_z[k] * pdf_76[k];

        t_117[k] = pb_x[k] * ppg0_72[k]
                   - f_7 * pc_x[k] * ppg1_72[k];

        t_118[k] = pb_x[k] * ppg0_73[k]
                   - f_7 * pc_x[k] * ppg1_73[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pa_y, pb_x, pb_z, pc_x, pc_y, pc_z, sdg0_30, \
                         sdg1_30, ppg0_46, ppg0_74, ppg1_46, ppg1_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_x[k] * ppg0_74[k]
                   - f_7 * pc_x[k] * ppg1_74[k];

        t_120[k] = pa_y[k] * sdg0_30[k]
                   - f_7 * pc_y[k] * sdg1_30[k];

        t_121[k] = pb_z[k] * ppg0_46[k]
                   - f_7 * pc_z[k] * ppg1_46[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pb_z, pc_y, pc_z, sdg0_35, sdg1_35, \
                         ppg0_48, ppf_30, ppf_31, ppg1_48, pdf_80, \
                         pdf_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * ppf_30[k]
                   + f_4 * pc_z[k] * pdf_80[k];

        t_123[k] = pb_z[k] * ppg0_48[k]
                   - f_7 * pc_z[k] * ppg1_48[k];

        t_124[k] = f_0 * ppf_31[k]
                   + f_4 * pc_z[k] * pdf_81[k];

        t_125[k] = pa_y[k] * sdg0_35[k]
                   - f_7 * pc_y[k] * sdg1_35[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, ppf_56, ppf_57, ppf_58, ppf_59, \
                         pdf_86, pdf_87, pdf_88, pdf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_0 * ppf_56[k]
                   + f_4 * pc_x[k] * pdf_86[k];

        t_127[k] = f_0 * ppf_57[k]
                   + f_4 * pc_x[k] * pdf_87[k];

        t_128[k] = f_0 * ppf_58[k]
                   + f_4 * pc_x[k] * pdf_88[k];

        t_129[k] = f_0 * ppf_59[k]
                   + f_4 * pc_x[k] * pdf_89[k];
    }
}

static auto
compute_prim_pdg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdg0, const size_t sdf,
                                                          const size_t sdg1, const size_t ppg0,
                                                          const size_t ppf, const size_t ppg1,
                                                          const size_t pdd0, const size_t pdd1,
                                                          const size_t pdf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 2.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);

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
    auto *t_256 = buffer.data(target + 256);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdg0_0 = buffer.data(sdg0 + 0);
    const auto *sdg0_2 = buffer.data(sdg0 + 2);
    const auto *sdg0_3 = buffer.data(sdg0 + 3);
    const auto *sdg0_5 = buffer.data(sdg0 + 5);
    const auto *sdg0_10 = buffer.data(sdg0 + 10);
    const auto *sdg0_14 = buffer.data(sdg0 + 14);
    const auto *sdg0_15 = buffer.data(sdg0 + 15);
    const auto *sdg0_18 = buffer.data(sdg0 + 18);
    const auto *sdg0_25 = buffer.data(sdg0 + 25);
    const auto *sdg0_44 = buffer.data(sdg0 + 44);
    const auto *sdg0_45 = buffer.data(sdg0 + 45);
    const auto *sdg0_48 = buffer.data(sdg0 + 48);
    const auto *sdg0_55 = buffer.data(sdg0 + 55);
    const auto *sdg0_56 = buffer.data(sdg0 + 56);
    const auto *sdg0_57 = buffer.data(sdg0 + 57);
    const auto *sdg0_59 = buffer.data(sdg0 + 59);
    const auto *sdg0_75 = buffer.data(sdg0 + 75);
    const auto *sdg0_80 = buffer.data(sdg0 + 80);
    const auto *sdg0_85 = buffer.data(sdg0 + 85);
    const auto *sdg0_87 = buffer.data(sdg0 + 87);
    const auto *sdg0_89 = buffer.data(sdg0 + 89);

    const auto *sdf_0 = buffer.data(sdf + 0);
    const auto *sdf_2 = buffer.data(sdf + 2);
    const auto *sdf_9 = buffer.data(sdf + 9);
    const auto *sdf_29 = buffer.data(sdf + 29);
    const auto *sdf_36 = buffer.data(sdf + 36);
    const auto *sdf_37 = buffer.data(sdf + 37);
    const auto *sdf_39 = buffer.data(sdf + 39);
    const auto *sdf_49 = buffer.data(sdf + 49);
    const auto *sdf_56 = buffer.data(sdf + 56);
    const auto *sdf_58 = buffer.data(sdf + 58);
    const auto *sdf_59 = buffer.data(sdf + 59);

    const auto *sdg1_0 = buffer.data(sdg1 + 0);
    const auto *sdg1_2 = buffer.data(sdg1 + 2);
    const auto *sdg1_3 = buffer.data(sdg1 + 3);
    const auto *sdg1_5 = buffer.data(sdg1 + 5);
    const auto *sdg1_10 = buffer.data(sdg1 + 10);
    const auto *sdg1_14 = buffer.data(sdg1 + 14);
    const auto *sdg1_15 = buffer.data(sdg1 + 15);
    const auto *sdg1_18 = buffer.data(sdg1 + 18);
    const auto *sdg1_25 = buffer.data(sdg1 + 25);
    const auto *sdg1_44 = buffer.data(sdg1 + 44);
    const auto *sdg1_45 = buffer.data(sdg1 + 45);
    const auto *sdg1_48 = buffer.data(sdg1 + 48);
    const auto *sdg1_55 = buffer.data(sdg1 + 55);
    const auto *sdg1_56 = buffer.data(sdg1 + 56);
    const auto *sdg1_57 = buffer.data(sdg1 + 57);
    const auto *sdg1_59 = buffer.data(sdg1 + 59);
    const auto *sdg1_75 = buffer.data(sdg1 + 75);
    const auto *sdg1_80 = buffer.data(sdg1 + 80);
    const auto *sdg1_85 = buffer.data(sdg1 + 85);
    const auto *sdg1_87 = buffer.data(sdg1 + 87);
    const auto *sdg1_89 = buffer.data(sdg1 + 89);

    const auto *ppg0_61 = buffer.data(ppg0 + 61);
    const auto *ppg0_63 = buffer.data(ppg0 + 63);
    const auto *ppg0_70 = buffer.data(ppg0 + 70);
    const auto *ppg0_85 = buffer.data(ppg0 + 85);
    const auto *ppg0_87 = buffer.data(ppg0 + 87);
    const auto *ppg0_92 = buffer.data(ppg0 + 92);
    const auto *ppg0_95 = buffer.data(ppg0 + 95);
    const auto *ppg0_116 = buffer.data(ppg0 + 116);
    const auto *ppg0_117 = buffer.data(ppg0 + 117);
    const auto *ppg0_119 = buffer.data(ppg0 + 119);
    const auto *ppg0_120 = buffer.data(ppg0 + 120);
    const auto *ppg0_122 = buffer.data(ppg0 + 122);
    const auto *ppg0_125 = buffer.data(ppg0 + 125);
    const auto *ppg0_130 = buffer.data(ppg0 + 130);
    const auto *ppg0_131 = buffer.data(ppg0 + 131);
    const auto *ppg0_132 = buffer.data(ppg0 + 132);
    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppf_36 = buffer.data(ppf + 36);
    const auto *ppf_40 = buffer.data(ppf + 40);
    const auto *ppf_41 = buffer.data(ppf + 41);
    const auto *ppf_46 = buffer.data(ppf + 46);
    const auto *ppf_47 = buffer.data(ppf + 47);
    const auto *ppf_49 = buffer.data(ppf + 49);
    const auto *ppf_50 = buffer.data(ppf + 50);
    const auto *ppf_51 = buffer.data(ppf + 51);
    const auto *ppf_56 = buffer.data(ppf + 56);
    const auto *ppf_59 = buffer.data(ppf + 59);
    const auto *ppf_60 = buffer.data(ppf + 60);
    const auto *ppf_62 = buffer.data(ppf + 62);
    const auto *ppf_66 = buffer.data(ppf + 66);
    const auto *ppf_67 = buffer.data(ppf + 67);
    const auto *ppf_68 = buffer.data(ppf + 68);
    const auto *ppf_69 = buffer.data(ppf + 69);
    const auto *ppf_70 = buffer.data(ppf + 70);
    const auto *ppf_72 = buffer.data(ppf + 72);
    const auto *ppf_76 = buffer.data(ppf + 76);
    const auto *ppf_77 = buffer.data(ppf + 77);
    const auto *ppf_78 = buffer.data(ppf + 78);
    const auto *ppf_79 = buffer.data(ppf + 79);
    const auto *ppf_80 = buffer.data(ppf + 80);
    const auto *ppf_82 = buffer.data(ppf + 82);
    const auto *ppf_83 = buffer.data(ppf + 83);
    const auto *ppf_85 = buffer.data(ppf + 85);
    const auto *ppf_86 = buffer.data(ppf + 86);
    const auto *ppf_87 = buffer.data(ppf + 87);
    const auto *ppf_88 = buffer.data(ppf + 88);
    const auto *ppf_89 = buffer.data(ppf + 89);

    const auto *ppg1_61 = buffer.data(ppg1 + 61);
    const auto *ppg1_63 = buffer.data(ppg1 + 63);
    const auto *ppg1_70 = buffer.data(ppg1 + 70);
    const auto *ppg1_85 = buffer.data(ppg1 + 85);
    const auto *ppg1_87 = buffer.data(ppg1 + 87);
    const auto *ppg1_92 = buffer.data(ppg1 + 92);
    const auto *ppg1_95 = buffer.data(ppg1 + 95);
    const auto *ppg1_116 = buffer.data(ppg1 + 116);
    const auto *ppg1_117 = buffer.data(ppg1 + 117);
    const auto *ppg1_119 = buffer.data(ppg1 + 119);
    const auto *ppg1_120 = buffer.data(ppg1 + 120);
    const auto *ppg1_122 = buffer.data(ppg1 + 122);
    const auto *ppg1_125 = buffer.data(ppg1 + 125);
    const auto *ppg1_130 = buffer.data(ppg1 + 130);
    const auto *ppg1_131 = buffer.data(ppg1 + 131);
    const auto *ppg1_132 = buffer.data(ppg1 + 132);
    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *pdd0_54 = buffer.data(pdd0 + 54);
    const auto *pdd0_55 = buffer.data(pdd0 + 55);
    const auto *pdd0_57 = buffer.data(pdd0 + 57);
    const auto *pdd0_59 = buffer.data(pdd0 + 59);
    const auto *pdd0_60 = buffer.data(pdd0 + 60);
    const auto *pdd0_63 = buffer.data(pdd0 + 63);
    const auto *pdd0_65 = buffer.data(pdd0 + 65);
    const auto *pdd0_67 = buffer.data(pdd0 + 67);
    const auto *pdd0_69 = buffer.data(pdd0 + 69);
    const auto *pdd0_76 = buffer.data(pdd0 + 76);
    const auto *pdd0_77 = buffer.data(pdd0 + 77);
    const auto *pdd0_84 = buffer.data(pdd0 + 84);
    const auto *pdd0_87 = buffer.data(pdd0 + 87);
    const auto *pdd0_92 = buffer.data(pdd0 + 92);
    const auto *pdd0_95 = buffer.data(pdd0 + 95);
    const auto *pdd0_99 = buffer.data(pdd0 + 99);
    const auto *pdd0_102 = buffer.data(pdd0 + 102);

    const auto *pdd1_54 = buffer.data(pdd1 + 54);
    const auto *pdd1_55 = buffer.data(pdd1 + 55);
    const auto *pdd1_57 = buffer.data(pdd1 + 57);
    const auto *pdd1_59 = buffer.data(pdd1 + 59);
    const auto *pdd1_60 = buffer.data(pdd1 + 60);
    const auto *pdd1_63 = buffer.data(pdd1 + 63);
    const auto *pdd1_65 = buffer.data(pdd1 + 65);
    const auto *pdd1_67 = buffer.data(pdd1 + 67);
    const auto *pdd1_69 = buffer.data(pdd1 + 69);
    const auto *pdd1_76 = buffer.data(pdd1 + 76);
    const auto *pdd1_77 = buffer.data(pdd1 + 77);
    const auto *pdd1_84 = buffer.data(pdd1 + 84);
    const auto *pdd1_87 = buffer.data(pdd1 + 87);
    const auto *pdd1_92 = buffer.data(pdd1 + 92);
    const auto *pdd1_95 = buffer.data(pdd1 + 95);
    const auto *pdd1_99 = buffer.data(pdd1 + 99);
    const auto *pdd1_102 = buffer.data(pdd1 + 102);

    const auto *pdf_86 = buffer.data(pdf + 86);
    const auto *pdf_89 = buffer.data(pdf + 89);
    const auto *pdf_90 = buffer.data(pdf + 90);
    const auto *pdf_91 = buffer.data(pdf + 91);
    const auto *pdf_93 = buffer.data(pdf + 93);
    const auto *pdf_95 = buffer.data(pdf + 95);
    const auto *pdf_96 = buffer.data(pdf + 96);
    const auto *pdf_97 = buffer.data(pdf + 97);
    const auto *pdf_98 = buffer.data(pdf + 98);
    const auto *pdf_99 = buffer.data(pdf + 99);
    const auto *pdf_100 = buffer.data(pdf + 100);
    const auto *pdf_101 = buffer.data(pdf + 101);
    const auto *pdf_105 = buffer.data(pdf + 105);
    const auto *pdf_106 = buffer.data(pdf + 106);
    const auto *pdf_107 = buffer.data(pdf + 107);
    const auto *pdf_108 = buffer.data(pdf + 108);
    const auto *pdf_109 = buffer.data(pdf + 109);
    const auto *pdf_110 = buffer.data(pdf + 110);
    const auto *pdf_111 = buffer.data(pdf + 111);
    const auto *pdf_113 = buffer.data(pdf + 113);
    const auto *pdf_116 = buffer.data(pdf + 116);
    const auto *pdf_117 = buffer.data(pdf + 117);
    const auto *pdf_118 = buffer.data(pdf + 118);
    const auto *pdf_119 = buffer.data(pdf + 119);
    const auto *pdf_120 = buffer.data(pdf + 120);
    const auto *pdf_122 = buffer.data(pdf + 122);
    const auto *pdf_126 = buffer.data(pdf + 126);
    const auto *pdf_127 = buffer.data(pdf + 127);
    const auto *pdf_128 = buffer.data(pdf + 128);
    const auto *pdf_129 = buffer.data(pdf + 129);
    const auto *pdf_130 = buffer.data(pdf + 130);
    const auto *pdf_132 = buffer.data(pdf + 132);
    const auto *pdf_136 = buffer.data(pdf + 136);
    const auto *pdf_137 = buffer.data(pdf + 137);
    const auto *pdf_138 = buffer.data(pdf + 138);
    const auto *pdf_139 = buffer.data(pdf + 139);
    const auto *pdf_140 = buffer.data(pdf + 140);
    const auto *pdf_142 = buffer.data(pdf + 142);
    const auto *pdf_143 = buffer.data(pdf + 143);
    const auto *pdf_146 = buffer.data(pdf + 146);
    const auto *pdf_147 = buffer.data(pdf + 147);
    const auto *pdf_148 = buffer.data(pdf + 148);
    const auto *pdf_149 = buffer.data(pdf + 149);
    const auto *pdf_150 = buffer.data(pdf + 150);
    const auto *pdf_152 = buffer.data(pdf + 152);
    const auto *pdf_155 = buffer.data(pdf + 155);
    const auto *pdf_156 = buffer.data(pdf + 156);
    const auto *pdf_157 = buffer.data(pdf + 157);
    const auto *pdf_158 = buffer.data(pdf + 158);
    const auto *pdf_159 = buffer.data(pdf + 159);
    const auto *pdf_160 = buffer.data(pdf + 160);
    const auto *pdf_162 = buffer.data(pdf + 162);
    const auto *pdf_163 = buffer.data(pdf + 163);
    const auto *pdf_166 = buffer.data(pdf + 166);
    const auto *pdf_167 = buffer.data(pdf + 167);
    const auto *pdf_168 = buffer.data(pdf + 168);
    const auto *pdf_169 = buffer.data(pdf + 169);
    const auto *pdf_170 = buffer.data(pdf + 170);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_x, pc_x, pc_y, pc_z, sdf_29, ppg0_85, \
                         ppg0_87, ppf_36, ppg1_85, ppg1_87, pdf_86, \
                         pdf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_x[k] * ppg0_85[k]
                   - f_7 * pc_x[k] * ppg1_85[k];

        t_131[k] = f_0 * ppf_36[k]
                   + f_4 * pc_z[k] * pdf_86[k];

        t_132[k] = pb_x[k] * ppg0_87[k]
                   - f_7 * pc_x[k] * ppg1_87[k];

        t_133[k] = f_0 * sdf_29[k]
                   + f_4 * pc_y[k] * pdf_89[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_y, pc_x, pc_y, pc_z, sdg0_44, sdg1_44, \
                         pdd0_54, pdd0_55, pdd1_54, pdd1_55, pdf_90, \
                         pdf_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_y[k] * sdg0_44[k]
                   - f_7 * pc_y[k] * sdg1_44[k];

        t_135[k] = f_2 * pdd0_54[k]
                   - f_3 * pdd1_54[k]
                   + f_4 * pc_x[k] * pdf_90[k];

        t_136[k] = f_10 * pdd0_55[k]
                   - f_11 * pdd1_55[k]
                   + f_4 * pc_x[k] * pdf_91[k];

        t_137[k] = f_4 * pc_z[k] * pdf_90[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, pc_x, pc_z, pdd0_57, pdd0_59, \
                         pdd1_57, pdd1_59, pdf_91, pdf_93, pdf_95, pdf_96, \
                         pdf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * pdd0_57[k]
                   - f_6 * pdd1_57[k]
                   + f_4 * pc_x[k] * pdf_93[k];

        t_139[k] = f_4 * pc_z[k] * pdf_91[k];

        t_140[k] = f_5 * pdd0_59[k]
                   - f_6 * pdd1_59[k]
                   + f_4 * pc_x[k] * pdf_95[k];

        t_141[k] = f_4 * pc_x[k] * pdf_96[k];

        t_142[k] = f_4 * pc_x[k] * pdf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, pc_x, pc_y, pc_z, sdf_36, ppf_46, \
                         pdd0_57, pdd1_57, pdf_96, pdf_97, pdf_98, \
                         pdf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * pc_x[k] * pdf_98[k];

        t_144[k] = f_4 * pc_x[k] * pdf_99[k];

        t_145[k] = f_0 * sdf_36[k]
                   + f_1 * ppf_46[k]
                   + f_2 * pdd0_57[k]
                   - f_3 * pdd1_57[k]
                   + f_4 * pc_y[k] * pdf_96[k];

        t_146[k] = f_4 * pc_z[k] * pdf_96[k];

        t_147[k] = f_5 * pdd0_57[k]
                   - f_6 * pdd1_57[k]
                   + f_4 * pc_z[k] * pdf_97[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pc_x, pc_y, pc_z, sdf_39, ppf_49, pdd0_59, \
                         pdd0_60, pdd1_59, pdd1_60, pdf_99, pdf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * sdf_39[k]
                   + f_1 * ppf_49[k]
                   + f_4 * pc_y[k] * pdf_99[k];

        t_149[k] = f_2 * pdd0_59[k]
                   - f_3 * pdd1_59[k]
                   + f_4 * pc_z[k] * pdf_99[k];

        t_150[k] = f_2 * pdd0_60[k]
                   - f_3 * pdd1_60[k]
                   + f_4 * pc_x[k] * pdf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pb_z, pc_z, ppg0_61, ppg0_63, ppf_40, \
                         ppf_41, ppg1_61, ppg1_63, pdf_100, pdf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = pb_z[k] * ppg0_61[k]
                   - f_7 * pc_z[k] * ppg1_61[k];

        t_152[k] = f_0 * ppf_40[k]
                   + f_4 * pc_z[k] * pdf_100[k];

        t_153[k] = pb_z[k] * ppg0_63[k]
                   - f_7 * pc_z[k] * ppg1_63[k];

        t_154[k] = f_0 * ppf_41[k]
                   + f_4 * pc_z[k] * pdf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pc_x, pdd0_65, pdd1_65, pdf_105, \
                         pdf_106, pdf_107, pdf_108, pdf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_5 * pdd0_65[k]
                   - f_6 * pdd1_65[k]
                   + f_4 * pc_x[k] * pdf_105[k];

        t_156[k] = f_4 * pc_x[k] * pdf_106[k];

        t_157[k] = f_4 * pc_x[k] * pdf_107[k];

        t_158[k] = f_4 * pc_x[k] * pdf_108[k];

        t_159[k] = f_4 * pc_x[k] * pdf_109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pb_z, pc_z, ppg0_70, ppf_46, ppf_47, ppg1_70, \
                         pdd0_63, pdd1_63, pdf_106, pdf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_z[k] * ppg0_70[k]
                   - f_7 * pc_z[k] * ppg1_70[k];

        t_161[k] = f_0 * ppf_46[k]
                   + f_4 * pc_z[k] * pdf_106[k];

        t_162[k] = f_0 * ppf_47[k]
                   + f_5 * pdd0_63[k]
                   - f_6 * pdd1_63[k]
                   + f_4 * pc_z[k] * pdf_107[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_y, pc_y, pc_z, sdg0_75, sdf_49, sdg1_75, \
                         ppf_49, ppf_59, pdd0_65, pdd1_65, pdf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_0 * sdf_49[k]
                   + f_0 * ppf_59[k]
                   + f_4 * pc_y[k] * pdf_109[k];

        t_164[k] = f_0 * ppf_49[k]
                   + f_2 * pdd0_65[k]
                   - f_3 * pdd1_65[k]
                   + f_4 * pc_z[k] * pdf_109[k];

        t_165[k] = pa_y[k] * sdg0_75[k]
                   - f_7 * pc_y[k] * sdg1_75[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pc_x, pc_z, ppf_50, ppf_51, pdd0_67, \
                         pdd0_69, pdd1_67, pdd1_69, pdf_110, pdf_111, \
                         pdf_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * pdd0_67[k]
                   - f_11 * pdd1_67[k]
                   + f_4 * pc_x[k] * pdf_111[k];

        t_167[k] = f_1 * ppf_50[k]
                   + f_4 * pc_z[k] * pdf_110[k];

        t_168[k] = f_5 * pdd0_69[k]
                   - f_6 * pdd1_69[k]
                   + f_4 * pc_x[k] * pdf_113[k];

        t_169[k] = f_1 * ppf_51[k]
                   + f_4 * pc_z[k] * pdf_111[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pa_y, pc_x, pc_y, sdg0_80, \
                         sdg1_80, pdf_116, pdf_117, pdf_118, pdf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pa_y[k] * sdg0_80[k]
                   - f_7 * pc_y[k] * sdg1_80[k];

        t_171[k] = f_4 * pc_x[k] * pdf_116[k];

        t_172[k] = f_4 * pc_x[k] * pdf_117[k];

        t_173[k] = f_4 * pc_x[k] * pdf_118[k];

        t_174[k] = f_4 * pc_x[k] * pdf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_y, pc_y, pc_z, sdg0_85, sdg0_87, sdf_56, \
                         sdf_58, sdg1_85, sdg1_87, ppf_56, pdf_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_y[k] * sdg0_85[k]
                   + f_8 * sdf_56[k]
                   - f_7 * pc_y[k] * sdg1_85[k];

        t_176[k] = f_1 * ppf_56[k]
                   + f_4 * pc_z[k] * pdf_116[k];

        t_177[k] = pa_y[k] * sdg0_87[k]
                   + f_1 * sdf_58[k]
                   - f_7 * pc_y[k] * sdg1_87[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_y, pa_z, pc_y, pc_z, sdg0_0, sdg0_89, \
                         sdf_59, sdg1_0, sdg1_89, pdf_119, pdf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_0 * sdf_59[k]
                   + f_4 * pc_y[k] * pdf_119[k];

        t_179[k] = pa_y[k] * sdg0_89[k]
                   - f_7 * pc_y[k] * sdg1_89[k];

        t_180[k] = pa_z[k] * sdg0_0[k]
                   - f_7 * pc_z[k] * sdg1_0[k];

        t_181[k] = f_4 * pc_y[k] * pdf_120[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pa_z, pc_y, pc_z, sdg0_2, sdg0_3, sdg0_5, \
                         sdf_0, sdf_2, sdg1_2, sdg1_3, sdg1_5, \
                         pdf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = pa_z[k] * sdg0_2[k]
                   + f_0 * sdf_0[k]
                   - f_7 * pc_z[k] * sdg1_2[k];

        t_183[k] = pa_z[k] * sdg0_3[k]
                   - f_7 * pc_z[k] * sdg1_3[k];

        t_184[k] = f_4 * pc_y[k] * pdf_122[k];

        t_185[k] = pa_z[k] * sdg0_5[k]
                   + f_1 * sdf_2[k]
                   - f_7 * pc_z[k] * sdg1_5[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, ppf_66, ppf_67, ppf_68, ppf_69, \
                         pdf_126, pdf_127, pdf_128, pdf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_1 * ppf_66[k]
                   + f_4 * pc_x[k] * pdf_126[k];

        t_187[k] = f_1 * ppf_67[k]
                   + f_4 * pc_x[k] * pdf_127[k];

        t_188[k] = f_1 * ppf_68[k]
                   + f_4 * pc_x[k] * pdf_128[k];

        t_189[k] = f_1 * ppf_69[k]
                   + f_4 * pc_x[k] * pdf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_z, pc_y, pc_z, sdg0_10, sdg1_10, \
                         pdd0_76, pdd0_77, pdd1_76, pdd1_77, pdf_127, pdf_128, \
                         pdf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_z[k] * sdg0_10[k]
                   - f_7 * pc_z[k] * sdg1_10[k];

        t_191[k] = f_10 * pdd0_76[k]
                   - f_11 * pdd1_76[k]
                   + f_4 * pc_y[k] * pdf_127[k];

        t_192[k] = f_5 * pdd0_77[k]
                   - f_6 * pdd1_77[k]
                   + f_4 * pc_y[k] * pdf_128[k];

        t_193[k] = f_4 * pc_y[k] * pdf_129[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_z, pc_y, pc_z, sdg0_14, sdg0_15, sdf_9, \
                         sdg1_14, sdg1_15, ppf_60, pdf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_z[k] * sdg0_14[k]
                   + f_8 * sdf_9[k]
                   - f_7 * pc_z[k] * sdg1_14[k];

        t_195[k] = pa_z[k] * sdg0_15[k]
                   - f_7 * pc_z[k] * sdg1_15[k];

        t_196[k] = f_0 * ppf_60[k]
                   + f_4 * pc_y[k] * pdf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pb_y, pc_y, pc_z, sdg0_18, sdg1_18, \
                         ppg0_92, ppg0_95, ppf_62, ppg1_92, ppg1_95, \
                         pdf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pb_y[k] * ppg0_92[k]
                   - f_7 * pc_y[k] * ppg1_92[k];

        t_198[k] = pa_z[k] * sdg0_18[k]
                   - f_7 * pc_z[k] * sdg1_18[k];

        t_199[k] = f_0 * ppf_62[k]
                   + f_4 * pc_y[k] * pdf_132[k];

        t_200[k] = pb_y[k] * ppg0_95[k]
                   - f_7 * pc_y[k] * ppg1_95[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, ppf_76, ppf_77, ppf_78, ppf_79, \
                         pdf_136, pdf_137, pdf_138, pdf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_0 * ppf_76[k]
                   + f_4 * pc_x[k] * pdf_136[k];

        t_202[k] = f_0 * ppf_77[k]
                   + f_4 * pc_x[k] * pdf_137[k];

        t_203[k] = f_0 * ppf_78[k]
                   + f_4 * pc_x[k] * pdf_138[k];

        t_204[k] = f_0 * ppf_79[k]
                   + f_4 * pc_x[k] * pdf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pa_z, pb_x, pc_x, pc_z, sdg0_25, sdg1_25, \
                         ppg0_116, ppg0_117, ppg1_116, ppg1_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = pa_z[k] * sdg0_25[k]
                   - f_7 * pc_z[k] * sdg1_25[k];

        t_206[k] = pb_x[k] * ppg0_116[k]
                   - f_7 * pc_x[k] * ppg1_116[k];

        t_207[k] = pb_x[k] * ppg0_117[k]
                   - f_7 * pc_x[k] * ppg1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_x, pc_x, pc_y, ppg0_119, ppf_69, \
                         ppf_80, ppg1_119, pdd0_84, pdd1_84, pdf_139, \
                         pdf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_0 * ppf_69[k]
                   + f_4 * pc_y[k] * pdf_139[k];

        t_209[k] = pb_x[k] * ppg0_119[k]
                   - f_7 * pc_x[k] * ppg1_119[k];

        t_210[k] = f_0 * ppf_80[k]
                   + f_2 * pdd0_84[k]
                   - f_3 * pdd1_84[k]
                   + f_4 * pc_x[k] * pdf_140[k];

        t_211[k] = f_4 * pc_y[k] * pdf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pb_x, pc_x, pc_y, ppg0_122, ppf_82, ppf_83, \
                         ppg1_122, pdd0_87, pdd1_87, pdf_142, pdf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_x[k] * ppg0_122[k]
                   + f_9 * ppf_82[k]
                   - f_7 * pc_x[k] * ppg1_122[k];

        t_213[k] = f_0 * ppf_83[k]
                   + f_5 * pdd0_87[k]
                   - f_6 * pdd1_87[k]
                   + f_4 * pc_x[k] * pdf_143[k];

        t_214[k] = f_4 * pc_y[k] * pdf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pb_x, pc_x, ppg0_125, ppf_85, ppf_86, \
                         ppf_87, ppf_88, ppg1_125, pdf_146, pdf_147, \
                         pdf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pb_x[k] * ppg0_125[k]
                   + f_1 * ppf_85[k]
                   - f_7 * pc_x[k] * ppg1_125[k];

        t_216[k] = f_0 * ppf_86[k]
                   + f_4 * pc_x[k] * pdf_146[k];

        t_217[k] = f_0 * ppf_87[k]
                   + f_4 * pc_x[k] * pdf_147[k];

        t_218[k] = f_0 * ppf_88[k]
                   + f_4 * pc_x[k] * pdf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pb_x, pc_x, pc_y, ppg0_130, \
                         ppg0_131, ppg0_132, ppf_89, ppg1_130, ppg1_131, ppg1_132, \
                         pdf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_0 * ppf_89[k]
                   + f_4 * pc_x[k] * pdf_149[k];

        t_220[k] = pb_x[k] * ppg0_130[k]
                   - f_7 * pc_x[k] * ppg1_130[k];

        t_221[k] = pb_x[k] * ppg0_131[k]
                   - f_7 * pc_x[k] * ppg1_131[k];

        t_222[k] = pb_x[k] * ppg0_132[k]
                   - f_7 * pc_x[k] * ppg1_132[k];

        t_223[k] = f_4 * pc_y[k] * pdf_149[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, pa_z, pb_x, pc_x, pc_y, pc_z, sdg0_45, sdg1_45, \
                         ppg0_134, ppf_70, ppg1_134, pdf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pb_x[k] * ppg0_134[k]
                   - f_7 * pc_x[k] * ppg1_134[k];

        t_225[k] = pa_z[k] * sdg0_45[k]
                   - f_7 * pc_z[k] * sdg1_45[k];

        t_226[k] = f_1 * ppf_70[k]
                   + f_4 * pc_y[k] * pdf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pa_z, pc_x, pc_y, pc_z, sdg0_48, sdg1_48, \
                         ppf_72, pdd0_92, pdd1_92, pdf_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_10 * pdd0_92[k]
                   - f_11 * pdd1_92[k]
                   + f_4 * pc_x[k] * pdf_152[k];

        t_228[k] = pa_z[k] * sdg0_48[k]
                   - f_7 * pc_z[k] * sdg1_48[k];

        t_229[k] = f_1 * ppf_72[k]
                   + f_4 * pc_y[k] * pdf_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, pc_x, pdd0_95, pdd1_95, pdf_155, \
                         pdf_156, pdf_157, pdf_158, pdf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_5 * pdd0_95[k]
                   - f_6 * pdd1_95[k]
                   + f_4 * pc_x[k] * pdf_155[k];

        t_231[k] = f_4 * pc_x[k] * pdf_156[k];

        t_232[k] = f_4 * pc_x[k] * pdf_157[k];

        t_233[k] = f_4 * pc_x[k] * pdf_158[k];

        t_234[k] = f_4 * pc_x[k] * pdf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_z, pc_z, sdg0_55, sdg0_56, sdg0_57, sdf_36, \
                         sdf_37, sdg1_55, sdg1_56, sdg1_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = pa_z[k] * sdg0_55[k]
                   - f_7 * pc_z[k] * sdg1_55[k];

        t_236[k] = pa_z[k] * sdg0_56[k]
                   + f_0 * sdf_36[k]
                   - f_7 * pc_z[k] * sdg1_56[k];

        t_237[k] = pa_z[k] * sdg0_57[k]
                   + f_1 * sdf_37[k]
                   - f_7 * pc_z[k] * sdg1_57[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_z, pb_y, pc_y, pc_z, sdg0_59, sdf_39, \
                         sdg1_59, ppg0_120, ppf_79, ppg1_120, pdf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_1 * ppf_79[k]
                   + f_4 * pc_y[k] * pdf_159[k];

        t_239[k] = pa_z[k] * sdg0_59[k]
                   + f_8 * sdf_39[k]
                   - f_7 * pc_z[k] * sdg1_59[k];

        t_240[k] = pb_y[k] * ppg0_120[k]
                   - f_7 * pc_y[k] * ppg1_120[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_y, pc_x, pc_y, ppg0_122, ppf_80, \
                         ppf_82, ppg1_122, pdd0_99, pdd1_99, pdf_160, pdf_162, \
                         pdf_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_0 * ppf_80[k]
                   + f_4 * pc_y[k] * pdf_160[k];

        t_242[k] = pb_y[k] * ppg0_122[k]
                   - f_7 * pc_y[k] * ppg1_122[k];

        t_243[k] = f_5 * pdd0_99[k]
                   - f_6 * pdd1_99[k]
                   + f_4 * pc_x[k] * pdf_163[k];

        t_244[k] = f_0 * ppf_82[k]
                   + f_4 * pc_y[k] * pdf_162[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, pb_y, pc_x, pc_y, ppg0_125, \
                         ppg1_125, pdf_166, pdf_167, pdf_168, pdf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pb_y[k] * ppg0_125[k]
                   - f_7 * pc_y[k] * ppg1_125[k];

        t_246[k] = f_4 * pc_x[k] * pdf_166[k];

        t_247[k] = f_4 * pc_x[k] * pdf_167[k];

        t_248[k] = f_4 * pc_x[k] * pdf_168[k];

        t_249[k] = f_4 * pc_x[k] * pdf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pb_y, pc_y, ppg0_131, ppg0_132, ppf_86, ppf_87, \
                         ppf_88, ppg1_131, ppg1_132, pdd0_99, pdd1_99, \
                         pdf_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * ppf_86[k]
                   + f_2 * pdd0_99[k]
                   - f_3 * pdd1_99[k]
                   + f_4 * pc_y[k] * pdf_166[k];

        t_251[k] = pb_y[k] * ppg0_131[k]
                   + f_9 * ppf_87[k]
                   - f_7 * pc_y[k] * ppg1_131[k];

        t_252[k] = pb_y[k] * ppg0_132[k]
                   + f_1 * ppf_88[k]
                   - f_7 * pc_y[k] * ppg1_132[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pb_y, pc_x, pc_y, ppg0_134, ppf_89, \
                         ppg1_134, pdd0_102, pdd1_102, pdf_169, \
                         pdf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_0 * ppf_89[k]
                   + f_4 * pc_y[k] * pdf_169[k];

        t_254[k] = pb_y[k] * ppg0_134[k]
                   - f_7 * pc_y[k] * ppg1_134[k];

        t_255[k] = f_2 * pdd0_102[k]
                   - f_3 * pdd1_102[k]
                   + f_4 * pc_x[k] * pdf_170[k];

        t_256[k] = f_4 * pc_y[k] * pdf_170[k];
    }
}

static auto
compute_prim_pdg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sdf, const size_t ppf,
                                                          const size_t pdd0, const size_t pdd1,
                                                          const size_t pdf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);

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

    const auto *sdf_59 = buffer.data(sdf + 59);

    const auto *ppf_89 = buffer.data(ppf + 89);

    const auto *pdd0_104 = buffer.data(pdd0 + 104);
    const auto *pdd0_105 = buffer.data(pdd0 + 105);
    const auto *pdd0_106 = buffer.data(pdd0 + 106);
    const auto *pdd0_107 = buffer.data(pdd0 + 107);

    const auto *pdd1_104 = buffer.data(pdd1 + 104);
    const auto *pdd1_105 = buffer.data(pdd1 + 105);
    const auto *pdd1_106 = buffer.data(pdd1 + 106);
    const auto *pdd1_107 = buffer.data(pdd1 + 107);

    const auto *pdf_172 = buffer.data(pdf + 172);
    const auto *pdf_173 = buffer.data(pdf + 173);
    const auto *pdf_175 = buffer.data(pdf + 175);
    const auto *pdf_176 = buffer.data(pdf + 176);
    const auto *pdf_177 = buffer.data(pdf + 177);
    const auto *pdf_178 = buffer.data(pdf + 178);
    const auto *pdf_179 = buffer.data(pdf + 179);

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pc_x, pc_y, pdd0_104, pdd0_105, pdd0_107, \
                         pdd1_104, pdd1_105, pdd1_107, pdf_172, pdf_173, \
                         pdf_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_10 * pdd0_104[k]
                   - f_11 * pdd1_104[k]
                   + f_4 * pc_x[k] * pdf_172[k];

        t_258[k] = f_5 * pdd0_105[k]
                   - f_6 * pdd1_105[k]
                   + f_4 * pc_x[k] * pdf_173[k];

        t_259[k] = f_4 * pc_y[k] * pdf_172[k];

        t_260[k] = f_5 * pdd0_107[k]
                   - f_6 * pdd1_107[k]
                   + f_4 * pc_x[k] * pdf_175[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, t_266, pc_x, pc_y, pdd0_105, \
                         pdd0_106, pdd1_105, pdd1_106, pdf_176, pdf_177, pdf_178, \
                         pdf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_4 * pc_x[k] * pdf_176[k];

        t_262[k] = f_4 * pc_x[k] * pdf_177[k];

        t_263[k] = f_4 * pc_x[k] * pdf_178[k];

        t_264[k] = f_4 * pc_x[k] * pdf_179[k];

        t_265[k] = f_2 * pdd0_105[k]
                   - f_3 * pdd1_105[k]
                   + f_4 * pc_y[k] * pdf_176[k];

        t_266[k] = f_10 * pdd0_106[k]
                   - f_11 * pdd1_106[k]
                   + f_4 * pc_y[k] * pdf_177[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, sdf_59, ppf_89, pdd0_107, pdd1_107, \
                         pdf_178, pdf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_5 * pdd0_107[k]
                   - f_6 * pdd1_107[k]
                   + f_4 * pc_y[k] * pdf_178[k];

        t_268[k] = f_4 * pc_y[k] * pdf_179[k];

        t_269[k] = f_0 * sdf_59[k]
                   + f_1 * ppf_89[k]
                   + f_2 * pdd0_107[k]
                   - f_3 * pdd1_107[k]
                   + f_4 * pc_z[k] * pdf_179[k];
    }
}

auto
compute_prim_pdg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sdg0,
                                                   const size_t sdf, const size_t sdg1,
                                                   const size_t ppg0, const size_t ppf,
                                                   const size_t ppg1, const size_t pdd0,
                                                   const size_t pdd1, const size_t pdf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pdg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sdg0,
                                                              sdf, sdg1, ppg0, ppf, ppg1, pdd0,
                                                              pdd1, pdf, ncols, gamma, p, q);

    compute_prim_pdg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sdg0,
                                                              sdf, sdg1, ppg0, ppf, ppg1, pdd0,
                                                              pdd1, pdf, ncols, gamma, p, q);

    compute_prim_pdg_three_center_electron_repulsion_0_piece2(buffer, target, pc, sdf, ppf,
                                                              pdd0, pdd1, pdf, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
