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


#include "SimdThreeCenterElectronRepulsionVrrRecPFI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfh,
                                                          const size_t pdi0, const size_t pdh,
                                                          const size_t pdi1, const size_t pfg0,
                                                          const size_t pfg1, const size_t pfh,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *sfh_0 = buffer.data(sfh + 0);
    const auto *sfh_15 = buffer.data(sfh + 15);
    const auto *sfh_17 = buffer.data(sfh + 17);
    const auto *sfh_18 = buffer.data(sfh + 18);
    const auto *sfh_20 = buffer.data(sfh + 20);
    const auto *sfh_36 = buffer.data(sfh + 36);
    const auto *sfh_38 = buffer.data(sfh + 38);
    const auto *sfh_39 = buffer.data(sfh + 39);
    const auto *sfh_59 = buffer.data(sfh + 59);
    const auto *sfh_60 = buffer.data(sfh + 60);
    const auto *sfh_62 = buffer.data(sfh + 62);
    const auto *sfh_63 = buffer.data(sfh + 63);
    const auto *sfh_78 = buffer.data(sfh + 78);
    const auto *sfh_80 = buffer.data(sfh + 80);
    const auto *sfh_81 = buffer.data(sfh + 81);
    const auto *sfh_83 = buffer.data(sfh + 83);

    const auto *pdi0_0 = buffer.data(pdi0 + 0);
    const auto *pdi0_3 = buffer.data(pdi0 + 3);
    const auto *pdi0_5 = buffer.data(pdi0 + 5);
    const auto *pdi0_6 = buffer.data(pdi0 + 6);
    const auto *pdi0_9 = buffer.data(pdi0 + 9);
    const auto *pdi0_10 = buffer.data(pdi0 + 10);
    const auto *pdi0_12 = buffer.data(pdi0 + 12);
    const auto *pdi0_14 = buffer.data(pdi0 + 14);
    const auto *pdi0_15 = buffer.data(pdi0 + 15);
    const auto *pdi0_20 = buffer.data(pdi0 + 20);
    const auto *pdi0_21 = buffer.data(pdi0 + 21);
    const auto *pdi0_27 = buffer.data(pdi0 + 27);
    const auto *pdi0_31 = buffer.data(pdi0 + 31);
    const auto *pdi0_34 = buffer.data(pdi0 + 34);
    const auto *pdi0_38 = buffer.data(pdi0 + 38);
    const auto *pdi0_56 = buffer.data(pdi0 + 56);
    const auto *pdi0_61 = buffer.data(pdi0 + 61);
    const auto *pdi0_65 = buffer.data(pdi0 + 65);

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
    const auto *pdh_35 = buffer.data(pdh + 35);
    const auto *pdh_36 = buffer.data(pdh + 36);
    const auto *pdh_38 = buffer.data(pdh + 38);
    const auto *pdh_39 = buffer.data(pdh + 39);
    const auto *pdh_40 = buffer.data(pdh + 40);
    const auto *pdh_41 = buffer.data(pdh + 41);
    const auto *pdh_42 = buffer.data(pdh + 42);
    const auto *pdh_44 = buffer.data(pdh + 44);
    const auto *pdh_47 = buffer.data(pdh + 47);
    const auto *pdh_50 = buffer.data(pdh + 50);
    const auto *pdh_59 = buffer.data(pdh + 59);
    const auto *pdh_60 = buffer.data(pdh + 60);
    const auto *pdh_62 = buffer.data(pdh + 62);
    const auto *pdh_63 = buffer.data(pdh + 63);
    const auto *pdh_78 = buffer.data(pdh + 78);
    const auto *pdh_80 = buffer.data(pdh + 80);
    const auto *pdh_81 = buffer.data(pdh + 81);
    const auto *pdh_83 = buffer.data(pdh + 83);

    const auto *pdi1_0 = buffer.data(pdi1 + 0);
    const auto *pdi1_3 = buffer.data(pdi1 + 3);
    const auto *pdi1_5 = buffer.data(pdi1 + 5);
    const auto *pdi1_6 = buffer.data(pdi1 + 6);
    const auto *pdi1_9 = buffer.data(pdi1 + 9);
    const auto *pdi1_10 = buffer.data(pdi1 + 10);
    const auto *pdi1_12 = buffer.data(pdi1 + 12);
    const auto *pdi1_14 = buffer.data(pdi1 + 14);
    const auto *pdi1_15 = buffer.data(pdi1 + 15);
    const auto *pdi1_20 = buffer.data(pdi1 + 20);
    const auto *pdi1_21 = buffer.data(pdi1 + 21);
    const auto *pdi1_27 = buffer.data(pdi1 + 27);
    const auto *pdi1_31 = buffer.data(pdi1 + 31);
    const auto *pdi1_34 = buffer.data(pdi1 + 34);
    const auto *pdi1_38 = buffer.data(pdi1 + 38);
    const auto *pdi1_56 = buffer.data(pdi1 + 56);
    const auto *pdi1_61 = buffer.data(pdi1 + 61);
    const auto *pdi1_65 = buffer.data(pdi1 + 65);

    const auto *pfg0_0 = buffer.data(pfg0 + 0);
    const auto *pfg0_1 = buffer.data(pfg0 + 1);
    const auto *pfg0_2 = buffer.data(pfg0 + 2);
    const auto *pfg0_3 = buffer.data(pfg0 + 3);
    const auto *pfg0_5 = buffer.data(pfg0 + 5);
    const auto *pfg0_10 = buffer.data(pfg0 + 10);
    const auto *pfg0_12 = buffer.data(pfg0 + 12);
    const auto *pfg0_13 = buffer.data(pfg0 + 13);
    const auto *pfg0_14 = buffer.data(pfg0 + 14);
    const auto *pfg0_25 = buffer.data(pfg0 + 25);
    const auto *pfg0_27 = buffer.data(pfg0 + 27);
    const auto *pfg0_28 = buffer.data(pfg0 + 28);
    const auto *pfg0_29 = buffer.data(pfg0 + 29);
    const auto *pfg0_35 = buffer.data(pfg0 + 35);
    const auto *pfg0_42 = buffer.data(pfg0 + 42);
    const auto *pfg0_43 = buffer.data(pfg0 + 43);
    const auto *pfg0_44 = buffer.data(pfg0 + 44);
    const auto *pfg0_45 = buffer.data(pfg0 + 45);
    const auto *pfg0_46 = buffer.data(pfg0 + 46);
    const auto *pfg0_47 = buffer.data(pfg0 + 47);
    const auto *pfg0_48 = buffer.data(pfg0 + 48);
    const auto *pfg0_50 = buffer.data(pfg0 + 50);
    const auto *pfg0_55 = buffer.data(pfg0 + 55);
    const auto *pfg0_57 = buffer.data(pfg0 + 57);
    const auto *pfg0_58 = buffer.data(pfg0 + 58);
    const auto *pfg0_59 = buffer.data(pfg0 + 59);
    const auto *pfg0_65 = buffer.data(pfg0 + 65);

    const auto *pfg1_0 = buffer.data(pfg1 + 0);
    const auto *pfg1_1 = buffer.data(pfg1 + 1);
    const auto *pfg1_2 = buffer.data(pfg1 + 2);
    const auto *pfg1_3 = buffer.data(pfg1 + 3);
    const auto *pfg1_5 = buffer.data(pfg1 + 5);
    const auto *pfg1_10 = buffer.data(pfg1 + 10);
    const auto *pfg1_12 = buffer.data(pfg1 + 12);
    const auto *pfg1_13 = buffer.data(pfg1 + 13);
    const auto *pfg1_14 = buffer.data(pfg1 + 14);
    const auto *pfg1_25 = buffer.data(pfg1 + 25);
    const auto *pfg1_27 = buffer.data(pfg1 + 27);
    const auto *pfg1_28 = buffer.data(pfg1 + 28);
    const auto *pfg1_29 = buffer.data(pfg1 + 29);
    const auto *pfg1_35 = buffer.data(pfg1 + 35);
    const auto *pfg1_42 = buffer.data(pfg1 + 42);
    const auto *pfg1_43 = buffer.data(pfg1 + 43);
    const auto *pfg1_44 = buffer.data(pfg1 + 44);
    const auto *pfg1_45 = buffer.data(pfg1 + 45);
    const auto *pfg1_46 = buffer.data(pfg1 + 46);
    const auto *pfg1_47 = buffer.data(pfg1 + 47);
    const auto *pfg1_48 = buffer.data(pfg1 + 48);
    const auto *pfg1_50 = buffer.data(pfg1 + 50);
    const auto *pfg1_55 = buffer.data(pfg1 + 55);
    const auto *pfg1_57 = buffer.data(pfg1 + 57);
    const auto *pfg1_58 = buffer.data(pfg1 + 58);
    const auto *pfg1_59 = buffer.data(pfg1 + 59);
    const auto *pfg1_65 = buffer.data(pfg1 + 65);

    const auto *pfh_0 = buffer.data(pfh + 0);
    const auto *pfh_1 = buffer.data(pfh + 1);
    const auto *pfh_2 = buffer.data(pfh + 2);
    const auto *pfh_3 = buffer.data(pfh + 3);
    const auto *pfh_5 = buffer.data(pfh + 5);
    const auto *pfh_6 = buffer.data(pfh + 6);
    const auto *pfh_8 = buffer.data(pfh + 8);
    const auto *pfh_9 = buffer.data(pfh + 9);
    const auto *pfh_10 = buffer.data(pfh + 10);
    const auto *pfh_14 = buffer.data(pfh + 14);
    const auto *pfh_15 = buffer.data(pfh + 15);
    const auto *pfh_17 = buffer.data(pfh + 17);
    const auto *pfh_18 = buffer.data(pfh + 18);
    const auto *pfh_19 = buffer.data(pfh + 19);
    const auto *pfh_20 = buffer.data(pfh + 20);
    const auto *pfh_21 = buffer.data(pfh + 21);
    const auto *pfh_23 = buffer.data(pfh + 23);
    const auto *pfh_24 = buffer.data(pfh + 24);
    const auto *pfh_26 = buffer.data(pfh + 26);
    const auto *pfh_27 = buffer.data(pfh + 27);
    const auto *pfh_30 = buffer.data(pfh + 30);
    const auto *pfh_31 = buffer.data(pfh + 31);
    const auto *pfh_35 = buffer.data(pfh + 35);
    const auto *pfh_36 = buffer.data(pfh + 36);
    const auto *pfh_38 = buffer.data(pfh + 38);
    const auto *pfh_39 = buffer.data(pfh + 39);
    const auto *pfh_40 = buffer.data(pfh + 40);
    const auto *pfh_41 = buffer.data(pfh + 41);
    const auto *pfh_42 = buffer.data(pfh + 42);
    const auto *pfh_44 = buffer.data(pfh + 44);
    const auto *pfh_45 = buffer.data(pfh + 45);
    const auto *pfh_47 = buffer.data(pfh + 47);
    const auto *pfh_48 = buffer.data(pfh + 48);
    const auto *pfh_50 = buffer.data(pfh + 50);
    const auto *pfh_51 = buffer.data(pfh + 51);
    const auto *pfh_52 = buffer.data(pfh + 52);
    const auto *pfh_56 = buffer.data(pfh + 56);
    const auto *pfh_57 = buffer.data(pfh + 57);
    const auto *pfh_59 = buffer.data(pfh + 59);
    const auto *pfh_60 = buffer.data(pfh + 60);
    const auto *pfh_61 = buffer.data(pfh + 61);
    const auto *pfh_62 = buffer.data(pfh + 62);
    const auto *pfh_63 = buffer.data(pfh + 63);
    const auto *pfh_64 = buffer.data(pfh + 64);
    const auto *pfh_65 = buffer.data(pfh + 65);
    const auto *pfh_66 = buffer.data(pfh + 66);
    const auto *pfh_68 = buffer.data(pfh + 68);
    const auto *pfh_69 = buffer.data(pfh + 69);
    const auto *pfh_71 = buffer.data(pfh + 71);
    const auto *pfh_72 = buffer.data(pfh + 72);
    const auto *pfh_73 = buffer.data(pfh + 73);
    const auto *pfh_77 = buffer.data(pfh + 77);
    const auto *pfh_78 = buffer.data(pfh + 78);
    const auto *pfh_80 = buffer.data(pfh + 80);
    const auto *pfh_81 = buffer.data(pfh + 81);
    const auto *pfh_82 = buffer.data(pfh + 82);
    const auto *pfh_83 = buffer.data(pfh + 83);
    const auto *pfh_84 = buffer.data(pfh + 84);
    const auto *pfh_86 = buffer.data(pfh + 86);
    const auto *pfh_87 = buffer.data(pfh + 87);
    const auto *pfh_89 = buffer.data(pfh + 89);
    const auto *pfh_90 = buffer.data(pfh + 90);
    const auto *pfh_92 = buffer.data(pfh + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sfh_0, pdh_0, pfg0_0, \
                         pfg1_0, pfh_0, pfh_1, pfh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfh_0[k]
                 + f_1 * pdh_0[k]
                 + f_2 * pfg0_0[k]
                 - f_3 * pfg1_0[k]
                 + f_4 * pc_x[k] * pfh_0[k];

        t_1[k] = f_4 * pc_y[k] * pfh_0[k];

        t_2[k] = f_4 * pc_z[k] * pfh_0[k];

        t_3[k] = f_5 * pfg0_0[k]
                 - f_6 * pfg1_0[k]
                 + f_4 * pc_y[k] * pfh_1[k];

        t_4[k] = f_4 * pc_y[k] * pfh_2[k];

        t_5[k] = f_5 * pfg0_0[k]
                 - f_6 * pfg1_0[k]
                 + f_4 * pc_z[k] * pfh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, pfg0_1, pfg0_2, pfg0_3, pfg1_1, \
                         pfg1_2, pfg1_3, pfh_3, pfh_5, pfh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * pfg0_1[k]
                 - f_8 * pfg1_1[k]
                 + f_4 * pc_y[k] * pfh_3[k];

        t_7[k] = f_4 * pc_z[k] * pfh_3[k];

        t_8[k] = f_4 * pc_y[k] * pfh_5[k];

        t_9[k] = f_7 * pfg0_2[k]
                 - f_8 * pfg1_2[k]
                 + f_4 * pc_z[k] * pfh_5[k];

        t_10[k] = f_9 * pfg0_3[k]
                  - f_10 * pfg1_3[k]
                  + f_4 * pc_y[k] * pfh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, sfh_15, pdh_15, \
                         pfg0_5, pfg1_5, pfh_6, pfh_8, pfh_9, pfh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * pc_z[k] * pfh_6[k];

        t_12[k] = f_5 * pfg0_5[k]
                  - f_6 * pfg1_5[k]
                  + f_4 * pc_y[k] * pfh_8[k];

        t_13[k] = f_4 * pc_y[k] * pfh_9[k];

        t_14[k] = f_9 * pfg0_5[k]
                  - f_10 * pfg1_5[k]
                  + f_4 * pc_z[k] * pfh_9[k];

        t_15[k] = f_0 * sfh_15[k]
                  + f_1 * pdh_15[k]
                  + f_4 * pc_x[k] * pfh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, sfh_17, sfh_18, pdh_17, \
                         pdh_18, pfh_10, pfh_14, pfh_17, pfh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * pc_z[k] * pfh_10[k];

        t_17[k] = f_0 * sfh_17[k]
                  + f_1 * pdh_17[k]
                  + f_4 * pc_x[k] * pfh_17[k];

        t_18[k] = f_0 * sfh_18[k]
                  + f_1 * pdh_18[k]
                  + f_4 * pc_x[k] * pfh_18[k];

        t_19[k] = f_4 * pc_y[k] * pfh_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sfh_20, pdh_20, pfg0_10, \
                         pfg0_12, pfg1_10, pfg1_12, pfh_15, pfh_17, \
                         pfh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * sfh_20[k]
                  + f_1 * pdh_20[k]
                  + f_4 * pc_x[k] * pfh_20[k];

        t_21[k] = f_2 * pfg0_10[k]
                  - f_3 * pfg1_10[k]
                  + f_4 * pc_y[k] * pfh_15[k];

        t_22[k] = f_4 * pc_z[k] * pfh_15[k];

        t_23[k] = f_9 * pfg0_12[k]
                  - f_10 * pfg1_12[k]
                  + f_4 * pc_y[k] * pfh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, pfg0_13, pfg0_14, pfg1_13, \
                         pfg1_14, pfh_18, pfh_19, pfh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * pfg0_13[k]
                  - f_8 * pfg1_13[k]
                  + f_4 * pc_y[k] * pfh_18[k];

        t_25[k] = f_5 * pfg0_14[k]
                  - f_6 * pfg1_14[k]
                  + f_4 * pc_y[k] * pfh_19[k];

        t_26[k] = f_4 * pc_y[k] * pfh_20[k];

        t_27[k] = f_2 * pfg0_14[k]
                  - f_3 * pfg1_14[k]
                  + f_4 * pc_z[k] * pfh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_y, pc_z, pdi0_0, pdi0_3, pdh_0, \
                         pdh_1, pdi1_0, pdi1_3, pfh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * pdi0_0[k]
                  - f_11 * pc_y[k] * pdi1_0[k];

        t_29[k] = f_0 * pdh_0[k]
                  + f_4 * pc_y[k] * pfh_21[k];

        t_30[k] = f_4 * pc_z[k] * pfh_21[k];

        t_31[k] = pb_y[k] * pdi0_3[k]
                  + f_12 * pdh_1[k]
                  - f_11 * pc_y[k] * pdi1_3[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pc_y, pc_z, pdi0_5, pdi0_6, pdh_2, \
                         pdh_3, pdi1_5, pdi1_6, pfh_23, pfh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * pdh_2[k]
                  + f_4 * pc_y[k] * pfh_23[k];

        t_33[k] = pb_y[k] * pdi0_5[k]
                  - f_11 * pc_y[k] * pdi1_5[k];

        t_34[k] = pb_y[k] * pdi0_6[k]
                  + f_1 * pdh_3[k]
                  - f_11 * pc_y[k] * pdi1_6[k];

        t_35[k] = f_4 * pc_z[k] * pfh_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, pdi0_9, pdi0_10, pdh_5, \
                         pdh_6, pdi1_9, pdi1_10, pfh_26, pfh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * pdh_5[k]
                  + f_4 * pc_y[k] * pfh_26[k];

        t_37[k] = pb_y[k] * pdi0_9[k]
                  - f_11 * pc_y[k] * pdi1_9[k];

        t_38[k] = pb_y[k] * pdi0_10[k]
                  + f_13 * pdh_6[k]
                  - f_11 * pc_y[k] * pdi1_10[k];

        t_39[k] = f_4 * pc_z[k] * pfh_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_y, pc_y, pdi0_12, pdi0_14, pdh_8, pdh_9, \
                         pdi1_12, pdi1_14, pfh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * pdi0_12[k]
                  + f_12 * pdh_8[k]
                  - f_11 * pc_y[k] * pdi1_12[k];

        t_41[k] = f_0 * pdh_9[k]
                  + f_4 * pc_y[k] * pfh_30[k];

        t_42[k] = pb_y[k] * pdi0_14[k]
                  - f_11 * pc_y[k] * pdi1_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pc_x, pc_z, sfh_36, sfh_38, sfh_39, pdh_36, \
                         pdh_38, pdh_39, pfh_31, pfh_36, pfh_38, \
                         pfh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * sfh_36[k]
                  + f_12 * pdh_36[k]
                  + f_4 * pc_x[k] * pfh_36[k];

        t_44[k] = f_4 * pc_z[k] * pfh_31[k];

        t_45[k] = f_0 * sfh_38[k]
                  + f_12 * pdh_38[k]
                  + f_4 * pc_x[k] * pfh_38[k];

        t_46[k] = f_0 * sfh_39[k]
                  + f_12 * pdh_39[k]
                  + f_4 * pc_x[k] * pfh_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pc_y, pc_z, pdi0_20, pdh_14, pdh_15, \
                         pdi1_20, pfg0_25, pfg1_25, pfh_35, pfh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_0 * pdh_14[k]
                  + f_4 * pc_y[k] * pfh_35[k];

        t_48[k] = pb_y[k] * pdi0_20[k]
                  - f_11 * pc_y[k] * pdi1_20[k];

        t_49[k] = f_0 * pdh_15[k]
                  + f_2 * pfg0_25[k]
                  - f_3 * pfg1_25[k]
                  + f_4 * pc_y[k] * pfh_36[k];

        t_50[k] = f_4 * pc_z[k] * pfh_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pc_y, pdh_17, pdh_18, pdh_19, pfg0_27, pfg0_28, \
                         pfg0_29, pfg1_27, pfg1_28, pfg1_29, pfh_38, pfh_39, \
                         pfh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * pdh_17[k]
                  + f_9 * pfg0_27[k]
                  - f_10 * pfg1_27[k]
                  + f_4 * pc_y[k] * pfh_38[k];

        t_52[k] = f_0 * pdh_18[k]
                  + f_7 * pfg0_28[k]
                  - f_8 * pfg1_28[k]
                  + f_4 * pc_y[k] * pfh_39[k];

        t_53[k] = f_0 * pdh_19[k]
                  + f_5 * pfg0_29[k]
                  - f_6 * pfg1_29[k]
                  + f_4 * pc_y[k] * pfh_40[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_y, pb_z, pc_y, pc_z, pdi0_0, pdi0_27, \
                         pdh_20, pdi1_0, pdi1_27, pfh_41, pfh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_0 * pdh_20[k]
                  + f_4 * pc_y[k] * pfh_41[k];

        t_55[k] = pb_y[k] * pdi0_27[k]
                  - f_11 * pc_y[k] * pdi1_27[k];

        t_56[k] = pb_z[k] * pdi0_0[k]
                  - f_11 * pc_z[k] * pdi1_0[k];

        t_57[k] = f_4 * pc_y[k] * pfh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_z, pc_y, pc_z, pdi0_3, pdi0_5, pdh_0, \
                         pdh_2, pdi1_3, pdi1_5, pfh_42, pfh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_0 * pdh_0[k]
                  + f_4 * pc_z[k] * pfh_42[k];

        t_59[k] = pb_z[k] * pdi0_3[k]
                  - f_11 * pc_z[k] * pdi1_3[k];

        t_60[k] = f_4 * pc_y[k] * pfh_44[k];

        t_61[k] = pb_z[k] * pdi0_5[k]
                  + f_12 * pdh_2[k]
                  - f_11 * pc_z[k] * pdi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_z, pc_y, pc_z, pdi0_6, pdi0_9, pdh_3, \
                         pdh_5, pdi1_6, pdi1_9, pfh_45, pfh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * pdi0_6[k]
                  - f_11 * pc_z[k] * pdi1_6[k];

        t_63[k] = f_0 * pdh_3[k]
                  + f_4 * pc_z[k] * pfh_45[k];

        t_64[k] = f_4 * pc_y[k] * pfh_47[k];

        t_65[k] = pb_z[k] * pdi0_9[k]
                  + f_1 * pdh_5[k]
                  - f_11 * pc_z[k] * pdi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_z, pc_y, pc_z, pdi0_10, pdh_6, pdi1_10, \
                         pfg0_35, pfg1_35, pfh_48, pfh_50, pfh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * pdi0_10[k]
                  - f_11 * pc_z[k] * pdi1_10[k];

        t_67[k] = f_0 * pdh_6[k]
                  + f_4 * pc_z[k] * pfh_48[k];

        t_68[k] = f_5 * pfg0_35[k]
                  - f_6 * pfg1_35[k]
                  + f_4 * pc_y[k] * pfh_50[k];

        t_69[k] = f_4 * pc_y[k] * pfh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_z, pc_z, pdi0_14, pdi0_15, pdh_9, pdh_10, \
                         pdi1_14, pdi1_15, pfh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * pdi0_14[k]
                  + f_13 * pdh_9[k]
                  - f_11 * pc_z[k] * pdi1_14[k];

        t_71[k] = pb_z[k] * pdi0_15[k]
                  - f_11 * pc_z[k] * pdi1_15[k];

        t_72[k] = f_0 * pdh_10[k]
                  + f_4 * pc_z[k] * pfh_52[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pc_x, pc_y, sfh_59, sfh_60, sfh_62, pdh_59, \
                         pdh_60, pdh_62, pfh_56, pfh_59, pfh_60, \
                         pfh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * sfh_59[k]
                  + f_12 * pdh_59[k]
                  + f_4 * pc_x[k] * pfh_59[k];

        t_74[k] = f_0 * sfh_60[k]
                  + f_12 * pdh_60[k]
                  + f_4 * pc_x[k] * pfh_60[k];

        t_75[k] = f_4 * pc_y[k] * pfh_56[k];

        t_76[k] = f_0 * sfh_62[k]
                  + f_12 * pdh_62[k]
                  + f_4 * pc_x[k] * pfh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_z, pc_y, pc_z, pdi0_21, pdh_15, pdi1_21, \
                         pfg0_42, pfg1_42, pfh_57, pfh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_z[k] * pdi0_21[k]
                  - f_11 * pc_z[k] * pdi1_21[k];

        t_78[k] = f_0 * pdh_15[k]
                  + f_4 * pc_z[k] * pfh_57[k];

        t_79[k] = f_9 * pfg0_42[k]
                  - f_10 * pfg1_42[k]
                  + f_4 * pc_y[k] * pfh_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, pdh_20, pfg0_43, pfg0_44, \
                         pfg1_43, pfg1_44, pfh_60, pfh_61, pfh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_7 * pfg0_43[k]
                  - f_8 * pfg1_43[k]
                  + f_4 * pc_y[k] * pfh_60[k];

        t_81[k] = f_5 * pfg0_44[k]
                  - f_6 * pfg1_44[k]
                  + f_4 * pc_y[k] * pfh_61[k];

        t_82[k] = f_4 * pc_y[k] * pfh_62[k];

        t_83[k] = f_0 * pdh_20[k]
                  + f_2 * pfg0_44[k]
                  - f_3 * pfg1_44[k]
                  + f_4 * pc_z[k] * pfh_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, sfh_63, pdh_21, pdh_22, \
                         pdh_63, pfg0_45, pfg1_45, pfh_63, pfh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_0 * sfh_63[k]
                  + f_0 * pdh_63[k]
                  + f_2 * pfg0_45[k]
                  - f_3 * pfg1_45[k]
                  + f_4 * pc_x[k] * pfh_63[k];

        t_85[k] = f_12 * pdh_21[k]
                  + f_4 * pc_y[k] * pfh_63[k];

        t_86[k] = f_4 * pc_z[k] * pfh_63[k];

        t_87[k] = f_12 * pdh_22[k]
                  + f_5 * pfg0_45[k]
                  - f_6 * pfg1_45[k]
                  + f_4 * pc_y[k] * pfh_64[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pc_y, pc_z, pdh_23, pdh_24, pfg0_45, pfg0_46, \
                         pfg1_45, pfg1_46, pfh_65, pfh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * pdh_23[k]
                  + f_4 * pc_y[k] * pfh_65[k];

        t_89[k] = f_5 * pfg0_45[k]
                  - f_6 * pfg1_45[k]
                  + f_4 * pc_z[k] * pfh_65[k];

        t_90[k] = f_12 * pdh_24[k]
                  + f_7 * pfg0_46[k]
                  - f_8 * pfg1_46[k]
                  + f_4 * pc_y[k] * pfh_66[k];

        t_91[k] = f_4 * pc_z[k] * pfh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_y, pc_z, pdh_26, pdh_27, pfg0_47, pfg0_48, \
                         pfg1_47, pfg1_48, pfh_68, pfh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * pdh_26[k]
                  + f_4 * pc_y[k] * pfh_68[k];

        t_93[k] = f_7 * pfg0_47[k]
                  - f_8 * pfg1_47[k]
                  + f_4 * pc_z[k] * pfh_68[k];

        t_94[k] = f_12 * pdh_27[k]
                  + f_9 * pfg0_48[k]
                  - f_10 * pfg1_48[k]
                  + f_4 * pc_y[k] * pfh_69[k];

        t_95[k] = f_4 * pc_z[k] * pfh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, sfh_78, pdh_29, pdh_30, \
                         pdh_78, pfg0_50, pfg1_50, pfh_71, pfh_72, \
                         pfh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_12 * pdh_29[k]
                  + f_5 * pfg0_50[k]
                  - f_6 * pfg1_50[k]
                  + f_4 * pc_y[k] * pfh_71[k];

        t_97[k] = f_12 * pdh_30[k]
                  + f_4 * pc_y[k] * pfh_72[k];

        t_98[k] = f_9 * pfg0_50[k]
                  - f_10 * pfg1_50[k]
                  + f_4 * pc_z[k] * pfh_72[k];

        t_99[k] = f_0 * sfh_78[k]
                  + f_0 * pdh_78[k]
                  + f_4 * pc_x[k] * pfh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, sfh_80, sfh_81, pdh_35, \
                         pdh_80, pdh_81, pfh_73, pfh_77, pfh_80, \
                         pfh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_4 * pc_z[k] * pfh_73[k];

        t_101[k] = f_0 * sfh_80[k]
                   + f_0 * pdh_80[k]
                   + f_4 * pc_x[k] * pfh_80[k];

        t_102[k] = f_0 * sfh_81[k]
                   + f_0 * pdh_81[k]
                   + f_4 * pc_x[k] * pfh_81[k];

        t_103[k] = f_12 * pdh_35[k]
                   + f_4 * pc_y[k] * pfh_77[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pc_x, pc_y, pc_z, sfh_83, pdh_36, pdh_83, \
                         pfg0_55, pfg1_55, pfh_78, pfh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_0 * sfh_83[k]
                   + f_0 * pdh_83[k]
                   + f_4 * pc_x[k] * pfh_83[k];

        t_105[k] = f_12 * pdh_36[k]
                   + f_2 * pfg0_55[k]
                   - f_3 * pfg1_55[k]
                   + f_4 * pc_y[k] * pfh_78[k];

        t_106[k] = f_4 * pc_z[k] * pfh_78[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pc_y, pdh_38, pdh_39, pdh_40, pfg0_57, pfg0_58, \
                         pfg0_59, pfg1_57, pfg1_58, pfg1_59, pfh_80, pfh_81, \
                         pfh_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_12 * pdh_38[k]
                   + f_9 * pfg0_57[k]
                   - f_10 * pfg1_57[k]
                   + f_4 * pc_y[k] * pfh_80[k];

        t_108[k] = f_12 * pdh_39[k]
                   + f_7 * pfg0_58[k]
                   - f_8 * pfg1_58[k]
                   + f_4 * pc_y[k] * pfh_81[k];

        t_109[k] = f_12 * pdh_40[k]
                   + f_5 * pfg0_59[k]
                   - f_6 * pfg1_59[k]
                   + f_4 * pc_y[k] * pfh_82[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_y, pc_y, pc_z, pdi0_56, pdh_41, \
                         pdh_42, pdi1_56, pfg0_59, pfg1_59, pfh_83, \
                         pfh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_12 * pdh_41[k]
                   + f_4 * pc_y[k] * pfh_83[k];

        t_111[k] = f_2 * pfg0_59[k]
                   - f_3 * pfg1_59[k]
                   + f_4 * pc_z[k] * pfh_83[k];

        t_112[k] = pb_y[k] * pdi0_56[k]
                   - f_11 * pc_y[k] * pdi1_56[k];

        t_113[k] = f_0 * pdh_42[k]
                   + f_4 * pc_y[k] * pfh_84[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_y, pb_z, pc_y, pc_z, pdi0_31, pdi0_61, \
                         pdh_21, pdh_44, pdi1_31, pdi1_61, pfh_84, \
                         pfh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * pdh_21[k]
                   + f_4 * pc_z[k] * pfh_84[k];

        t_115[k] = pb_z[k] * pdi0_31[k]
                   - f_11 * pc_z[k] * pdi1_31[k];

        t_116[k] = f_0 * pdh_44[k]
                   + f_4 * pc_y[k] * pfh_86[k];

        t_117[k] = pb_y[k] * pdi0_61[k]
                   - f_11 * pc_y[k] * pdi1_61[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pb_y, pb_z, pc_y, pc_z, pdi0_34, pdi0_65, \
                         pdh_24, pdh_47, pdi1_34, pdi1_65, pfh_87, \
                         pfh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pb_z[k] * pdi0_34[k]
                   - f_11 * pc_z[k] * pdi1_34[k];

        t_119[k] = f_0 * pdh_24[k]
                   + f_4 * pc_z[k] * pfh_87[k];

        t_120[k] = f_0 * pdh_47[k]
                   + f_4 * pc_y[k] * pfh_89[k];

        t_121[k] = pb_y[k] * pdi0_65[k]
                   - f_11 * pc_y[k] * pdi1_65[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pb_z, pc_y, pc_z, pdi0_38, pdh_27, pdh_50, \
                         pdi1_38, pfg0_65, pfg1_65, pfh_90, pfh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = pb_z[k] * pdi0_38[k]
                   - f_11 * pc_z[k] * pdi1_38[k];

        t_123[k] = f_0 * pdh_27[k]
                   + f_4 * pc_z[k] * pfh_90[k];

        t_124[k] = f_0 * pdh_50[k]
                   + f_5 * pfg0_65[k]
                   - f_6 * pfg1_65[k]
                   + f_4 * pc_y[k] * pfh_92[k];
    }
}

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfi0, const size_t sfh,
                                                          const size_t sfi1, const size_t pdi0,
                                                          const size_t pdh, const size_t pdi1,
                                                          const size_t pfg0, const size_t pfg1,
                                                          const size_t pfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_168 = buffer.data(sfi0 + 168);
    const auto *sfi0_171 = buffer.data(sfi0 + 171);
    const auto *sfi0_174 = buffer.data(sfi0 + 174);
    const auto *sfi0_178 = buffer.data(sfi0 + 178);
    const auto *sfi0_180 = buffer.data(sfi0 + 180);
    const auto *sfi0_189 = buffer.data(sfi0 + 189);
    const auto *sfi0_191 = buffer.data(sfi0 + 191);
    const auto *sfi0_192 = buffer.data(sfi0 + 192);
    const auto *sfi0_193 = buffer.data(sfi0 + 193);
    const auto *sfi0_195 = buffer.data(sfi0 + 195);
    const auto *sfi0_201 = buffer.data(sfi0 + 201);
    const auto *sfi0_205 = buffer.data(sfi0 + 205);
    const auto *sfi0_208 = buffer.data(sfi0 + 208);
    const auto *sfi0_210 = buffer.data(sfi0 + 210);
    const auto *sfi0_217 = buffer.data(sfi0 + 217);
    const auto *sfi0_219 = buffer.data(sfi0 + 219);
    const auto *sfi0_220 = buffer.data(sfi0 + 220);
    const auto *sfi0_221 = buffer.data(sfi0 + 221);
    const auto *sfi0_223 = buffer.data(sfi0 + 223);
    const auto *sfi0_227 = buffer.data(sfi0 + 227);
    const auto *sfi0_230 = buffer.data(sfi0 + 230);
    const auto *sfi0_234 = buffer.data(sfi0 + 234);
    const auto *sfi0_236 = buffer.data(sfi0 + 236);

    const auto *sfh_101 = buffer.data(sfh + 101);
    const auto *sfh_102 = buffer.data(sfh + 102);
    const auto *sfh_105 = buffer.data(sfh + 105);
    const auto *sfh_120 = buffer.data(sfh + 120);
    const auto *sfh_122 = buffer.data(sfh + 122);
    const auto *sfh_123 = buffer.data(sfh + 123);
    const auto *sfh_125 = buffer.data(sfh + 125);
    const auto *sfh_126 = buffer.data(sfh + 126);
    const auto *sfh_129 = buffer.data(sfh + 129);
    const auto *sfh_132 = buffer.data(sfh + 132);
    const auto *sfh_136 = buffer.data(sfh + 136);
    const auto *sfh_138 = buffer.data(sfh + 138);
    const auto *sfh_141 = buffer.data(sfh + 141);
    const auto *sfh_143 = buffer.data(sfh + 143);
    const auto *sfh_144 = buffer.data(sfh + 144);
    const auto *sfh_146 = buffer.data(sfh + 146);
    const auto *sfh_152 = buffer.data(sfh + 152);
    const auto *sfh_156 = buffer.data(sfh + 156);
    const auto *sfh_159 = buffer.data(sfh + 159);
    const auto *sfh_161 = buffer.data(sfh + 161);
    const auto *sfh_162 = buffer.data(sfh + 162);
    const auto *sfh_164 = buffer.data(sfh + 164);
    const auto *sfh_165 = buffer.data(sfh + 165);
    const auto *sfh_167 = buffer.data(sfh + 167);
    const auto *sfh_171 = buffer.data(sfh + 171);
    const auto *sfh_174 = buffer.data(sfh + 174);
    const auto *sfh_178 = buffer.data(sfh + 178);
    const auto *sfh_180 = buffer.data(sfh + 180);
    const auto *sfh_183 = buffer.data(sfh + 183);
    const auto *sfh_185 = buffer.data(sfh + 185);
    const auto *sfh_186 = buffer.data(sfh + 186);

    const auto *sfi1_168 = buffer.data(sfi1 + 168);
    const auto *sfi1_171 = buffer.data(sfi1 + 171);
    const auto *sfi1_174 = buffer.data(sfi1 + 174);
    const auto *sfi1_178 = buffer.data(sfi1 + 178);
    const auto *sfi1_180 = buffer.data(sfi1 + 180);
    const auto *sfi1_189 = buffer.data(sfi1 + 189);
    const auto *sfi1_191 = buffer.data(sfi1 + 191);
    const auto *sfi1_192 = buffer.data(sfi1 + 192);
    const auto *sfi1_193 = buffer.data(sfi1 + 193);
    const auto *sfi1_195 = buffer.data(sfi1 + 195);
    const auto *sfi1_201 = buffer.data(sfi1 + 201);
    const auto *sfi1_205 = buffer.data(sfi1 + 205);
    const auto *sfi1_208 = buffer.data(sfi1 + 208);
    const auto *sfi1_210 = buffer.data(sfi1 + 210);
    const auto *sfi1_217 = buffer.data(sfi1 + 217);
    const auto *sfi1_219 = buffer.data(sfi1 + 219);
    const auto *sfi1_220 = buffer.data(sfi1 + 220);
    const auto *sfi1_221 = buffer.data(sfi1 + 221);
    const auto *sfi1_223 = buffer.data(sfi1 + 223);
    const auto *sfi1_227 = buffer.data(sfi1 + 227);
    const auto *sfi1_230 = buffer.data(sfi1 + 230);
    const auto *sfi1_234 = buffer.data(sfi1 + 234);
    const auto *sfi1_236 = buffer.data(sfi1 + 236);

    const auto *pdi0_43 = buffer.data(pdi0 + 43);
    const auto *pdi0_70 = buffer.data(pdi0 + 70);
    const auto *pdi0_76 = buffer.data(pdi0 + 76);
    const auto *pdi0_84 = buffer.data(pdi0 + 84);
    const auto *pdi0_87 = buffer.data(pdi0 + 87);
    const auto *pdi0_90 = buffer.data(pdi0 + 90);
    const auto *pdi0_94 = buffer.data(pdi0 + 94);
    const auto *pdi0_140 = buffer.data(pdi0 + 140);
    const auto *pdi0_145 = buffer.data(pdi0 + 145);
    const auto *pdi0_149 = buffer.data(pdi0 + 149);
    const auto *pdi0_154 = buffer.data(pdi0 + 154);

    const auto *pdh_31 = buffer.data(pdh + 31);
    const auto *pdh_36 = buffer.data(pdh + 36);
    const auto *pdh_41 = buffer.data(pdh + 41);
    const auto *pdh_42 = buffer.data(pdh + 42);
    const auto *pdh_44 = buffer.data(pdh + 44);
    const auto *pdh_45 = buffer.data(pdh + 45);
    const auto *pdh_47 = buffer.data(pdh + 47);
    const auto *pdh_48 = buffer.data(pdh + 48);
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
    const auto *pdh_83 = buffer.data(pdh + 83);
    const auto *pdh_84 = buffer.data(pdh + 84);
    const auto *pdh_86 = buffer.data(pdh + 86);
    const auto *pdh_87 = buffer.data(pdh + 87);
    const auto *pdh_89 = buffer.data(pdh + 89);
    const auto *pdh_90 = buffer.data(pdh + 90);
    const auto *pdh_93 = buffer.data(pdh + 93);
    const auto *pdh_94 = buffer.data(pdh + 94);
    const auto *pdh_98 = buffer.data(pdh + 98);
    const auto *pdh_101 = buffer.data(pdh + 101);
    const auto *pdh_102 = buffer.data(pdh + 102);
    const auto *pdh_104 = buffer.data(pdh + 104);
    const auto *pdh_105 = buffer.data(pdh + 105);
    const auto *pdh_107 = buffer.data(pdh + 107);
    const auto *pdh_110 = buffer.data(pdh + 110);
    const auto *pdh_114 = buffer.data(pdh + 114);
    const auto *pdh_120 = buffer.data(pdh + 120);
    const auto *pdh_122 = buffer.data(pdh + 122);
    const auto *pdh_123 = buffer.data(pdh + 123);
    const auto *pdh_125 = buffer.data(pdh + 125);

    const auto *pdi1_43 = buffer.data(pdi1 + 43);
    const auto *pdi1_70 = buffer.data(pdi1 + 70);
    const auto *pdi1_76 = buffer.data(pdi1 + 76);
    const auto *pdi1_84 = buffer.data(pdi1 + 84);
    const auto *pdi1_87 = buffer.data(pdi1 + 87);
    const auto *pdi1_90 = buffer.data(pdi1 + 90);
    const auto *pdi1_94 = buffer.data(pdi1 + 94);
    const auto *pdi1_140 = buffer.data(pdi1 + 140);
    const auto *pdi1_145 = buffer.data(pdi1 + 145);
    const auto *pdi1_149 = buffer.data(pdi1 + 149);
    const auto *pdi1_154 = buffer.data(pdi1 + 154);

    const auto *pfg0_70 = buffer.data(pfg0 + 70);
    const auto *pfg0_72 = buffer.data(pfg0 + 72);
    const auto *pfg0_73 = buffer.data(pfg0 + 73);
    const auto *pfg0_74 = buffer.data(pfg0 + 74);
    const auto *pfg0_75 = buffer.data(pfg0 + 75);
    const auto *pfg0_76 = buffer.data(pfg0 + 76);
    const auto *pfg0_77 = buffer.data(pfg0 + 77);
    const auto *pfg0_78 = buffer.data(pfg0 + 78);
    const auto *pfg0_80 = buffer.data(pfg0 + 80);
    const auto *pfg0_85 = buffer.data(pfg0 + 85);
    const auto *pfg0_87 = buffer.data(pfg0 + 87);
    const auto *pfg0_88 = buffer.data(pfg0 + 88);
    const auto *pfg0_89 = buffer.data(pfg0 + 89);
    const auto *pfg0_90 = buffer.data(pfg0 + 90);
    const auto *pfg0_92 = buffer.data(pfg0 + 92);
    const auto *pfg0_95 = buffer.data(pfg0 + 95);

    const auto *pfg1_70 = buffer.data(pfg1 + 70);
    const auto *pfg1_72 = buffer.data(pfg1 + 72);
    const auto *pfg1_73 = buffer.data(pfg1 + 73);
    const auto *pfg1_74 = buffer.data(pfg1 + 74);
    const auto *pfg1_75 = buffer.data(pfg1 + 75);
    const auto *pfg1_76 = buffer.data(pfg1 + 76);
    const auto *pfg1_77 = buffer.data(pfg1 + 77);
    const auto *pfg1_78 = buffer.data(pfg1 + 78);
    const auto *pfg1_80 = buffer.data(pfg1 + 80);
    const auto *pfg1_85 = buffer.data(pfg1 + 85);
    const auto *pfg1_87 = buffer.data(pfg1 + 87);
    const auto *pfg1_88 = buffer.data(pfg1 + 88);
    const auto *pfg1_89 = buffer.data(pfg1 + 89);
    const auto *pfg1_90 = buffer.data(pfg1 + 90);
    const auto *pfg1_92 = buffer.data(pfg1 + 92);
    const auto *pfg1_95 = buffer.data(pfg1 + 95);

    const auto *pfh_93 = buffer.data(pfh + 93);
    const auto *pfh_94 = buffer.data(pfh + 94);
    const auto *pfh_98 = buffer.data(pfh + 98);
    const auto *pfh_99 = buffer.data(pfh + 99);
    const auto *pfh_101 = buffer.data(pfh + 101);
    const auto *pfh_102 = buffer.data(pfh + 102);
    const auto *pfh_103 = buffer.data(pfh + 103);
    const auto *pfh_104 = buffer.data(pfh + 104);
    const auto *pfh_105 = buffer.data(pfh + 105);
    const auto *pfh_106 = buffer.data(pfh + 106);
    const auto *pfh_107 = buffer.data(pfh + 107);
    const auto *pfh_108 = buffer.data(pfh + 108);
    const auto *pfh_110 = buffer.data(pfh + 110);
    const auto *pfh_111 = buffer.data(pfh + 111);
    const auto *pfh_113 = buffer.data(pfh + 113);
    const auto *pfh_114 = buffer.data(pfh + 114);
    const auto *pfh_115 = buffer.data(pfh + 115);
    const auto *pfh_119 = buffer.data(pfh + 119);
    const auto *pfh_120 = buffer.data(pfh + 120);
    const auto *pfh_122 = buffer.data(pfh + 122);
    const auto *pfh_123 = buffer.data(pfh + 123);
    const auto *pfh_124 = buffer.data(pfh + 124);
    const auto *pfh_125 = buffer.data(pfh + 125);
    const auto *pfh_126 = buffer.data(pfh + 126);
    const auto *pfh_128 = buffer.data(pfh + 128);
    const auto *pfh_129 = buffer.data(pfh + 129);
    const auto *pfh_131 = buffer.data(pfh + 131);
    const auto *pfh_132 = buffer.data(pfh + 132);
    const auto *pfh_135 = buffer.data(pfh + 135);
    const auto *pfh_136 = buffer.data(pfh + 136);
    const auto *pfh_140 = buffer.data(pfh + 140);
    const auto *pfh_141 = buffer.data(pfh + 141);
    const auto *pfh_143 = buffer.data(pfh + 143);
    const auto *pfh_144 = buffer.data(pfh + 144);
    const auto *pfh_146 = buffer.data(pfh + 146);
    const auto *pfh_147 = buffer.data(pfh + 147);
    const auto *pfh_149 = buffer.data(pfh + 149);
    const auto *pfh_150 = buffer.data(pfh + 150);
    const auto *pfh_152 = buffer.data(pfh + 152);
    const auto *pfh_153 = buffer.data(pfh + 153);
    const auto *pfh_156 = buffer.data(pfh + 156);
    const auto *pfh_157 = buffer.data(pfh + 157);
    const auto *pfh_161 = buffer.data(pfh + 161);
    const auto *pfh_162 = buffer.data(pfh + 162);
    const auto *pfh_164 = buffer.data(pfh + 164);
    const auto *pfh_165 = buffer.data(pfh + 165);
    const auto *pfh_167 = buffer.data(pfh + 167);
    const auto *pfh_168 = buffer.data(pfh + 168);
    const auto *pfh_170 = buffer.data(pfh + 170);
    const auto *pfh_171 = buffer.data(pfh + 171);
    const auto *pfh_173 = buffer.data(pfh + 173);
    const auto *pfh_174 = buffer.data(pfh + 174);
    const auto *pfh_177 = buffer.data(pfh + 177);
    const auto *pfh_178 = buffer.data(pfh + 178);
    const auto *pfh_183 = buffer.data(pfh + 183);
    const auto *pfh_185 = buffer.data(pfh + 185);
    const auto *pfh_186 = buffer.data(pfh + 186);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, pc_y, pc_z, pdi0_43, pdi0_70, \
                         pdh_31, pdh_51, pdi1_43, pdi1_70, pfh_93, \
                         pfh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * pdh_51[k]
                   + f_4 * pc_y[k] * pfh_93[k];

        t_126[k] = pb_y[k] * pdi0_70[k]
                   - f_11 * pc_y[k] * pdi1_70[k];

        t_127[k] = pb_z[k] * pdi0_43[k]
                   - f_11 * pc_z[k] * pdi1_43[k];

        t_128[k] = f_0 * pdh_31[k]
                   + f_4 * pc_z[k] * pfh_94[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, sfh_101, sfh_102, pdh_56, pdh_101, \
                         pdh_102, pfh_98, pfh_101, pfh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_0 * sfh_101[k]
                   + f_0 * pdh_101[k]
                   + f_4 * pc_x[k] * pfh_101[k];

        t_130[k] = f_0 * sfh_102[k]
                   + f_0 * pdh_102[k]
                   + f_4 * pc_x[k] * pfh_102[k];

        t_131[k] = f_0 * pdh_56[k]
                   + f_4 * pc_y[k] * pfh_98[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, pc_z, pdi0_76, pdh_36, pdh_57, \
                         pdi1_76, pfg0_70, pfg1_70, pfh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_y[k] * pdi0_76[k]
                   - f_11 * pc_y[k] * pdi1_76[k];

        t_133[k] = f_0 * pdh_57[k]
                   + f_2 * pfg0_70[k]
                   - f_3 * pfg1_70[k]
                   + f_4 * pc_y[k] * pfh_99[k];

        t_134[k] = f_0 * pdh_36[k]
                   + f_4 * pc_z[k] * pfh_99[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pc_y, pdh_59, pdh_60, pdh_61, pfg0_72, pfg0_73, \
                         pfg0_74, pfg1_72, pfg1_73, pfg1_74, pfh_101, pfh_102, \
                         pfh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * pdh_59[k]
                   + f_9 * pfg0_72[k]
                   - f_10 * pfg1_72[k]
                   + f_4 * pc_y[k] * pfh_101[k];

        t_136[k] = f_0 * pdh_60[k]
                   + f_7 * pfg0_73[k]
                   - f_8 * pfg1_73[k]
                   + f_4 * pc_y[k] * pfh_102[k];

        t_137[k] = f_0 * pdh_61[k]
                   + f_5 * pfg0_74[k]
                   - f_6 * pfg1_74[k]
                   + f_4 * pc_y[k] * pfh_103[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, pc_z, sfh_105, pdh_41, pdh_62, \
                         pdh_105, pfg0_74, pfg0_75, pfg1_74, pfg1_75, pfh_104, \
                         pfh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_0 * pdh_62[k]
                   + f_4 * pc_y[k] * pfh_104[k];

        t_139[k] = f_0 * pdh_41[k]
                   + f_2 * pfg0_74[k]
                   - f_3 * pfg1_74[k]
                   + f_4 * pc_z[k] * pfh_104[k];

        t_140[k] = f_0 * sfh_105[k]
                   + f_0 * pdh_105[k]
                   + f_2 * pfg0_75[k]
                   - f_3 * pfg1_75[k]
                   + f_4 * pc_x[k] * pfh_105[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pc_y, pc_z, pdh_42, pdh_44, \
                         pfg0_75, pfg1_75, pfh_105, pfh_106, pfh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_4 * pc_y[k] * pfh_105[k];

        t_142[k] = f_12 * pdh_42[k]
                   + f_4 * pc_z[k] * pfh_105[k];

        t_143[k] = f_5 * pfg0_75[k]
                   - f_6 * pfg1_75[k]
                   + f_4 * pc_y[k] * pfh_106[k];

        t_144[k] = f_4 * pc_y[k] * pfh_107[k];

        t_145[k] = f_12 * pdh_44[k]
                   + f_5 * pfg0_75[k]
                   - f_6 * pfg1_75[k]
                   + f_4 * pc_z[k] * pfh_107[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, pdh_45, pdh_47, pfg0_76, \
                         pfg0_77, pfg1_76, pfg1_77, pfh_108, pfh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * pfg0_76[k]
                   - f_8 * pfg1_76[k]
                   + f_4 * pc_y[k] * pfh_108[k];

        t_147[k] = f_12 * pdh_45[k]
                   + f_4 * pc_z[k] * pfh_108[k];

        t_148[k] = f_4 * pc_y[k] * pfh_110[k];

        t_149[k] = f_12 * pdh_47[k]
                   + f_7 * pfg0_77[k]
                   - f_8 * pfg1_77[k]
                   + f_4 * pc_z[k] * pfh_110[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pc_y, pc_z, pdh_48, pdh_51, \
                         pfg0_78, pfg0_80, pfg1_78, pfg1_80, pfh_111, pfh_113, \
                         pfh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_9 * pfg0_78[k]
                   - f_10 * pfg1_78[k]
                   + f_4 * pc_y[k] * pfh_111[k];

        t_151[k] = f_12 * pdh_48[k]
                   + f_4 * pc_z[k] * pfh_111[k];

        t_152[k] = f_5 * pfg0_80[k]
                   - f_6 * pfg1_80[k]
                   + f_4 * pc_y[k] * pfh_113[k];

        t_153[k] = f_4 * pc_y[k] * pfh_114[k];

        t_154[k] = f_12 * pdh_51[k]
                   + f_9 * pfg0_80[k]
                   - f_10 * pfg1_80[k]
                   + f_4 * pc_z[k] * pfh_114[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pc_x, pc_z, sfh_120, sfh_122, pdh_52, pdh_120, \
                         pdh_122, pfh_115, pfh_120, pfh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_0 * sfh_120[k]
                   + f_0 * pdh_120[k]
                   + f_4 * pc_x[k] * pfh_120[k];

        t_156[k] = f_12 * pdh_52[k]
                   + f_4 * pc_z[k] * pfh_115[k];

        t_157[k] = f_0 * sfh_122[k]
                   + f_0 * pdh_122[k]
                   + f_4 * pc_x[k] * pfh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pc_x, pc_y, sfh_123, sfh_125, pdh_123, \
                         pdh_125, pfg0_85, pfg1_85, pfh_119, pfh_120, pfh_123, \
                         pfh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * sfh_123[k]
                   + f_0 * pdh_123[k]
                   + f_4 * pc_x[k] * pfh_123[k];

        t_159[k] = f_4 * pc_y[k] * pfh_119[k];

        t_160[k] = f_0 * sfh_125[k]
                   + f_0 * pdh_125[k]
                   + f_4 * pc_x[k] * pfh_125[k];

        t_161[k] = f_2 * pfg0_85[k]
                   - f_3 * pfg1_85[k]
                   + f_4 * pc_y[k] * pfh_120[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pc_y, pc_z, pdh_57, pfg0_87, pfg0_88, pfg1_87, \
                         pfg1_88, pfh_120, pfh_122, pfh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_12 * pdh_57[k]
                   + f_4 * pc_z[k] * pfh_120[k];

        t_163[k] = f_9 * pfg0_87[k]
                   - f_10 * pfg1_87[k]
                   + f_4 * pc_y[k] * pfh_122[k];

        t_164[k] = f_7 * pfg0_88[k]
                   - f_8 * pfg1_88[k]
                   + f_4 * pc_y[k] * pfh_123[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pa_x, pc_x, pc_y, pc_z, sfi0_168, \
                         sfh_126, sfi1_168, pdh_62, pfg0_89, pfg1_89, pfh_124, \
                         pfh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_5 * pfg0_89[k]
                   - f_6 * pfg1_89[k]
                   + f_4 * pc_y[k] * pfh_124[k];

        t_166[k] = f_4 * pc_y[k] * pfh_125[k];

        t_167[k] = f_12 * pdh_62[k]
                   + f_2 * pfg0_89[k]
                   - f_3 * pfg1_89[k]
                   + f_4 * pc_z[k] * pfh_125[k];

        t_168[k] = pa_x[k] * sfi0_168[k]
                   + f_14 * sfh_126[k]
                   - f_11 * pc_x[k] * sfi1_168[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_x, pc_x, pc_y, pc_z, sfi0_171, \
                         sfh_129, sfi1_171, pdh_63, pdh_65, pfh_126, \
                         pfh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_1 * pdh_63[k]
                   + f_4 * pc_y[k] * pfh_126[k];

        t_170[k] = f_4 * pc_z[k] * pfh_126[k];

        t_171[k] = pa_x[k] * sfi0_171[k]
                   + f_13 * sfh_129[k]
                   - f_11 * pc_x[k] * sfi1_171[k];

        t_172[k] = f_1 * pdh_65[k]
                   + f_4 * pc_y[k] * pfh_128[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pa_x, pc_x, pc_z, sfi0_174, sfh_132, sfi1_174, \
                         pfg0_90, pfg1_90, pfh_128, pfh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_5 * pfg0_90[k]
                   - f_6 * pfg1_90[k]
                   + f_4 * pc_z[k] * pfh_128[k];

        t_174[k] = pa_x[k] * sfi0_174[k]
                   + f_1 * sfh_132[k]
                   - f_11 * pc_x[k] * sfi1_174[k];

        t_175[k] = f_4 * pc_z[k] * pfh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pc_x, pc_y, pc_z, sfi0_178, \
                         sfh_136, sfi1_178, pdh_68, pfg0_92, pfg1_92, pfh_131, \
                         pfh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_1 * pdh_68[k]
                   + f_4 * pc_y[k] * pfh_131[k];

        t_177[k] = f_7 * pfg0_92[k]
                   - f_8 * pfg1_92[k]
                   + f_4 * pc_z[k] * pfh_131[k];

        t_178[k] = pa_x[k] * sfi0_178[k]
                   + f_12 * sfh_136[k]
                   - f_11 * pc_x[k] * sfi1_178[k];

        t_179[k] = f_4 * pc_z[k] * pfh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_x, pc_x, pc_y, pc_z, sfi0_180, sfh_138, \
                         sfi1_180, pdh_72, pfg0_95, pfg1_95, pfh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_x[k] * sfi0_180[k]
                   + f_12 * sfh_138[k]
                   - f_11 * pc_x[k] * sfi1_180[k];

        t_181[k] = f_1 * pdh_72[k]
                   + f_4 * pc_y[k] * pfh_135[k];

        t_182[k] = f_9 * pfg0_95[k]
                   - f_10 * pfg1_95[k]
                   + f_4 * pc_z[k] * pfh_135[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pc_x, pc_z, sfh_141, sfh_143, sfh_144, \
                         pfh_136, pfh_141, pfh_143, pfh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_0 * sfh_141[k]
                   + f_4 * pc_x[k] * pfh_141[k];

        t_184[k] = f_4 * pc_z[k] * pfh_136[k];

        t_185[k] = f_0 * sfh_143[k]
                   + f_4 * pc_x[k] * pfh_143[k];

        t_186[k] = f_0 * sfh_144[k]
                   + f_4 * pc_x[k] * pfh_144[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pa_x, pc_x, pc_y, pc_z, sfi0_189, \
                         sfh_146, sfi1_189, pdh_77, pfh_140, pfh_141, \
                         pfh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_1 * pdh_77[k]
                   + f_4 * pc_y[k] * pfh_140[k];

        t_188[k] = f_0 * sfh_146[k]
                   + f_4 * pc_x[k] * pfh_146[k];

        t_189[k] = pa_x[k] * sfi0_189[k]
                   - f_11 * pc_x[k] * sfi1_189[k];

        t_190[k] = f_4 * pc_z[k] * pfh_141[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_x, pc_x, pc_y, sfi0_191, sfi0_192, \
                         sfi0_193, sfi1_191, sfi1_192, sfi1_193, pdh_83, \
                         pfh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = pa_x[k] * sfi0_191[k]
                   - f_11 * pc_x[k] * sfi1_191[k];

        t_192[k] = pa_x[k] * sfi0_192[k]
                   - f_11 * pc_x[k] * sfi1_192[k];

        t_193[k] = pa_x[k] * sfi0_193[k]
                   - f_11 * pc_x[k] * sfi1_193[k];

        t_194[k] = f_1 * pdh_83[k]
                   + f_4 * pc_y[k] * pfh_146[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_x, pb_z, pc_x, pc_y, pc_z, sfi0_195, \
                         sfi1_195, pdi0_84, pdh_63, pdh_84, pdi1_84, \
                         pfh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_x[k] * sfi0_195[k]
                   - f_11 * pc_x[k] * sfi1_195[k];

        t_196[k] = pb_z[k] * pdi0_84[k]
                   - f_11 * pc_z[k] * pdi1_84[k];

        t_197[k] = f_12 * pdh_84[k]
                   + f_4 * pc_y[k] * pfh_147[k];

        t_198[k] = f_0 * pdh_63[k]
                   + f_4 * pc_z[k] * pfh_147[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, pa_x, pb_z, pc_x, pc_y, pc_z, sfi0_201, sfh_152, \
                         sfi1_201, pdi0_87, pdh_86, pdi1_87, pfh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_z[k] * pdi0_87[k]
                   - f_11 * pc_z[k] * pdi1_87[k];

        t_200[k] = f_12 * pdh_86[k]
                   + f_4 * pc_y[k] * pfh_149[k];

        t_201[k] = pa_x[k] * sfi0_201[k]
                   + f_13 * sfh_152[k]
                   - f_11 * pc_x[k] * sfi1_201[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pb_z, pc_y, pc_z, pdi0_90, pdh_66, pdh_89, \
                         pdi1_90, pfh_150, pfh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = pb_z[k] * pdi0_90[k]
                   - f_11 * pc_z[k] * pdi1_90[k];

        t_203[k] = f_0 * pdh_66[k]
                   + f_4 * pc_z[k] * pfh_150[k];

        t_204[k] = f_12 * pdh_89[k]
                   + f_4 * pc_y[k] * pfh_152[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pa_x, pb_z, pc_x, pc_z, sfi0_205, sfh_156, \
                         sfi1_205, pdi0_94, pdh_69, pdi1_94, pfh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = pa_x[k] * sfi0_205[k]
                   + f_1 * sfh_156[k]
                   - f_11 * pc_x[k] * sfi1_205[k];

        t_206[k] = pb_z[k] * pdi0_94[k]
                   - f_11 * pc_z[k] * pdi1_94[k];

        t_207[k] = f_0 * pdh_69[k]
                   + f_4 * pc_z[k] * pfh_153[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_x, pc_x, pc_y, sfi0_208, sfi0_210, sfh_159, \
                         sfh_161, sfi1_208, sfi1_210, pdh_93, pfh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_x[k] * sfi0_208[k]
                   + f_12 * sfh_159[k]
                   - f_11 * pc_x[k] * sfi1_208[k];

        t_209[k] = f_12 * pdh_93[k]
                   + f_4 * pc_y[k] * pfh_156[k];

        t_210[k] = pa_x[k] * sfi0_210[k]
                   + f_12 * sfh_161[k]
                   - f_11 * pc_x[k] * sfi1_210[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pc_x, pc_z, sfh_162, sfh_164, sfh_165, \
                         pdh_73, pfh_157, pfh_162, pfh_164, pfh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_0 * sfh_162[k]
                   + f_4 * pc_x[k] * pfh_162[k];

        t_212[k] = f_0 * pdh_73[k]
                   + f_4 * pc_z[k] * pfh_157[k];

        t_213[k] = f_0 * sfh_164[k]
                   + f_4 * pc_x[k] * pfh_164[k];

        t_214[k] = f_0 * sfh_165[k]
                   + f_4 * pc_x[k] * pfh_165[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pa_x, pc_x, pc_y, pc_z, sfi0_217, \
                         sfh_167, sfi1_217, pdh_78, pdh_98, pfh_161, pfh_162, \
                         pfh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_12 * pdh_98[k]
                   + f_4 * pc_y[k] * pfh_161[k];

        t_216[k] = f_0 * sfh_167[k]
                   + f_4 * pc_x[k] * pfh_167[k];

        t_217[k] = pa_x[k] * sfi0_217[k]
                   - f_11 * pc_x[k] * sfi1_217[k];

        t_218[k] = f_0 * pdh_78[k]
                   + f_4 * pc_z[k] * pfh_162[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pa_x, pc_x, pc_y, sfi0_219, sfi0_220, \
                         sfi0_221, sfi1_219, sfi1_220, sfi1_221, pdh_104, \
                         pfh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = pa_x[k] * sfi0_219[k]
                   - f_11 * pc_x[k] * sfi1_219[k];

        t_220[k] = pa_x[k] * sfi0_220[k]
                   - f_11 * pc_x[k] * sfi1_220[k];

        t_221[k] = pa_x[k] * sfi0_221[k]
                   - f_11 * pc_x[k] * sfi1_221[k];

        t_222[k] = f_12 * pdh_104[k]
                   + f_4 * pc_y[k] * pfh_167[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_x, pb_y, pc_x, pc_y, pc_z, sfi0_223, \
                         sfi1_223, pdi0_140, pdh_84, pdh_105, pdi1_140, \
                         pfh_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pa_x[k] * sfi0_223[k]
                   - f_11 * pc_x[k] * sfi1_223[k];

        t_224[k] = pb_y[k] * pdi0_140[k]
                   - f_11 * pc_y[k] * pdi1_140[k];

        t_225[k] = f_0 * pdh_105[k]
                   + f_4 * pc_y[k] * pfh_168[k];

        t_226[k] = f_12 * pdh_84[k]
                   + f_4 * pc_z[k] * pfh_168[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pa_x, pb_y, pc_x, pc_y, sfi0_227, sfh_171, \
                         sfi1_227, pdi0_145, pdh_107, pdi1_145, \
                         pfh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pa_x[k] * sfi0_227[k]
                   + f_13 * sfh_171[k]
                   - f_11 * pc_x[k] * sfi1_227[k];

        t_228[k] = f_0 * pdh_107[k]
                   + f_4 * pc_y[k] * pfh_170[k];

        t_229[k] = pb_y[k] * pdi0_145[k]
                   - f_11 * pc_y[k] * pdi1_145[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pa_x, pc_x, pc_y, pc_z, sfi0_230, sfh_174, \
                         sfi1_230, pdh_87, pdh_110, pfh_171, pfh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pa_x[k] * sfi0_230[k]
                   + f_1 * sfh_174[k]
                   - f_11 * pc_x[k] * sfi1_230[k];

        t_231[k] = f_12 * pdh_87[k]
                   + f_4 * pc_z[k] * pfh_171[k];

        t_232[k] = f_0 * pdh_110[k]
                   + f_4 * pc_y[k] * pfh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_x, pb_y, pc_x, pc_y, pc_z, sfi0_234, sfh_178, \
                         sfi1_234, pdi0_149, pdh_90, pdi1_149, \
                         pfh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * pdi0_149[k]
                   - f_11 * pc_y[k] * pdi1_149[k];

        t_234[k] = pa_x[k] * sfi0_234[k]
                   + f_12 * sfh_178[k]
                   - f_11 * pc_x[k] * sfi1_234[k];

        t_235[k] = f_12 * pdh_90[k]
                   + f_4 * pc_z[k] * pfh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_x, pb_y, pc_x, pc_y, sfi0_236, sfh_180, \
                         sfi1_236, pdi0_154, pdh_114, pdi1_154, \
                         pfh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_x[k] * sfi0_236[k]
                   + f_12 * sfh_180[k]
                   - f_11 * pc_x[k] * sfi1_236[k];

        t_237[k] = f_0 * pdh_114[k]
                   + f_4 * pc_y[k] * pfh_177[k];

        t_238[k] = pb_y[k] * pdi0_154[k]
                   - f_11 * pc_y[k] * pdi1_154[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pc_x, pc_z, sfh_183, sfh_185, sfh_186, \
                         pdh_94, pfh_178, pfh_183, pfh_185, pfh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_0 * sfh_183[k]
                   + f_4 * pc_x[k] * pfh_183[k];

        t_240[k] = f_12 * pdh_94[k]
                   + f_4 * pc_z[k] * pfh_178[k];

        t_241[k] = f_0 * sfh_185[k]
                   + f_4 * pc_x[k] * pfh_185[k];

        t_242[k] = f_0 * sfh_186[k]
                   + f_4 * pc_x[k] * pfh_186[k];
    }
}

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfi0, const size_t sfh,
                                                          const size_t sfi1, const size_t ppi0,
                                                          const size_t ppi1, const size_t pdi0,
                                                          const size_t pdh, const size_t pdi1,
                                                          const size_t pfg0, const size_t pfg1,
                                                          const size_t pfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_0 = buffer.data(sfi0 + 0);
    const auto *sfi0_1 = buffer.data(sfi0 + 1);
    const auto *sfi0_3 = buffer.data(sfi0 + 3);
    const auto *sfi0_5 = buffer.data(sfi0 + 5);
    const auto *sfi0_6 = buffer.data(sfi0 + 6);
    const auto *sfi0_8 = buffer.data(sfi0 + 8);
    const auto *sfi0_9 = buffer.data(sfi0 + 9);
    const auto *sfi0_10 = buffer.data(sfi0 + 10);
    const auto *sfi0_12 = buffer.data(sfi0 + 12);
    const auto *sfi0_13 = buffer.data(sfi0 + 13);
    const auto *sfi0_14 = buffer.data(sfi0 + 14);
    const auto *sfi0_21 = buffer.data(sfi0 + 21);
    const auto *sfi0_27 = buffer.data(sfi0 + 27);
    const auto *sfi0_56 = buffer.data(sfi0 + 56);
    const auto *sfi0_61 = buffer.data(sfi0 + 61);
    const auto *sfi0_64 = buffer.data(sfi0 + 64);
    const auto *sfi0_65 = buffer.data(sfi0 + 65);
    const auto *sfi0_68 = buffer.data(sfi0 + 68);
    const auto *sfi0_69 = buffer.data(sfi0 + 69);
    const auto *sfi0_70 = buffer.data(sfi0 + 70);
    const auto *sfi0_79 = buffer.data(sfi0 + 79);
    const auto *sfi0_80 = buffer.data(sfi0 + 80);
    const auto *sfi0_81 = buffer.data(sfi0 + 81);
    const auto *sfi0_83 = buffer.data(sfi0 + 83);
    const auto *sfi0_245 = buffer.data(sfi0 + 245);
    const auto *sfi0_247 = buffer.data(sfi0 + 247);
    const auto *sfi0_248 = buffer.data(sfi0 + 248);
    const auto *sfi0_249 = buffer.data(sfi0 + 249);
    const auto *sfi0_251 = buffer.data(sfi0 + 251);
    const auto *sfi0_252 = buffer.data(sfi0 + 252);
    const auto *sfi0_257 = buffer.data(sfi0 + 257);
    const auto *sfi0_261 = buffer.data(sfi0 + 261);
    const auto *sfi0_266 = buffer.data(sfi0 + 266);
    const auto *sfi0_273 = buffer.data(sfi0 + 273);
    const auto *sfi0_275 = buffer.data(sfi0 + 275);
    const auto *sfi0_276 = buffer.data(sfi0 + 276);
    const auto *sfi0_277 = buffer.data(sfi0 + 277);
    const auto *sfi0_279 = buffer.data(sfi0 + 279);

    const auto *sfh_0 = buffer.data(sfh + 0);
    const auto *sfh_1 = buffer.data(sfh + 1);
    const auto *sfh_3 = buffer.data(sfh + 3);
    const auto *sfh_5 = buffer.data(sfh + 5);
    const auto *sfh_6 = buffer.data(sfh + 6);
    const auto *sfh_8 = buffer.data(sfh + 8);
    const auto *sfh_9 = buffer.data(sfh + 9);
    const auto *sfh_15 = buffer.data(sfh + 15);
    const auto *sfh_20 = buffer.data(sfh + 20);
    const auto *sfh_41 = buffer.data(sfh + 41);
    const auto *sfh_47 = buffer.data(sfh + 47);
    const auto *sfh_50 = buffer.data(sfh + 50);
    const auto *sfh_51 = buffer.data(sfh + 51);
    const auto *sfh_59 = buffer.data(sfh + 59);
    const auto *sfh_60 = buffer.data(sfh + 60);
    const auto *sfh_61 = buffer.data(sfh + 61);
    const auto *sfh_62 = buffer.data(sfh + 62);
    const auto *sfh_188 = buffer.data(sfh + 188);
    const auto *sfh_189 = buffer.data(sfh + 189);
    const auto *sfh_194 = buffer.data(sfh + 194);
    const auto *sfh_198 = buffer.data(sfh + 198);
    const auto *sfh_203 = buffer.data(sfh + 203);
    const auto *sfh_204 = buffer.data(sfh + 204);
    const auto *sfh_206 = buffer.data(sfh + 206);
    const auto *sfh_207 = buffer.data(sfh + 207);
    const auto *sfh_209 = buffer.data(sfh + 209);

    const auto *sfi1_0 = buffer.data(sfi1 + 0);
    const auto *sfi1_1 = buffer.data(sfi1 + 1);
    const auto *sfi1_3 = buffer.data(sfi1 + 3);
    const auto *sfi1_5 = buffer.data(sfi1 + 5);
    const auto *sfi1_6 = buffer.data(sfi1 + 6);
    const auto *sfi1_8 = buffer.data(sfi1 + 8);
    const auto *sfi1_9 = buffer.data(sfi1 + 9);
    const auto *sfi1_10 = buffer.data(sfi1 + 10);
    const auto *sfi1_12 = buffer.data(sfi1 + 12);
    const auto *sfi1_13 = buffer.data(sfi1 + 13);
    const auto *sfi1_14 = buffer.data(sfi1 + 14);
    const auto *sfi1_21 = buffer.data(sfi1 + 21);
    const auto *sfi1_27 = buffer.data(sfi1 + 27);
    const auto *sfi1_56 = buffer.data(sfi1 + 56);
    const auto *sfi1_61 = buffer.data(sfi1 + 61);
    const auto *sfi1_64 = buffer.data(sfi1 + 64);
    const auto *sfi1_65 = buffer.data(sfi1 + 65);
    const auto *sfi1_68 = buffer.data(sfi1 + 68);
    const auto *sfi1_69 = buffer.data(sfi1 + 69);
    const auto *sfi1_70 = buffer.data(sfi1 + 70);
    const auto *sfi1_79 = buffer.data(sfi1 + 79);
    const auto *sfi1_80 = buffer.data(sfi1 + 80);
    const auto *sfi1_81 = buffer.data(sfi1 + 81);
    const auto *sfi1_83 = buffer.data(sfi1 + 83);
    const auto *sfi1_245 = buffer.data(sfi1 + 245);
    const auto *sfi1_247 = buffer.data(sfi1 + 247);
    const auto *sfi1_248 = buffer.data(sfi1 + 248);
    const auto *sfi1_249 = buffer.data(sfi1 + 249);
    const auto *sfi1_251 = buffer.data(sfi1 + 251);
    const auto *sfi1_252 = buffer.data(sfi1 + 252);
    const auto *sfi1_257 = buffer.data(sfi1 + 257);
    const auto *sfi1_261 = buffer.data(sfi1 + 261);
    const auto *sfi1_266 = buffer.data(sfi1 + 266);
    const auto *sfi1_273 = buffer.data(sfi1 + 273);
    const auto *sfi1_275 = buffer.data(sfi1 + 275);
    const auto *sfi1_276 = buffer.data(sfi1 + 276);
    const auto *sfi1_277 = buffer.data(sfi1 + 277);
    const auto *sfi1_279 = buffer.data(sfi1 + 279);

    const auto *ppi0_133 = buffer.data(ppi0 + 133);

    const auto *ppi1_133 = buffer.data(ppi1 + 133);

    const auto *pdi0_169 = buffer.data(pdi0 + 169);
    const auto *pdi0_171 = buffer.data(pdi0 + 171);
    const auto *pdi0_174 = buffer.data(pdi0 + 174);
    const auto *pdi0_178 = buffer.data(pdi0 + 178);
    const auto *pdi0_189 = buffer.data(pdi0 + 189);
    const auto *pdi0_217 = buffer.data(pdi0 + 217);

    const auto *pdh_99 = buffer.data(pdh + 99);
    const auto *pdh_105 = buffer.data(pdh + 105);
    const auto *pdh_108 = buffer.data(pdh + 108);
    const auto *pdh_111 = buffer.data(pdh + 111);
    const auto *pdh_115 = buffer.data(pdh + 115);
    const auto *pdh_119 = buffer.data(pdh + 119);
    const auto *pdh_120 = buffer.data(pdh + 120);
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
    const auto *pdh_155 = buffer.data(pdh + 155);
    const auto *pdh_156 = buffer.data(pdh + 156);
    const auto *pdh_157 = buffer.data(pdh + 157);
    const auto *pdh_159 = buffer.data(pdh + 159);
    const auto *pdh_160 = buffer.data(pdh + 160);
    const auto *pdh_161 = buffer.data(pdh + 161);
    const auto *pdh_162 = buffer.data(pdh + 162);
    const auto *pdh_163 = buffer.data(pdh + 163);
    const auto *pdh_164 = buffer.data(pdh + 164);
    const auto *pdh_165 = buffer.data(pdh + 165);
    const auto *pdh_166 = buffer.data(pdh + 166);
    const auto *pdh_167 = buffer.data(pdh + 167);
    const auto *pdh_183 = buffer.data(pdh + 183);
    const auto *pdh_184 = buffer.data(pdh + 184);
    const auto *pdh_185 = buffer.data(pdh + 185);
    const auto *pdh_186 = buffer.data(pdh + 186);
    const auto *pdh_187 = buffer.data(pdh + 187);
    const auto *pdh_188 = buffer.data(pdh + 188);

    const auto *pdi1_169 = buffer.data(pdi1 + 169);
    const auto *pdi1_171 = buffer.data(pdi1 + 171);
    const auto *pdi1_174 = buffer.data(pdi1 + 174);
    const auto *pdi1_178 = buffer.data(pdi1 + 178);
    const auto *pdi1_189 = buffer.data(pdi1 + 189);
    const auto *pdi1_217 = buffer.data(pdi1 + 217);

    const auto *pfg0_135 = buffer.data(pfg0 + 135);
    const auto *pfg0_136 = buffer.data(pfg0 + 136);
    const auto *pfg0_138 = buffer.data(pfg0 + 138);
    const auto *pfg0_140 = buffer.data(pfg0 + 140);
    const auto *pfg0_160 = buffer.data(pfg0 + 160);
    const auto *pfg0_161 = buffer.data(pfg0 + 161);
    const auto *pfg0_162 = buffer.data(pfg0 + 162);
    const auto *pfg0_165 = buffer.data(pfg0 + 165);
    const auto *pfg0_166 = buffer.data(pfg0 + 166);
    const auto *pfg0_168 = buffer.data(pfg0 + 168);
    const auto *pfg0_170 = buffer.data(pfg0 + 170);
    const auto *pfg0_171 = buffer.data(pfg0 + 171);
    const auto *pfg0_173 = buffer.data(pfg0 + 173);
    const auto *pfg0_174 = buffer.data(pfg0 + 174);
    const auto *pfg0_175 = buffer.data(pfg0 + 175);
    const auto *pfg0_176 = buffer.data(pfg0 + 176);
    const auto *pfg0_177 = buffer.data(pfg0 + 177);
    const auto *pfg0_178 = buffer.data(pfg0 + 178);
    const auto *pfg0_179 = buffer.data(pfg0 + 179);

    const auto *pfg1_135 = buffer.data(pfg1 + 135);
    const auto *pfg1_136 = buffer.data(pfg1 + 136);
    const auto *pfg1_138 = buffer.data(pfg1 + 138);
    const auto *pfg1_140 = buffer.data(pfg1 + 140);
    const auto *pfg1_160 = buffer.data(pfg1 + 160);
    const auto *pfg1_161 = buffer.data(pfg1 + 161);
    const auto *pfg1_162 = buffer.data(pfg1 + 162);
    const auto *pfg1_165 = buffer.data(pfg1 + 165);
    const auto *pfg1_166 = buffer.data(pfg1 + 166);
    const auto *pfg1_168 = buffer.data(pfg1 + 168);
    const auto *pfg1_170 = buffer.data(pfg1 + 170);
    const auto *pfg1_171 = buffer.data(pfg1 + 171);
    const auto *pfg1_173 = buffer.data(pfg1 + 173);
    const auto *pfg1_174 = buffer.data(pfg1 + 174);
    const auto *pfg1_175 = buffer.data(pfg1 + 175);
    const auto *pfg1_176 = buffer.data(pfg1 + 176);
    const auto *pfg1_177 = buffer.data(pfg1 + 177);
    const auto *pfg1_178 = buffer.data(pfg1 + 178);
    const auto *pfg1_179 = buffer.data(pfg1 + 179);

    const auto *pfh_182 = buffer.data(pfh + 182);
    const auto *pfh_183 = buffer.data(pfh + 183);
    const auto *pfh_188 = buffer.data(pfh + 188);
    const auto *pfh_189 = buffer.data(pfh + 189);
    const auto *pfh_190 = buffer.data(pfh + 190);
    const auto *pfh_191 = buffer.data(pfh + 191);
    const auto *pfh_192 = buffer.data(pfh + 192);
    const auto *pfh_194 = buffer.data(pfh + 194);
    const auto *pfh_195 = buffer.data(pfh + 195);
    const auto *pfh_197 = buffer.data(pfh + 197);
    const auto *pfh_198 = buffer.data(pfh + 198);
    const auto *pfh_199 = buffer.data(pfh + 199);
    const auto *pfh_203 = buffer.data(pfh + 203);
    const auto *pfh_204 = buffer.data(pfh + 204);
    const auto *pfh_206 = buffer.data(pfh + 206);
    const auto *pfh_207 = buffer.data(pfh + 207);
    const auto *pfh_209 = buffer.data(pfh + 209);
    const auto *pfh_210 = buffer.data(pfh + 210);
    const auto *pfh_211 = buffer.data(pfh + 211);
    const auto *pfh_213 = buffer.data(pfh + 213);
    const auto *pfh_216 = buffer.data(pfh + 216);
    const auto *pfh_225 = buffer.data(pfh + 225);
    const auto *pfh_226 = buffer.data(pfh + 226);
    const auto *pfh_227 = buffer.data(pfh + 227);
    const auto *pfh_228 = buffer.data(pfh + 228);
    const auto *pfh_229 = buffer.data(pfh + 229);
    const auto *pfh_230 = buffer.data(pfh + 230);
    const auto *pfh_231 = buffer.data(pfh + 231);
    const auto *pfh_232 = buffer.data(pfh + 232);
    const auto *pfh_234 = buffer.data(pfh + 234);
    const auto *pfh_236 = buffer.data(pfh + 236);
    const auto *pfh_237 = buffer.data(pfh + 237);
    const auto *pfh_239 = buffer.data(pfh + 239);
    const auto *pfh_240 = buffer.data(pfh + 240);
    const auto *pfh_241 = buffer.data(pfh + 241);
    const auto *pfh_243 = buffer.data(pfh + 243);
    const auto *pfh_244 = buffer.data(pfh + 244);
    const auto *pfh_245 = buffer.data(pfh + 245);
    const auto *pfh_246 = buffer.data(pfh + 246);
    const auto *pfh_247 = buffer.data(pfh + 247);
    const auto *pfh_248 = buffer.data(pfh + 248);
    const auto *pfh_249 = buffer.data(pfh + 249);
    const auto *pfh_250 = buffer.data(pfh + 250);
    const auto *pfh_251 = buffer.data(pfh + 251);
    const auto *pfh_252 = buffer.data(pfh + 252);
    const auto *pfh_253 = buffer.data(pfh + 253);
    const auto *pfh_255 = buffer.data(pfh + 255);
    const auto *pfh_258 = buffer.data(pfh + 258);
    const auto *pfh_267 = buffer.data(pfh + 267);
    const auto *pfh_268 = buffer.data(pfh + 268);
    const auto *pfh_269 = buffer.data(pfh + 269);
    const auto *pfh_270 = buffer.data(pfh + 270);
    const auto *pfh_271 = buffer.data(pfh + 271);
    const auto *pfh_272 = buffer.data(pfh + 272);

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_x, pc_x, pc_y, pc_z, sfi0_245, \
                         sfh_188, sfi1_245, pdh_99, pdh_119, pfh_182, pfh_183, \
                         pfh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_0 * pdh_119[k]
                   + f_4 * pc_y[k] * pfh_182[k];

        t_244[k] = f_0 * sfh_188[k]
                   + f_4 * pc_x[k] * pfh_188[k];

        t_245[k] = pa_x[k] * sfi0_245[k]
                   - f_11 * pc_x[k] * sfi1_245[k];

        t_246[k] = f_12 * pdh_99[k]
                   + f_4 * pc_z[k] * pfh_183[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_x, pc_x, pc_y, sfi0_247, sfi0_248, \
                         sfi0_249, sfi1_247, sfi1_248, sfi1_249, pdh_125, \
                         pfh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = pa_x[k] * sfi0_247[k]
                   - f_11 * pc_x[k] * sfi1_247[k];

        t_248[k] = pa_x[k] * sfi0_248[k]
                   - f_11 * pc_x[k] * sfi1_248[k];

        t_249[k] = pa_x[k] * sfi0_249[k]
                   - f_11 * pc_x[k] * sfi1_249[k];

        t_250[k] = f_0 * pdh_125[k]
                   + f_4 * pc_y[k] * pfh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_x, pc_x, pc_y, pc_z, sfi0_251, \
                         sfi0_252, sfh_189, sfi1_251, sfi1_252, pdh_105, \
                         pfh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_x[k] * sfi0_251[k]
                   - f_11 * pc_x[k] * sfi1_251[k];

        t_252[k] = pa_x[k] * sfi0_252[k]
                   + f_14 * sfh_189[k]
                   - f_11 * pc_x[k] * sfi1_252[k];

        t_253[k] = f_4 * pc_y[k] * pfh_189[k];

        t_254[k] = f_1 * pdh_105[k]
                   + f_4 * pc_z[k] * pfh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pa_x, pc_x, pc_y, sfi0_257, sfh_194, sfi1_257, \
                         pfg0_135, pfg1_135, pfh_190, pfh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_5 * pfg0_135[k]
                   - f_6 * pfg1_135[k]
                   + f_4 * pc_y[k] * pfh_190[k];

        t_256[k] = f_4 * pc_y[k] * pfh_191[k];

        t_257[k] = pa_x[k] * sfi0_257[k]
                   + f_13 * sfh_194[k]
                   - f_11 * pc_x[k] * sfi1_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pa_x, pc_x, pc_y, pc_z, sfi0_261, \
                         sfh_198, sfi1_261, pdh_108, pfg0_136, pfg1_136, pfh_192, \
                         pfh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_7 * pfg0_136[k]
                   - f_8 * pfg1_136[k]
                   + f_4 * pc_y[k] * pfh_192[k];

        t_259[k] = f_1 * pdh_108[k]
                   + f_4 * pc_z[k] * pfh_192[k];

        t_260[k] = f_4 * pc_y[k] * pfh_194[k];

        t_261[k] = pa_x[k] * sfi0_261[k]
                   + f_1 * sfh_198[k]
                   - f_11 * pc_x[k] * sfi1_261[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, pc_y, pc_z, pdh_111, pfg0_138, pfg0_140, \
                         pfg1_138, pfg1_140, pfh_195, pfh_197, \
                         pfh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_9 * pfg0_138[k]
                   - f_10 * pfg1_138[k]
                   + f_4 * pc_y[k] * pfh_195[k];

        t_263[k] = f_1 * pdh_111[k]
                   + f_4 * pc_z[k] * pfh_195[k];

        t_264[k] = f_5 * pfg0_140[k]
                   - f_6 * pfg1_140[k]
                   + f_4 * pc_y[k] * pfh_197[k];

        t_265[k] = f_4 * pc_y[k] * pfh_198[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_x, pc_x, pc_z, sfi0_266, sfh_203, \
                         sfh_204, sfh_206, sfi1_266, pdh_115, pfh_199, pfh_204, \
                         pfh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_x[k] * sfi0_266[k]
                   + f_12 * sfh_203[k]
                   - f_11 * pc_x[k] * sfi1_266[k];

        t_267[k] = f_0 * sfh_204[k]
                   + f_4 * pc_x[k] * pfh_204[k];

        t_268[k] = f_1 * pdh_115[k]
                   + f_4 * pc_z[k] * pfh_199[k];

        t_269[k] = f_0 * sfh_206[k]
                   + f_4 * pc_x[k] * pfh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pa_x, pc_x, pc_y, sfi0_273, sfh_207, \
                         sfh_209, sfi1_273, pfh_203, pfh_207, pfh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_0 * sfh_207[k]
                   + f_4 * pc_x[k] * pfh_207[k];

        t_271[k] = f_4 * pc_y[k] * pfh_203[k];

        t_272[k] = f_0 * sfh_209[k]
                   + f_4 * pc_x[k] * pfh_209[k];

        t_273[k] = pa_x[k] * sfi0_273[k]
                   - f_11 * pc_x[k] * sfi1_273[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pa_x, pc_x, pc_z, sfi0_275, sfi0_276, \
                         sfi0_277, sfi1_275, sfi1_276, sfi1_277, pdh_120, \
                         pfh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_1 * pdh_120[k]
                   + f_4 * pc_z[k] * pfh_204[k];

        t_275[k] = pa_x[k] * sfi0_275[k]
                   - f_11 * pc_x[k] * sfi1_275[k];

        t_276[k] = pa_x[k] * sfi0_276[k]
                   - f_11 * pc_x[k] * sfi1_276[k];

        t_277[k] = pa_x[k] * sfi0_277[k]
                   - f_11 * pc_x[k] * sfi1_277[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, pa_x, pa_y, pc_x, pc_y, sfi0_0, sfi0_1, \
                         sfi0_279, sfh_0, sfi1_0, sfi1_1, sfi1_279, \
                         pfh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_4 * pc_y[k] * pfh_209[k];

        t_279[k] = pa_x[k] * sfi0_279[k]
                   - f_11 * pc_x[k] * sfi1_279[k];

        t_280[k] = pa_y[k] * sfi0_0[k]
                   - f_11 * pc_y[k] * sfi1_0[k];

        t_281[k] = pa_y[k] * sfi0_1[k]
                   + f_0 * sfh_0[k]
                   - f_11 * pc_y[k] * sfi1_1[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_y, pc_y, pc_z, sfi0_3, sfi0_5, sfh_1, \
                         sfi1_3, sfi1_5, pfh_210, pfh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_4 * pc_z[k] * pfh_210[k];

        t_283[k] = pa_y[k] * sfi0_3[k]
                   + f_12 * sfh_1[k]
                   - f_11 * pc_y[k] * sfi1_3[k];

        t_284[k] = f_4 * pc_z[k] * pfh_211[k];

        t_285[k] = pa_y[k] * sfi0_5[k]
                   - f_11 * pc_y[k] * sfi1_5[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pc_y, pc_z, sfi0_6, sfi0_8, sfi0_9, \
                         sfh_3, sfh_5, sfi1_6, sfi1_8, sfi1_9, \
                         pfh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = pa_y[k] * sfi0_6[k]
                   + f_1 * sfh_3[k]
                   - f_11 * pc_y[k] * sfi1_6[k];

        t_287[k] = f_4 * pc_z[k] * pfh_213[k];

        t_288[k] = pa_y[k] * sfi0_8[k]
                   + f_0 * sfh_5[k]
                   - f_11 * pc_y[k] * sfi1_8[k];

        t_289[k] = pa_y[k] * sfi0_9[k]
                   - f_11 * pc_y[k] * sfi1_9[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pa_y, pc_y, pc_z, sfi0_10, sfi0_12, sfh_6, \
                         sfh_8, sfi1_10, sfi1_12, pfh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_y[k] * sfi0_10[k]
                   + f_13 * sfh_6[k]
                   - f_11 * pc_y[k] * sfi1_10[k];

        t_291[k] = f_4 * pc_z[k] * pfh_216[k];

        t_292[k] = pa_y[k] * sfi0_12[k]
                   + f_12 * sfh_8[k]
                   - f_11 * pc_y[k] * sfi1_12[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_y, pc_x, pc_y, sfi0_13, sfi0_14, \
                         sfh_9, sfi1_13, sfi1_14, pdh_141, pdh_142, pfh_225, \
                         pfh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_y[k] * sfi0_13[k]
                   + f_0 * sfh_9[k]
                   - f_11 * pc_y[k] * sfi1_13[k];

        t_294[k] = pa_y[k] * sfi0_14[k]
                   - f_11 * pc_y[k] * sfi1_14[k];

        t_295[k] = f_1 * pdh_141[k]
                   + f_4 * pc_x[k] * pfh_225[k];

        t_296[k] = f_1 * pdh_142[k]
                   + f_4 * pc_x[k] * pfh_226[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pc_x, pdh_143, pdh_144, pdh_145, pdh_146, \
                         pfh_227, pfh_228, pfh_229, pfh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_1 * pdh_143[k]
                   + f_4 * pc_x[k] * pfh_227[k];

        t_298[k] = f_1 * pdh_144[k]
                   + f_4 * pc_x[k] * pfh_228[k];

        t_299[k] = f_1 * pdh_145[k]
                   + f_4 * pc_x[k] * pfh_229[k];

        t_300[k] = f_1 * pdh_146[k]
                   + f_4 * pc_x[k] * pfh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pa_y, pc_y, pc_z, sfi0_21, sfh_15, sfi1_21, \
                         pfg0_160, pfg1_160, pfh_225, pfh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = pa_y[k] * sfi0_21[k]
                   + f_14 * sfh_15[k]
                   - f_11 * pc_y[k] * sfi1_21[k];

        t_302[k] = f_4 * pc_z[k] * pfh_225[k];

        t_303[k] = f_5 * pfg0_160[k]
                   - f_6 * pfg1_160[k]
                   + f_4 * pc_z[k] * pfh_226[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, pc_y, pc_z, sfh_20, pfg0_161, pfg0_162, \
                         pfg1_161, pfg1_162, pfh_227, pfh_228, \
                         pfh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_7 * pfg0_161[k]
                   - f_8 * pfg1_161[k]
                   + f_4 * pc_z[k] * pfh_227[k];

        t_305[k] = f_9 * pfg0_162[k]
                   - f_10 * pfg1_162[k]
                   + f_4 * pc_z[k] * pfh_228[k];

        t_306[k] = f_0 * sfh_20[k]
                   + f_4 * pc_y[k] * pfh_230[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, pa_y, pc_x, pc_y, sfi0_27, sfi1_27, pdh_147, \
                         pdh_148, pfg0_165, pfg0_166, pfg1_165, pfg1_166, pfh_231, \
                         pfh_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = pa_y[k] * sfi0_27[k]
                   - f_11 * pc_y[k] * sfi1_27[k];

        t_308[k] = f_12 * pdh_147[k]
                   + f_2 * pfg0_165[k]
                   - f_3 * pfg1_165[k]
                   + f_4 * pc_x[k] * pfh_231[k];

        t_309[k] = f_12 * pdh_148[k]
                   + f_15 * pfg0_166[k]
                   - f_16 * pfg1_166[k]
                   + f_4 * pc_x[k] * pfh_232[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pc_x, pc_z, pdh_150, pdh_152, pfg0_168, \
                         pfg0_170, pfg1_168, pfg1_170, pfh_231, pfh_232, pfh_234, \
                         pfh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_4 * pc_z[k] * pfh_231[k];

        t_311[k] = f_12 * pdh_150[k]
                   + f_9 * pfg0_168[k]
                   - f_10 * pfg1_168[k]
                   + f_4 * pc_x[k] * pfh_234[k];

        t_312[k] = f_4 * pc_z[k] * pfh_232[k];

        t_313[k] = f_12 * pdh_152[k]
                   + f_9 * pfg0_170[k]
                   - f_10 * pfg1_170[k]
                   + f_4 * pc_x[k] * pfh_236[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, pc_x, pc_z, pdh_153, pdh_155, pfg0_171, \
                         pfg0_173, pfg1_171, pfg1_173, pfh_234, pfh_237, \
                         pfh_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_12 * pdh_153[k]
                   + f_7 * pfg0_171[k]
                   - f_8 * pfg1_171[k]
                   + f_4 * pc_x[k] * pfh_237[k];

        t_315[k] = f_4 * pc_z[k] * pfh_234[k];

        t_316[k] = f_12 * pdh_155[k]
                   + f_7 * pfg0_173[k]
                   - f_8 * pfg1_173[k]
                   + f_4 * pc_x[k] * pfh_239[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pc_x, pc_z, pdh_156, pdh_157, pfg0_174, \
                         pfg0_175, pfg1_174, pfg1_175, pfh_237, pfh_240, \
                         pfh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_12 * pdh_156[k]
                   + f_7 * pfg0_174[k]
                   - f_8 * pfg1_174[k]
                   + f_4 * pc_x[k] * pfh_240[k];

        t_318[k] = f_12 * pdh_157[k]
                   + f_5 * pfg0_175[k]
                   - f_6 * pfg1_175[k]
                   + f_4 * pc_x[k] * pfh_241[k];

        t_319[k] = f_4 * pc_z[k] * pfh_237[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pdh_159, pdh_160, pdh_161, pfg0_177, \
                         pfg0_178, pfg0_179, pfg1_177, pfg1_178, pfg1_179, pfh_243, pfh_244, \
                         pfh_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_12 * pdh_159[k]
                   + f_5 * pfg0_177[k]
                   - f_6 * pfg1_177[k]
                   + f_4 * pc_x[k] * pfh_243[k];

        t_321[k] = f_12 * pdh_160[k]
                   + f_5 * pfg0_178[k]
                   - f_6 * pfg1_178[k]
                   + f_4 * pc_x[k] * pfh_244[k];

        t_322[k] = f_12 * pdh_161[k]
                   + f_5 * pfg0_179[k]
                   - f_6 * pfg1_179[k]
                   + f_4 * pc_x[k] * pfh_245[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, pc_x, pdh_162, pdh_163, pdh_164, \
                         pdh_165, pdh_166, pfh_246, pfh_247, pfh_248, pfh_249, \
                         pfh_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_12 * pdh_162[k]
                   + f_4 * pc_x[k] * pfh_246[k];

        t_324[k] = f_12 * pdh_163[k]
                   + f_4 * pc_x[k] * pfh_247[k];

        t_325[k] = f_12 * pdh_164[k]
                   + f_4 * pc_x[k] * pfh_248[k];

        t_326[k] = f_12 * pdh_165[k]
                   + f_4 * pc_x[k] * pfh_249[k];

        t_327[k] = f_12 * pdh_166[k]
                   + f_4 * pc_x[k] * pfh_250[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pb_x, pc_x, pc_z, ppi0_133, ppi1_133, pdi0_217, \
                         pdh_167, pdi1_217, pfh_246, pfh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_12 * pdh_167[k]
                   + f_4 * pc_x[k] * pfh_251[k];

        t_329[k] = f_17 * ppi0_133[k]
                   - f_18 * ppi1_133[k]
                   + pb_x[k] * pdi0_217[k]
                   - f_11 * pc_x[k] * pdi1_217[k];

        t_330[k] = f_4 * pc_z[k] * pfh_246[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_z, pfg0_175, pfg0_176, pfg0_177, pfg1_175, \
                         pfg1_176, pfg1_177, pfh_247, pfh_248, \
                         pfh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_5 * pfg0_175[k]
                   - f_6 * pfg1_175[k]
                   + f_4 * pc_z[k] * pfh_247[k];

        t_332[k] = f_7 * pfg0_176[k]
                   - f_8 * pfg1_176[k]
                   + f_4 * pc_z[k] * pfh_248[k];

        t_333[k] = f_9 * pfg0_177[k]
                   - f_10 * pfg1_177[k]
                   + f_4 * pc_z[k] * pfh_249[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pa_y, pc_y, pc_z, sfi0_56, sfh_41, sfi1_56, \
                         pdh_146, pfg0_179, pfg1_179, pfh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * sfh_41[k]
                   + f_0 * pdh_146[k]
                   + f_4 * pc_y[k] * pfh_251[k];

        t_335[k] = f_2 * pfg0_179[k]
                   - f_3 * pfg1_179[k]
                   + f_4 * pc_z[k] * pfh_251[k];

        t_336[k] = pa_y[k] * sfi0_56[k]
                   - f_11 * pc_y[k] * sfi1_56[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pb_z, pc_z, pdi0_169, pdi0_171, pdh_126, \
                         pdh_127, pdi1_169, pdi1_171, pfh_252, \
                         pfh_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = pb_z[k] * pdi0_169[k]
                   - f_11 * pc_z[k] * pdi1_169[k];

        t_338[k] = f_0 * pdh_126[k]
                   + f_4 * pc_z[k] * pfh_252[k];

        t_339[k] = pb_z[k] * pdi0_171[k]
                   - f_11 * pc_z[k] * pdi1_171[k];

        t_340[k] = f_0 * pdh_127[k]
                   + f_4 * pc_z[k] * pfh_253[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pa_y, pb_z, pc_y, pc_z, sfi0_61, sfi1_61, \
                         pdi0_174, pdh_129, pdi1_174, pfh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pa_y[k] * sfi0_61[k]
                   - f_11 * pc_y[k] * sfi1_61[k];

        t_342[k] = pb_z[k] * pdi0_174[k]
                   - f_11 * pc_z[k] * pdi1_174[k];

        t_343[k] = f_0 * pdh_129[k]
                   + f_4 * pc_z[k] * pfh_255[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, pa_y, pb_z, pc_y, pc_z, sfi0_64, sfi0_65, \
                         sfh_47, sfi1_64, sfi1_65, pdi0_178, pdi1_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pa_y[k] * sfi0_64[k]
                   + f_0 * sfh_47[k]
                   - f_11 * pc_y[k] * sfi1_64[k];

        t_345[k] = pa_y[k] * sfi0_65[k]
                   - f_11 * pc_y[k] * sfi1_65[k];

        t_346[k] = pb_z[k] * pdi0_178[k]
                   - f_11 * pc_z[k] * pdi1_178[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pa_y, pc_y, pc_z, sfi0_68, sfi0_69, sfh_50, \
                         sfh_51, sfi1_68, sfi1_69, pdh_132, pfh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_0 * pdh_132[k]
                   + f_4 * pc_z[k] * pfh_258[k];

        t_348[k] = pa_y[k] * sfi0_68[k]
                   + f_12 * sfh_50[k]
                   - f_11 * pc_y[k] * sfi1_68[k];

        t_349[k] = pa_y[k] * sfi0_69[k]
                   + f_0 * sfh_51[k]
                   - f_11 * pc_y[k] * sfi1_69[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_y, pc_x, pc_y, sfi0_70, sfi1_70, \
                         pdh_183, pdh_184, pdh_185, pfh_267, pfh_268, \
                         pfh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pa_y[k] * sfi0_70[k]
                   - f_11 * pc_y[k] * sfi1_70[k];

        t_351[k] = f_12 * pdh_183[k]
                   + f_4 * pc_x[k] * pfh_267[k];

        t_352[k] = f_12 * pdh_184[k]
                   + f_4 * pc_x[k] * pfh_268[k];

        t_353[k] = f_12 * pdh_185[k]
                   + f_4 * pc_x[k] * pfh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pb_z, pc_x, pc_z, pdi0_189, pdh_186, \
                         pdh_187, pdh_188, pdi1_189, pfh_270, pfh_271, \
                         pfh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_12 * pdh_186[k]
                   + f_4 * pc_x[k] * pfh_270[k];

        t_355[k] = f_12 * pdh_187[k]
                   + f_4 * pc_x[k] * pfh_271[k];

        t_356[k] = f_12 * pdh_188[k]
                   + f_4 * pc_x[k] * pfh_272[k];

        t_357[k] = pb_z[k] * pdi0_189[k]
                   - f_11 * pc_z[k] * pdi1_189[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, pa_y, pc_y, pc_z, sfi0_79, sfi0_80, sfh_59, \
                         sfh_60, sfi1_79, sfi1_80, pdh_141, pfh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_0 * pdh_141[k]
                   + f_4 * pc_z[k] * pfh_267[k];

        t_359[k] = pa_y[k] * sfi0_79[k]
                   + f_13 * sfh_59[k]
                   - f_11 * pc_y[k] * sfi1_79[k];

        t_360[k] = pa_y[k] * sfi0_80[k]
                   + f_1 * sfh_60[k]
                   - f_11 * pc_y[k] * sfi1_80[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pa_y, pc_y, sfi0_81, sfi0_83, sfh_61, sfh_62, \
                         sfi1_81, sfi1_83, pfh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = pa_y[k] * sfi0_81[k]
                   + f_12 * sfh_61[k]
                   - f_11 * pc_y[k] * sfi1_81[k];

        t_362[k] = f_0 * sfh_62[k]
                   + f_4 * pc_y[k] * pfh_272[k];

        t_363[k] = pa_y[k] * sfi0_83[k]
                   - f_11 * pc_y[k] * sfi1_83[k];
    }
}

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfi0, const size_t sfh,
                                                          const size_t sfi1, const size_t pdi0,
                                                          const size_t pdh, const size_t pdi1,
                                                          const size_t pfg0, const size_t pfg1,
                                                          const size_t pfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_140 = buffer.data(sfi0 + 140);
    const auto *sfi0_141 = buffer.data(sfi0 + 141);
    const auto *sfi0_143 = buffer.data(sfi0 + 143);
    const auto *sfi0_145 = buffer.data(sfi0 + 145);
    const auto *sfi0_146 = buffer.data(sfi0 + 146);
    const auto *sfi0_148 = buffer.data(sfi0 + 148);
    const auto *sfi0_149 = buffer.data(sfi0 + 149);
    const auto *sfi0_150 = buffer.data(sfi0 + 150);
    const auto *sfi0_152 = buffer.data(sfi0 + 152);
    const auto *sfi0_153 = buffer.data(sfi0 + 153);
    const auto *sfi0_154 = buffer.data(sfi0 + 154);
    const auto *sfi0_167 = buffer.data(sfi0 + 167);

    const auto *sfh_105 = buffer.data(sfh + 105);
    const auto *sfh_106 = buffer.data(sfh + 106);
    const auto *sfh_108 = buffer.data(sfh + 108);
    const auto *sfh_110 = buffer.data(sfh + 110);
    const auto *sfh_111 = buffer.data(sfh + 111);
    const auto *sfh_113 = buffer.data(sfh + 113);
    const auto *sfh_114 = buffer.data(sfh + 114);
    const auto *sfh_125 = buffer.data(sfh + 125);
    const auto *sfh_141 = buffer.data(sfh + 141);
    const auto *sfh_146 = buffer.data(sfh + 146);

    const auto *sfi1_140 = buffer.data(sfi1 + 140);
    const auto *sfi1_141 = buffer.data(sfi1 + 141);
    const auto *sfi1_143 = buffer.data(sfi1 + 143);
    const auto *sfi1_145 = buffer.data(sfi1 + 145);
    const auto *sfi1_146 = buffer.data(sfi1 + 146);
    const auto *sfi1_148 = buffer.data(sfi1 + 148);
    const auto *sfi1_149 = buffer.data(sfi1 + 149);
    const auto *sfi1_150 = buffer.data(sfi1 + 150);
    const auto *sfi1_152 = buffer.data(sfi1 + 152);
    const auto *sfi1_153 = buffer.data(sfi1 + 153);
    const auto *sfi1_154 = buffer.data(sfi1 + 154);
    const auto *sfi1_167 = buffer.data(sfi1 + 167);

    const auto *pdi0_197 = buffer.data(pdi0 + 197);
    const auto *pdi0_199 = buffer.data(pdi0 + 199);
    const auto *pdi0_202 = buffer.data(pdi0 + 202);
    const auto *pdi0_206 = buffer.data(pdi0 + 206);
    const auto *pdi0_252 = buffer.data(pdi0 + 252);
    const auto *pdi0_253 = buffer.data(pdi0 + 253);
    const auto *pdi0_255 = buffer.data(pdi0 + 255);
    const auto *pdi0_257 = buffer.data(pdi0 + 257);
    const auto *pdi0_258 = buffer.data(pdi0 + 258);
    const auto *pdi0_260 = buffer.data(pdi0 + 260);
    const auto *pdi0_261 = buffer.data(pdi0 + 261);
    const auto *pdi0_262 = buffer.data(pdi0 + 262);
    const auto *pdi0_264 = buffer.data(pdi0 + 264);
    const auto *pdi0_265 = buffer.data(pdi0 + 265);
    const auto *pdi0_266 = buffer.data(pdi0 + 266);
    const auto *pdi0_273 = buffer.data(pdi0 + 273);
    const auto *pdi0_275 = buffer.data(pdi0 + 275);
    const auto *pdi0_276 = buffer.data(pdi0 + 276);
    const auto *pdi0_277 = buffer.data(pdi0 + 277);
    const auto *pdi0_278 = buffer.data(pdi0 + 278);
    const auto *pdi0_279 = buffer.data(pdi0 + 279);
    const auto *pdi0_288 = buffer.data(pdi0 + 288);
    const auto *pdi0_292 = buffer.data(pdi0 + 292);
    const auto *pdi0_293 = buffer.data(pdi0 + 293);
    const auto *pdi0_301 = buffer.data(pdi0 + 301);
    const auto *pdi0_303 = buffer.data(pdi0 + 303);
    const auto *pdi0_304 = buffer.data(pdi0 + 304);
    const auto *pdi0_305 = buffer.data(pdi0 + 305);
    const auto *pdi0_306 = buffer.data(pdi0 + 306);
    const auto *pdi0_307 = buffer.data(pdi0 + 307);
    const auto *pdi0_329 = buffer.data(pdi0 + 329);
    const auto *pdi0_331 = buffer.data(pdi0 + 331);
    const auto *pdi0_332 = buffer.data(pdi0 + 332);
    const auto *pdi0_333 = buffer.data(pdi0 + 333);

    const auto *pdh_147 = buffer.data(pdh + 147);
    const auto *pdh_148 = buffer.data(pdh + 148);
    const auto *pdh_150 = buffer.data(pdh + 150);
    const auto *pdh_153 = buffer.data(pdh + 153);
    const auto *pdh_162 = buffer.data(pdh + 162);
    const auto *pdh_168 = buffer.data(pdh + 168);
    const auto *pdh_169 = buffer.data(pdh + 169);
    const auto *pdh_171 = buffer.data(pdh + 171);
    const auto *pdh_174 = buffer.data(pdh + 174);
    const auto *pdh_183 = buffer.data(pdh + 183);
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
    const auto *pdh_215 = buffer.data(pdh + 215);
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
    const auto *pdh_246 = buffer.data(pdh + 246);
    const auto *pdh_247 = buffer.data(pdh + 247);
    const auto *pdh_248 = buffer.data(pdh + 248);
    const auto *pdh_249 = buffer.data(pdh + 249);
    const auto *pdh_250 = buffer.data(pdh + 250);
    const auto *pdh_251 = buffer.data(pdh + 251);

    const auto *pdi1_197 = buffer.data(pdi1 + 197);
    const auto *pdi1_199 = buffer.data(pdi1 + 199);
    const auto *pdi1_202 = buffer.data(pdi1 + 202);
    const auto *pdi1_206 = buffer.data(pdi1 + 206);
    const auto *pdi1_252 = buffer.data(pdi1 + 252);
    const auto *pdi1_253 = buffer.data(pdi1 + 253);
    const auto *pdi1_255 = buffer.data(pdi1 + 255);
    const auto *pdi1_257 = buffer.data(pdi1 + 257);
    const auto *pdi1_258 = buffer.data(pdi1 + 258);
    const auto *pdi1_260 = buffer.data(pdi1 + 260);
    const auto *pdi1_261 = buffer.data(pdi1 + 261);
    const auto *pdi1_262 = buffer.data(pdi1 + 262);
    const auto *pdi1_264 = buffer.data(pdi1 + 264);
    const auto *pdi1_265 = buffer.data(pdi1 + 265);
    const auto *pdi1_266 = buffer.data(pdi1 + 266);
    const auto *pdi1_273 = buffer.data(pdi1 + 273);
    const auto *pdi1_275 = buffer.data(pdi1 + 275);
    const auto *pdi1_276 = buffer.data(pdi1 + 276);
    const auto *pdi1_277 = buffer.data(pdi1 + 277);
    const auto *pdi1_278 = buffer.data(pdi1 + 278);
    const auto *pdi1_279 = buffer.data(pdi1 + 279);
    const auto *pdi1_288 = buffer.data(pdi1 + 288);
    const auto *pdi1_292 = buffer.data(pdi1 + 292);
    const auto *pdi1_293 = buffer.data(pdi1 + 293);
    const auto *pdi1_301 = buffer.data(pdi1 + 301);
    const auto *pdi1_303 = buffer.data(pdi1 + 303);
    const auto *pdi1_304 = buffer.data(pdi1 + 304);
    const auto *pdi1_305 = buffer.data(pdi1 + 305);
    const auto *pdi1_306 = buffer.data(pdi1 + 306);
    const auto *pdi1_307 = buffer.data(pdi1 + 307);
    const auto *pdi1_329 = buffer.data(pdi1 + 329);
    const auto *pdi1_331 = buffer.data(pdi1 + 331);
    const auto *pdi1_332 = buffer.data(pdi1 + 332);
    const auto *pdi1_333 = buffer.data(pdi1 + 333);

    const auto *pfg0_210 = buffer.data(pfg0 + 210);
    const auto *pfg0_215 = buffer.data(pfg0 + 215);
    const auto *pfg0_219 = buffer.data(pfg0 + 219);
    const auto *pfg0_224 = buffer.data(pfg0 + 224);
    const auto *pfg0_240 = buffer.data(pfg0 + 240);
    const auto *pfg0_241 = buffer.data(pfg0 + 241);
    const auto *pfg0_243 = buffer.data(pfg0 + 243);
    const auto *pfg0_245 = buffer.data(pfg0 + 245);
    const auto *pfg0_246 = buffer.data(pfg0 + 246);
    const auto *pfg0_248 = buffer.data(pfg0 + 248);
    const auto *pfg0_249 = buffer.data(pfg0 + 249);
    const auto *pfg0_250 = buffer.data(pfg0 + 250);
    const auto *pfg0_251 = buffer.data(pfg0 + 251);
    const auto *pfg0_252 = buffer.data(pfg0 + 252);
    const auto *pfg0_253 = buffer.data(pfg0 + 253);
    const auto *pfg0_254 = buffer.data(pfg0 + 254);
    const auto *pfg0_260 = buffer.data(pfg0 + 260);
    const auto *pfg0_263 = buffer.data(pfg0 + 263);
    const auto *pfg0_264 = buffer.data(pfg0 + 264);

    const auto *pfg1_210 = buffer.data(pfg1 + 210);
    const auto *pfg1_215 = buffer.data(pfg1 + 215);
    const auto *pfg1_219 = buffer.data(pfg1 + 219);
    const auto *pfg1_224 = buffer.data(pfg1 + 224);
    const auto *pfg1_240 = buffer.data(pfg1 + 240);
    const auto *pfg1_241 = buffer.data(pfg1 + 241);
    const auto *pfg1_243 = buffer.data(pfg1 + 243);
    const auto *pfg1_245 = buffer.data(pfg1 + 245);
    const auto *pfg1_246 = buffer.data(pfg1 + 246);
    const auto *pfg1_248 = buffer.data(pfg1 + 248);
    const auto *pfg1_249 = buffer.data(pfg1 + 249);
    const auto *pfg1_250 = buffer.data(pfg1 + 250);
    const auto *pfg1_251 = buffer.data(pfg1 + 251);
    const auto *pfg1_252 = buffer.data(pfg1 + 252);
    const auto *pfg1_253 = buffer.data(pfg1 + 253);
    const auto *pfg1_254 = buffer.data(pfg1 + 254);
    const auto *pfg1_260 = buffer.data(pfg1 + 260);
    const auto *pfg1_263 = buffer.data(pfg1 + 263);
    const auto *pfg1_264 = buffer.data(pfg1 + 264);

    const auto *pfh_273 = buffer.data(pfh + 273);
    const auto *pfh_274 = buffer.data(pfh + 274);
    const auto *pfh_276 = buffer.data(pfh + 276);
    const auto *pfh_279 = buffer.data(pfh + 279);
    const auto *pfh_288 = buffer.data(pfh + 288);
    const auto *pfh_289 = buffer.data(pfh + 289);
    const auto *pfh_290 = buffer.data(pfh + 290);
    const auto *pfh_291 = buffer.data(pfh + 291);
    const auto *pfh_292 = buffer.data(pfh + 292);
    const auto *pfh_293 = buffer.data(pfh + 293);
    const auto *pfh_294 = buffer.data(pfh + 294);
    const auto *pfh_295 = buffer.data(pfh + 295);
    const auto *pfh_297 = buffer.data(pfh + 297);
    const auto *pfh_299 = buffer.data(pfh + 299);
    const auto *pfh_300 = buffer.data(pfh + 300);
    const auto *pfh_303 = buffer.data(pfh + 303);
    const auto *pfh_308 = buffer.data(pfh + 308);
    const auto *pfh_309 = buffer.data(pfh + 309);
    const auto *pfh_310 = buffer.data(pfh + 310);
    const auto *pfh_311 = buffer.data(pfh + 311);
    const auto *pfh_312 = buffer.data(pfh + 312);
    const auto *pfh_313 = buffer.data(pfh + 313);
    const auto *pfh_314 = buffer.data(pfh + 314);
    const auto *pfh_315 = buffer.data(pfh + 315);
    const auto *pfh_316 = buffer.data(pfh + 316);
    const auto *pfh_318 = buffer.data(pfh + 318);
    const auto *pfh_321 = buffer.data(pfh + 321);
    const auto *pfh_330 = buffer.data(pfh + 330);
    const auto *pfh_331 = buffer.data(pfh + 331);
    const auto *pfh_332 = buffer.data(pfh + 332);
    const auto *pfh_333 = buffer.data(pfh + 333);
    const auto *pfh_334 = buffer.data(pfh + 334);
    const auto *pfh_335 = buffer.data(pfh + 335);
    const auto *pfh_336 = buffer.data(pfh + 336);
    const auto *pfh_337 = buffer.data(pfh + 337);
    const auto *pfh_339 = buffer.data(pfh + 339);
    const auto *pfh_341 = buffer.data(pfh + 341);
    const auto *pfh_342 = buffer.data(pfh + 342);
    const auto *pfh_344 = buffer.data(pfh + 344);
    const auto *pfh_345 = buffer.data(pfh + 345);
    const auto *pfh_346 = buffer.data(pfh + 346);
    const auto *pfh_348 = buffer.data(pfh + 348);
    const auto *pfh_349 = buffer.data(pfh + 349);
    const auto *pfh_350 = buffer.data(pfh + 350);
    const auto *pfh_351 = buffer.data(pfh + 351);
    const auto *pfh_352 = buffer.data(pfh + 352);
    const auto *pfh_353 = buffer.data(pfh + 353);
    const auto *pfh_354 = buffer.data(pfh + 354);
    const auto *pfh_355 = buffer.data(pfh + 355);
    const auto *pfh_356 = buffer.data(pfh + 356);
    const auto *pfh_357 = buffer.data(pfh + 357);
    const auto *pfh_358 = buffer.data(pfh + 358);
    const auto *pfh_360 = buffer.data(pfh + 360);
    const auto *pfh_362 = buffer.data(pfh + 362);
    const auto *pfh_363 = buffer.data(pfh + 363);
    const auto *pfh_365 = buffer.data(pfh + 365);
    const auto *pfh_366 = buffer.data(pfh + 366);

#pragma omp simd aligned(t_364, t_365, t_366, pb_x, pc_x, pc_z, pdi0_252, pdi0_253, pdh_189, \
                         pdh_190, pdi1_252, pdi1_253, pfh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = pb_x[k] * pdi0_252[k]
                   + f_14 * pdh_189[k]
                   - f_11 * pc_x[k] * pdi1_252[k];

        t_365[k] = pb_x[k] * pdi0_253[k]
                   + f_19 * pdh_190[k]
                   - f_11 * pc_x[k] * pdi1_253[k];

        t_366[k] = f_4 * pc_z[k] * pfh_273[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pb_x, pc_x, pc_z, pdi0_255, pdi0_257, pdh_192, \
                         pdh_194, pdi1_255, pdi1_257, pfh_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pb_x[k] * pdi0_255[k]
                   + f_13 * pdh_192[k]
                   - f_11 * pc_x[k] * pdi1_255[k];

        t_368[k] = f_4 * pc_z[k] * pfh_274[k];

        t_369[k] = pb_x[k] * pdi0_257[k]
                   + f_13 * pdh_194[k]
                   - f_11 * pc_x[k] * pdi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pb_x, pc_x, pc_z, pdi0_258, pdi0_260, pdh_195, \
                         pdh_197, pdi1_258, pdi1_260, pfh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_x[k] * pdi0_258[k]
                   + f_1 * pdh_195[k]
                   - f_11 * pc_x[k] * pdi1_258[k];

        t_371[k] = f_4 * pc_z[k] * pfh_276[k];

        t_372[k] = pb_x[k] * pdi0_260[k]
                   + f_1 * pdh_197[k]
                   - f_11 * pc_x[k] * pdi1_260[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pb_x, pc_x, pc_z, pdi0_261, pdi0_262, pdh_198, \
                         pdh_199, pdi1_261, pdi1_262, pfh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = pb_x[k] * pdi0_261[k]
                   + f_1 * pdh_198[k]
                   - f_11 * pc_x[k] * pdi1_261[k];

        t_374[k] = pb_x[k] * pdi0_262[k]
                   + f_12 * pdh_199[k]
                   - f_11 * pc_x[k] * pdi1_262[k];

        t_375[k] = f_4 * pc_z[k] * pfh_279[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pb_x, pc_x, pdi0_264, pdi0_265, pdi0_266, \
                         pdh_201, pdh_202, pdh_203, pdi1_264, pdi1_265, \
                         pdi1_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pb_x[k] * pdi0_264[k]
                   + f_12 * pdh_201[k]
                   - f_11 * pc_x[k] * pdi1_264[k];

        t_377[k] = pb_x[k] * pdi0_265[k]
                   + f_12 * pdh_202[k]
                   - f_11 * pc_x[k] * pdi1_265[k];

        t_378[k] = pb_x[k] * pdi0_266[k]
                   + f_12 * pdh_203[k]
                   - f_11 * pc_x[k] * pdi1_266[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, pc_x, pdh_204, pdh_205, pdh_206, \
                         pdh_207, pdh_208, pfh_288, pfh_289, pfh_290, pfh_291, \
                         pfh_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_0 * pdh_204[k]
                   + f_4 * pc_x[k] * pfh_288[k];

        t_380[k] = f_0 * pdh_205[k]
                   + f_4 * pc_x[k] * pfh_289[k];

        t_381[k] = f_0 * pdh_206[k]
                   + f_4 * pc_x[k] * pfh_290[k];

        t_382[k] = f_0 * pdh_207[k]
                   + f_4 * pc_x[k] * pfh_291[k];

        t_383[k] = f_0 * pdh_208[k]
                   + f_4 * pc_x[k] * pfh_292[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pb_x, pc_x, pc_z, pdi0_273, pdi0_275, \
                         pdh_209, pdi1_273, pdi1_275, pfh_288, \
                         pfh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_0 * pdh_209[k]
                   + f_4 * pc_x[k] * pfh_293[k];

        t_385[k] = pb_x[k] * pdi0_273[k]
                   - f_11 * pc_x[k] * pdi1_273[k];

        t_386[k] = f_4 * pc_z[k] * pfh_288[k];

        t_387[k] = pb_x[k] * pdi0_275[k]
                   - f_11 * pc_x[k] * pdi1_275[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_x, pc_x, pdi0_276, pdi0_277, pdi0_278, \
                         pdi0_279, pdi1_276, pdi1_277, pdi1_278, \
                         pdi1_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pb_x[k] * pdi0_276[k]
                   - f_11 * pc_x[k] * pdi1_276[k];

        t_389[k] = pb_x[k] * pdi0_277[k]
                   - f_11 * pc_x[k] * pdi1_277[k];

        t_390[k] = pb_x[k] * pdi0_278[k]
                   - f_11 * pc_x[k] * pdi1_278[k];

        t_391[k] = pb_x[k] * pdi0_279[k]
                   - f_11 * pc_x[k] * pdi1_279[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pb_z, pc_x, pc_z, pdi0_197, pdi0_199, \
                         pdh_147, pdh_210, pdi1_197, pdi1_199, pfg0_210, pfg1_210, \
                         pfh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_0 * pdh_210[k]
                   + f_2 * pfg0_210[k]
                   - f_3 * pfg1_210[k]
                   + f_4 * pc_x[k] * pfh_294[k];

        t_393[k] = pb_z[k] * pdi0_197[k]
                   - f_11 * pc_z[k] * pdi1_197[k];

        t_394[k] = f_0 * pdh_147[k]
                   + f_4 * pc_z[k] * pfh_294[k];

        t_395[k] = pb_z[k] * pdi0_199[k]
                   - f_11 * pc_z[k] * pdi1_199[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, pb_z, pc_x, pc_z, pdi0_202, pdh_148, pdh_215, \
                         pdi1_202, pfg0_215, pfg1_215, pfh_295, \
                         pfh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_0 * pdh_148[k]
                   + f_4 * pc_z[k] * pfh_295[k];

        t_397[k] = f_0 * pdh_215[k]
                   + f_9 * pfg0_215[k]
                   - f_10 * pfg1_215[k]
                   + f_4 * pc_x[k] * pfh_299[k];

        t_398[k] = pb_z[k] * pdi0_202[k]
                   - f_11 * pc_z[k] * pdi1_202[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pb_x, pc_x, pc_z, pdi0_288, pdh_150, pdh_218, \
                         pdh_219, pdi1_288, pfg0_219, pfg1_219, pfh_297, \
                         pfh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_0 * pdh_150[k]
                   + f_4 * pc_z[k] * pfh_297[k];

        t_400[k] = pb_x[k] * pdi0_288[k]
                   + f_1 * pdh_218[k]
                   - f_11 * pc_x[k] * pdi1_288[k];

        t_401[k] = f_0 * pdh_219[k]
                   + f_7 * pfg0_219[k]
                   - f_8 * pfg1_219[k]
                   + f_4 * pc_x[k] * pfh_303[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pb_x, pb_z, pc_x, pc_z, pdi0_206, pdi0_292, \
                         pdh_153, pdh_222, pdi1_206, pdi1_292, \
                         pfh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pb_z[k] * pdi0_206[k]
                   - f_11 * pc_z[k] * pdi1_206[k];

        t_403[k] = f_0 * pdh_153[k]
                   + f_4 * pc_z[k] * pfh_300[k];

        t_404[k] = pb_x[k] * pdi0_292[k]
                   + f_12 * pdh_222[k]
                   - f_11 * pc_x[k] * pdi1_292[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, pb_x, pc_x, pdi0_293, pdh_223, pdh_224, pdh_225, \
                         pdi1_293, pfg0_224, pfg1_224, pfh_308, \
                         pfh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pb_x[k] * pdi0_293[k]
                   + f_12 * pdh_223[k]
                   - f_11 * pc_x[k] * pdi1_293[k];

        t_406[k] = f_0 * pdh_224[k]
                   + f_5 * pfg0_224[k]
                   - f_6 * pfg1_224[k]
                   + f_4 * pc_x[k] * pfh_308[k];

        t_407[k] = f_0 * pdh_225[k]
                   + f_4 * pc_x[k] * pfh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pdh_226, pdh_227, pdh_228, \
                         pdh_229, pdh_230, pfh_310, pfh_311, pfh_312, pfh_313, \
                         pfh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_0 * pdh_226[k]
                   + f_4 * pc_x[k] * pfh_310[k];

        t_409[k] = f_0 * pdh_227[k]
                   + f_4 * pc_x[k] * pfh_311[k];

        t_410[k] = f_0 * pdh_228[k]
                   + f_4 * pc_x[k] * pfh_312[k];

        t_411[k] = f_0 * pdh_229[k]
                   + f_4 * pc_x[k] * pfh_313[k];

        t_412[k] = f_0 * pdh_230[k]
                   + f_4 * pc_x[k] * pfh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, pb_x, pc_x, pc_z, pdi0_301, pdi0_303, \
                         pdi0_304, pdh_162, pdi1_301, pdi1_303, pdi1_304, \
                         pfh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_x[k] * pdi0_301[k]
                   - f_11 * pc_x[k] * pdi1_301[k];

        t_414[k] = f_0 * pdh_162[k]
                   + f_4 * pc_z[k] * pfh_309[k];

        t_415[k] = pb_x[k] * pdi0_303[k]
                   - f_11 * pc_x[k] * pdi1_303[k];

        t_416[k] = pb_x[k] * pdi0_304[k]
                   - f_11 * pc_x[k] * pdi1_304[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pa_y, pb_x, pc_x, pc_y, sfi0_140, \
                         sfi1_140, pdi0_305, pdi0_306, pdi0_307, pdi1_305, pdi1_306, \
                         pdi1_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = pb_x[k] * pdi0_305[k]
                   - f_11 * pc_x[k] * pdi1_305[k];

        t_418[k] = pb_x[k] * pdi0_306[k]
                   - f_11 * pc_x[k] * pdi1_306[k];

        t_419[k] = pb_x[k] * pdi0_307[k]
                   - f_11 * pc_x[k] * pdi1_307[k];

        t_420[k] = pa_y[k] * sfi0_140[k]
                   - f_11 * pc_y[k] * sfi1_140[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, pa_y, pc_y, pc_z, sfi0_141, sfi0_143, sfh_105, \
                         sfh_106, sfi1_141, sfi1_143, pdh_168, \
                         pfh_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = pa_y[k] * sfi0_141[k]
                   + f_0 * sfh_105[k]
                   - f_11 * pc_y[k] * sfi1_141[k];

        t_422[k] = f_12 * pdh_168[k]
                   + f_4 * pc_z[k] * pfh_315[k];

        t_423[k] = pa_y[k] * sfi0_143[k]
                   + f_12 * sfh_106[k]
                   - f_11 * pc_y[k] * sfi1_143[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pa_y, pc_y, pc_z, sfi0_145, sfi0_146, \
                         sfh_108, sfi1_145, sfi1_146, pdh_169, pdh_171, pfh_316, \
                         pfh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_12 * pdh_169[k]
                   + f_4 * pc_z[k] * pfh_316[k];

        t_425[k] = pa_y[k] * sfi0_145[k]
                   - f_11 * pc_y[k] * sfi1_145[k];

        t_426[k] = pa_y[k] * sfi0_146[k]
                   + f_1 * sfh_108[k]
                   - f_11 * pc_y[k] * sfi1_146[k];

        t_427[k] = f_12 * pdh_171[k]
                   + f_4 * pc_z[k] * pfh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, pa_y, pc_y, sfi0_148, sfi0_149, sfi0_150, \
                         sfh_110, sfh_111, sfi1_148, sfi1_149, \
                         sfi1_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = pa_y[k] * sfi0_148[k]
                   + f_0 * sfh_110[k]
                   - f_11 * pc_y[k] * sfi1_148[k];

        t_429[k] = pa_y[k] * sfi0_149[k]
                   - f_11 * pc_y[k] * sfi1_149[k];

        t_430[k] = pa_y[k] * sfi0_150[k]
                   + f_13 * sfh_111[k]
                   - f_11 * pc_y[k] * sfi1_150[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, pa_y, pc_y, pc_z, sfi0_152, sfi0_153, sfh_113, \
                         sfh_114, sfi1_152, sfi1_153, pdh_174, \
                         pfh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_12 * pdh_174[k]
                   + f_4 * pc_z[k] * pfh_321[k];

        t_432[k] = pa_y[k] * sfi0_152[k]
                   + f_12 * sfh_113[k]
                   - f_11 * pc_y[k] * sfi1_152[k];

        t_433[k] = pa_y[k] * sfi0_153[k]
                   + f_0 * sfh_114[k]
                   - f_11 * pc_y[k] * sfi1_153[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_y, pc_x, pc_y, sfi0_154, sfi1_154, \
                         pdh_246, pdh_247, pdh_248, pfh_330, pfh_331, \
                         pfh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = pa_y[k] * sfi0_154[k]
                   - f_11 * pc_y[k] * sfi1_154[k];

        t_435[k] = f_0 * pdh_246[k]
                   + f_4 * pc_x[k] * pfh_330[k];

        t_436[k] = f_0 * pdh_247[k]
                   + f_4 * pc_x[k] * pfh_331[k];

        t_437[k] = f_0 * pdh_248[k]
                   + f_4 * pc_x[k] * pfh_332[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, pb_x, pc_x, pdi0_329, pdh_249, pdh_250, \
                         pdh_251, pdi1_329, pfh_333, pfh_334, pfh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_0 * pdh_249[k]
                   + f_4 * pc_x[k] * pfh_333[k];

        t_439[k] = f_0 * pdh_250[k]
                   + f_4 * pc_x[k] * pfh_334[k];

        t_440[k] = f_0 * pdh_251[k]
                   + f_4 * pc_x[k] * pfh_335[k];

        t_441[k] = pb_x[k] * pdi0_329[k]
                   - f_11 * pc_x[k] * pdi1_329[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pb_x, pc_x, pc_z, pdi0_331, pdi0_332, \
                         pdi0_333, pdh_183, pdi1_331, pdi1_332, pdi1_333, \
                         pfh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_12 * pdh_183[k]
                   + f_4 * pc_z[k] * pfh_330[k];

        t_443[k] = pb_x[k] * pdi0_331[k]
                   - f_11 * pc_x[k] * pdi1_331[k];

        t_444[k] = pb_x[k] * pdi0_332[k]
                   - f_11 * pc_x[k] * pdi1_332[k];

        t_445[k] = pb_x[k] * pdi0_333[k]
                   - f_11 * pc_x[k] * pdi1_333[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, pa_y, pc_x, pc_y, sfi0_167, sfh_125, sfi1_167, \
                         pfg0_240, pfg1_240, pfh_335, pfh_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_0 * sfh_125[k]
                   + f_4 * pc_y[k] * pfh_335[k];

        t_447[k] = pa_y[k] * sfi0_167[k]
                   - f_11 * pc_y[k] * sfi1_167[k];

        t_448[k] = f_2 * pfg0_240[k]
                   - f_3 * pfg1_240[k]
                   + f_4 * pc_x[k] * pfh_336[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_z, pfg0_241, pfg0_243, pfg1_241, \
                         pfg1_243, pfh_336, pfh_337, pfh_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_15 * pfg0_241[k]
                   - f_16 * pfg1_241[k]
                   + f_4 * pc_x[k] * pfh_337[k];

        t_450[k] = f_4 * pc_z[k] * pfh_336[k];

        t_451[k] = f_9 * pfg0_243[k]
                   - f_10 * pfg1_243[k]
                   + f_4 * pc_x[k] * pfh_339[k];

        t_452[k] = f_4 * pc_z[k] * pfh_337[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_z, pfg0_245, pfg0_246, pfg0_248, \
                         pfg1_245, pfg1_246, pfg1_248, pfh_339, pfh_341, pfh_342, \
                         pfh_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_9 * pfg0_245[k]
                   - f_10 * pfg1_245[k]
                   + f_4 * pc_x[k] * pfh_341[k];

        t_454[k] = f_7 * pfg0_246[k]
                   - f_8 * pfg1_246[k]
                   + f_4 * pc_x[k] * pfh_342[k];

        t_455[k] = f_4 * pc_z[k] * pfh_339[k];

        t_456[k] = f_7 * pfg0_248[k]
                   - f_8 * pfg1_248[k]
                   + f_4 * pc_x[k] * pfh_344[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_x, pc_z, pfg0_249, pfg0_250, pfg0_252, \
                         pfg1_249, pfg1_250, pfg1_252, pfh_342, pfh_345, pfh_346, \
                         pfh_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_7 * pfg0_249[k]
                   - f_8 * pfg1_249[k]
                   + f_4 * pc_x[k] * pfh_345[k];

        t_458[k] = f_5 * pfg0_250[k]
                   - f_6 * pfg1_250[k]
                   + f_4 * pc_x[k] * pfh_346[k];

        t_459[k] = f_4 * pc_z[k] * pfh_342[k];

        t_460[k] = f_5 * pfg0_252[k]
                   - f_6 * pfg1_252[k]
                   + f_4 * pc_x[k] * pfh_348[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, pc_x, pfg0_253, pfg0_254, \
                         pfg1_253, pfg1_254, pfh_349, pfh_350, pfh_351, pfh_352, \
                         pfh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_5 * pfg0_253[k]
                   - f_6 * pfg1_253[k]
                   + f_4 * pc_x[k] * pfh_349[k];

        t_462[k] = f_5 * pfg0_254[k]
                   - f_6 * pfg1_254[k]
                   + f_4 * pc_x[k] * pfh_350[k];

        t_463[k] = f_4 * pc_x[k] * pfh_351[k];

        t_464[k] = f_4 * pc_x[k] * pfh_352[k];

        t_465[k] = f_4 * pc_x[k] * pfh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, pc_x, pc_y, pc_z, sfh_141, \
                         pdh_204, pfg0_250, pfg1_250, pfh_351, pfh_354, pfh_355, \
                         pfh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_4 * pc_x[k] * pfh_354[k];

        t_467[k] = f_4 * pc_x[k] * pfh_355[k];

        t_468[k] = f_4 * pc_x[k] * pfh_356[k];

        t_469[k] = f_0 * sfh_141[k]
                   + f_1 * pdh_204[k]
                   + f_2 * pfg0_250[k]
                   - f_3 * pfg1_250[k]
                   + f_4 * pc_y[k] * pfh_351[k];

        t_470[k] = f_4 * pc_z[k] * pfh_351[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_z, pfg0_250, pfg0_251, pfg0_252, pfg1_250, \
                         pfg1_251, pfg1_252, pfh_352, pfh_353, \
                         pfh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_5 * pfg0_250[k]
                   - f_6 * pfg1_250[k]
                   + f_4 * pc_z[k] * pfh_352[k];

        t_472[k] = f_7 * pfg0_251[k]
                   - f_8 * pfg1_251[k]
                   + f_4 * pc_z[k] * pfh_353[k];

        t_473[k] = f_9 * pfg0_252[k]
                   - f_10 * pfg1_252[k]
                   + f_4 * pc_z[k] * pfh_354[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_z, pc_y, pc_z, sfh_146, pdi0_252, \
                         pdi0_253, pdh_209, pdi1_252, pdi1_253, pfg0_254, pfg1_254, \
                         pfh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_0 * sfh_146[k]
                   + f_1 * pdh_209[k]
                   + f_4 * pc_y[k] * pfh_356[k];

        t_475[k] = f_2 * pfg0_254[k]
                   - f_3 * pfg1_254[k]
                   + f_4 * pc_z[k] * pfh_356[k];

        t_476[k] = pb_z[k] * pdi0_252[k]
                   - f_11 * pc_z[k] * pdi1_252[k];

        t_477[k] = pb_z[k] * pdi0_253[k]
                   - f_11 * pc_z[k] * pdi1_253[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pb_z, pc_x, pc_z, pdi0_255, pdh_189, \
                         pdh_190, pdi1_255, pfg0_260, pfg1_260, pfh_357, pfh_358, \
                         pfh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_0 * pdh_189[k]
                   + f_4 * pc_z[k] * pfh_357[k];

        t_479[k] = pb_z[k] * pdi0_255[k]
                   - f_11 * pc_z[k] * pdi1_255[k];

        t_480[k] = f_0 * pdh_190[k]
                   + f_4 * pc_z[k] * pfh_358[k];

        t_481[k] = f_9 * pfg0_260[k]
                   - f_10 * pfg1_260[k]
                   + f_4 * pc_x[k] * pfh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pb_z, pc_x, pc_z, pdi0_258, pdh_192, pdi1_258, \
                         pfg0_263, pfg1_263, pfh_360, pfh_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = pb_z[k] * pdi0_258[k]
                   - f_11 * pc_z[k] * pdi1_258[k];

        t_483[k] = f_0 * pdh_192[k]
                   + f_4 * pc_z[k] * pfh_360[k];

        t_484[k] = f_7 * pfg0_263[k]
                   - f_8 * pfg1_263[k]
                   + f_4 * pc_x[k] * pfh_365[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pb_z, pc_x, pc_z, pdi0_262, pdh_195, pdi1_262, \
                         pfg0_264, pfg1_264, pfh_363, pfh_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_7 * pfg0_264[k]
                   - f_8 * pfg1_264[k]
                   + f_4 * pc_x[k] * pfh_366[k];

        t_486[k] = pb_z[k] * pdi0_262[k]
                   - f_11 * pc_z[k] * pdi1_262[k];

        t_487[k] = f_0 * pdh_195[k]
                   + f_4 * pc_z[k] * pfh_363[k];
    }
}

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfi0, const size_t sfh,
                                                          const size_t sfi1, const size_t pdi0,
                                                          const size_t pdh, const size_t pdi1,
                                                          const size_t pfg0, const size_t pfg1,
                                                          const size_t pfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_0 = buffer.data(sfi0 + 0);
    const auto *sfi0_2 = buffer.data(sfi0 + 2);
    const auto *sfi0_3 = buffer.data(sfi0 + 3);
    const auto *sfi0_5 = buffer.data(sfi0 + 5);
    const auto *sfi0_6 = buffer.data(sfi0 + 6);
    const auto *sfi0_7 = buffer.data(sfi0 + 7);
    const auto *sfi0_9 = buffer.data(sfi0 + 9);
    const auto *sfi0_10 = buffer.data(sfi0 + 10);
    const auto *sfi0_11 = buffer.data(sfi0 + 11);
    const auto *sfi0_12 = buffer.data(sfi0 + 12);
    const auto *sfi0_14 = buffer.data(sfi0 + 14);
    const auto *sfi0_21 = buffer.data(sfi0 + 21);
    const auto *sfi0_27 = buffer.data(sfi0 + 27);
    const auto *sfi0_28 = buffer.data(sfi0 + 28);
    const auto *sfi0_31 = buffer.data(sfi0 + 31);
    const auto *sfi0_34 = buffer.data(sfi0 + 34);
    const auto *sfi0_35 = buffer.data(sfi0 + 35);
    const auto *sfi0_38 = buffer.data(sfi0 + 38);
    const auto *sfi0_39 = buffer.data(sfi0 + 39);
    const auto *sfi0_40 = buffer.data(sfi0 + 40);
    const auto *sfi0_49 = buffer.data(sfi0 + 49);
    const auto *sfi0_252 = buffer.data(sfi0 + 252);
    const auto *sfi0_257 = buffer.data(sfi0 + 257);
    const auto *sfi0_261 = buffer.data(sfi0 + 261);
    const auto *sfi0_266 = buffer.data(sfi0 + 266);
    const auto *sfi0_273 = buffer.data(sfi0 + 273);
    const auto *sfi0_275 = buffer.data(sfi0 + 275);
    const auto *sfi0_276 = buffer.data(sfi0 + 276);
    const auto *sfi0_277 = buffer.data(sfi0 + 277);
    const auto *sfi0_279 = buffer.data(sfi0 + 279);

    const auto *sfh_0 = buffer.data(sfh + 0);
    const auto *sfh_2 = buffer.data(sfh + 2);
    const auto *sfh_3 = buffer.data(sfh + 3);
    const auto *sfh_5 = buffer.data(sfh + 5);
    const auto *sfh_6 = buffer.data(sfh + 6);
    const auto *sfh_7 = buffer.data(sfh + 7);
    const auto *sfh_9 = buffer.data(sfh + 9);
    const auto *sfh_20 = buffer.data(sfh + 20);
    const auto *sfh_24 = buffer.data(sfh + 24);
    const auto *sfh_27 = buffer.data(sfh + 27);
    const auto *sfh_28 = buffer.data(sfh + 28);
    const auto *sfh_167 = buffer.data(sfh + 167);
    const auto *sfh_183 = buffer.data(sfh + 183);
    const auto *sfh_188 = buffer.data(sfh + 188);
    const auto *sfh_204 = buffer.data(sfh + 204);
    const auto *sfh_206 = buffer.data(sfh + 206);
    const auto *sfh_207 = buffer.data(sfh + 207);
    const auto *sfh_208 = buffer.data(sfh + 208);
    const auto *sfh_209 = buffer.data(sfh + 209);

    const auto *sfi1_0 = buffer.data(sfi1 + 0);
    const auto *sfi1_2 = buffer.data(sfi1 + 2);
    const auto *sfi1_3 = buffer.data(sfi1 + 3);
    const auto *sfi1_5 = buffer.data(sfi1 + 5);
    const auto *sfi1_6 = buffer.data(sfi1 + 6);
    const auto *sfi1_7 = buffer.data(sfi1 + 7);
    const auto *sfi1_9 = buffer.data(sfi1 + 9);
    const auto *sfi1_10 = buffer.data(sfi1 + 10);
    const auto *sfi1_11 = buffer.data(sfi1 + 11);
    const auto *sfi1_12 = buffer.data(sfi1 + 12);
    const auto *sfi1_14 = buffer.data(sfi1 + 14);
    const auto *sfi1_21 = buffer.data(sfi1 + 21);
    const auto *sfi1_27 = buffer.data(sfi1 + 27);
    const auto *sfi1_28 = buffer.data(sfi1 + 28);
    const auto *sfi1_31 = buffer.data(sfi1 + 31);
    const auto *sfi1_34 = buffer.data(sfi1 + 34);
    const auto *sfi1_35 = buffer.data(sfi1 + 35);
    const auto *sfi1_38 = buffer.data(sfi1 + 38);
    const auto *sfi1_39 = buffer.data(sfi1 + 39);
    const auto *sfi1_40 = buffer.data(sfi1 + 40);
    const auto *sfi1_49 = buffer.data(sfi1 + 49);
    const auto *sfi1_252 = buffer.data(sfi1 + 252);
    const auto *sfi1_257 = buffer.data(sfi1 + 257);
    const auto *sfi1_261 = buffer.data(sfi1 + 261);
    const auto *sfi1_266 = buffer.data(sfi1 + 266);
    const auto *sfi1_273 = buffer.data(sfi1 + 273);
    const auto *sfi1_275 = buffer.data(sfi1 + 275);
    const auto *sfi1_276 = buffer.data(sfi1 + 276);
    const auto *sfi1_277 = buffer.data(sfi1 + 277);
    const auto *sfi1_279 = buffer.data(sfi1 + 279);

    const auto *pdi0_273 = buffer.data(pdi0 + 273);
    const auto *pdi0_275 = buffer.data(pdi0 + 275);
    const auto *pdi0_276 = buffer.data(pdi0 + 276);
    const auto *pdi0_277 = buffer.data(pdi0 + 277);
    const auto *pdi0_338 = buffer.data(pdi0 + 338);
    const auto *pdi0_341 = buffer.data(pdi0 + 341);
    const auto *pdi0_345 = buffer.data(pdi0 + 345);
    const auto *pdi0_350 = buffer.data(pdi0 + 350);

    const auto *pdh_204 = buffer.data(pdh + 204);
    const auto *pdh_205 = buffer.data(pdh + 205);
    const auto *pdh_206 = buffer.data(pdh + 206);
    const auto *pdh_207 = buffer.data(pdh + 207);
    const auto *pdh_209 = buffer.data(pdh + 209);
    const auto *pdh_210 = buffer.data(pdh + 210);
    const auto *pdh_211 = buffer.data(pdh + 211);
    const auto *pdh_213 = buffer.data(pdh + 213);
    const auto *pdh_216 = buffer.data(pdh + 216);
    const auto *pdh_225 = buffer.data(pdh + 225);
    const auto *pdh_226 = buffer.data(pdh + 226);
    const auto *pdh_227 = buffer.data(pdh + 227);
    const auto *pdh_228 = buffer.data(pdh + 228);
    const auto *pdh_230 = buffer.data(pdh + 230);
    const auto *pdh_231 = buffer.data(pdh + 231);
    const auto *pdh_232 = buffer.data(pdh + 232);
    const auto *pdh_234 = buffer.data(pdh + 234);
    const auto *pdh_237 = buffer.data(pdh + 237);
    const auto *pdh_246 = buffer.data(pdh + 246);
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
    const auto *pdh_288 = buffer.data(pdh + 288);
    const auto *pdh_289 = buffer.data(pdh + 289);
    const auto *pdh_290 = buffer.data(pdh + 290);
    const auto *pdh_291 = buffer.data(pdh + 291);
    const auto *pdh_292 = buffer.data(pdh + 292);
    const auto *pdh_293 = buffer.data(pdh + 293);

    const auto *pdi1_273 = buffer.data(pdi1 + 273);
    const auto *pdi1_275 = buffer.data(pdi1 + 275);
    const auto *pdi1_276 = buffer.data(pdi1 + 276);
    const auto *pdi1_277 = buffer.data(pdi1 + 277);
    const auto *pdi1_338 = buffer.data(pdi1 + 338);
    const auto *pdi1_341 = buffer.data(pdi1 + 341);
    const auto *pdi1_345 = buffer.data(pdi1 + 345);
    const auto *pdi1_350 = buffer.data(pdi1 + 350);

    const auto *pfg0_267 = buffer.data(pfg0 + 267);
    const auto *pfg0_268 = buffer.data(pfg0 + 268);
    const auto *pfg0_269 = buffer.data(pfg0 + 269);
    const auto *pfg0_270 = buffer.data(pfg0 + 270);
    const auto *pfg0_271 = buffer.data(pfg0 + 271);
    const auto *pfg0_273 = buffer.data(pfg0 + 273);
    const auto *pfg0_275 = buffer.data(pfg0 + 275);
    const auto *pfg0_276 = buffer.data(pfg0 + 276);
    const auto *pfg0_278 = buffer.data(pfg0 + 278);
    const auto *pfg0_279 = buffer.data(pfg0 + 279);
    const auto *pfg0_280 = buffer.data(pfg0 + 280);
    const auto *pfg0_281 = buffer.data(pfg0 + 281);
    const auto *pfg0_282 = buffer.data(pfg0 + 282);
    const auto *pfg0_283 = buffer.data(pfg0 + 283);
    const auto *pfg0_284 = buffer.data(pfg0 + 284);
    const auto *pfg0_286 = buffer.data(pfg0 + 286);
    const auto *pfg0_288 = buffer.data(pfg0 + 288);
    const auto *pfg0_291 = buffer.data(pfg0 + 291);
    const auto *pfg0_293 = buffer.data(pfg0 + 293);
    const auto *pfg0_295 = buffer.data(pfg0 + 295);
    const auto *pfg0_297 = buffer.data(pfg0 + 297);
    const auto *pfg0_298 = buffer.data(pfg0 + 298);
    const auto *pfg0_311 = buffer.data(pfg0 + 311);
    const auto *pfg0_312 = buffer.data(pfg0 + 312);
    const auto *pfg0_313 = buffer.data(pfg0 + 313);
    const auto *pfg0_314 = buffer.data(pfg0 + 314);

    const auto *pfg1_267 = buffer.data(pfg1 + 267);
    const auto *pfg1_268 = buffer.data(pfg1 + 268);
    const auto *pfg1_269 = buffer.data(pfg1 + 269);
    const auto *pfg1_270 = buffer.data(pfg1 + 270);
    const auto *pfg1_271 = buffer.data(pfg1 + 271);
    const auto *pfg1_273 = buffer.data(pfg1 + 273);
    const auto *pfg1_275 = buffer.data(pfg1 + 275);
    const auto *pfg1_276 = buffer.data(pfg1 + 276);
    const auto *pfg1_278 = buffer.data(pfg1 + 278);
    const auto *pfg1_279 = buffer.data(pfg1 + 279);
    const auto *pfg1_280 = buffer.data(pfg1 + 280);
    const auto *pfg1_281 = buffer.data(pfg1 + 281);
    const auto *pfg1_282 = buffer.data(pfg1 + 282);
    const auto *pfg1_283 = buffer.data(pfg1 + 283);
    const auto *pfg1_284 = buffer.data(pfg1 + 284);
    const auto *pfg1_286 = buffer.data(pfg1 + 286);
    const auto *pfg1_288 = buffer.data(pfg1 + 288);
    const auto *pfg1_291 = buffer.data(pfg1 + 291);
    const auto *pfg1_293 = buffer.data(pfg1 + 293);
    const auto *pfg1_295 = buffer.data(pfg1 + 295);
    const auto *pfg1_297 = buffer.data(pfg1 + 297);
    const auto *pfg1_298 = buffer.data(pfg1 + 298);
    const auto *pfg1_311 = buffer.data(pfg1 + 311);
    const auto *pfg1_312 = buffer.data(pfg1 + 312);
    const auto *pfg1_313 = buffer.data(pfg1 + 313);
    const auto *pfg1_314 = buffer.data(pfg1 + 314);

    const auto *pfh_369 = buffer.data(pfh + 369);
    const auto *pfh_370 = buffer.data(pfh + 370);
    const auto *pfh_371 = buffer.data(pfh + 371);
    const auto *pfh_372 = buffer.data(pfh + 372);
    const auto *pfh_373 = buffer.data(pfh + 373);
    const auto *pfh_374 = buffer.data(pfh + 374);
    const auto *pfh_375 = buffer.data(pfh + 375);
    const auto *pfh_376 = buffer.data(pfh + 376);
    const auto *pfh_377 = buffer.data(pfh + 377);
    const auto *pfh_378 = buffer.data(pfh + 378);
    const auto *pfh_379 = buffer.data(pfh + 379);
    const auto *pfh_381 = buffer.data(pfh + 381);
    const auto *pfh_383 = buffer.data(pfh + 383);
    const auto *pfh_384 = buffer.data(pfh + 384);
    const auto *pfh_386 = buffer.data(pfh + 386);
    const auto *pfh_387 = buffer.data(pfh + 387);
    const auto *pfh_388 = buffer.data(pfh + 388);
    const auto *pfh_390 = buffer.data(pfh + 390);
    const auto *pfh_391 = buffer.data(pfh + 391);
    const auto *pfh_392 = buffer.data(pfh + 392);
    const auto *pfh_393 = buffer.data(pfh + 393);
    const auto *pfh_394 = buffer.data(pfh + 394);
    const auto *pfh_395 = buffer.data(pfh + 395);
    const auto *pfh_396 = buffer.data(pfh + 396);
    const auto *pfh_397 = buffer.data(pfh + 397);
    const auto *pfh_398 = buffer.data(pfh + 398);
    const auto *pfh_399 = buffer.data(pfh + 399);
    const auto *pfh_400 = buffer.data(pfh + 400);
    const auto *pfh_402 = buffer.data(pfh + 402);
    const auto *pfh_405 = buffer.data(pfh + 405);
    const auto *pfh_407 = buffer.data(pfh + 407);
    const auto *pfh_409 = buffer.data(pfh + 409);
    const auto *pfh_411 = buffer.data(pfh + 411);
    const auto *pfh_412 = buffer.data(pfh + 412);
    const auto *pfh_414 = buffer.data(pfh + 414);
    const auto *pfh_415 = buffer.data(pfh + 415);
    const auto *pfh_416 = buffer.data(pfh + 416);
    const auto *pfh_417 = buffer.data(pfh + 417);
    const auto *pfh_418 = buffer.data(pfh + 418);
    const auto *pfh_419 = buffer.data(pfh + 419);
    const auto *pfh_420 = buffer.data(pfh + 420);
    const auto *pfh_422 = buffer.data(pfh + 422);
    const auto *pfh_425 = buffer.data(pfh + 425);
    const auto *pfh_429 = buffer.data(pfh + 429);
    const auto *pfh_435 = buffer.data(pfh + 435);
    const auto *pfh_436 = buffer.data(pfh + 436);
    const auto *pfh_437 = buffer.data(pfh + 437);
    const auto *pfh_438 = buffer.data(pfh + 438);
    const auto *pfh_439 = buffer.data(pfh + 439);
    const auto *pfh_440 = buffer.data(pfh + 440);
    const auto *pfh_441 = buffer.data(pfh + 441);
    const auto *pfh_443 = buffer.data(pfh + 443);
    const auto *pfh_446 = buffer.data(pfh + 446);
    const auto *pfh_450 = buffer.data(pfh + 450);
    const auto *pfh_456 = buffer.data(pfh + 456);
    const auto *pfh_457 = buffer.data(pfh + 457);
    const auto *pfh_458 = buffer.data(pfh + 458);
    const auto *pfh_459 = buffer.data(pfh + 459);
    const auto *pfh_460 = buffer.data(pfh + 460);
    const auto *pfh_461 = buffer.data(pfh + 461);

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pc_x, pfg0_267, pfg0_268, pfg0_269, \
                         pfg1_267, pfg1_268, pfg1_269, pfh_369, pfh_370, pfh_371, \
                         pfh_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * pfg0_267[k]
                   - f_6 * pfg1_267[k]
                   + f_4 * pc_x[k] * pfh_369[k];

        t_489[k] = f_5 * pfg0_268[k]
                   - f_6 * pfg1_268[k]
                   + f_4 * pc_x[k] * pfh_370[k];

        t_490[k] = f_5 * pfg0_269[k]
                   - f_6 * pfg1_269[k]
                   + f_4 * pc_x[k] * pfh_371[k];

        t_491[k] = f_4 * pc_x[k] * pfh_372[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, t_497, pb_z, pc_x, pc_z, pdi0_273, \
                         pdi1_273, pfh_373, pfh_374, pfh_375, pfh_376, \
                         pfh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_4 * pc_x[k] * pfh_373[k];

        t_493[k] = f_4 * pc_x[k] * pfh_374[k];

        t_494[k] = f_4 * pc_x[k] * pfh_375[k];

        t_495[k] = f_4 * pc_x[k] * pfh_376[k];

        t_496[k] = f_4 * pc_x[k] * pfh_377[k];

        t_497[k] = pb_z[k] * pdi0_273[k]
                   - f_11 * pc_z[k] * pdi1_273[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pb_z, pc_z, pdi0_275, pdi0_276, pdh_204, \
                         pdh_205, pdh_206, pdi1_275, pdi1_276, \
                         pfh_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_0 * pdh_204[k]
                   + f_4 * pc_z[k] * pfh_372[k];

        t_499[k] = pb_z[k] * pdi0_275[k]
                   + f_12 * pdh_205[k]
                   - f_11 * pc_z[k] * pdi1_275[k];

        t_500[k] = pb_z[k] * pdi0_276[k]
                   + f_1 * pdh_206[k]
                   - f_11 * pc_z[k] * pdi1_276[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pb_z, pc_y, pc_z, sfh_167, pdi0_277, pdh_207, \
                         pdh_209, pdh_230, pdi1_277, pfg0_269, pfg1_269, \
                         pfh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = pb_z[k] * pdi0_277[k]
                   + f_13 * pdh_207[k]
                   - f_11 * pc_z[k] * pdi1_277[k];

        t_502[k] = f_0 * sfh_167[k]
                   + f_12 * pdh_230[k]
                   + f_4 * pc_y[k] * pfh_377[k];

        t_503[k] = f_0 * pdh_209[k]
                   + f_2 * pfg0_269[k]
                   - f_3 * pfg1_269[k]
                   + f_4 * pc_z[k] * pfh_377[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_z, pdh_210, pfg0_270, pfg0_271, \
                         pfg0_273, pfg1_270, pfg1_271, pfg1_273, pfh_378, pfh_379, \
                         pfh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_2 * pfg0_270[k]
                   - f_3 * pfg1_270[k]
                   + f_4 * pc_x[k] * pfh_378[k];

        t_505[k] = f_15 * pfg0_271[k]
                   - f_16 * pfg1_271[k]
                   + f_4 * pc_x[k] * pfh_379[k];

        t_506[k] = f_12 * pdh_210[k]
                   + f_4 * pc_z[k] * pfh_378[k];

        t_507[k] = f_9 * pfg0_273[k]
                   - f_10 * pfg1_273[k]
                   + f_4 * pc_x[k] * pfh_381[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pc_x, pc_z, pdh_211, pdh_213, pfg0_275, \
                         pfg0_276, pfg1_275, pfg1_276, pfh_379, pfh_381, pfh_383, \
                         pfh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_12 * pdh_211[k]
                   + f_4 * pc_z[k] * pfh_379[k];

        t_509[k] = f_9 * pfg0_275[k]
                   - f_10 * pfg1_275[k]
                   + f_4 * pc_x[k] * pfh_383[k];

        t_510[k] = f_7 * pfg0_276[k]
                   - f_8 * pfg1_276[k]
                   + f_4 * pc_x[k] * pfh_384[k];

        t_511[k] = f_12 * pdh_213[k]
                   + f_4 * pc_z[k] * pfh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pfg0_278, pfg0_279, pfg0_280, pfg1_278, \
                         pfg1_279, pfg1_280, pfh_386, pfh_387, \
                         pfh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_7 * pfg0_278[k]
                   - f_8 * pfg1_278[k]
                   + f_4 * pc_x[k] * pfh_386[k];

        t_513[k] = f_7 * pfg0_279[k]
                   - f_8 * pfg1_279[k]
                   + f_4 * pc_x[k] * pfh_387[k];

        t_514[k] = f_5 * pfg0_280[k]
                   - f_6 * pfg1_280[k]
                   + f_4 * pc_x[k] * pfh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_z, pdh_216, pfg0_282, pfg0_283, \
                         pfg1_282, pfg1_283, pfh_384, pfh_390, \
                         pfh_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_12 * pdh_216[k]
                   + f_4 * pc_z[k] * pfh_384[k];

        t_516[k] = f_5 * pfg0_282[k]
                   - f_6 * pfg1_282[k]
                   + f_4 * pc_x[k] * pfh_390[k];

        t_517[k] = f_5 * pfg0_283[k]
                   - f_6 * pfg1_283[k]
                   + f_4 * pc_x[k] * pfh_391[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, t_522, t_523, pc_x, pfg0_284, pfg1_284, \
                         pfh_392, pfh_393, pfh_394, pfh_395, pfh_396, \
                         pfh_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_5 * pfg0_284[k]
                   - f_6 * pfg1_284[k]
                   + f_4 * pc_x[k] * pfh_392[k];

        t_519[k] = f_4 * pc_x[k] * pfh_393[k];

        t_520[k] = f_4 * pc_x[k] * pfh_394[k];

        t_521[k] = f_4 * pc_x[k] * pfh_395[k];

        t_522[k] = f_4 * pc_x[k] * pfh_396[k];

        t_523[k] = f_4 * pc_x[k] * pfh_397[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, pc_x, pc_y, pc_z, sfh_183, pdh_225, \
                         pdh_226, pdh_246, pfg0_280, pfg1_280, pfh_393, pfh_394, \
                         pfh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_4 * pc_x[k] * pfh_398[k];

        t_525[k] = f_0 * sfh_183[k]
                   + f_0 * pdh_246[k]
                   + f_2 * pfg0_280[k]
                   - f_3 * pfg1_280[k]
                   + f_4 * pc_y[k] * pfh_393[k];

        t_526[k] = f_12 * pdh_225[k]
                   + f_4 * pc_z[k] * pfh_393[k];

        t_527[k] = f_12 * pdh_226[k]
                   + f_5 * pfg0_280[k]
                   - f_6 * pfg1_280[k]
                   + f_4 * pc_z[k] * pfh_394[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pc_y, pc_z, sfh_188, pdh_227, pdh_228, pdh_251, \
                         pfg0_281, pfg0_282, pfg1_281, pfg1_282, pfh_395, pfh_396, \
                         pfh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_12 * pdh_227[k]
                   + f_7 * pfg0_281[k]
                   - f_8 * pfg1_281[k]
                   + f_4 * pc_z[k] * pfh_395[k];

        t_529[k] = f_12 * pdh_228[k]
                   + f_9 * pfg0_282[k]
                   - f_10 * pfg1_282[k]
                   + f_4 * pc_z[k] * pfh_396[k];

        t_530[k] = f_0 * sfh_188[k]
                   + f_0 * pdh_251[k]
                   + f_4 * pc_y[k] * pfh_398[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pa_y, pc_x, pc_y, pc_z, sfi0_252, sfi1_252, \
                         pdh_230, pfg0_284, pfg0_286, pfg1_284, pfg1_286, pfh_398, \
                         pfh_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_12 * pdh_230[k]
                   + f_2 * pfg0_284[k]
                   - f_3 * pfg1_284[k]
                   + f_4 * pc_z[k] * pfh_398[k];

        t_532[k] = pa_y[k] * sfi0_252[k]
                   - f_11 * pc_y[k] * sfi1_252[k];

        t_533[k] = f_15 * pfg0_286[k]
                   - f_16 * pfg1_286[k]
                   + f_4 * pc_x[k] * pfh_400[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_x, pc_z, pdh_231, pdh_232, pfg0_288, \
                         pfg1_288, pfh_399, pfh_400, pfh_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_1 * pdh_231[k]
                   + f_4 * pc_z[k] * pfh_399[k];

        t_535[k] = f_9 * pfg0_288[k]
                   - f_10 * pfg1_288[k]
                   + f_4 * pc_x[k] * pfh_402[k];

        t_536[k] = f_1 * pdh_232[k]
                   + f_4 * pc_z[k] * pfh_400[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pa_y, pc_x, pc_y, pc_z, sfi0_257, sfi1_257, \
                         pdh_234, pfg0_291, pfg1_291, pfh_402, \
                         pfh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * sfi0_257[k]
                   - f_11 * pc_y[k] * sfi1_257[k];

        t_538[k] = f_7 * pfg0_291[k]
                   - f_8 * pfg1_291[k]
                   + f_4 * pc_x[k] * pfh_405[k];

        t_539[k] = f_1 * pdh_234[k]
                   + f_4 * pc_z[k] * pfh_402[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_y, pc_x, pc_y, sfi0_261, sfi1_261, pfg0_293, \
                         pfg0_295, pfg1_293, pfg1_295, pfh_407, \
                         pfh_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_7 * pfg0_293[k]
                   - f_8 * pfg1_293[k]
                   + f_4 * pc_x[k] * pfh_407[k];

        t_541[k] = pa_y[k] * sfi0_261[k]
                   - f_11 * pc_y[k] * sfi1_261[k];

        t_542[k] = f_5 * pfg0_295[k]
                   - f_6 * pfg1_295[k]
                   + f_4 * pc_x[k] * pfh_409[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pc_x, pc_z, pdh_237, pfg0_297, pfg0_298, \
                         pfg1_297, pfg1_298, pfh_405, pfh_411, \
                         pfh_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_1 * pdh_237[k]
                   + f_4 * pc_z[k] * pfh_405[k];

        t_544[k] = f_5 * pfg0_297[k]
                   - f_6 * pfg1_297[k]
                   + f_4 * pc_x[k] * pfh_411[k];

        t_545[k] = f_5 * pfg0_298[k]
                   - f_6 * pfg1_298[k]
                   + f_4 * pc_x[k] * pfh_412[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, t_551, pa_y, pc_x, pc_y, sfi0_266, \
                         sfi1_266, pfh_414, pfh_415, pfh_416, pfh_417, \
                         pfh_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * sfi0_266[k]
                   - f_11 * pc_y[k] * sfi1_266[k];

        t_547[k] = f_4 * pc_x[k] * pfh_414[k];

        t_548[k] = f_4 * pc_x[k] * pfh_415[k];

        t_549[k] = f_4 * pc_x[k] * pfh_416[k];

        t_550[k] = f_4 * pc_x[k] * pfh_417[k];

        t_551[k] = f_4 * pc_x[k] * pfh_418[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_y, pc_x, pc_y, pc_z, sfi0_273, sfh_204, \
                         sfi1_273, pdh_246, pfh_414, pfh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_4 * pc_x[k] * pfh_419[k];

        t_553[k] = pa_y[k] * sfi0_273[k]
                   + f_14 * sfh_204[k]
                   - f_11 * pc_y[k] * sfi1_273[k];

        t_554[k] = f_1 * pdh_246[k]
                   + f_4 * pc_z[k] * pfh_414[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pa_y, pc_y, sfi0_275, sfi0_276, sfi0_277, \
                         sfh_206, sfh_207, sfh_208, sfi1_275, sfi1_276, \
                         sfi1_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_y[k] * sfi0_275[k]
                   + f_13 * sfh_206[k]
                   - f_11 * pc_y[k] * sfi1_275[k];

        t_556[k] = pa_y[k] * sfi0_276[k]
                   + f_1 * sfh_207[k]
                   - f_11 * pc_y[k] * sfi1_276[k];

        t_557[k] = pa_y[k] * sfi0_277[k]
                   + f_12 * sfh_208[k]
                   - f_11 * pc_y[k] * sfi1_277[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pa_y, pa_z, pc_y, pc_z, sfi0_0, sfi0_279, \
                         sfh_209, sfi1_0, sfi1_279, pfh_419, pfh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_0 * sfh_209[k]
                   + f_4 * pc_y[k] * pfh_419[k];

        t_559[k] = pa_y[k] * sfi0_279[k]
                   - f_11 * pc_y[k] * sfi1_279[k];

        t_560[k] = pa_z[k] * sfi0_0[k]
                   - f_11 * pc_z[k] * sfi1_0[k];

        t_561[k] = f_4 * pc_y[k] * pfh_420[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, pa_z, pc_y, pc_z, sfi0_2, sfi0_3, sfi0_5, \
                         sfh_0, sfh_2, sfi1_2, sfi1_3, sfi1_5, \
                         pfh_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = pa_z[k] * sfi0_2[k]
                   + f_0 * sfh_0[k]
                   - f_11 * pc_z[k] * sfi1_2[k];

        t_563[k] = pa_z[k] * sfi0_3[k]
                   - f_11 * pc_z[k] * sfi1_3[k];

        t_564[k] = f_4 * pc_y[k] * pfh_422[k];

        t_565[k] = pa_z[k] * sfi0_5[k]
                   + f_12 * sfh_2[k]
                   - f_11 * pc_z[k] * sfi1_5[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_z, pc_y, pc_z, sfi0_6, sfi0_7, sfi0_9, \
                         sfh_3, sfh_5, sfi1_6, sfi1_7, sfi1_9, \
                         pfh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pa_z[k] * sfi0_6[k]
                   - f_11 * pc_z[k] * sfi1_6[k];

        t_567[k] = pa_z[k] * sfi0_7[k]
                   + f_0 * sfh_3[k]
                   - f_11 * pc_z[k] * sfi1_7[k];

        t_568[k] = f_4 * pc_y[k] * pfh_425[k];

        t_569[k] = pa_z[k] * sfi0_9[k]
                   + f_1 * sfh_5[k]
                   - f_11 * pc_z[k] * sfi1_9[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pa_z, pc_y, pc_z, sfi0_10, sfi0_11, \
                         sfi0_12, sfh_6, sfh_7, sfi1_10, sfi1_11, sfi1_12, \
                         pfh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = pa_z[k] * sfi0_10[k]
                   - f_11 * pc_z[k] * sfi1_10[k];

        t_571[k] = pa_z[k] * sfi0_11[k]
                   + f_0 * sfh_6[k]
                   - f_11 * pc_z[k] * sfi1_11[k];

        t_572[k] = pa_z[k] * sfi0_12[k]
                   + f_12 * sfh_7[k]
                   - f_11 * pc_z[k] * sfi1_12[k];

        t_573[k] = f_4 * pc_y[k] * pfh_429[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_z, pc_x, pc_z, sfi0_14, sfh_9, \
                         sfi1_14, pdh_267, pdh_268, pdh_269, pfh_435, pfh_436, \
                         pfh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = pa_z[k] * sfi0_14[k]
                   + f_13 * sfh_9[k]
                   - f_11 * pc_z[k] * sfi1_14[k];

        t_575[k] = f_1 * pdh_267[k]
                   + f_4 * pc_x[k] * pfh_435[k];

        t_576[k] = f_1 * pdh_268[k]
                   + f_4 * pc_x[k] * pfh_436[k];

        t_577[k] = f_1 * pdh_269[k]
                   + f_4 * pc_x[k] * pfh_437[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pa_z, pc_x, pc_z, sfi0_21, sfi1_21, \
                         pdh_270, pdh_271, pdh_272, pfh_438, pfh_439, \
                         pfh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_1 * pdh_270[k]
                   + f_4 * pc_x[k] * pfh_438[k];

        t_579[k] = f_1 * pdh_271[k]
                   + f_4 * pc_x[k] * pfh_439[k];

        t_580[k] = f_1 * pdh_272[k]
                   + f_4 * pc_x[k] * pfh_440[k];

        t_581[k] = pa_z[k] * sfi0_21[k]
                   - f_11 * pc_z[k] * sfi1_21[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, pfg0_311, pfg0_312, pfg0_313, pfg1_311, \
                         pfg1_312, pfg1_313, pfh_436, pfh_437, \
                         pfh_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_15 * pfg0_311[k]
                   - f_16 * pfg1_311[k]
                   + f_4 * pc_y[k] * pfh_436[k];

        t_583[k] = f_9 * pfg0_312[k]
                   - f_10 * pfg1_312[k]
                   + f_4 * pc_y[k] * pfh_437[k];

        t_584[k] = f_7 * pfg0_313[k]
                   - f_8 * pfg1_313[k]
                   + f_4 * pc_y[k] * pfh_438[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pa_z, pc_y, pc_z, sfi0_27, sfi0_28, \
                         sfh_20, sfi1_27, sfi1_28, pfg0_314, pfg1_314, pfh_439, \
                         pfh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_5 * pfg0_314[k]
                   - f_6 * pfg1_314[k]
                   + f_4 * pc_y[k] * pfh_439[k];

        t_586[k] = f_4 * pc_y[k] * pfh_440[k];

        t_587[k] = pa_z[k] * sfi0_27[k]
                   + f_14 * sfh_20[k]
                   - f_11 * pc_z[k] * sfi1_27[k];

        t_588[k] = pa_z[k] * sfi0_28[k]
                   - f_11 * pc_z[k] * sfi1_28[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, pa_z, pb_y, pc_y, pc_z, sfi0_31, sfi1_31, \
                         pdi0_338, pdh_252, pdh_254, pdi1_338, pfh_441, \
                         pfh_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_0 * pdh_252[k]
                   + f_4 * pc_y[k] * pfh_441[k];

        t_590[k] = pb_y[k] * pdi0_338[k]
                   - f_11 * pc_y[k] * pdi1_338[k];

        t_591[k] = pa_z[k] * sfi0_31[k]
                   - f_11 * pc_z[k] * sfi1_31[k];

        t_592[k] = f_0 * pdh_254[k]
                   + f_4 * pc_y[k] * pfh_443[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, pa_z, pb_y, pc_y, pc_z, sfi0_34, sfi0_35, \
                         sfh_24, sfi1_34, sfi1_35, pdi0_341, pdi1_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_y[k] * pdi0_341[k]
                   - f_11 * pc_y[k] * pdi1_341[k];

        t_594[k] = pa_z[k] * sfi0_34[k]
                   - f_11 * pc_z[k] * sfi1_34[k];

        t_595[k] = pa_z[k] * sfi0_35[k]
                   + f_0 * sfh_24[k]
                   - f_11 * pc_z[k] * sfi1_35[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, pa_z, pb_y, pc_y, pc_z, sfi0_38, sfi1_38, \
                         pdi0_345, pdh_257, pdi1_345, pfh_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_0 * pdh_257[k]
                   + f_4 * pc_y[k] * pfh_446[k];

        t_597[k] = pb_y[k] * pdi0_345[k]
                   - f_11 * pc_y[k] * pdi1_345[k];

        t_598[k] = pa_z[k] * sfi0_38[k]
                   - f_11 * pc_z[k] * sfi1_38[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, pa_z, pc_y, pc_z, sfi0_39, sfi0_40, sfh_27, \
                         sfh_28, sfi1_39, sfi1_40, pdh_261, pfh_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pa_z[k] * sfi0_39[k]
                   + f_0 * sfh_27[k]
                   - f_11 * pc_z[k] * sfi1_39[k];

        t_600[k] = pa_z[k] * sfi0_40[k]
                   + f_12 * sfh_28[k]
                   - f_11 * pc_z[k] * sfi1_40[k];

        t_601[k] = f_0 * pdh_261[k]
                   + f_4 * pc_y[k] * pfh_450[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, pb_y, pc_x, pc_y, pdi0_350, pdh_288, \
                         pdh_289, pdh_290, pdi1_350, pfh_456, pfh_457, \
                         pfh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = pb_y[k] * pdi0_350[k]
                   - f_11 * pc_y[k] * pdi1_350[k];

        t_603[k] = f_12 * pdh_288[k]
                   + f_4 * pc_x[k] * pfh_456[k];

        t_604[k] = f_12 * pdh_289[k]
                   + f_4 * pc_x[k] * pfh_457[k];

        t_605[k] = f_12 * pdh_290[k]
                   + f_4 * pc_x[k] * pfh_458[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, pa_z, pc_x, pc_z, sfi0_49, sfi1_49, \
                         pdh_291, pdh_292, pdh_293, pfh_459, pfh_460, \
                         pfh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = f_12 * pdh_291[k]
                   + f_4 * pc_x[k] * pfh_459[k];

        t_607[k] = f_12 * pdh_292[k]
                   + f_4 * pc_x[k] * pfh_460[k];

        t_608[k] = f_12 * pdh_293[k]
                   + f_4 * pc_x[k] * pfh_461[k];

        t_609[k] = pa_z[k] * sfi0_49[k]
                   - f_11 * pc_z[k] * sfi1_49[k];
    }
}

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfi0, const size_t sfh,
                                                          const size_t sfi1, const size_t ppi0,
                                                          const size_t ppi1, const size_t pdi0,
                                                          const size_t pdh, const size_t pdi1,
                                                          const size_t pfg0, const size_t pfg1,
                                                          const size_t pfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_50 = buffer.data(sfi0 + 50);
    const auto *sfi0_51 = buffer.data(sfi0 + 51);
    const auto *sfi0_52 = buffer.data(sfi0 + 52);
    const auto *sfi0_53 = buffer.data(sfi0 + 53);
    const auto *sfi0_84 = buffer.data(sfi0 + 84);
    const auto *sfi0_86 = buffer.data(sfi0 + 86);
    const auto *sfi0_87 = buffer.data(sfi0 + 87);
    const auto *sfi0_89 = buffer.data(sfi0 + 89);
    const auto *sfi0_90 = buffer.data(sfi0 + 90);
    const auto *sfi0_91 = buffer.data(sfi0 + 91);
    const auto *sfi0_93 = buffer.data(sfi0 + 93);
    const auto *sfi0_94 = buffer.data(sfi0 + 94);
    const auto *sfi0_95 = buffer.data(sfi0 + 95);
    const auto *sfi0_96 = buffer.data(sfi0 + 96);
    const auto *sfi0_98 = buffer.data(sfi0 + 98);
    const auto *sfi0_105 = buffer.data(sfi0 + 105);

    const auto *sfh_36 = buffer.data(sfh + 36);
    const auto *sfh_37 = buffer.data(sfh + 37);
    const auto *sfh_38 = buffer.data(sfh + 38);
    const auto *sfh_39 = buffer.data(sfh + 39);
    const auto *sfh_63 = buffer.data(sfh + 63);
    const auto *sfh_65 = buffer.data(sfh + 65);
    const auto *sfh_66 = buffer.data(sfh + 66);
    const auto *sfh_68 = buffer.data(sfh + 68);
    const auto *sfh_69 = buffer.data(sfh + 69);
    const auto *sfh_70 = buffer.data(sfh + 70);
    const auto *sfh_72 = buffer.data(sfh + 72);

    const auto *sfi1_50 = buffer.data(sfi1 + 50);
    const auto *sfi1_51 = buffer.data(sfi1 + 51);
    const auto *sfi1_52 = buffer.data(sfi1 + 52);
    const auto *sfi1_53 = buffer.data(sfi1 + 53);
    const auto *sfi1_84 = buffer.data(sfi1 + 84);
    const auto *sfi1_86 = buffer.data(sfi1 + 86);
    const auto *sfi1_87 = buffer.data(sfi1 + 87);
    const auto *sfi1_89 = buffer.data(sfi1 + 89);
    const auto *sfi1_90 = buffer.data(sfi1 + 90);
    const auto *sfi1_91 = buffer.data(sfi1 + 91);
    const auto *sfi1_93 = buffer.data(sfi1 + 93);
    const auto *sfi1_94 = buffer.data(sfi1 + 94);
    const auto *sfi1_95 = buffer.data(sfi1 + 95);
    const auto *sfi1_96 = buffer.data(sfi1 + 96);
    const auto *sfi1_98 = buffer.data(sfi1 + 98);
    const auto *sfi1_105 = buffer.data(sfi1 + 105);

    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *pdi0_363 = buffer.data(pdi0 + 363);
    const auto *pdi0_392 = buffer.data(pdi0 + 392);
    const auto *pdi0_394 = buffer.data(pdi0 + 394);
    const auto *pdi0_397 = buffer.data(pdi0 + 397);
    const auto *pdi0_401 = buffer.data(pdi0 + 401);
    const auto *pdi0_406 = buffer.data(pdi0 + 406);
    const auto *pdi0_419 = buffer.data(pdi0 + 419);
    const auto *pdi0_442 = buffer.data(pdi0 + 442);
    const auto *pdi0_443 = buffer.data(pdi0 + 443);
    const auto *pdi0_444 = buffer.data(pdi0 + 444);
    const auto *pdi0_445 = buffer.data(pdi0 + 445);
    const auto *pdi0_447 = buffer.data(pdi0 + 447);
    const auto *pdi0_455 = buffer.data(pdi0 + 455);
    const auto *pdi0_459 = buffer.data(pdi0 + 459);
    const auto *pdi0_460 = buffer.data(pdi0 + 460);
    const auto *pdi0_469 = buffer.data(pdi0 + 469);
    const auto *pdi0_470 = buffer.data(pdi0 + 470);
    const auto *pdi0_471 = buffer.data(pdi0 + 471);
    const auto *pdi0_472 = buffer.data(pdi0 + 472);
    const auto *pdi0_473 = buffer.data(pdi0 + 473);
    const auto *pdi0_475 = buffer.data(pdi0 + 475);
    const auto *pdi0_476 = buffer.data(pdi0 + 476);
    const auto *pdi0_478 = buffer.data(pdi0 + 478);
    const auto *pdi0_479 = buffer.data(pdi0 + 479);
    const auto *pdi0_481 = buffer.data(pdi0 + 481);
    const auto *pdi0_482 = buffer.data(pdi0 + 482);
    const auto *pdi0_483 = buffer.data(pdi0 + 483);
    const auto *pdi0_485 = buffer.data(pdi0 + 485);
    const auto *pdi0_486 = buffer.data(pdi0 + 486);
    const auto *pdi0_487 = buffer.data(pdi0 + 487);
    const auto *pdi0_488 = buffer.data(pdi0 + 488);
    const auto *pdi0_490 = buffer.data(pdi0 + 490);
    const auto *pdi0_497 = buffer.data(pdi0 + 497);
    const auto *pdi0_498 = buffer.data(pdi0 + 498);
    const auto *pdi0_499 = buffer.data(pdi0 + 499);
    const auto *pdi0_500 = buffer.data(pdi0 + 500);
    const auto *pdi0_501 = buffer.data(pdi0 + 501);

    const auto *pdh_272 = buffer.data(pdh + 272);
    const auto *pdh_273 = buffer.data(pdh + 273);
    const auto *pdh_275 = buffer.data(pdh + 275);
    const auto *pdh_278 = buffer.data(pdh + 278);
    const auto *pdh_282 = buffer.data(pdh + 282);
    const auto *pdh_293 = buffer.data(pdh + 293);
    const auto *pdh_294 = buffer.data(pdh + 294);
    const auto *pdh_296 = buffer.data(pdh + 296);
    const auto *pdh_297 = buffer.data(pdh + 297);
    const auto *pdh_299 = buffer.data(pdh + 299);
    const auto *pdh_300 = buffer.data(pdh + 300);
    const auto *pdh_301 = buffer.data(pdh + 301);
    const auto *pdh_303 = buffer.data(pdh + 303);
    const auto *pdh_304 = buffer.data(pdh + 304);
    const auto *pdh_305 = buffer.data(pdh + 305);
    const auto *pdh_306 = buffer.data(pdh + 306);
    const auto *pdh_308 = buffer.data(pdh + 308);
    const auto *pdh_309 = buffer.data(pdh + 309);
    const auto *pdh_310 = buffer.data(pdh + 310);
    const auto *pdh_311 = buffer.data(pdh + 311);
    const auto *pdh_312 = buffer.data(pdh + 312);
    const auto *pdh_313 = buffer.data(pdh + 313);
    const auto *pdh_314 = buffer.data(pdh + 314);
    const auto *pdh_330 = buffer.data(pdh + 330);
    const auto *pdh_331 = buffer.data(pdh + 331);
    const auto *pdh_332 = buffer.data(pdh + 332);
    const auto *pdh_333 = buffer.data(pdh + 333);
    const auto *pdh_334 = buffer.data(pdh + 334);
    const auto *pdh_335 = buffer.data(pdh + 335);
    const auto *pdh_339 = buffer.data(pdh + 339);
    const auto *pdh_342 = buffer.data(pdh + 342);
    const auto *pdh_343 = buffer.data(pdh + 343);
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

    const auto *pdi1_363 = buffer.data(pdi1 + 363);
    const auto *pdi1_392 = buffer.data(pdi1 + 392);
    const auto *pdi1_394 = buffer.data(pdi1 + 394);
    const auto *pdi1_397 = buffer.data(pdi1 + 397);
    const auto *pdi1_401 = buffer.data(pdi1 + 401);
    const auto *pdi1_406 = buffer.data(pdi1 + 406);
    const auto *pdi1_419 = buffer.data(pdi1 + 419);
    const auto *pdi1_442 = buffer.data(pdi1 + 442);
    const auto *pdi1_443 = buffer.data(pdi1 + 443);
    const auto *pdi1_444 = buffer.data(pdi1 + 444);
    const auto *pdi1_445 = buffer.data(pdi1 + 445);
    const auto *pdi1_447 = buffer.data(pdi1 + 447);
    const auto *pdi1_455 = buffer.data(pdi1 + 455);
    const auto *pdi1_459 = buffer.data(pdi1 + 459);
    const auto *pdi1_460 = buffer.data(pdi1 + 460);
    const auto *pdi1_469 = buffer.data(pdi1 + 469);
    const auto *pdi1_470 = buffer.data(pdi1 + 470);
    const auto *pdi1_471 = buffer.data(pdi1 + 471);
    const auto *pdi1_472 = buffer.data(pdi1 + 472);
    const auto *pdi1_473 = buffer.data(pdi1 + 473);
    const auto *pdi1_475 = buffer.data(pdi1 + 475);
    const auto *pdi1_476 = buffer.data(pdi1 + 476);
    const auto *pdi1_478 = buffer.data(pdi1 + 478);
    const auto *pdi1_479 = buffer.data(pdi1 + 479);
    const auto *pdi1_481 = buffer.data(pdi1 + 481);
    const auto *pdi1_482 = buffer.data(pdi1 + 482);
    const auto *pdi1_483 = buffer.data(pdi1 + 483);
    const auto *pdi1_485 = buffer.data(pdi1 + 485);
    const auto *pdi1_486 = buffer.data(pdi1 + 486);
    const auto *pdi1_487 = buffer.data(pdi1 + 487);
    const auto *pdi1_488 = buffer.data(pdi1 + 488);
    const auto *pdi1_490 = buffer.data(pdi1 + 490);
    const auto *pdi1_497 = buffer.data(pdi1 + 497);
    const auto *pdi1_498 = buffer.data(pdi1 + 498);
    const auto *pdi1_499 = buffer.data(pdi1 + 499);
    const auto *pdi1_500 = buffer.data(pdi1 + 500);
    const auto *pdi1_501 = buffer.data(pdi1 + 501);

    const auto *pfg0_330 = buffer.data(pfg0 + 330);
    const auto *pfg0_332 = buffer.data(pfg0 + 332);
    const auto *pfg0_333 = buffer.data(pfg0 + 333);
    const auto *pfg0_335 = buffer.data(pfg0 + 335);
    const auto *pfg0_336 = buffer.data(pfg0 + 336);
    const auto *pfg0_337 = buffer.data(pfg0 + 337);
    const auto *pfg0_339 = buffer.data(pfg0 + 339);
    const auto *pfg0_340 = buffer.data(pfg0 + 340);
    const auto *pfg0_341 = buffer.data(pfg0 + 341);
    const auto *pfg0_342 = buffer.data(pfg0 + 342);
    const auto *pfg0_343 = buffer.data(pfg0 + 343);
    const auto *pfg0_344 = buffer.data(pfg0 + 344);
    const auto *pfg0_363 = buffer.data(pfg0 + 363);
    const auto *pfg0_366 = buffer.data(pfg0 + 366);
    const auto *pfg0_370 = buffer.data(pfg0 + 370);

    const auto *pfg1_330 = buffer.data(pfg1 + 330);
    const auto *pfg1_332 = buffer.data(pfg1 + 332);
    const auto *pfg1_333 = buffer.data(pfg1 + 333);
    const auto *pfg1_335 = buffer.data(pfg1 + 335);
    const auto *pfg1_336 = buffer.data(pfg1 + 336);
    const auto *pfg1_337 = buffer.data(pfg1 + 337);
    const auto *pfg1_339 = buffer.data(pfg1 + 339);
    const auto *pfg1_340 = buffer.data(pfg1 + 340);
    const auto *pfg1_341 = buffer.data(pfg1 + 341);
    const auto *pfg1_342 = buffer.data(pfg1 + 342);
    const auto *pfg1_343 = buffer.data(pfg1 + 343);
    const auto *pfg1_344 = buffer.data(pfg1 + 344);
    const auto *pfg1_363 = buffer.data(pfg1 + 363);
    const auto *pfg1_366 = buffer.data(pfg1 + 366);
    const auto *pfg1_370 = buffer.data(pfg1 + 370);

    const auto *pfh_461 = buffer.data(pfh + 461);
    const auto *pfh_462 = buffer.data(pfh + 462);
    const auto *pfh_464 = buffer.data(pfh + 464);
    const auto *pfh_465 = buffer.data(pfh + 465);
    const auto *pfh_467 = buffer.data(pfh + 467);
    const auto *pfh_468 = buffer.data(pfh + 468);
    const auto *pfh_469 = buffer.data(pfh + 469);
    const auto *pfh_471 = buffer.data(pfh + 471);
    const auto *pfh_472 = buffer.data(pfh + 472);
    const auto *pfh_473 = buffer.data(pfh + 473);
    const auto *pfh_474 = buffer.data(pfh + 474);
    const auto *pfh_476 = buffer.data(pfh + 476);
    const auto *pfh_477 = buffer.data(pfh + 477);
    const auto *pfh_478 = buffer.data(pfh + 478);
    const auto *pfh_479 = buffer.data(pfh + 479);
    const auto *pfh_480 = buffer.data(pfh + 480);
    const auto *pfh_481 = buffer.data(pfh + 481);
    const auto *pfh_482 = buffer.data(pfh + 482);
    const auto *pfh_483 = buffer.data(pfh + 483);
    const auto *pfh_485 = buffer.data(pfh + 485);
    const auto *pfh_488 = buffer.data(pfh + 488);
    const auto *pfh_492 = buffer.data(pfh + 492);
    const auto *pfh_498 = buffer.data(pfh + 498);
    const auto *pfh_499 = buffer.data(pfh + 499);
    const auto *pfh_500 = buffer.data(pfh + 500);
    const auto *pfh_501 = buffer.data(pfh + 501);
    const auto *pfh_502 = buffer.data(pfh + 502);
    const auto *pfh_503 = buffer.data(pfh + 503);
    const auto *pfh_504 = buffer.data(pfh + 504);
    const auto *pfh_506 = buffer.data(pfh + 506);
    const auto *pfh_507 = buffer.data(pfh + 507);
    const auto *pfh_509 = buffer.data(pfh + 509);
    const auto *pfh_510 = buffer.data(pfh + 510);
    const auto *pfh_513 = buffer.data(pfh + 513);
    const auto *pfh_514 = buffer.data(pfh + 514);
    const auto *pfh_519 = buffer.data(pfh + 519);
    const auto *pfh_520 = buffer.data(pfh + 520);
    const auto *pfh_521 = buffer.data(pfh + 521);
    const auto *pfh_522 = buffer.data(pfh + 522);
    const auto *pfh_523 = buffer.data(pfh + 523);
    const auto *pfh_524 = buffer.data(pfh + 524);
    const auto *pfh_525 = buffer.data(pfh + 525);
    const auto *pfh_527 = buffer.data(pfh + 527);
    const auto *pfh_530 = buffer.data(pfh + 530);
    const auto *pfh_534 = buffer.data(pfh + 534);
    const auto *pfh_540 = buffer.data(pfh + 540);
    const auto *pfh_541 = buffer.data(pfh + 541);
    const auto *pfh_542 = buffer.data(pfh + 542);
    const auto *pfh_543 = buffer.data(pfh + 543);
    const auto *pfh_544 = buffer.data(pfh + 544);
    const auto *pfh_545 = buffer.data(pfh + 545);

#pragma omp simd aligned(t_610, t_611, t_612, pa_z, pc_z, sfi0_50, sfi0_51, sfi0_52, sfh_36, \
                         sfh_37, sfh_38, sfi1_50, sfi1_51, sfi1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_z[k] * sfi0_50[k]
                   + f_0 * sfh_36[k]
                   - f_11 * pc_z[k] * sfi1_50[k];

        t_611[k] = pa_z[k] * sfi0_51[k]
                   + f_12 * sfh_37[k]
                   - f_11 * pc_z[k] * sfi1_51[k];

        t_612[k] = pa_z[k] * sfi0_52[k]
                   + f_1 * sfh_38[k]
                   - f_11 * pc_z[k] * sfi1_52[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pa_z, pb_y, pc_y, pc_z, sfi0_53, sfh_39, \
                         sfi1_53, pdi0_363, pdh_272, pdi1_363, \
                         pfh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = pa_z[k] * sfi0_53[k]
                   + f_13 * sfh_39[k]
                   - f_11 * pc_z[k] * sfi1_53[k];

        t_614[k] = f_0 * pdh_272[k]
                   + f_4 * pc_y[k] * pfh_461[k];

        t_615[k] = pb_y[k] * pdi0_363[k]
                   - f_11 * pc_y[k] * pdi1_363[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, pc_x, pc_y, pdh_294, pdh_296, pfg0_330, \
                         pfg0_332, pfg1_330, pfg1_332, pfh_462, \
                         pfh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = f_12 * pdh_294[k]
                   + f_2 * pfg0_330[k]
                   - f_3 * pfg1_330[k]
                   + f_4 * pc_x[k] * pfh_462[k];

        t_617[k] = f_4 * pc_y[k] * pfh_462[k];

        t_618[k] = f_12 * pdh_296[k]
                   + f_15 * pfg0_332[k]
                   - f_16 * pfg1_332[k]
                   + f_4 * pc_x[k] * pfh_464[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, pc_x, pc_y, pdh_297, pdh_299, pfg0_333, \
                         pfg0_335, pfg1_333, pfg1_335, pfh_464, pfh_465, \
                         pfh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_12 * pdh_297[k]
                   + f_9 * pfg0_333[k]
                   - f_10 * pfg1_333[k]
                   + f_4 * pc_x[k] * pfh_465[k];

        t_620[k] = f_4 * pc_y[k] * pfh_464[k];

        t_621[k] = f_12 * pdh_299[k]
                   + f_9 * pfg0_335[k]
                   - f_10 * pfg1_335[k]
                   + f_4 * pc_x[k] * pfh_467[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, pc_x, pc_y, pdh_300, pdh_301, pfg0_336, \
                         pfg0_337, pfg1_336, pfg1_337, pfh_467, pfh_468, \
                         pfh_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_12 * pdh_300[k]
                   + f_7 * pfg0_336[k]
                   - f_8 * pfg1_336[k]
                   + f_4 * pc_x[k] * pfh_468[k];

        t_623[k] = f_12 * pdh_301[k]
                   + f_7 * pfg0_337[k]
                   - f_8 * pfg1_337[k]
                   + f_4 * pc_x[k] * pfh_469[k];

        t_624[k] = f_4 * pc_y[k] * pfh_467[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_x, pdh_303, pdh_304, pdh_305, pfg0_339, \
                         pfg0_340, pfg0_341, pfg1_339, pfg1_340, pfg1_341, pfh_471, pfh_472, \
                         pfh_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_12 * pdh_303[k]
                   + f_7 * pfg0_339[k]
                   - f_8 * pfg1_339[k]
                   + f_4 * pc_x[k] * pfh_471[k];

        t_626[k] = f_12 * pdh_304[k]
                   + f_5 * pfg0_340[k]
                   - f_6 * pfg1_340[k]
                   + f_4 * pc_x[k] * pfh_472[k];

        t_627[k] = f_12 * pdh_305[k]
                   + f_5 * pfg0_341[k]
                   - f_6 * pfg1_341[k]
                   + f_4 * pc_x[k] * pfh_473[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pdh_306, pdh_308, pfg0_342, \
                         pfg0_344, pfg1_342, pfg1_344, pfh_471, pfh_474, \
                         pfh_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_12 * pdh_306[k]
                   + f_5 * pfg0_342[k]
                   - f_6 * pfg1_342[k]
                   + f_4 * pc_x[k] * pfh_474[k];

        t_629[k] = f_4 * pc_y[k] * pfh_471[k];

        t_630[k] = f_12 * pdh_308[k]
                   + f_5 * pfg0_344[k]
                   - f_6 * pfg1_344[k]
                   + f_4 * pc_x[k] * pfh_476[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, t_635, pc_x, pdh_309, pdh_310, pdh_311, \
                         pdh_312, pdh_313, pfh_477, pfh_478, pfh_479, pfh_480, \
                         pfh_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_12 * pdh_309[k]
                   + f_4 * pc_x[k] * pfh_477[k];

        t_632[k] = f_12 * pdh_310[k]
                   + f_4 * pc_x[k] * pfh_478[k];

        t_633[k] = f_12 * pdh_311[k]
                   + f_4 * pc_x[k] * pfh_479[k];

        t_634[k] = f_12 * pdh_312[k]
                   + f_4 * pc_x[k] * pfh_480[k];

        t_635[k] = f_12 * pdh_313[k]
                   + f_4 * pc_x[k] * pfh_481[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_x, pc_y, pdh_314, pfg0_340, pfg0_341, \
                         pfg1_340, pfg1_341, pfh_477, pfh_478, \
                         pfh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_12 * pdh_314[k]
                   + f_4 * pc_x[k] * pfh_482[k];

        t_637[k] = f_2 * pfg0_340[k]
                   - f_3 * pfg1_340[k]
                   + f_4 * pc_y[k] * pfh_477[k];

        t_638[k] = f_15 * pfg0_341[k]
                   - f_16 * pfg1_341[k]
                   + f_4 * pc_y[k] * pfh_478[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pc_y, pfg0_342, pfg0_343, pfg0_344, \
                         pfg1_342, pfg1_343, pfg1_344, pfh_479, pfh_480, pfh_481, \
                         pfh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_9 * pfg0_342[k]
                   - f_10 * pfg1_342[k]
                   + f_4 * pc_y[k] * pfh_479[k];

        t_640[k] = f_7 * pfg0_343[k]
                   - f_8 * pfg1_343[k]
                   + f_4 * pc_y[k] * pfh_480[k];

        t_641[k] = f_5 * pfg0_344[k]
                   - f_6 * pfg1_344[k]
                   + f_4 * pc_y[k] * pfh_481[k];

        t_642[k] = f_4 * pc_y[k] * pfh_482[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pa_z, pb_x, pc_x, pc_y, pc_z, sfi0_84, sfi1_84, \
                         ppi0_251, ppi1_251, pdi0_419, pdh_273, pdi1_419, \
                         pfh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_17 * ppi0_251[k]
                   - f_18 * ppi1_251[k]
                   + pb_x[k] * pdi0_419[k]
                   - f_11 * pc_x[k] * pdi1_419[k];

        t_644[k] = pa_z[k] * sfi0_84[k]
                   - f_11 * pc_z[k] * sfi1_84[k];

        t_645[k] = f_12 * pdh_273[k]
                   + f_4 * pc_y[k] * pfh_483[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pa_z, pc_y, pc_z, sfi0_86, sfi0_87, sfh_63, \
                         sfi1_86, sfi1_87, pdh_275, pfh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = pa_z[k] * sfi0_86[k]
                   + f_0 * sfh_63[k]
                   - f_11 * pc_z[k] * sfi1_86[k];

        t_647[k] = pa_z[k] * sfi0_87[k]
                   - f_11 * pc_z[k] * sfi1_87[k];

        t_648[k] = f_12 * pdh_275[k]
                   + f_4 * pc_y[k] * pfh_485[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pa_z, pc_z, sfi0_89, sfi0_90, sfi0_91, sfh_65, \
                         sfh_66, sfi1_89, sfi1_90, sfi1_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pa_z[k] * sfi0_89[k]
                   + f_12 * sfh_65[k]
                   - f_11 * pc_z[k] * sfi1_89[k];

        t_650[k] = pa_z[k] * sfi0_90[k]
                   - f_11 * pc_z[k] * sfi1_90[k];

        t_651[k] = pa_z[k] * sfi0_91[k]
                   + f_0 * sfh_66[k]
                   - f_11 * pc_z[k] * sfi1_91[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pa_z, pc_y, pc_z, sfi0_93, sfi0_94, sfh_68, \
                         sfi1_93, sfi1_94, pdh_278, pfh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_12 * pdh_278[k]
                   + f_4 * pc_y[k] * pfh_488[k];

        t_653[k] = pa_z[k] * sfi0_93[k]
                   + f_1 * sfh_68[k]
                   - f_11 * pc_z[k] * sfi1_93[k];

        t_654[k] = pa_z[k] * sfi0_94[k]
                   - f_11 * pc_z[k] * sfi1_94[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pa_z, pc_y, pc_z, sfi0_95, sfi0_96, sfh_69, \
                         sfh_70, sfi1_95, sfi1_96, pdh_282, pfh_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = pa_z[k] * sfi0_95[k]
                   + f_0 * sfh_69[k]
                   - f_11 * pc_z[k] * sfi1_95[k];

        t_656[k] = pa_z[k] * sfi0_96[k]
                   + f_12 * sfh_70[k]
                   - f_11 * pc_z[k] * sfi1_96[k];

        t_657[k] = f_12 * pdh_282[k]
                   + f_4 * pc_y[k] * pfh_492[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pa_z, pc_x, pc_z, sfi0_98, sfh_72, \
                         sfi1_98, pdh_330, pdh_331, pdh_332, pfh_498, pfh_499, \
                         pfh_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = pa_z[k] * sfi0_98[k]
                   + f_13 * sfh_72[k]
                   - f_11 * pc_z[k] * sfi1_98[k];

        t_659[k] = f_0 * pdh_330[k]
                   + f_4 * pc_x[k] * pfh_498[k];

        t_660[k] = f_0 * pdh_331[k]
                   + f_4 * pc_x[k] * pfh_499[k];

        t_661[k] = f_0 * pdh_332[k]
                   + f_4 * pc_x[k] * pfh_500[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pa_z, pc_x, pc_z, sfi0_105, sfi1_105, \
                         pdh_333, pdh_334, pdh_335, pfh_501, pfh_502, \
                         pfh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_0 * pdh_333[k]
                   + f_4 * pc_x[k] * pfh_501[k];

        t_663[k] = f_0 * pdh_334[k]
                   + f_4 * pc_x[k] * pfh_502[k];

        t_664[k] = f_0 * pdh_335[k]
                   + f_4 * pc_x[k] * pfh_503[k];

        t_665[k] = pa_z[k] * sfi0_105[k]
                   - f_11 * pc_z[k] * sfi1_105[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_x, pc_x, pdi0_442, pdi0_443, pdi0_444, \
                         pdi0_445, pdi1_442, pdi1_443, pdi1_444, \
                         pdi1_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = pb_x[k] * pdi0_442[k]
                   - f_11 * pc_x[k] * pdi1_442[k];

        t_667[k] = pb_x[k] * pdi0_443[k]
                   - f_11 * pc_x[k] * pdi1_443[k];

        t_668[k] = pb_x[k] * pdi0_444[k]
                   - f_11 * pc_x[k] * pdi1_444[k];

        t_669[k] = pb_x[k] * pdi0_445[k]
                   - f_11 * pc_x[k] * pdi1_445[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_x, pb_y, pc_x, pc_y, pdi0_392, \
                         pdi0_447, pdh_293, pdh_294, pdi1_392, pdi1_447, pfh_503, \
                         pfh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_12 * pdh_293[k]
                   + f_4 * pc_y[k] * pfh_503[k];

        t_671[k] = pb_x[k] * pdi0_447[k]
                   - f_11 * pc_x[k] * pdi1_447[k];

        t_672[k] = pb_y[k] * pdi0_392[k]
                   - f_11 * pc_y[k] * pdi1_392[k];

        t_673[k] = f_0 * pdh_294[k]
                   + f_4 * pc_y[k] * pfh_504[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pb_y, pc_x, pc_y, pdi0_394, pdh_296, pdh_339, \
                         pdi1_394, pfg0_363, pfg1_363, pfh_506, \
                         pfh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = pb_y[k] * pdi0_394[k]
                   - f_11 * pc_y[k] * pdi1_394[k];

        t_675[k] = f_0 * pdh_339[k]
                   + f_9 * pfg0_363[k]
                   - f_10 * pfg1_363[k]
                   + f_4 * pc_x[k] * pfh_507[k];

        t_676[k] = f_0 * pdh_296[k]
                   + f_4 * pc_y[k] * pfh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pb_x, pb_y, pc_x, pc_y, pdi0_397, pdi0_455, \
                         pdh_342, pdh_343, pdi1_397, pdi1_455, pfg0_366, pfg1_366, \
                         pfh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = pb_y[k] * pdi0_397[k]
                   - f_11 * pc_y[k] * pdi1_397[k];

        t_678[k] = f_0 * pdh_342[k]
                   + f_7 * pfg0_366[k]
                   - f_8 * pfg1_366[k]
                   + f_4 * pc_x[k] * pfh_510[k];

        t_679[k] = pb_x[k] * pdi0_455[k]
                   + f_1 * pdh_343[k]
                   - f_11 * pc_x[k] * pdi1_455[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pb_y, pc_x, pc_y, pdi0_401, pdh_299, pdh_346, \
                         pdi1_401, pfg0_370, pfg1_370, pfh_509, \
                         pfh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_0 * pdh_299[k]
                   + f_4 * pc_y[k] * pfh_509[k];

        t_681[k] = pb_y[k] * pdi0_401[k]
                   - f_11 * pc_y[k] * pdi1_401[k];

        t_682[k] = f_0 * pdh_346[k]
                   + f_5 * pfg0_370[k]
                   - f_6 * pfg1_370[k]
                   + f_4 * pc_x[k] * pfh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pb_x, pc_x, pc_y, pdi0_459, pdi0_460, pdh_303, \
                         pdh_347, pdh_348, pdi1_459, pdi1_460, \
                         pfh_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = pb_x[k] * pdi0_459[k]
                   + f_12 * pdh_347[k]
                   - f_11 * pc_x[k] * pdi1_459[k];

        t_684[k] = pb_x[k] * pdi0_460[k]
                   + f_12 * pdh_348[k]
                   - f_11 * pc_x[k] * pdi1_460[k];

        t_685[k] = f_0 * pdh_303[k]
                   + f_4 * pc_y[k] * pfh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pb_y, pc_x, pc_y, pdi0_406, pdh_351, \
                         pdh_352, pdh_353, pdi1_406, pfh_519, pfh_520, \
                         pfh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pb_y[k] * pdi0_406[k]
                   - f_11 * pc_y[k] * pdi1_406[k];

        t_687[k] = f_0 * pdh_351[k]
                   + f_4 * pc_x[k] * pfh_519[k];

        t_688[k] = f_0 * pdh_352[k]
                   + f_4 * pc_x[k] * pfh_520[k];

        t_689[k] = f_0 * pdh_353[k]
                   + f_4 * pc_x[k] * pfh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pb_x, pc_x, pdi0_469, pdh_354, pdh_355, \
                         pdh_356, pdi1_469, pfh_522, pfh_523, pfh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_0 * pdh_354[k]
                   + f_4 * pc_x[k] * pfh_522[k];

        t_691[k] = f_0 * pdh_355[k]
                   + f_4 * pc_x[k] * pfh_523[k];

        t_692[k] = f_0 * pdh_356[k]
                   + f_4 * pc_x[k] * pfh_524[k];

        t_693[k] = pb_x[k] * pdi0_469[k]
                   - f_11 * pc_x[k] * pdi1_469[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pb_x, pc_x, pdi0_470, pdi0_471, pdi0_472, \
                         pdi0_473, pdi1_470, pdi1_471, pdi1_472, \
                         pdi1_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pb_x[k] * pdi0_470[k]
                   - f_11 * pc_x[k] * pdi1_470[k];

        t_695[k] = pb_x[k] * pdi0_471[k]
                   - f_11 * pc_x[k] * pdi1_471[k];

        t_696[k] = pb_x[k] * pdi0_472[k]
                   - f_11 * pc_x[k] * pdi1_472[k];

        t_697[k] = pb_x[k] * pdi0_473[k]
                   - f_11 * pc_x[k] * pdi1_473[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, t_701, pb_x, pc_x, pc_y, pdi0_475, pdi0_476, \
                         pdh_314, pdh_357, pdi1_475, pdi1_476, pfh_524, \
                         pfh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_0 * pdh_314[k]
                   + f_4 * pc_y[k] * pfh_524[k];

        t_699[k] = pb_x[k] * pdi0_475[k]
                   - f_11 * pc_x[k] * pdi1_475[k];

        t_700[k] = pb_x[k] * pdi0_476[k]
                   + f_14 * pdh_357[k]
                   - f_11 * pc_x[k] * pdi1_476[k];

        t_701[k] = f_4 * pc_y[k] * pfh_525[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, pb_x, pc_x, pc_y, pdi0_478, pdi0_479, pdh_359, \
                         pdh_360, pdi1_478, pdi1_479, pfh_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = pb_x[k] * pdi0_478[k]
                   + f_19 * pdh_359[k]
                   - f_11 * pc_x[k] * pdi1_478[k];

        t_703[k] = pb_x[k] * pdi0_479[k]
                   + f_13 * pdh_360[k]
                   - f_11 * pc_x[k] * pdi1_479[k];

        t_704[k] = f_4 * pc_y[k] * pfh_527[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, pb_x, pc_x, pdi0_481, pdi0_482, pdi0_483, \
                         pdh_362, pdh_363, pdh_364, pdi1_481, pdi1_482, \
                         pdi1_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = pb_x[k] * pdi0_481[k]
                   + f_13 * pdh_362[k]
                   - f_11 * pc_x[k] * pdi1_481[k];

        t_706[k] = pb_x[k] * pdi0_482[k]
                   + f_1 * pdh_363[k]
                   - f_11 * pc_x[k] * pdi1_482[k];

        t_707[k] = pb_x[k] * pdi0_483[k]
                   + f_1 * pdh_364[k]
                   - f_11 * pc_x[k] * pdi1_483[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, pb_x, pc_x, pc_y, pdi0_485, pdi0_486, pdh_366, \
                         pdh_367, pdi1_485, pdi1_486, pfh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_4 * pc_y[k] * pfh_530[k];

        t_709[k] = pb_x[k] * pdi0_485[k]
                   + f_1 * pdh_366[k]
                   - f_11 * pc_x[k] * pdi1_485[k];

        t_710[k] = pb_x[k] * pdi0_486[k]
                   + f_12 * pdh_367[k]
                   - f_11 * pc_x[k] * pdi1_486[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pb_x, pc_x, pc_y, pdi0_487, pdi0_488, pdh_368, \
                         pdh_369, pdi1_487, pdi1_488, pfh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = pb_x[k] * pdi0_487[k]
                   + f_12 * pdh_368[k]
                   - f_11 * pc_x[k] * pdi1_487[k];

        t_712[k] = pb_x[k] * pdi0_488[k]
                   + f_12 * pdh_369[k]
                   - f_11 * pc_x[k] * pdi1_488[k];

        t_713[k] = f_4 * pc_y[k] * pfh_534[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pb_x, pc_x, pdi0_490, pdh_371, pdh_372, \
                         pdh_373, pdh_374, pdi1_490, pfh_540, pfh_541, \
                         pfh_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_x[k] * pdi0_490[k]
                   + f_12 * pdh_371[k]
                   - f_11 * pc_x[k] * pdi1_490[k];

        t_715[k] = f_0 * pdh_372[k]
                   + f_4 * pc_x[k] * pfh_540[k];

        t_716[k] = f_0 * pdh_373[k]
                   + f_4 * pc_x[k] * pfh_541[k];

        t_717[k] = f_0 * pdh_374[k]
                   + f_4 * pc_x[k] * pfh_542[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pb_x, pc_x, pdi0_497, pdh_375, pdh_376, \
                         pdh_377, pdi1_497, pfh_543, pfh_544, pfh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_0 * pdh_375[k]
                   + f_4 * pc_x[k] * pfh_543[k];

        t_719[k] = f_0 * pdh_376[k]
                   + f_4 * pc_x[k] * pfh_544[k];

        t_720[k] = f_0 * pdh_377[k]
                   + f_4 * pc_x[k] * pfh_545[k];

        t_721[k] = pb_x[k] * pdi0_497[k]
                   - f_11 * pc_x[k] * pdi1_497[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, pb_x, pc_x, pdi0_498, pdi0_499, pdi0_500, \
                         pdi0_501, pdi1_498, pdi1_499, pdi1_500, \
                         pdi1_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = pb_x[k] * pdi0_498[k]
                   - f_11 * pc_x[k] * pdi1_498[k];

        t_723[k] = pb_x[k] * pdi0_499[k]
                   - f_11 * pc_x[k] * pdi1_499[k];

        t_724[k] = pb_x[k] * pdi0_500[k]
                   - f_11 * pc_x[k] * pdi1_500[k];

        t_725[k] = pb_x[k] * pdi0_501[k]
                   - f_11 * pc_x[k] * pdi1_501[k];
    }
}

static auto
compute_prim_pfi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfi0, const size_t sfh,
                                                          const size_t sfi1, const size_t ppi0,
                                                          const size_t ppi1, const size_t pdi0,
                                                          const size_t pdh, const size_t pdi1,
                                                          const size_t pfg0, const size_t pfg1,
                                                          const size_t pfh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfi0_168 = buffer.data(sfi0 + 168);
    const auto *sfi0_171 = buffer.data(sfi0 + 171);
    const auto *sfi0_174 = buffer.data(sfi0 + 174);
    const auto *sfi0_178 = buffer.data(sfi0 + 178);
    const auto *sfi0_189 = buffer.data(sfi0 + 189);
    const auto *sfi0_190 = buffer.data(sfi0 + 190);
    const auto *sfi0_191 = buffer.data(sfi0 + 191);
    const auto *sfi0_192 = buffer.data(sfi0 + 192);
    const auto *sfi0_193 = buffer.data(sfi0 + 193);
    const auto *sfi0_195 = buffer.data(sfi0 + 195);

    const auto *sfh_141 = buffer.data(sfh + 141);
    const auto *sfh_142 = buffer.data(sfh + 142);
    const auto *sfh_143 = buffer.data(sfh + 143);
    const auto *sfh_144 = buffer.data(sfh + 144);
    const auto *sfh_146 = buffer.data(sfh + 146);
    const auto *sfh_209 = buffer.data(sfh + 209);

    const auto *sfi1_168 = buffer.data(sfi1 + 168);
    const auto *sfi1_171 = buffer.data(sfi1 + 171);
    const auto *sfi1_174 = buffer.data(sfi1 + 174);
    const auto *sfi1_178 = buffer.data(sfi1 + 178);
    const auto *sfi1_189 = buffer.data(sfi1 + 189);
    const auto *sfi1_190 = buffer.data(sfi1 + 190);
    const auto *sfi1_191 = buffer.data(sfi1 + 191);
    const auto *sfi1_192 = buffer.data(sfi1 + 192);
    const auto *sfi1_193 = buffer.data(sfi1 + 193);
    const auto *sfi1_195 = buffer.data(sfi1 + 195);

    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *pdi0_475 = buffer.data(pdi0 + 475);
    const auto *pdi0_476 = buffer.data(pdi0 + 476);
    const auto *pdi0_478 = buffer.data(pdi0 + 478);
    const auto *pdi0_481 = buffer.data(pdi0 + 481);
    const auto *pdi0_485 = buffer.data(pdi0 + 485);
    const auto *pdi0_490 = buffer.data(pdi0 + 490);
    const auto *pdi0_497 = buffer.data(pdi0 + 497);
    const auto *pdi0_498 = buffer.data(pdi0 + 498);
    const auto *pdi0_499 = buffer.data(pdi0 + 499);
    const auto *pdi0_500 = buffer.data(pdi0 + 500);
    const auto *pdi0_501 = buffer.data(pdi0 + 501);
    const auto *pdi0_503 = buffer.data(pdi0 + 503);

    const auto *pdh_315 = buffer.data(pdh + 315);
    const auto *pdh_317 = buffer.data(pdh + 317);
    const auto *pdh_320 = buffer.data(pdh + 320);
    const auto *pdh_324 = buffer.data(pdh + 324);
    const auto *pdh_335 = buffer.data(pdh + 335);
    const auto *pdh_336 = buffer.data(pdh + 336);
    const auto *pdh_338 = buffer.data(pdh + 338);
    const auto *pdh_341 = buffer.data(pdh + 341);
    const auto *pdh_345 = buffer.data(pdh + 345);
    const auto *pdh_351 = buffer.data(pdh + 351);
    const auto *pdh_352 = buffer.data(pdh + 352);
    const auto *pdh_353 = buffer.data(pdh + 353);
    const auto *pdh_354 = buffer.data(pdh + 354);
    const auto *pdh_355 = buffer.data(pdh + 355);
    const auto *pdh_356 = buffer.data(pdh + 356);
    const auto *pdh_357 = buffer.data(pdh + 357);
    const auto *pdh_359 = buffer.data(pdh + 359);
    const auto *pdh_362 = buffer.data(pdh + 362);
    const auto *pdh_366 = buffer.data(pdh + 366);
    const auto *pdh_372 = buffer.data(pdh + 372);
    const auto *pdh_373 = buffer.data(pdh + 373);
    const auto *pdh_374 = buffer.data(pdh + 374);
    const auto *pdh_375 = buffer.data(pdh + 375);
    const auto *pdh_376 = buffer.data(pdh + 376);
    const auto *pdh_377 = buffer.data(pdh + 377);

    const auto *pdi1_475 = buffer.data(pdi1 + 475);
    const auto *pdi1_476 = buffer.data(pdi1 + 476);
    const auto *pdi1_478 = buffer.data(pdi1 + 478);
    const auto *pdi1_481 = buffer.data(pdi1 + 481);
    const auto *pdi1_485 = buffer.data(pdi1 + 485);
    const auto *pdi1_490 = buffer.data(pdi1 + 490);
    const auto *pdi1_497 = buffer.data(pdi1 + 497);
    const auto *pdi1_498 = buffer.data(pdi1 + 498);
    const auto *pdi1_499 = buffer.data(pdi1 + 499);
    const auto *pdi1_500 = buffer.data(pdi1 + 500);
    const auto *pdi1_501 = buffer.data(pdi1 + 501);
    const auto *pdi1_503 = buffer.data(pdi1 + 503);

    const auto *pfg0_392 = buffer.data(pfg0 + 392);
    const auto *pfg0_395 = buffer.data(pfg0 + 395);
    const auto *pfg0_397 = buffer.data(pfg0 + 397);
    const auto *pfg0_399 = buffer.data(pfg0 + 399);
    const auto *pfg0_401 = buffer.data(pfg0 + 401);
    const auto *pfg0_402 = buffer.data(pfg0 + 402);
    const auto *pfg0_404 = buffer.data(pfg0 + 404);
    const auto *pfg0_405 = buffer.data(pfg0 + 405);
    const auto *pfg0_407 = buffer.data(pfg0 + 407);
    const auto *pfg0_408 = buffer.data(pfg0 + 408);
    const auto *pfg0_410 = buffer.data(pfg0 + 410);
    const auto *pfg0_411 = buffer.data(pfg0 + 411);
    const auto *pfg0_412 = buffer.data(pfg0 + 412);
    const auto *pfg0_414 = buffer.data(pfg0 + 414);
    const auto *pfg0_415 = buffer.data(pfg0 + 415);
    const auto *pfg0_416 = buffer.data(pfg0 + 416);
    const auto *pfg0_417 = buffer.data(pfg0 + 417);
    const auto *pfg0_418 = buffer.data(pfg0 + 418);
    const auto *pfg0_419 = buffer.data(pfg0 + 419);
    const auto *pfg0_423 = buffer.data(pfg0 + 423);
    const auto *pfg0_426 = buffer.data(pfg0 + 426);
    const auto *pfg0_427 = buffer.data(pfg0 + 427);
    const auto *pfg0_430 = buffer.data(pfg0 + 430);
    const auto *pfg0_431 = buffer.data(pfg0 + 431);
    const auto *pfg0_432 = buffer.data(pfg0 + 432);
    const auto *pfg0_435 = buffer.data(pfg0 + 435);
    const auto *pfg0_437 = buffer.data(pfg0 + 437);
    const auto *pfg0_438 = buffer.data(pfg0 + 438);
    const auto *pfg0_440 = buffer.data(pfg0 + 440);
    const auto *pfg0_441 = buffer.data(pfg0 + 441);
    const auto *pfg0_442 = buffer.data(pfg0 + 442);
    const auto *pfg0_444 = buffer.data(pfg0 + 444);
    const auto *pfg0_445 = buffer.data(pfg0 + 445);
    const auto *pfg0_446 = buffer.data(pfg0 + 446);
    const auto *pfg0_447 = buffer.data(pfg0 + 447);
    const auto *pfg0_448 = buffer.data(pfg0 + 448);
    const auto *pfg0_449 = buffer.data(pfg0 + 449);

    const auto *pfg1_392 = buffer.data(pfg1 + 392);
    const auto *pfg1_395 = buffer.data(pfg1 + 395);
    const auto *pfg1_397 = buffer.data(pfg1 + 397);
    const auto *pfg1_399 = buffer.data(pfg1 + 399);
    const auto *pfg1_401 = buffer.data(pfg1 + 401);
    const auto *pfg1_402 = buffer.data(pfg1 + 402);
    const auto *pfg1_404 = buffer.data(pfg1 + 404);
    const auto *pfg1_405 = buffer.data(pfg1 + 405);
    const auto *pfg1_407 = buffer.data(pfg1 + 407);
    const auto *pfg1_408 = buffer.data(pfg1 + 408);
    const auto *pfg1_410 = buffer.data(pfg1 + 410);
    const auto *pfg1_411 = buffer.data(pfg1 + 411);
    const auto *pfg1_412 = buffer.data(pfg1 + 412);
    const auto *pfg1_414 = buffer.data(pfg1 + 414);
    const auto *pfg1_415 = buffer.data(pfg1 + 415);
    const auto *pfg1_416 = buffer.data(pfg1 + 416);
    const auto *pfg1_417 = buffer.data(pfg1 + 417);
    const auto *pfg1_418 = buffer.data(pfg1 + 418);
    const auto *pfg1_419 = buffer.data(pfg1 + 419);
    const auto *pfg1_423 = buffer.data(pfg1 + 423);
    const auto *pfg1_426 = buffer.data(pfg1 + 426);
    const auto *pfg1_427 = buffer.data(pfg1 + 427);
    const auto *pfg1_430 = buffer.data(pfg1 + 430);
    const auto *pfg1_431 = buffer.data(pfg1 + 431);
    const auto *pfg1_432 = buffer.data(pfg1 + 432);
    const auto *pfg1_435 = buffer.data(pfg1 + 435);
    const auto *pfg1_437 = buffer.data(pfg1 + 437);
    const auto *pfg1_438 = buffer.data(pfg1 + 438);
    const auto *pfg1_440 = buffer.data(pfg1 + 440);
    const auto *pfg1_441 = buffer.data(pfg1 + 441);
    const auto *pfg1_442 = buffer.data(pfg1 + 442);
    const auto *pfg1_444 = buffer.data(pfg1 + 444);
    const auto *pfg1_445 = buffer.data(pfg1 + 445);
    const auto *pfg1_446 = buffer.data(pfg1 + 446);
    const auto *pfg1_447 = buffer.data(pfg1 + 447);
    const auto *pfg1_448 = buffer.data(pfg1 + 448);
    const auto *pfg1_449 = buffer.data(pfg1 + 449);

    const auto *pfh_545 = buffer.data(pfh + 545);
    const auto *pfh_546 = buffer.data(pfh + 546);
    const auto *pfh_548 = buffer.data(pfh + 548);
    const auto *pfh_551 = buffer.data(pfh + 551);
    const auto *pfh_553 = buffer.data(pfh + 553);
    const auto *pfh_555 = buffer.data(pfh + 555);
    const auto *pfh_557 = buffer.data(pfh + 557);
    const auto *pfh_558 = buffer.data(pfh + 558);
    const auto *pfh_560 = buffer.data(pfh + 560);
    const auto *pfh_561 = buffer.data(pfh + 561);
    const auto *pfh_562 = buffer.data(pfh + 562);
    const auto *pfh_563 = buffer.data(pfh + 563);
    const auto *pfh_564 = buffer.data(pfh + 564);
    const auto *pfh_565 = buffer.data(pfh + 565);
    const auto *pfh_566 = buffer.data(pfh + 566);
    const auto *pfh_567 = buffer.data(pfh + 567);
    const auto *pfh_569 = buffer.data(pfh + 569);
    const auto *pfh_570 = buffer.data(pfh + 570);
    const auto *pfh_572 = buffer.data(pfh + 572);
    const auto *pfh_573 = buffer.data(pfh + 573);
    const auto *pfh_574 = buffer.data(pfh + 574);
    const auto *pfh_576 = buffer.data(pfh + 576);
    const auto *pfh_577 = buffer.data(pfh + 577);
    const auto *pfh_578 = buffer.data(pfh + 578);
    const auto *pfh_579 = buffer.data(pfh + 579);
    const auto *pfh_581 = buffer.data(pfh + 581);
    const auto *pfh_582 = buffer.data(pfh + 582);
    const auto *pfh_583 = buffer.data(pfh + 583);
    const auto *pfh_584 = buffer.data(pfh + 584);
    const auto *pfh_585 = buffer.data(pfh + 585);
    const auto *pfh_586 = buffer.data(pfh + 586);
    const auto *pfh_587 = buffer.data(pfh + 587);
    const auto *pfh_588 = buffer.data(pfh + 588);
    const auto *pfh_590 = buffer.data(pfh + 590);
    const auto *pfh_591 = buffer.data(pfh + 591);
    const auto *pfh_593 = buffer.data(pfh + 593);
    const auto *pfh_594 = buffer.data(pfh + 594);
    const auto *pfh_595 = buffer.data(pfh + 595);
    const auto *pfh_597 = buffer.data(pfh + 597);
    const auto *pfh_598 = buffer.data(pfh + 598);
    const auto *pfh_599 = buffer.data(pfh + 599);
    const auto *pfh_600 = buffer.data(pfh + 600);
    const auto *pfh_603 = buffer.data(pfh + 603);
    const auto *pfh_604 = buffer.data(pfh + 604);
    const auto *pfh_605 = buffer.data(pfh + 605);
    const auto *pfh_606 = buffer.data(pfh + 606);
    const auto *pfh_607 = buffer.data(pfh + 607);
    const auto *pfh_608 = buffer.data(pfh + 608);
    const auto *pfh_609 = buffer.data(pfh + 609);
    const auto *pfh_611 = buffer.data(pfh + 611);
    const auto *pfh_612 = buffer.data(pfh + 612);
    const auto *pfh_614 = buffer.data(pfh + 614);
    const auto *pfh_615 = buffer.data(pfh + 615);
    const auto *pfh_616 = buffer.data(pfh + 616);
    const auto *pfh_618 = buffer.data(pfh + 618);
    const auto *pfh_619 = buffer.data(pfh + 619);
    const auto *pfh_620 = buffer.data(pfh + 620);
    const auto *pfh_621 = buffer.data(pfh + 621);
    const auto *pfh_623 = buffer.data(pfh + 623);
    const auto *pfh_624 = buffer.data(pfh + 624);
    const auto *pfh_625 = buffer.data(pfh + 625);
    const auto *pfh_626 = buffer.data(pfh + 626);
    const auto *pfh_627 = buffer.data(pfh + 627);
    const auto *pfh_628 = buffer.data(pfh + 628);
    const auto *pfh_629 = buffer.data(pfh + 629);

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_z, pb_x, pc_x, pc_y, pc_z, sfi0_168, \
                         sfi1_168, pdi0_503, pdh_315, pdi1_503, pfh_545, \
                         pfh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_4 * pc_y[k] * pfh_545[k];

        t_727[k] = pb_x[k] * pdi0_503[k]
                   - f_11 * pc_x[k] * pdi1_503[k];

        t_728[k] = pa_z[k] * sfi0_168[k]
                   - f_11 * pc_z[k] * sfi1_168[k];

        t_729[k] = f_1 * pdh_315[k]
                   + f_4 * pc_y[k] * pfh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pa_z, pc_x, pc_y, pc_z, sfi0_171, sfi1_171, \
                         pdh_317, pfg0_392, pfg1_392, pfh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_15 * pfg0_392[k]
                   - f_16 * pfg1_392[k]
                   + f_4 * pc_x[k] * pfh_548[k];

        t_731[k] = pa_z[k] * sfi0_171[k]
                   - f_11 * pc_z[k] * sfi1_171[k];

        t_732[k] = f_1 * pdh_317[k]
                   + f_4 * pc_y[k] * pfh_548[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pa_z, pc_x, pc_z, sfi0_174, sfi1_174, pfg0_395, \
                         pfg0_397, pfg1_395, pfg1_397, pfh_551, \
                         pfh_553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_9 * pfg0_395[k]
                   - f_10 * pfg1_395[k]
                   + f_4 * pc_x[k] * pfh_551[k];

        t_734[k] = pa_z[k] * sfi0_174[k]
                   - f_11 * pc_z[k] * sfi1_174[k];

        t_735[k] = f_7 * pfg0_397[k]
                   - f_8 * pfg1_397[k]
                   + f_4 * pc_x[k] * pfh_553[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pa_z, pc_x, pc_y, pc_z, sfi0_178, sfi1_178, \
                         pdh_320, pfg0_399, pfg1_399, pfh_551, \
                         pfh_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_1 * pdh_320[k]
                   + f_4 * pc_y[k] * pfh_551[k];

        t_737[k] = f_7 * pfg0_399[k]
                   - f_8 * pfg1_399[k]
                   + f_4 * pc_x[k] * pfh_555[k];

        t_738[k] = pa_z[k] * sfi0_178[k]
                   - f_11 * pc_z[k] * sfi1_178[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pc_x, pc_y, pdh_324, pfg0_401, pfg0_402, \
                         pfg1_401, pfg1_402, pfh_555, pfh_557, \
                         pfh_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_5 * pfg0_401[k]
                   - f_6 * pfg1_401[k]
                   + f_4 * pc_x[k] * pfh_557[k];

        t_740[k] = f_5 * pfg0_402[k]
                   - f_6 * pfg1_402[k]
                   + f_4 * pc_x[k] * pfh_558[k];

        t_741[k] = f_1 * pdh_324[k]
                   + f_4 * pc_y[k] * pfh_555[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, t_746, t_747, pc_x, pfg0_404, pfg1_404, \
                         pfh_560, pfh_561, pfh_562, pfh_563, pfh_564, \
                         pfh_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_5 * pfg0_404[k]
                   - f_6 * pfg1_404[k]
                   + f_4 * pc_x[k] * pfh_560[k];

        t_743[k] = f_4 * pc_x[k] * pfh_561[k];

        t_744[k] = f_4 * pc_x[k] * pfh_562[k];

        t_745[k] = f_4 * pc_x[k] * pfh_563[k];

        t_746[k] = f_4 * pc_x[k] * pfh_564[k];

        t_747[k] = f_4 * pc_x[k] * pfh_565[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pa_z, pc_x, pc_z, sfi0_189, sfi0_190, \
                         sfi0_191, sfh_141, sfh_142, sfi1_189, sfi1_190, sfi1_191, \
                         pfh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_4 * pc_x[k] * pfh_566[k];

        t_749[k] = pa_z[k] * sfi0_189[k]
                   - f_11 * pc_z[k] * sfi1_189[k];

        t_750[k] = pa_z[k] * sfi0_190[k]
                   + f_0 * sfh_141[k]
                   - f_11 * pc_z[k] * sfi1_190[k];

        t_751[k] = pa_z[k] * sfi0_191[k]
                   + f_12 * sfh_142[k]
                   - f_11 * pc_z[k] * sfi1_191[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pa_z, pc_y, pc_z, sfi0_192, sfi0_193, sfh_143, \
                         sfh_144, sfi1_192, sfi1_193, pdh_335, \
                         pfh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = pa_z[k] * sfi0_192[k]
                   + f_1 * sfh_143[k]
                   - f_11 * pc_z[k] * sfi1_192[k];

        t_753[k] = pa_z[k] * sfi0_193[k]
                   + f_13 * sfh_144[k]
                   - f_11 * pc_z[k] * sfi1_193[k];

        t_754[k] = f_1 * pdh_335[k]
                   + f_4 * pc_y[k] * pfh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pa_z, pc_x, pc_y, pc_z, sfi0_195, sfh_146, \
                         sfi1_195, pdh_336, pfg0_405, pfg1_405, \
                         pfh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_z[k] * sfi0_195[k]
                   + f_14 * sfh_146[k]
                   - f_11 * pc_z[k] * sfi1_195[k];

        t_756[k] = f_2 * pfg0_405[k]
                   - f_3 * pfg1_405[k]
                   + f_4 * pc_x[k] * pfh_567[k];

        t_757[k] = f_12 * pdh_336[k]
                   + f_4 * pc_y[k] * pfh_567[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, t_761, pc_x, pc_y, pdh_338, pfg0_407, pfg0_408, \
                         pfg0_410, pfg1_407, pfg1_408, pfg1_410, pfh_569, pfh_570, \
                         pfh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_15 * pfg0_407[k]
                   - f_16 * pfg1_407[k]
                   + f_4 * pc_x[k] * pfh_569[k];

        t_759[k] = f_9 * pfg0_408[k]
                   - f_10 * pfg1_408[k]
                   + f_4 * pc_x[k] * pfh_570[k];

        t_760[k] = f_12 * pdh_338[k]
                   + f_4 * pc_y[k] * pfh_569[k];

        t_761[k] = f_9 * pfg0_410[k]
                   - f_10 * pfg1_410[k]
                   + f_4 * pc_x[k] * pfh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_x, pc_y, pdh_341, pfg0_411, pfg0_412, \
                         pfg1_411, pfg1_412, pfh_572, pfh_573, \
                         pfh_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_7 * pfg0_411[k]
                   - f_8 * pfg1_411[k]
                   + f_4 * pc_x[k] * pfh_573[k];

        t_763[k] = f_7 * pfg0_412[k]
                   - f_8 * pfg1_412[k]
                   + f_4 * pc_x[k] * pfh_574[k];

        t_764[k] = f_12 * pdh_341[k]
                   + f_4 * pc_y[k] * pfh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pfg0_414, pfg0_415, pfg0_416, pfg1_414, \
                         pfg1_415, pfg1_416, pfh_576, pfh_577, \
                         pfh_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_7 * pfg0_414[k]
                   - f_8 * pfg1_414[k]
                   + f_4 * pc_x[k] * pfh_576[k];

        t_766[k] = f_5 * pfg0_415[k]
                   - f_6 * pfg1_415[k]
                   + f_4 * pc_x[k] * pfh_577[k];

        t_767[k] = f_5 * pfg0_416[k]
                   - f_6 * pfg1_416[k]
                   + f_4 * pc_x[k] * pfh_578[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, pdh_345, pfg0_417, pfg0_419, \
                         pfg1_417, pfg1_419, pfh_576, pfh_579, pfh_581, \
                         pfh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_5 * pfg0_417[k]
                   - f_6 * pfg1_417[k]
                   + f_4 * pc_x[k] * pfh_579[k];

        t_769[k] = f_12 * pdh_345[k]
                   + f_4 * pc_y[k] * pfh_576[k];

        t_770[k] = f_5 * pfg0_419[k]
                   - f_6 * pfg1_419[k]
                   + f_4 * pc_x[k] * pfh_581[k];

        t_771[k] = f_4 * pc_x[k] * pfh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pfh_583, pfh_584, pfh_585, \
                         pfh_586, pfh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_4 * pc_x[k] * pfh_583[k];

        t_773[k] = f_4 * pc_x[k] * pfh_584[k];

        t_774[k] = f_4 * pc_x[k] * pfh_585[k];

        t_775[k] = f_4 * pc_x[k] * pfh_586[k];

        t_776[k] = f_4 * pc_x[k] * pfh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, pdh_351, pdh_352, pdh_353, pfg0_415, \
                         pfg0_416, pfg0_417, pfg1_415, pfg1_416, pfg1_417, pfh_582, pfh_583, \
                         pfh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_12 * pdh_351[k]
                   + f_2 * pfg0_415[k]
                   - f_3 * pfg1_415[k]
                   + f_4 * pc_y[k] * pfh_582[k];

        t_778[k] = f_12 * pdh_352[k]
                   + f_15 * pfg0_416[k]
                   - f_16 * pfg1_416[k]
                   + f_4 * pc_y[k] * pfh_583[k];

        t_779[k] = f_12 * pdh_353[k]
                   + f_9 * pfg0_417[k]
                   - f_10 * pfg1_417[k]
                   + f_4 * pc_y[k] * pfh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pc_y, pdh_354, pdh_355, pdh_356, pfg0_418, \
                         pfg0_419, pfg1_418, pfg1_419, pfh_585, pfh_586, \
                         pfh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_12 * pdh_354[k]
                   + f_7 * pfg0_418[k]
                   - f_8 * pfg1_418[k]
                   + f_4 * pc_y[k] * pfh_585[k];

        t_781[k] = f_12 * pdh_355[k]
                   + f_5 * pfg0_419[k]
                   - f_6 * pfg1_419[k]
                   + f_4 * pc_y[k] * pfh_586[k];

        t_782[k] = f_12 * pdh_356[k]
                   + f_4 * pc_y[k] * pfh_587[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, t_786, pb_y, pc_y, ppi0_251, ppi1_251, pdi0_475, \
                         pdi0_476, pdi0_478, pdh_357, pdi1_475, pdi1_476, pdi1_478, \
                         pfh_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_17 * ppi0_251[k]
                   - f_18 * ppi1_251[k]
                   + pb_y[k] * pdi0_475[k]
                   - f_11 * pc_y[k] * pdi1_475[k];

        t_784[k] = pb_y[k] * pdi0_476[k]
                   - f_11 * pc_y[k] * pdi1_476[k];

        t_785[k] = f_0 * pdh_357[k]
                   + f_4 * pc_y[k] * pfh_588[k];

        t_786[k] = pb_y[k] * pdi0_478[k]
                   - f_11 * pc_y[k] * pdi1_478[k];
    }

#pragma omp simd aligned(t_787, t_788, t_789, pb_y, pc_x, pc_y, pdi0_481, pdh_359, pdi1_481, \
                         pfg0_423, pfg1_423, pfh_590, pfh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_787[k] = f_9 * pfg0_423[k]
                   - f_10 * pfg1_423[k]
                   + f_4 * pc_x[k] * pfh_591[k];

        t_788[k] = f_0 * pdh_359[k]
                   + f_4 * pc_y[k] * pfh_590[k];

        t_789[k] = pb_y[k] * pdi0_481[k]
                   - f_11 * pc_y[k] * pdi1_481[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, pc_x, pc_y, pdh_362, pfg0_426, pfg0_427, \
                         pfg1_426, pfg1_427, pfh_593, pfh_594, \
                         pfh_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_7 * pfg0_426[k]
                   - f_8 * pfg1_426[k]
                   + f_4 * pc_x[k] * pfh_594[k];

        t_791[k] = f_7 * pfg0_427[k]
                   - f_8 * pfg1_427[k]
                   + f_4 * pc_x[k] * pfh_595[k];

        t_792[k] = f_0 * pdh_362[k]
                   + f_4 * pc_y[k] * pfh_593[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, pb_y, pc_x, pc_y, pdi0_485, pdi1_485, pfg0_430, \
                         pfg0_431, pfg1_430, pfg1_431, pfh_598, \
                         pfh_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = pb_y[k] * pdi0_485[k]
                   - f_11 * pc_y[k] * pdi1_485[k];

        t_794[k] = f_5 * pfg0_430[k]
                   - f_6 * pfg1_430[k]
                   + f_4 * pc_x[k] * pfh_598[k];

        t_795[k] = f_5 * pfg0_431[k]
                   - f_6 * pfg1_431[k]
                   + f_4 * pc_x[k] * pfh_599[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pb_y, pc_x, pc_y, pdi0_490, pdh_366, \
                         pdi1_490, pfg0_432, pfg1_432, pfh_597, pfh_600, \
                         pfh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_5 * pfg0_432[k]
                   - f_6 * pfg1_432[k]
                   + f_4 * pc_x[k] * pfh_600[k];

        t_797[k] = f_0 * pdh_366[k]
                   + f_4 * pc_y[k] * pfh_597[k];

        t_798[k] = pb_y[k] * pdi0_490[k]
                   - f_11 * pc_y[k] * pdi1_490[k];

        t_799[k] = f_4 * pc_x[k] * pfh_603[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pc_x, pfh_604, pfh_605, pfh_606, \
                         pfh_607, pfh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_4 * pc_x[k] * pfh_604[k];

        t_801[k] = f_4 * pc_x[k] * pfh_605[k];

        t_802[k] = f_4 * pc_x[k] * pfh_606[k];

        t_803[k] = f_4 * pc_x[k] * pfh_607[k];

        t_804[k] = f_4 * pc_x[k] * pfh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pb_y, pc_y, pdi0_497, pdi0_498, pdi0_499, \
                         pdh_372, pdh_373, pdh_374, pdi1_497, pdi1_498, \
                         pdi1_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = pb_y[k] * pdi0_497[k]
                   + f_14 * pdh_372[k]
                   - f_11 * pc_y[k] * pdi1_497[k];

        t_806[k] = pb_y[k] * pdi0_498[k]
                   + f_19 * pdh_373[k]
                   - f_11 * pc_y[k] * pdi1_498[k];

        t_807[k] = pb_y[k] * pdi0_499[k]
                   + f_13 * pdh_374[k]
                   - f_11 * pc_y[k] * pdi1_499[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pb_y, pc_y, pdi0_500, pdi0_501, pdi0_503, \
                         pdh_375, pdh_376, pdh_377, pdi1_500, pdi1_501, pdi1_503, \
                         pfh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = pb_y[k] * pdi0_500[k]
                   + f_1 * pdh_375[k]
                   - f_11 * pc_y[k] * pdi1_500[k];

        t_809[k] = pb_y[k] * pdi0_501[k]
                   + f_12 * pdh_376[k]
                   - f_11 * pc_y[k] * pdi1_501[k];

        t_810[k] = f_0 * pdh_377[k]
                   + f_4 * pc_y[k] * pfh_608[k];

        t_811[k] = pb_y[k] * pdi0_503[k]
                   - f_11 * pc_y[k] * pdi1_503[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, t_816, pc_x, pc_y, pfg0_435, pfg0_437, \
                         pfg0_438, pfg1_435, pfg1_437, pfg1_438, pfh_609, pfh_611, \
                         pfh_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_2 * pfg0_435[k]
                   - f_3 * pfg1_435[k]
                   + f_4 * pc_x[k] * pfh_609[k];

        t_813[k] = f_4 * pc_y[k] * pfh_609[k];

        t_814[k] = f_15 * pfg0_437[k]
                   - f_16 * pfg1_437[k]
                   + f_4 * pc_x[k] * pfh_611[k];

        t_815[k] = f_9 * pfg0_438[k]
                   - f_10 * pfg1_438[k]
                   + f_4 * pc_x[k] * pfh_612[k];

        t_816[k] = f_4 * pc_y[k] * pfh_611[k];
    }

#pragma omp simd aligned(t_817, t_818, t_819, t_820, pc_x, pc_y, pfg0_440, pfg0_441, pfg0_442, \
                         pfg1_440, pfg1_441, pfg1_442, pfh_614, pfh_615, \
                         pfh_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_817[k] = f_9 * pfg0_440[k]
                   - f_10 * pfg1_440[k]
                   + f_4 * pc_x[k] * pfh_614[k];

        t_818[k] = f_7 * pfg0_441[k]
                   - f_8 * pfg1_441[k]
                   + f_4 * pc_x[k] * pfh_615[k];

        t_819[k] = f_7 * pfg0_442[k]
                   - f_8 * pfg1_442[k]
                   + f_4 * pc_x[k] * pfh_616[k];

        t_820[k] = f_4 * pc_y[k] * pfh_614[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pfg0_444, pfg0_445, pfg0_446, pfg1_444, \
                         pfg1_445, pfg1_446, pfh_618, pfh_619, \
                         pfh_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_7 * pfg0_444[k]
                   - f_8 * pfg1_444[k]
                   + f_4 * pc_x[k] * pfh_618[k];

        t_822[k] = f_5 * pfg0_445[k]
                   - f_6 * pfg1_445[k]
                   + f_4 * pc_x[k] * pfh_619[k];

        t_823[k] = f_5 * pfg0_446[k]
                   - f_6 * pfg1_446[k]
                   + f_4 * pc_x[k] * pfh_620[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, t_827, t_828, pc_x, pc_y, pfg0_447, pfg0_449, \
                         pfg1_447, pfg1_449, pfh_618, pfh_621, pfh_623, pfh_624, \
                         pfh_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_5 * pfg0_447[k]
                   - f_6 * pfg1_447[k]
                   + f_4 * pc_x[k] * pfh_621[k];

        t_825[k] = f_4 * pc_y[k] * pfh_618[k];

        t_826[k] = f_5 * pfg0_449[k]
                   - f_6 * pfg1_449[k]
                   + f_4 * pc_x[k] * pfh_623[k];

        t_827[k] = f_4 * pc_x[k] * pfh_624[k];

        t_828[k] = f_4 * pc_x[k] * pfh_625[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, t_833, pc_x, pc_y, pfg0_445, pfg1_445, \
                         pfh_624, pfh_626, pfh_627, pfh_628, pfh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_4 * pc_x[k] * pfh_626[k];

        t_830[k] = f_4 * pc_x[k] * pfh_627[k];

        t_831[k] = f_4 * pc_x[k] * pfh_628[k];

        t_832[k] = f_4 * pc_x[k] * pfh_629[k];

        t_833[k] = f_2 * pfg0_445[k]
                   - f_3 * pfg1_445[k]
                   + f_4 * pc_y[k] * pfh_624[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, pc_y, pfg0_446, pfg0_447, pfg0_448, pfg1_446, \
                         pfg1_447, pfg1_448, pfh_625, pfh_626, \
                         pfh_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_15 * pfg0_446[k]
                   - f_16 * pfg1_446[k]
                   + f_4 * pc_y[k] * pfh_625[k];

        t_835[k] = f_9 * pfg0_447[k]
                   - f_10 * pfg1_447[k]
                   + f_4 * pc_y[k] * pfh_626[k];

        t_836[k] = f_7 * pfg0_448[k]
                   - f_8 * pfg1_448[k]
                   + f_4 * pc_y[k] * pfh_627[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, pc_y, pc_z, sfh_209, pdh_377, pfg0_449, \
                         pfg1_449, pfh_628, pfh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = f_5 * pfg0_449[k]
                   - f_6 * pfg1_449[k]
                   + f_4 * pc_y[k] * pfh_628[k];

        t_838[k] = f_4 * pc_y[k] * pfh_629[k];

        t_839[k] = f_0 * sfh_209[k]
                   + f_1 * pdh_377[k]
                   + f_2 * pfg0_449[k]
                   - f_3 * pfg1_449[k]
                   + f_4 * pc_z[k] * pfh_629[k];
    }
}

auto
compute_prim_pfi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sfi0,
                                                   const size_t sfh, const size_t sfi1,
                                                   const size_t ppi0, const size_t ppi1,
                                                   const size_t pdi0, const size_t pdh,
                                                   const size_t pdi1, const size_t pfg0,
                                                   const size_t pfg1, const size_t pfh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pfi_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfh, pdi0,
                                                              pdh, pdi1, pfg0, pfg1, pfh, ncols,
                                                              gamma, p, q);

    compute_prim_pfi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sfi0,
                                                              sfh, sfi1, pdi0, pdh, pdi1, pfg0,
                                                              pfg1, pfh, ncols, gamma, p, q);

    compute_prim_pfi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sfi0,
                                                              sfh, sfi1, ppi0, ppi1, pdi0, pdh,
                                                              pdi1, pfg0, pfg1, pfh, ncols,
                                                              gamma, p, q);

    compute_prim_pfi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sfi0,
                                                              sfh, sfi1, pdi0, pdh, pdi1, pfg0,
                                                              pfg1, pfh, ncols, gamma, p, q);

    compute_prim_pfi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, sfi0,
                                                              sfh, sfi1, pdi0, pdh, pdi1, pfg0,
                                                              pfg1, pfh, ncols, gamma, p, q);

    compute_prim_pfi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pb, pc, sfi0,
                                                              sfh, sfi1, ppi0, ppi1, pdi0, pdh,
                                                              pdi1, pfg0, pfg1, pfh, ncols,
                                                              gamma, p, q);

    compute_prim_pfi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pb, pc, sfi0,
                                                              sfh, sfi1, ppi0, ppi1, pdi0, pdh,
                                                              pdi1, pfg0, pfg1, pfh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
