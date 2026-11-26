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


#include "SimdThreeCenterElectronRepulsionVrrRecPFH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pfh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfg,
                                                          const size_t pdh0, const size_t pdg,
                                                          const size_t pdh1, const size_t pff0,
                                                          const size_t pff1, const size_t pfg,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;

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

    const auto *sfg_0 = buffer.data(sfg + 0);
    const auto *sfg_10 = buffer.data(sfg + 10);
    const auto *sfg_12 = buffer.data(sfg + 12);
    const auto *sfg_14 = buffer.data(sfg + 14);
    const auto *sfg_25 = buffer.data(sfg + 25);
    const auto *sfg_27 = buffer.data(sfg + 27);
    const auto *sfg_42 = buffer.data(sfg + 42);
    const auto *sfg_44 = buffer.data(sfg + 44);
    const auto *sfg_45 = buffer.data(sfg + 45);
    const auto *sfg_55 = buffer.data(sfg + 55);
    const auto *sfg_57 = buffer.data(sfg + 57);
    const auto *sfg_59 = buffer.data(sfg + 59);
    const auto *sfg_72 = buffer.data(sfg + 72);
    const auto *sfg_75 = buffer.data(sfg + 75);
    const auto *sfg_85 = buffer.data(sfg + 85);
    const auto *sfg_87 = buffer.data(sfg + 87);
    const auto *sfg_89 = buffer.data(sfg + 89);

    const auto *pdh0_0 = buffer.data(pdh0 + 0);
    const auto *pdh0_3 = buffer.data(pdh0 + 3);
    const auto *pdh0_5 = buffer.data(pdh0 + 5);
    const auto *pdh0_6 = buffer.data(pdh0 + 6);
    const auto *pdh0_9 = buffer.data(pdh0 + 9);
    const auto *pdh0_10 = buffer.data(pdh0 + 10);
    const auto *pdh0_14 = buffer.data(pdh0 + 14);
    const auto *pdh0_15 = buffer.data(pdh0 + 15);
    const auto *pdh0_20 = buffer.data(pdh0 + 20);
    const auto *pdh0_24 = buffer.data(pdh0 + 24);
    const auto *pdh0_27 = buffer.data(pdh0 + 27);
    const auto *pdh0_31 = buffer.data(pdh0 + 31);
    const auto *pdh0_42 = buffer.data(pdh0 + 42);
    const auto *pdh0_47 = buffer.data(pdh0 + 47);
    const auto *pdh0_51 = buffer.data(pdh0 + 51);
    const auto *pdh0_56 = buffer.data(pdh0 + 56);

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
    const auto *pdg_55 = buffer.data(pdg + 55);
    const auto *pdg_57 = buffer.data(pdg + 57);
    const auto *pdg_59 = buffer.data(pdg + 59);
    const auto *pdg_72 = buffer.data(pdg + 72);
    const auto *pdg_75 = buffer.data(pdg + 75);
    const auto *pdg_85 = buffer.data(pdg + 85);
    const auto *pdg_87 = buffer.data(pdg + 87);
    const auto *pdg_89 = buffer.data(pdg + 89);

    const auto *pdh1_0 = buffer.data(pdh1 + 0);
    const auto *pdh1_3 = buffer.data(pdh1 + 3);
    const auto *pdh1_5 = buffer.data(pdh1 + 5);
    const auto *pdh1_6 = buffer.data(pdh1 + 6);
    const auto *pdh1_9 = buffer.data(pdh1 + 9);
    const auto *pdh1_10 = buffer.data(pdh1 + 10);
    const auto *pdh1_14 = buffer.data(pdh1 + 14);
    const auto *pdh1_15 = buffer.data(pdh1 + 15);
    const auto *pdh1_20 = buffer.data(pdh1 + 20);
    const auto *pdh1_24 = buffer.data(pdh1 + 24);
    const auto *pdh1_27 = buffer.data(pdh1 + 27);
    const auto *pdh1_31 = buffer.data(pdh1 + 31);
    const auto *pdh1_42 = buffer.data(pdh1 + 42);
    const auto *pdh1_47 = buffer.data(pdh1 + 47);
    const auto *pdh1_51 = buffer.data(pdh1 + 51);
    const auto *pdh1_56 = buffer.data(pdh1 + 56);

    const auto *pff0_0 = buffer.data(pff0 + 0);
    const auto *pff0_1 = buffer.data(pff0 + 1);
    const auto *pff0_2 = buffer.data(pff0 + 2);
    const auto *pff0_6 = buffer.data(pff0 + 6);
    const auto *pff0_8 = buffer.data(pff0 + 8);
    const auto *pff0_9 = buffer.data(pff0 + 9);
    const auto *pff0_16 = buffer.data(pff0 + 16);
    const auto *pff0_18 = buffer.data(pff0 + 18);
    const auto *pff0_19 = buffer.data(pff0 + 19);
    const auto *pff0_28 = buffer.data(pff0 + 28);
    const auto *pff0_29 = buffer.data(pff0 + 29);
    const auto *pff0_30 = buffer.data(pff0 + 30);
    const auto *pff0_31 = buffer.data(pff0 + 31);
    const auto *pff0_32 = buffer.data(pff0 + 32);
    const auto *pff0_36 = buffer.data(pff0 + 36);
    const auto *pff0_38 = buffer.data(pff0 + 38);
    const auto *pff0_39 = buffer.data(pff0 + 39);
    const auto *pff0_46 = buffer.data(pff0 + 46);
    const auto *pff0_48 = buffer.data(pff0 + 48);
    const auto *pff0_49 = buffer.data(pff0 + 49);
    const auto *pff0_50 = buffer.data(pff0 + 50);
    const auto *pff0_51 = buffer.data(pff0 + 51);
    const auto *pff0_52 = buffer.data(pff0 + 52);
    const auto *pff0_56 = buffer.data(pff0 + 56);
    const auto *pff0_58 = buffer.data(pff0 + 58);
    const auto *pff0_59 = buffer.data(pff0 + 59);

    const auto *pff1_0 = buffer.data(pff1 + 0);
    const auto *pff1_1 = buffer.data(pff1 + 1);
    const auto *pff1_2 = buffer.data(pff1 + 2);
    const auto *pff1_6 = buffer.data(pff1 + 6);
    const auto *pff1_8 = buffer.data(pff1 + 8);
    const auto *pff1_9 = buffer.data(pff1 + 9);
    const auto *pff1_16 = buffer.data(pff1 + 16);
    const auto *pff1_18 = buffer.data(pff1 + 18);
    const auto *pff1_19 = buffer.data(pff1 + 19);
    const auto *pff1_28 = buffer.data(pff1 + 28);
    const auto *pff1_29 = buffer.data(pff1 + 29);
    const auto *pff1_30 = buffer.data(pff1 + 30);
    const auto *pff1_31 = buffer.data(pff1 + 31);
    const auto *pff1_32 = buffer.data(pff1 + 32);
    const auto *pff1_36 = buffer.data(pff1 + 36);
    const auto *pff1_38 = buffer.data(pff1 + 38);
    const auto *pff1_39 = buffer.data(pff1 + 39);
    const auto *pff1_46 = buffer.data(pff1 + 46);
    const auto *pff1_48 = buffer.data(pff1 + 48);
    const auto *pff1_49 = buffer.data(pff1 + 49);
    const auto *pff1_50 = buffer.data(pff1 + 50);
    const auto *pff1_51 = buffer.data(pff1 + 51);
    const auto *pff1_52 = buffer.data(pff1 + 52);
    const auto *pff1_56 = buffer.data(pff1 + 56);
    const auto *pff1_58 = buffer.data(pff1 + 58);
    const auto *pff1_59 = buffer.data(pff1 + 59);

    const auto *pfg_0 = buffer.data(pfg + 0);
    const auto *pfg_1 = buffer.data(pfg + 1);
    const auto *pfg_2 = buffer.data(pfg + 2);
    const auto *pfg_3 = buffer.data(pfg + 3);
    const auto *pfg_5 = buffer.data(pfg + 5);
    const auto *pfg_6 = buffer.data(pfg + 6);
    const auto *pfg_9 = buffer.data(pfg + 9);
    const auto *pfg_10 = buffer.data(pfg + 10);
    const auto *pfg_12 = buffer.data(pfg + 12);
    const auto *pfg_13 = buffer.data(pfg + 13);
    const auto *pfg_14 = buffer.data(pfg + 14);
    const auto *pfg_15 = buffer.data(pfg + 15);
    const auto *pfg_17 = buffer.data(pfg + 17);
    const auto *pfg_18 = buffer.data(pfg + 18);
    const auto *pfg_20 = buffer.data(pfg + 20);
    const auto *pfg_21 = buffer.data(pfg + 21);
    const auto *pfg_24 = buffer.data(pfg + 24);
    const auto *pfg_25 = buffer.data(pfg + 25);
    const auto *pfg_27 = buffer.data(pfg + 27);
    const auto *pfg_28 = buffer.data(pfg + 28);
    const auto *pfg_29 = buffer.data(pfg + 29);
    const auto *pfg_30 = buffer.data(pfg + 30);
    const auto *pfg_32 = buffer.data(pfg + 32);
    const auto *pfg_33 = buffer.data(pfg + 33);
    const auto *pfg_35 = buffer.data(pfg + 35);
    const auto *pfg_36 = buffer.data(pfg + 36);
    const auto *pfg_39 = buffer.data(pfg + 39);
    const auto *pfg_40 = buffer.data(pfg + 40);
    const auto *pfg_42 = buffer.data(pfg + 42);
    const auto *pfg_43 = buffer.data(pfg + 43);
    const auto *pfg_44 = buffer.data(pfg + 44);
    const auto *pfg_45 = buffer.data(pfg + 45);
    const auto *pfg_46 = buffer.data(pfg + 46);
    const auto *pfg_47 = buffer.data(pfg + 47);
    const auto *pfg_48 = buffer.data(pfg + 48);
    const auto *pfg_50 = buffer.data(pfg + 50);
    const auto *pfg_51 = buffer.data(pfg + 51);
    const auto *pfg_54 = buffer.data(pfg + 54);
    const auto *pfg_55 = buffer.data(pfg + 55);
    const auto *pfg_57 = buffer.data(pfg + 57);
    const auto *pfg_58 = buffer.data(pfg + 58);
    const auto *pfg_59 = buffer.data(pfg + 59);
    const auto *pfg_60 = buffer.data(pfg + 60);
    const auto *pfg_62 = buffer.data(pfg + 62);
    const auto *pfg_63 = buffer.data(pfg + 63);
    const auto *pfg_65 = buffer.data(pfg + 65);
    const auto *pfg_66 = buffer.data(pfg + 66);
    const auto *pfg_69 = buffer.data(pfg + 69);
    const auto *pfg_70 = buffer.data(pfg + 70);
    const auto *pfg_72 = buffer.data(pfg + 72);
    const auto *pfg_73 = buffer.data(pfg + 73);
    const auto *pfg_74 = buffer.data(pfg + 74);
    const auto *pfg_75 = buffer.data(pfg + 75);
    const auto *pfg_76 = buffer.data(pfg + 76);
    const auto *pfg_77 = buffer.data(pfg + 77);
    const auto *pfg_78 = buffer.data(pfg + 78);
    const auto *pfg_80 = buffer.data(pfg + 80);
    const auto *pfg_81 = buffer.data(pfg + 81);
    const auto *pfg_84 = buffer.data(pfg + 84);
    const auto *pfg_85 = buffer.data(pfg + 85);
    const auto *pfg_87 = buffer.data(pfg + 87);
    const auto *pfg_88 = buffer.data(pfg + 88);
    const auto *pfg_89 = buffer.data(pfg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sfg_0, pdg_0, pff0_0, \
                         pff1_0, pfg_0, pfg_1, pfg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfg_0[k]
                 + f_1 * pdg_0[k]
                 + f_2 * pff0_0[k]
                 - f_3 * pff1_0[k]
                 + f_4 * pc_x[k] * pfg_0[k];

        t_1[k] = f_4 * pc_y[k] * pfg_0[k];

        t_2[k] = f_4 * pc_z[k] * pfg_0[k];

        t_3[k] = f_5 * pff0_0[k]
                 - f_6 * pff1_0[k]
                 + f_4 * pc_y[k] * pfg_1[k];

        t_4[k] = f_4 * pc_y[k] * pfg_2[k];

        t_5[k] = f_5 * pff0_0[k]
                 - f_6 * pff1_0[k]
                 + f_4 * pc_z[k] * pfg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, pff0_1, pff0_2, pff1_1, pff1_2, \
                         pfg_3, pfg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * pff0_1[k]
                 - f_8 * pff1_1[k]
                 + f_4 * pc_y[k] * pfg_3[k];

        t_7[k] = f_4 * pc_z[k] * pfg_3[k];

        t_8[k] = f_4 * pc_y[k] * pfg_5[k];

        t_9[k] = f_7 * pff0_2[k]
                 - f_8 * pff1_2[k]
                 + f_4 * pc_z[k] * pfg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, sfg_10, sfg_12, pdg_10, \
                         pdg_12, pfg_6, pfg_9, pfg_10, pfg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sfg_10[k]
                  + f_1 * pdg_10[k]
                  + f_4 * pc_x[k] * pfg_10[k];

        t_11[k] = f_4 * pc_z[k] * pfg_6[k];

        t_12[k] = f_0 * sfg_12[k]
                  + f_1 * pdg_12[k]
                  + f_4 * pc_x[k] * pfg_12[k];

        t_13[k] = f_4 * pc_y[k] * pfg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sfg_14, pdg_14, pff0_6, \
                         pff0_8, pff1_6, pff1_8, pfg_10, pfg_12, \
                         pfg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * sfg_14[k]
                  + f_1 * pdg_14[k]
                  + f_4 * pc_x[k] * pfg_14[k];

        t_15[k] = f_2 * pff0_6[k]
                  - f_3 * pff1_6[k]
                  + f_4 * pc_y[k] * pfg_10[k];

        t_16[k] = f_4 * pc_z[k] * pfg_10[k];

        t_17[k] = f_7 * pff0_8[k]
                  - f_8 * pff1_8[k]
                  + f_4 * pc_y[k] * pfg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, pdh0_0, pdg_0, \
                         pdh1_0, pff0_9, pff1_9, pfg_13, pfg_14, \
                         pfg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * pff0_9[k]
                  - f_6 * pff1_9[k]
                  + f_4 * pc_y[k] * pfg_13[k];

        t_19[k] = f_4 * pc_y[k] * pfg_14[k];

        t_20[k] = f_2 * pff0_9[k]
                  - f_3 * pff1_9[k]
                  + f_4 * pc_z[k] * pfg_14[k];

        t_21[k] = pb_y[k] * pdh0_0[k]
                  - f_9 * pc_y[k] * pdh1_0[k];

        t_22[k] = f_0 * pdg_0[k]
                  + f_4 * pc_y[k] * pfg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, pdh0_3, pdh0_5, pdg_1, \
                         pdg_2, pdh1_3, pdh1_5, pfg_15, pfg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * pc_z[k] * pfg_15[k];

        t_24[k] = pb_y[k] * pdh0_3[k]
                  + f_10 * pdg_1[k]
                  - f_9 * pc_y[k] * pdh1_3[k];

        t_25[k] = f_0 * pdg_2[k]
                  + f_4 * pc_y[k] * pfg_17[k];

        t_26[k] = pb_y[k] * pdh0_5[k]
                  - f_9 * pc_y[k] * pdh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_y, pc_z, pdh0_6, pdh0_9, pdg_3, \
                         pdg_5, pdh1_6, pdh1_9, pfg_18, pfg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * pdh0_6[k]
                  + f_1 * pdg_3[k]
                  - f_9 * pc_y[k] * pdh1_6[k];

        t_28[k] = f_4 * pc_z[k] * pfg_18[k];

        t_29[k] = f_0 * pdg_5[k]
                  + f_4 * pc_y[k] * pfg_20[k];

        t_30[k] = pb_y[k] * pdh0_9[k]
                  - f_9 * pc_y[k] * pdh1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, sfg_25, sfg_27, pdg_9, \
                         pdg_25, pdg_27, pfg_21, pfg_24, pfg_25, \
                         pfg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * sfg_25[k]
                  + f_10 * pdg_25[k]
                  + f_4 * pc_x[k] * pfg_25[k];

        t_32[k] = f_4 * pc_z[k] * pfg_21[k];

        t_33[k] = f_0 * sfg_27[k]
                  + f_10 * pdg_27[k]
                  + f_4 * pc_x[k] * pfg_27[k];

        t_34[k] = f_0 * pdg_9[k]
                  + f_4 * pc_y[k] * pfg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_y, pc_z, pdh0_14, pdg_10, pdh1_14, \
                         pff0_16, pff1_16, pfg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_y[k] * pdh0_14[k]
                  - f_9 * pc_y[k] * pdh1_14[k];

        t_36[k] = f_0 * pdg_10[k]
                  + f_2 * pff0_16[k]
                  - f_3 * pff1_16[k]
                  + f_4 * pc_y[k] * pfg_25[k];

        t_37[k] = f_4 * pc_z[k] * pfg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pc_y, pdg_12, pdg_13, pdg_14, pff0_18, pff0_19, \
                         pff1_18, pff1_19, pfg_27, pfg_28, pfg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_0 * pdg_12[k]
                  + f_7 * pff0_18[k]
                  - f_8 * pff1_18[k]
                  + f_4 * pc_y[k] * pfg_27[k];

        t_39[k] = f_0 * pdg_13[k]
                  + f_5 * pff0_19[k]
                  - f_6 * pff1_19[k]
                  + f_4 * pc_y[k] * pfg_28[k];

        t_40[k] = f_0 * pdg_14[k]
                  + f_4 * pc_y[k] * pfg_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pb_y, pb_z, pc_y, pc_z, pdh0_0, pdh0_20, \
                         pdg_0, pdh1_0, pdh1_20, pfg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pb_y[k] * pdh0_20[k]
                  - f_9 * pc_y[k] * pdh1_20[k];

        t_42[k] = pb_z[k] * pdh0_0[k]
                  - f_9 * pc_z[k] * pdh1_0[k];

        t_43[k] = f_4 * pc_y[k] * pfg_30[k];

        t_44[k] = f_0 * pdg_0[k]
                  + f_4 * pc_z[k] * pfg_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_z, pc_y, pc_z, pdh0_3, pdh0_5, pdh0_6, \
                         pdg_2, pdh1_3, pdh1_5, pdh1_6, pfg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_z[k] * pdh0_3[k]
                  - f_9 * pc_z[k] * pdh1_3[k];

        t_46[k] = f_4 * pc_y[k] * pfg_32[k];

        t_47[k] = pb_z[k] * pdh0_5[k]
                  + f_10 * pdg_2[k]
                  - f_9 * pc_z[k] * pdh1_5[k];

        t_48[k] = pb_z[k] * pdh0_6[k]
                  - f_9 * pc_z[k] * pdh1_6[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_z, pc_y, pc_z, pdh0_9, pdh0_10, pdg_3, \
                         pdg_5, pdh1_9, pdh1_10, pfg_33, pfg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_0 * pdg_3[k]
                  + f_4 * pc_z[k] * pfg_33[k];

        t_50[k] = f_4 * pc_y[k] * pfg_35[k];

        t_51[k] = pb_z[k] * pdh0_9[k]
                  + f_1 * pdg_5[k]
                  - f_9 * pc_z[k] * pdh1_9[k];

        t_52[k] = pb_z[k] * pdh0_10[k]
                  - f_9 * pc_z[k] * pdh1_10[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, sfg_42, sfg_44, pdg_6, \
                         pdg_42, pdg_44, pfg_36, pfg_39, pfg_42, \
                         pfg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * pdg_6[k]
                  + f_4 * pc_z[k] * pfg_36[k];

        t_54[k] = f_0 * sfg_42[k]
                  + f_10 * pdg_42[k]
                  + f_4 * pc_x[k] * pfg_42[k];

        t_55[k] = f_4 * pc_y[k] * pfg_39[k];

        t_56[k] = f_0 * sfg_44[k]
                  + f_10 * pdg_44[k]
                  + f_4 * pc_x[k] * pfg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_z, pc_y, pc_z, pdh0_15, pdg_10, pdh1_15, \
                         pff0_28, pff1_28, pfg_40, pfg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_z[k] * pdh0_15[k]
                  - f_9 * pc_z[k] * pdh1_15[k];

        t_58[k] = f_0 * pdg_10[k]
                  + f_4 * pc_z[k] * pfg_40[k];

        t_59[k] = f_7 * pff0_28[k]
                  - f_8 * pff1_28[k]
                  + f_4 * pc_y[k] * pfg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pc_y, pc_z, pdg_14, pff0_29, pff1_29, pfg_43, \
                         pfg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_5 * pff0_29[k]
                  - f_6 * pff1_29[k]
                  + f_4 * pc_y[k] * pfg_43[k];

        t_61[k] = f_4 * pc_y[k] * pfg_44[k];

        t_62[k] = f_0 * pdg_14[k]
                  + f_2 * pff0_29[k]
                  - f_3 * pff1_29[k]
                  + f_4 * pc_z[k] * pfg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, sfg_45, pdg_15, pdg_16, \
                         pdg_45, pff0_30, pff1_30, pfg_45, pfg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * sfg_45[k]
                  + f_0 * pdg_45[k]
                  + f_2 * pff0_30[k]
                  - f_3 * pff1_30[k]
                  + f_4 * pc_x[k] * pfg_45[k];

        t_64[k] = f_10 * pdg_15[k]
                  + f_4 * pc_y[k] * pfg_45[k];

        t_65[k] = f_4 * pc_z[k] * pfg_45[k];

        t_66[k] = f_10 * pdg_16[k]
                  + f_5 * pff0_30[k]
                  - f_6 * pff1_30[k]
                  + f_4 * pc_y[k] * pfg_46[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pc_y, pc_z, pdg_17, pdg_18, pff0_30, pff0_31, \
                         pff1_30, pff1_31, pfg_47, pfg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * pdg_17[k]
                  + f_4 * pc_y[k] * pfg_47[k];

        t_68[k] = f_5 * pff0_30[k]
                  - f_6 * pff1_30[k]
                  + f_4 * pc_z[k] * pfg_47[k];

        t_69[k] = f_10 * pdg_18[k]
                  + f_7 * pff0_31[k]
                  - f_8 * pff1_31[k]
                  + f_4 * pc_y[k] * pfg_48[k];

        t_70[k] = f_4 * pc_z[k] * pfg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, sfg_55, pdg_20, pdg_55, \
                         pff0_32, pff1_32, pfg_50, pfg_51, pfg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * pdg_20[k]
                  + f_4 * pc_y[k] * pfg_50[k];

        t_72[k] = f_7 * pff0_32[k]
                  - f_8 * pff1_32[k]
                  + f_4 * pc_z[k] * pfg_50[k];

        t_73[k] = f_0 * sfg_55[k]
                  + f_0 * pdg_55[k]
                  + f_4 * pc_x[k] * pfg_55[k];

        t_74[k] = f_4 * pc_z[k] * pfg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pc_x, pc_y, sfg_57, sfg_59, pdg_24, pdg_57, pdg_59, \
                         pfg_54, pfg_57, pfg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * sfg_57[k]
                  + f_0 * pdg_57[k]
                  + f_4 * pc_x[k] * pfg_57[k];

        t_76[k] = f_10 * pdg_24[k]
                  + f_4 * pc_y[k] * pfg_54[k];

        t_77[k] = f_0 * sfg_59[k]
                  + f_0 * pdg_59[k]
                  + f_4 * pc_x[k] * pfg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, pdg_25, pdg_27, pff0_36, pff0_38, \
                         pff1_36, pff1_38, pfg_55, pfg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * pdg_25[k]
                  + f_2 * pff0_36[k]
                  - f_3 * pff1_36[k]
                  + f_4 * pc_y[k] * pfg_55[k];

        t_79[k] = f_4 * pc_z[k] * pfg_55[k];

        t_80[k] = f_10 * pdg_27[k]
                  + f_7 * pff0_38[k]
                  - f_8 * pff1_38[k]
                  + f_4 * pc_y[k] * pfg_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, pdh0_42, pdg_28, pdg_29, \
                         pdh1_42, pff0_39, pff1_39, pfg_58, pfg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * pdg_28[k]
                  + f_5 * pff0_39[k]
                  - f_6 * pff1_39[k]
                  + f_4 * pc_y[k] * pfg_58[k];

        t_82[k] = f_10 * pdg_29[k]
                  + f_4 * pc_y[k] * pfg_59[k];

        t_83[k] = f_2 * pff0_39[k]
                  - f_3 * pff1_39[k]
                  + f_4 * pc_z[k] * pfg_59[k];

        t_84[k] = pb_y[k] * pdh0_42[k]
                  - f_9 * pc_y[k] * pdh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, pdh0_24, pdg_15, pdg_30, \
                         pdg_32, pdh1_24, pfg_60, pfg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_0 * pdg_30[k]
                  + f_4 * pc_y[k] * pfg_60[k];

        t_86[k] = f_0 * pdg_15[k]
                  + f_4 * pc_z[k] * pfg_60[k];

        t_87[k] = pb_z[k] * pdh0_24[k]
                  - f_9 * pc_z[k] * pdh1_24[k];

        t_88[k] = f_0 * pdg_32[k]
                  + f_4 * pc_y[k] * pfg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, pdh0_27, pdh0_47, \
                         pdg_18, pdg_35, pdh1_27, pdh1_47, pfg_63, \
                         pfg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * pdh0_47[k]
                  - f_9 * pc_y[k] * pdh1_47[k];

        t_90[k] = pb_z[k] * pdh0_27[k]
                  - f_9 * pc_z[k] * pdh1_27[k];

        t_91[k] = f_0 * pdg_18[k]
                  + f_4 * pc_z[k] * pfg_63[k];

        t_92[k] = f_0 * pdg_35[k]
                  + f_4 * pc_y[k] * pfg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_y, pb_z, pc_y, pc_z, pdh0_31, pdh0_51, pdg_21, \
                         pdh1_31, pdh1_51, pfg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * pdh0_51[k]
                  - f_9 * pc_y[k] * pdh1_51[k];

        t_94[k] = pb_z[k] * pdh0_31[k]
                  - f_9 * pc_z[k] * pdh1_31[k];

        t_95[k] = f_0 * pdg_21[k]
                  + f_4 * pc_z[k] * pfg_66[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pc_x, pc_y, sfg_72, pdh0_56, pdg_39, pdg_72, \
                         pdh1_56, pfg_69, pfg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_0 * sfg_72[k]
                  + f_0 * pdg_72[k]
                  + f_4 * pc_x[k] * pfg_72[k];

        t_97[k] = f_0 * pdg_39[k]
                  + f_4 * pc_y[k] * pfg_69[k];

        t_98[k] = pb_y[k] * pdh0_56[k]
                  - f_9 * pc_y[k] * pdh1_56[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_y, pc_z, pdg_25, pdg_40, pdg_42, pff0_46, \
                         pff0_48, pff1_46, pff1_48, pfg_70, pfg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_0 * pdg_40[k]
                  + f_2 * pff0_46[k]
                  - f_3 * pff1_46[k]
                  + f_4 * pc_y[k] * pfg_70[k];

        t_100[k] = f_0 * pdg_25[k]
                   + f_4 * pc_z[k] * pfg_70[k];

        t_101[k] = f_0 * pdg_42[k]
                   + f_7 * pff0_48[k]
                   - f_8 * pff1_48[k]
                   + f_4 * pc_y[k] * pfg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pc_y, pc_z, pdg_29, pdg_43, pdg_44, pff0_49, \
                         pff1_49, pfg_73, pfg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_0 * pdg_43[k]
                   + f_5 * pff0_49[k]
                   - f_6 * pff1_49[k]
                   + f_4 * pc_y[k] * pfg_73[k];

        t_103[k] = f_0 * pdg_44[k]
                   + f_4 * pc_y[k] * pfg_74[k];

        t_104[k] = f_0 * pdg_29[k]
                   + f_2 * pff0_49[k]
                   - f_3 * pff1_49[k]
                   + f_4 * pc_z[k] * pfg_74[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, sfg_75, pdg_30, \
                         pdg_75, pff0_50, pff1_50, pfg_75, pfg_76, \
                         pfg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_0 * sfg_75[k]
                   + f_0 * pdg_75[k]
                   + f_2 * pff0_50[k]
                   - f_3 * pff1_50[k]
                   + f_4 * pc_x[k] * pfg_75[k];

        t_106[k] = f_4 * pc_y[k] * pfg_75[k];

        t_107[k] = f_10 * pdg_30[k]
                   + f_4 * pc_z[k] * pfg_75[k];

        t_108[k] = f_5 * pff0_50[k]
                   - f_6 * pff1_50[k]
                   + f_4 * pc_y[k] * pfg_76[k];

        t_109[k] = f_4 * pc_y[k] * pfg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_y, pc_z, pdg_32, pdg_33, pff0_50, \
                         pff0_51, pff1_50, pff1_51, pfg_77, pfg_78, \
                         pfg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_10 * pdg_32[k]
                   + f_5 * pff0_50[k]
                   - f_6 * pff1_50[k]
                   + f_4 * pc_z[k] * pfg_77[k];

        t_111[k] = f_7 * pff0_51[k]
                   - f_8 * pff1_51[k]
                   + f_4 * pc_y[k] * pfg_78[k];

        t_112[k] = f_10 * pdg_33[k]
                   + f_4 * pc_z[k] * pfg_78[k];

        t_113[k] = f_4 * pc_y[k] * pfg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pc_x, pc_z, sfg_85, pdg_35, pdg_36, pdg_85, \
                         pff0_52, pff1_52, pfg_80, pfg_81, pfg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_10 * pdg_35[k]
                   + f_7 * pff0_52[k]
                   - f_8 * pff1_52[k]
                   + f_4 * pc_z[k] * pfg_80[k];

        t_115[k] = f_0 * sfg_85[k]
                   + f_0 * pdg_85[k]
                   + f_4 * pc_x[k] * pfg_85[k];

        t_116[k] = f_10 * pdg_36[k]
                   + f_4 * pc_z[k] * pfg_81[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pc_x, pc_y, sfg_87, sfg_89, pdg_87, \
                         pdg_89, pff0_56, pff1_56, pfg_84, pfg_85, pfg_87, \
                         pfg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * sfg_87[k]
                   + f_0 * pdg_87[k]
                   + f_4 * pc_x[k] * pfg_87[k];

        t_118[k] = f_4 * pc_y[k] * pfg_84[k];

        t_119[k] = f_0 * sfg_89[k]
                   + f_0 * pdg_89[k]
                   + f_4 * pc_x[k] * pfg_89[k];

        t_120[k] = f_2 * pff0_56[k]
                   - f_3 * pff1_56[k]
                   + f_4 * pc_y[k] * pfg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_y, pc_z, pdg_40, pff0_58, pff0_59, \
                         pff1_58, pff1_59, pfg_85, pfg_87, pfg_88, \
                         pfg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * pdg_40[k]
                   + f_4 * pc_z[k] * pfg_85[k];

        t_122[k] = f_7 * pff0_58[k]
                   - f_8 * pff1_58[k]
                   + f_4 * pc_y[k] * pfg_87[k];

        t_123[k] = f_5 * pff0_59[k]
                   - f_6 * pff1_59[k]
                   + f_4 * pc_y[k] * pfg_88[k];

        t_124[k] = f_4 * pc_y[k] * pfg_89[k];
    }
}

static auto
compute_prim_pfh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfh0, const size_t sfg,
                                                          const size_t sfh1, const size_t pph0,
                                                          const size_t pph1, const size_t pdh0,
                                                          const size_t pdg, const size_t pdh1,
                                                          const size_t pff0, const size_t pff1,
                                                          const size_t pfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_0 = buffer.data(sfh0 + 0);
    const auto *sfh0_1 = buffer.data(sfh0 + 1);
    const auto *sfh0_3 = buffer.data(sfh0 + 3);
    const auto *sfh0_5 = buffer.data(sfh0 + 5);
    const auto *sfh0_6 = buffer.data(sfh0 + 6);
    const auto *sfh0_8 = buffer.data(sfh0 + 8);
    const auto *sfh0_9 = buffer.data(sfh0 + 9);
    const auto *sfh0_15 = buffer.data(sfh0 + 15);
    const auto *sfh0_20 = buffer.data(sfh0 + 20);
    const auto *sfh0_126 = buffer.data(sfh0 + 126);
    const auto *sfh0_129 = buffer.data(sfh0 + 129);
    const auto *sfh0_132 = buffer.data(sfh0 + 132);
    const auto *sfh0_141 = buffer.data(sfh0 + 141);
    const auto *sfh0_143 = buffer.data(sfh0 + 143);
    const auto *sfh0_144 = buffer.data(sfh0 + 144);
    const auto *sfh0_146 = buffer.data(sfh0 + 146);
    const auto *sfh0_152 = buffer.data(sfh0 + 152);
    const auto *sfh0_156 = buffer.data(sfh0 + 156);
    const auto *sfh0_162 = buffer.data(sfh0 + 162);
    const auto *sfh0_164 = buffer.data(sfh0 + 164);
    const auto *sfh0_165 = buffer.data(sfh0 + 165);
    const auto *sfh0_167 = buffer.data(sfh0 + 167);
    const auto *sfh0_171 = buffer.data(sfh0 + 171);
    const auto *sfh0_174 = buffer.data(sfh0 + 174);
    const auto *sfh0_183 = buffer.data(sfh0 + 183);
    const auto *sfh0_185 = buffer.data(sfh0 + 185);
    const auto *sfh0_186 = buffer.data(sfh0 + 186);
    const auto *sfh0_188 = buffer.data(sfh0 + 188);
    const auto *sfh0_189 = buffer.data(sfh0 + 189);
    const auto *sfh0_194 = buffer.data(sfh0 + 194);
    const auto *sfh0_198 = buffer.data(sfh0 + 198);
    const auto *sfh0_204 = buffer.data(sfh0 + 204);
    const auto *sfh0_206 = buffer.data(sfh0 + 206);
    const auto *sfh0_207 = buffer.data(sfh0 + 207);
    const auto *sfh0_209 = buffer.data(sfh0 + 209);

    const auto *sfg_0 = buffer.data(sfg + 0);
    const auto *sfg_1 = buffer.data(sfg + 1);
    const auto *sfg_3 = buffer.data(sfg + 3);
    const auto *sfg_5 = buffer.data(sfg + 5);
    const auto *sfg_10 = buffer.data(sfg + 10);
    const auto *sfg_14 = buffer.data(sfg + 14);
    const auto *sfg_90 = buffer.data(sfg + 90);
    const auto *sfg_93 = buffer.data(sfg + 93);
    const auto *sfg_96 = buffer.data(sfg + 96);
    const auto *sfg_100 = buffer.data(sfg + 100);
    const auto *sfg_102 = buffer.data(sfg + 102);
    const auto *sfg_104 = buffer.data(sfg + 104);
    const auto *sfg_110 = buffer.data(sfg + 110);
    const auto *sfg_114 = buffer.data(sfg + 114);
    const auto *sfg_115 = buffer.data(sfg + 115);
    const auto *sfg_117 = buffer.data(sfg + 117);
    const auto *sfg_119 = buffer.data(sfg + 119);
    const auto *sfg_123 = buffer.data(sfg + 123);
    const auto *sfg_126 = buffer.data(sfg + 126);
    const auto *sfg_130 = buffer.data(sfg + 130);
    const auto *sfg_132 = buffer.data(sfg + 132);
    const auto *sfg_134 = buffer.data(sfg + 134);
    const auto *sfg_135 = buffer.data(sfg + 135);
    const auto *sfg_140 = buffer.data(sfg + 140);
    const auto *sfg_144 = buffer.data(sfg + 144);
    const auto *sfg_145 = buffer.data(sfg + 145);
    const auto *sfg_147 = buffer.data(sfg + 147);
    const auto *sfg_149 = buffer.data(sfg + 149);

    const auto *sfh1_0 = buffer.data(sfh1 + 0);
    const auto *sfh1_1 = buffer.data(sfh1 + 1);
    const auto *sfh1_3 = buffer.data(sfh1 + 3);
    const auto *sfh1_5 = buffer.data(sfh1 + 5);
    const auto *sfh1_6 = buffer.data(sfh1 + 6);
    const auto *sfh1_8 = buffer.data(sfh1 + 8);
    const auto *sfh1_9 = buffer.data(sfh1 + 9);
    const auto *sfh1_15 = buffer.data(sfh1 + 15);
    const auto *sfh1_20 = buffer.data(sfh1 + 20);
    const auto *sfh1_126 = buffer.data(sfh1 + 126);
    const auto *sfh1_129 = buffer.data(sfh1 + 129);
    const auto *sfh1_132 = buffer.data(sfh1 + 132);
    const auto *sfh1_141 = buffer.data(sfh1 + 141);
    const auto *sfh1_143 = buffer.data(sfh1 + 143);
    const auto *sfh1_144 = buffer.data(sfh1 + 144);
    const auto *sfh1_146 = buffer.data(sfh1 + 146);
    const auto *sfh1_152 = buffer.data(sfh1 + 152);
    const auto *sfh1_156 = buffer.data(sfh1 + 156);
    const auto *sfh1_162 = buffer.data(sfh1 + 162);
    const auto *sfh1_164 = buffer.data(sfh1 + 164);
    const auto *sfh1_165 = buffer.data(sfh1 + 165);
    const auto *sfh1_167 = buffer.data(sfh1 + 167);
    const auto *sfh1_171 = buffer.data(sfh1 + 171);
    const auto *sfh1_174 = buffer.data(sfh1 + 174);
    const auto *sfh1_183 = buffer.data(sfh1 + 183);
    const auto *sfh1_185 = buffer.data(sfh1 + 185);
    const auto *sfh1_186 = buffer.data(sfh1 + 186);
    const auto *sfh1_188 = buffer.data(sfh1 + 188);
    const auto *sfh1_189 = buffer.data(sfh1 + 189);
    const auto *sfh1_194 = buffer.data(sfh1 + 194);
    const auto *sfh1_198 = buffer.data(sfh1 + 198);
    const auto *sfh1_204 = buffer.data(sfh1 + 204);
    const auto *sfh1_206 = buffer.data(sfh1 + 206);
    const auto *sfh1_207 = buffer.data(sfh1 + 207);
    const auto *sfh1_209 = buffer.data(sfh1 + 209);

    const auto *pph0_99 = buffer.data(pph0 + 99);

    const auto *pph1_99 = buffer.data(pph1 + 99);

    const auto *pdh0_63 = buffer.data(pdh0 + 63);
    const auto *pdh0_66 = buffer.data(pdh0 + 66);
    const auto *pdh0_69 = buffer.data(pdh0 + 69);
    const auto *pdh0_105 = buffer.data(pdh0 + 105);
    const auto *pdh0_110 = buffer.data(pdh0 + 110);
    const auto *pdh0_114 = buffer.data(pdh0 + 114);
    const auto *pdh0_162 = buffer.data(pdh0 + 162);

    const auto *pdg_44 = buffer.data(pdg + 44);
    const auto *pdg_45 = buffer.data(pdg + 45);
    const auto *pdg_47 = buffer.data(pdg + 47);
    const auto *pdg_48 = buffer.data(pdg + 48);
    const auto *pdg_50 = buffer.data(pdg + 50);
    const auto *pdg_51 = buffer.data(pdg + 51);
    const auto *pdg_54 = buffer.data(pdg + 54);
    const auto *pdg_55 = buffer.data(pdg + 55);
    const auto *pdg_59 = buffer.data(pdg + 59);
    const auto *pdg_60 = buffer.data(pdg + 60);
    const auto *pdg_62 = buffer.data(pdg + 62);
    const auto *pdg_63 = buffer.data(pdg + 63);
    const auto *pdg_65 = buffer.data(pdg + 65);
    const auto *pdg_66 = buffer.data(pdg + 66);
    const auto *pdg_69 = buffer.data(pdg + 69);
    const auto *pdg_70 = buffer.data(pdg + 70);
    const auto *pdg_74 = buffer.data(pdg + 74);
    const auto *pdg_75 = buffer.data(pdg + 75);
    const auto *pdg_77 = buffer.data(pdg + 77);
    const auto *pdg_78 = buffer.data(pdg + 78);
    const auto *pdg_80 = buffer.data(pdg + 80);
    const auto *pdg_81 = buffer.data(pdg + 81);
    const auto *pdg_84 = buffer.data(pdg + 84);
    const auto *pdg_85 = buffer.data(pdg + 85);
    const auto *pdg_89 = buffer.data(pdg + 89);
    const auto *pdg_100 = buffer.data(pdg + 100);
    const auto *pdg_101 = buffer.data(pdg + 101);
    const auto *pdg_102 = buffer.data(pdg + 102);
    const auto *pdg_103 = buffer.data(pdg + 103);
    const auto *pdg_104 = buffer.data(pdg + 104);
    const auto *pdg_105 = buffer.data(pdg + 105);
    const auto *pdg_106 = buffer.data(pdg + 106);
    const auto *pdg_108 = buffer.data(pdg + 108);
    const auto *pdg_110 = buffer.data(pdg + 110);
    const auto *pdg_111 = buffer.data(pdg + 111);
    const auto *pdg_113 = buffer.data(pdg + 113);
    const auto *pdg_114 = buffer.data(pdg + 114);
    const auto *pdg_115 = buffer.data(pdg + 115);
    const auto *pdg_116 = buffer.data(pdg + 116);
    const auto *pdg_117 = buffer.data(pdg + 117);
    const auto *pdg_118 = buffer.data(pdg + 118);
    const auto *pdg_119 = buffer.data(pdg + 119);

    const auto *pdh1_63 = buffer.data(pdh1 + 63);
    const auto *pdh1_66 = buffer.data(pdh1 + 66);
    const auto *pdh1_69 = buffer.data(pdh1 + 69);
    const auto *pdh1_105 = buffer.data(pdh1 + 105);
    const auto *pdh1_110 = buffer.data(pdh1 + 110);
    const auto *pdh1_114 = buffer.data(pdh1 + 114);
    const auto *pdh1_162 = buffer.data(pdh1 + 162);

    const auto *pff0_59 = buffer.data(pff0 + 59);
    const auto *pff0_60 = buffer.data(pff0 + 60);
    const auto *pff0_62 = buffer.data(pff0 + 62);
    const auto *pff0_90 = buffer.data(pff0 + 90);
    const auto *pff0_91 = buffer.data(pff0 + 91);
    const auto *pff0_106 = buffer.data(pff0 + 106);
    const auto *pff0_107 = buffer.data(pff0 + 107);
    const auto *pff0_110 = buffer.data(pff0 + 110);
    const auto *pff0_111 = buffer.data(pff0 + 111);
    const auto *pff0_113 = buffer.data(pff0 + 113);
    const auto *pff0_115 = buffer.data(pff0 + 115);
    const auto *pff0_116 = buffer.data(pff0 + 116);
    const auto *pff0_118 = buffer.data(pff0 + 118);
    const auto *pff0_119 = buffer.data(pff0 + 119);

    const auto *pff1_59 = buffer.data(pff1 + 59);
    const auto *pff1_60 = buffer.data(pff1 + 60);
    const auto *pff1_62 = buffer.data(pff1 + 62);
    const auto *pff1_90 = buffer.data(pff1 + 90);
    const auto *pff1_91 = buffer.data(pff1 + 91);
    const auto *pff1_106 = buffer.data(pff1 + 106);
    const auto *pff1_107 = buffer.data(pff1 + 107);
    const auto *pff1_110 = buffer.data(pff1 + 110);
    const auto *pff1_111 = buffer.data(pff1 + 111);
    const auto *pff1_113 = buffer.data(pff1 + 113);
    const auto *pff1_115 = buffer.data(pff1 + 115);
    const auto *pff1_116 = buffer.data(pff1 + 116);
    const auto *pff1_118 = buffer.data(pff1 + 118);
    const auto *pff1_119 = buffer.data(pff1 + 119);

    const auto *pfg_89 = buffer.data(pfg + 89);
    const auto *pfg_90 = buffer.data(pfg + 90);
    const auto *pfg_92 = buffer.data(pfg + 92);
    const auto *pfg_93 = buffer.data(pfg + 93);
    const auto *pfg_95 = buffer.data(pfg + 95);
    const auto *pfg_96 = buffer.data(pfg + 96);
    const auto *pfg_99 = buffer.data(pfg + 99);
    const auto *pfg_100 = buffer.data(pfg + 100);
    const auto *pfg_102 = buffer.data(pfg + 102);
    const auto *pfg_104 = buffer.data(pfg + 104);
    const auto *pfg_105 = buffer.data(pfg + 105);
    const auto *pfg_107 = buffer.data(pfg + 107);
    const auto *pfg_108 = buffer.data(pfg + 108);
    const auto *pfg_110 = buffer.data(pfg + 110);
    const auto *pfg_111 = buffer.data(pfg + 111);
    const auto *pfg_114 = buffer.data(pfg + 114);
    const auto *pfg_115 = buffer.data(pfg + 115);
    const auto *pfg_117 = buffer.data(pfg + 117);
    const auto *pfg_119 = buffer.data(pfg + 119);
    const auto *pfg_120 = buffer.data(pfg + 120);
    const auto *pfg_122 = buffer.data(pfg + 122);
    const auto *pfg_123 = buffer.data(pfg + 123);
    const auto *pfg_125 = buffer.data(pfg + 125);
    const auto *pfg_126 = buffer.data(pfg + 126);
    const auto *pfg_129 = buffer.data(pfg + 129);
    const auto *pfg_130 = buffer.data(pfg + 130);
    const auto *pfg_132 = buffer.data(pfg + 132);
    const auto *pfg_134 = buffer.data(pfg + 134);
    const auto *pfg_135 = buffer.data(pfg + 135);
    const auto *pfg_136 = buffer.data(pfg + 136);
    const auto *pfg_137 = buffer.data(pfg + 137);
    const auto *pfg_138 = buffer.data(pfg + 138);
    const auto *pfg_140 = buffer.data(pfg + 140);
    const auto *pfg_141 = buffer.data(pfg + 141);
    const auto *pfg_144 = buffer.data(pfg + 144);
    const auto *pfg_145 = buffer.data(pfg + 145);
    const auto *pfg_147 = buffer.data(pfg + 147);
    const auto *pfg_149 = buffer.data(pfg + 149);
    const auto *pfg_150 = buffer.data(pfg + 150);
    const auto *pfg_151 = buffer.data(pfg + 151);
    const auto *pfg_153 = buffer.data(pfg + 153);
    const auto *pfg_160 = buffer.data(pfg + 160);
    const auto *pfg_161 = buffer.data(pfg + 161);
    const auto *pfg_162 = buffer.data(pfg + 162);
    const auto *pfg_163 = buffer.data(pfg + 163);
    const auto *pfg_164 = buffer.data(pfg + 164);
    const auto *pfg_165 = buffer.data(pfg + 165);
    const auto *pfg_166 = buffer.data(pfg + 166);
    const auto *pfg_168 = buffer.data(pfg + 168);
    const auto *pfg_170 = buffer.data(pfg + 170);
    const auto *pfg_171 = buffer.data(pfg + 171);
    const auto *pfg_173 = buffer.data(pfg + 173);
    const auto *pfg_174 = buffer.data(pfg + 174);
    const auto *pfg_175 = buffer.data(pfg + 175);
    const auto *pfg_176 = buffer.data(pfg + 176);
    const auto *pfg_177 = buffer.data(pfg + 177);
    const auto *pfg_178 = buffer.data(pfg + 178);
    const auto *pfg_179 = buffer.data(pfg + 179);

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pc_x, pc_y, pc_z, sfh0_126, sfg_90, \
                         sfh1_126, pdg_44, pdg_45, pff0_59, pff1_59, pfg_89, \
                         pfg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_10 * pdg_44[k]
                   + f_2 * pff0_59[k]
                   - f_3 * pff1_59[k]
                   + f_4 * pc_z[k] * pfg_89[k];

        t_126[k] = pa_x[k] * sfh0_126[k]
                   + f_11 * sfg_90[k]
                   - f_9 * pc_x[k] * sfh1_126[k];

        t_127[k] = f_1 * pdg_45[k]
                   + f_4 * pc_y[k] * pfg_90[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_x, pc_x, pc_y, pc_z, sfh0_129, sfg_93, \
                         sfh1_129, pdg_47, pff0_60, pff1_60, pfg_90, \
                         pfg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_4 * pc_z[k] * pfg_90[k];

        t_129[k] = pa_x[k] * sfh0_129[k]
                   + f_1 * sfg_93[k]
                   - f_9 * pc_x[k] * sfh1_129[k];

        t_130[k] = f_1 * pdg_47[k]
                   + f_4 * pc_y[k] * pfg_92[k];

        t_131[k] = f_5 * pff0_60[k]
                   - f_6 * pff1_60[k]
                   + f_4 * pc_z[k] * pfg_92[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pc_x, pc_y, pc_z, sfh0_132, sfg_96, \
                         sfh1_132, pdg_50, pff0_62, pff1_62, pfg_93, \
                         pfg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_x[k] * sfh0_132[k]
                   + f_10 * sfg_96[k]
                   - f_9 * pc_x[k] * sfh1_132[k];

        t_133[k] = f_4 * pc_z[k] * pfg_93[k];

        t_134[k] = f_1 * pdg_50[k]
                   + f_4 * pc_y[k] * pfg_95[k];

        t_135[k] = f_7 * pff0_62[k]
                   - f_8 * pff1_62[k]
                   + f_4 * pc_z[k] * pfg_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pc_x, pc_y, pc_z, sfg_100, sfg_102, \
                         pdg_54, pfg_96, pfg_99, pfg_100, pfg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_0 * sfg_100[k]
                   + f_4 * pc_x[k] * pfg_100[k];

        t_137[k] = f_4 * pc_z[k] * pfg_96[k];

        t_138[k] = f_0 * sfg_102[k]
                   + f_4 * pc_x[k] * pfg_102[k];

        t_139[k] = f_1 * pdg_54[k]
                   + f_4 * pc_y[k] * pfg_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pc_x, pc_z, sfh0_141, sfh0_143, \
                         sfg_104, sfh1_141, sfh1_143, pfg_100, \
                         pfg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * sfg_104[k]
                   + f_4 * pc_x[k] * pfg_104[k];

        t_141[k] = pa_x[k] * sfh0_141[k]
                   - f_9 * pc_x[k] * sfh1_141[k];

        t_142[k] = f_4 * pc_z[k] * pfg_100[k];

        t_143[k] = pa_x[k] * sfh0_143[k]
                   - f_9 * pc_x[k] * sfh1_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pa_x, pc_x, pc_y, sfh0_144, sfh0_146, sfh1_144, \
                         sfh1_146, pdg_59, pfg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pa_x[k] * sfh0_144[k]
                   - f_9 * pc_x[k] * sfh1_144[k];

        t_145[k] = f_1 * pdg_59[k]
                   + f_4 * pc_y[k] * pfg_104[k];

        t_146[k] = pa_x[k] * sfh0_146[k]
                   - f_9 * pc_x[k] * sfh1_146[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pb_z, pc_y, pc_z, pdh0_63, pdh0_66, \
                         pdg_45, pdg_60, pdh1_63, pdh1_66, pfg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pb_z[k] * pdh0_63[k]
                   - f_9 * pc_z[k] * pdh1_63[k];

        t_148[k] = f_10 * pdg_60[k]
                   + f_4 * pc_y[k] * pfg_105[k];

        t_149[k] = f_0 * pdg_45[k]
                   + f_4 * pc_z[k] * pfg_105[k];

        t_150[k] = pb_z[k] * pdh0_66[k]
                   - f_9 * pc_z[k] * pdh1_66[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_x, pb_z, pc_x, pc_y, pc_z, sfh0_152, sfg_110, \
                         sfh1_152, pdh0_69, pdg_62, pdh1_69, pfg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_10 * pdg_62[k]
                   + f_4 * pc_y[k] * pfg_107[k];

        t_152[k] = pa_x[k] * sfh0_152[k]
                   + f_1 * sfg_110[k]
                   - f_9 * pc_x[k] * sfh1_152[k];

        t_153[k] = pb_z[k] * pdh0_69[k]
                   - f_9 * pc_z[k] * pdh1_69[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_x, pc_x, pc_y, pc_z, sfh0_156, sfg_114, \
                         sfh1_156, pdg_48, pdg_65, pfg_108, pfg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_0 * pdg_48[k]
                   + f_4 * pc_z[k] * pfg_108[k];

        t_155[k] = f_10 * pdg_65[k]
                   + f_4 * pc_y[k] * pfg_110[k];

        t_156[k] = pa_x[k] * sfh0_156[k]
                   + f_10 * sfg_114[k]
                   - f_9 * pc_x[k] * sfh1_156[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, pc_z, sfg_115, sfg_117, \
                         pdg_51, pdg_69, pfg_111, pfg_114, pfg_115, \
                         pfg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_0 * sfg_115[k]
                   + f_4 * pc_x[k] * pfg_115[k];

        t_158[k] = f_0 * pdg_51[k]
                   + f_4 * pc_z[k] * pfg_111[k];

        t_159[k] = f_0 * sfg_117[k]
                   + f_4 * pc_x[k] * pfg_117[k];

        t_160[k] = f_10 * pdg_69[k]
                   + f_4 * pc_y[k] * pfg_114[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pa_x, pc_x, pc_z, sfh0_162, sfh0_164, \
                         sfg_119, sfh1_162, sfh1_164, pdg_55, pfg_115, \
                         pfg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_0 * sfg_119[k]
                   + f_4 * pc_x[k] * pfg_119[k];

        t_162[k] = pa_x[k] * sfh0_162[k]
                   - f_9 * pc_x[k] * sfh1_162[k];

        t_163[k] = f_0 * pdg_55[k]
                   + f_4 * pc_z[k] * pfg_115[k];

        t_164[k] = pa_x[k] * sfh0_164[k]
                   - f_9 * pc_x[k] * sfh1_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pa_x, pb_y, pc_x, pc_y, sfh0_165, \
                         sfh0_167, sfh1_165, sfh1_167, pdh0_105, pdg_74, pdh1_105, \
                         pfg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pa_x[k] * sfh0_165[k]
                   - f_9 * pc_x[k] * sfh1_165[k];

        t_166[k] = f_10 * pdg_74[k]
                   + f_4 * pc_y[k] * pfg_119[k];

        t_167[k] = pa_x[k] * sfh0_167[k]
                   - f_9 * pc_x[k] * sfh1_167[k];

        t_168[k] = pb_y[k] * pdh0_105[k]
                   - f_9 * pc_y[k] * pdh1_105[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_x, pc_x, pc_y, pc_z, sfh0_171, \
                         sfg_123, sfh1_171, pdg_60, pdg_75, pdg_77, pfg_120, \
                         pfg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_0 * pdg_75[k]
                   + f_4 * pc_y[k] * pfg_120[k];

        t_170[k] = f_10 * pdg_60[k]
                   + f_4 * pc_z[k] * pfg_120[k];

        t_171[k] = pa_x[k] * sfh0_171[k]
                   + f_1 * sfg_123[k]
                   - f_9 * pc_x[k] * sfh1_171[k];

        t_172[k] = f_0 * pdg_77[k]
                   + f_4 * pc_y[k] * pfg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pa_x, pb_y, pc_x, pc_y, pc_z, sfh0_174, sfg_126, \
                         sfh1_174, pdh0_110, pdg_63, pdh1_110, \
                         pfg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pb_y[k] * pdh0_110[k]
                   - f_9 * pc_y[k] * pdh1_110[k];

        t_174[k] = pa_x[k] * sfh0_174[k]
                   + f_10 * sfg_126[k]
                   - f_9 * pc_x[k] * sfh1_174[k];

        t_175[k] = f_10 * pdg_63[k]
                   + f_4 * pc_z[k] * pfg_123[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pc_x, pc_y, pc_z, sfg_130, \
                         pdh0_114, pdg_66, pdg_80, pdh1_114, pfg_125, pfg_126, \
                         pfg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_0 * pdg_80[k]
                   + f_4 * pc_y[k] * pfg_125[k];

        t_177[k] = pb_y[k] * pdh0_114[k]
                   - f_9 * pc_y[k] * pdh1_114[k];

        t_178[k] = f_0 * sfg_130[k]
                   + f_4 * pc_x[k] * pfg_130[k];

        t_179[k] = f_10 * pdg_66[k]
                   + f_4 * pc_z[k] * pfg_126[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pc_x, pc_y, sfh0_183, sfg_132, \
                         sfg_134, sfh1_183, pdg_84, pfg_129, pfg_132, \
                         pfg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_0 * sfg_132[k]
                   + f_4 * pc_x[k] * pfg_132[k];

        t_181[k] = f_0 * pdg_84[k]
                   + f_4 * pc_y[k] * pfg_129[k];

        t_182[k] = f_0 * sfg_134[k]
                   + f_4 * pc_x[k] * pfg_134[k];

        t_183[k] = pa_x[k] * sfh0_183[k]
                   - f_9 * pc_x[k] * sfh1_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pc_x, pc_y, pc_z, sfh0_185, \
                         sfh0_186, sfh1_185, sfh1_186, pdg_70, pdg_89, pfg_130, \
                         pfg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_10 * pdg_70[k]
                   + f_4 * pc_z[k] * pfg_130[k];

        t_185[k] = pa_x[k] * sfh0_185[k]
                   - f_9 * pc_x[k] * sfh1_185[k];

        t_186[k] = pa_x[k] * sfh0_186[k]
                   - f_9 * pc_x[k] * sfh1_186[k];

        t_187[k] = f_0 * pdg_89[k]
                   + f_4 * pc_y[k] * pfg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_x, pc_x, pc_y, pc_z, sfh0_188, \
                         sfh0_189, sfg_135, sfh1_188, sfh1_189, pdg_75, \
                         pfg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_x[k] * sfh0_188[k]
                   - f_9 * pc_x[k] * sfh1_188[k];

        t_189[k] = pa_x[k] * sfh0_189[k]
                   + f_11 * sfg_135[k]
                   - f_9 * pc_x[k] * sfh1_189[k];

        t_190[k] = f_4 * pc_y[k] * pfg_135[k];

        t_191[k] = f_1 * pdg_75[k]
                   + f_4 * pc_z[k] * pfg_135[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_x, pc_x, pc_y, sfh0_194, sfg_140, sfh1_194, \
                         pff0_90, pff1_90, pfg_136, pfg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_5 * pff0_90[k]
                   - f_6 * pff1_90[k]
                   + f_4 * pc_y[k] * pfg_136[k];

        t_193[k] = f_4 * pc_y[k] * pfg_137[k];

        t_194[k] = pa_x[k] * sfh0_194[k]
                   + f_1 * sfg_140[k]
                   - f_9 * pc_x[k] * sfh1_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_x, pc_x, pc_y, pc_z, sfh0_198, \
                         sfg_144, sfh1_198, pdg_78, pff0_91, pff1_91, pfg_138, \
                         pfg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_7 * pff0_91[k]
                   - f_8 * pff1_91[k]
                   + f_4 * pc_y[k] * pfg_138[k];

        t_196[k] = f_1 * pdg_78[k]
                   + f_4 * pc_z[k] * pfg_138[k];

        t_197[k] = f_4 * pc_y[k] * pfg_140[k];

        t_198[k] = pa_x[k] * sfh0_198[k]
                   + f_10 * sfg_144[k]
                   - f_9 * pc_x[k] * sfh1_198[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pc_x, pc_y, pc_z, sfg_145, sfg_147, \
                         pdg_81, pfg_141, pfg_144, pfg_145, pfg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_0 * sfg_145[k]
                   + f_4 * pc_x[k] * pfg_145[k];

        t_200[k] = f_1 * pdg_81[k]
                   + f_4 * pc_z[k] * pfg_141[k];

        t_201[k] = f_0 * sfg_147[k]
                   + f_4 * pc_x[k] * pfg_147[k];

        t_202[k] = f_4 * pc_y[k] * pfg_144[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pa_x, pc_x, pc_z, sfh0_204, sfh0_206, \
                         sfg_149, sfh1_204, sfh1_206, pdg_85, pfg_145, \
                         pfg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_0 * sfg_149[k]
                   + f_4 * pc_x[k] * pfg_149[k];

        t_204[k] = pa_x[k] * sfh0_204[k]
                   - f_9 * pc_x[k] * sfh1_204[k];

        t_205[k] = f_1 * pdg_85[k]
                   + f_4 * pc_z[k] * pfg_145[k];

        t_206[k] = pa_x[k] * sfh0_206[k]
                   - f_9 * pc_x[k] * sfh1_206[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pa_x, pa_y, pc_x, pc_y, sfh0_0, sfh0_207, \
                         sfh0_209, sfh1_0, sfh1_207, sfh1_209, \
                         pfg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_x[k] * sfh0_207[k]
                   - f_9 * pc_x[k] * sfh1_207[k];

        t_208[k] = f_4 * pc_y[k] * pfg_149[k];

        t_209[k] = pa_x[k] * sfh0_209[k]
                   - f_9 * pc_x[k] * sfh1_209[k];

        t_210[k] = pa_y[k] * sfh0_0[k]
                   - f_9 * pc_y[k] * sfh1_0[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pa_y, pc_y, pc_z, sfh0_1, sfh0_3, sfg_0, \
                         sfg_1, sfh1_1, sfh1_3, pfg_150, pfg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = pa_y[k] * sfh0_1[k]
                   + f_0 * sfg_0[k]
                   - f_9 * pc_y[k] * sfh1_1[k];

        t_212[k] = f_4 * pc_z[k] * pfg_150[k];

        t_213[k] = pa_y[k] * sfh0_3[k]
                   + f_10 * sfg_1[k]
                   - f_9 * pc_y[k] * sfh1_3[k];

        t_214[k] = f_4 * pc_z[k] * pfg_151[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pa_y, pc_y, pc_z, sfh0_5, sfh0_6, sfh0_8, \
                         sfg_3, sfg_5, sfh1_5, sfh1_6, sfh1_8, \
                         pfg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pa_y[k] * sfh0_5[k]
                   - f_9 * pc_y[k] * sfh1_5[k];

        t_216[k] = pa_y[k] * sfh0_6[k]
                   + f_1 * sfg_3[k]
                   - f_9 * pc_y[k] * sfh1_6[k];

        t_217[k] = f_4 * pc_z[k] * pfg_153[k];

        t_218[k] = pa_y[k] * sfh0_8[k]
                   + f_0 * sfg_5[k]
                   - f_9 * pc_y[k] * sfh1_8[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pa_y, pc_x, pc_y, sfh0_9, sfh1_9, \
                         pdg_100, pdg_101, pdg_102, pfg_160, pfg_161, \
                         pfg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = pa_y[k] * sfh0_9[k]
                   - f_9 * pc_y[k] * sfh1_9[k];

        t_220[k] = f_1 * pdg_100[k]
                   + f_4 * pc_x[k] * pfg_160[k];

        t_221[k] = f_1 * pdg_101[k]
                   + f_4 * pc_x[k] * pfg_161[k];

        t_222[k] = f_1 * pdg_102[k]
                   + f_4 * pc_x[k] * pfg_162[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pa_y, pc_x, pc_y, pc_z, sfh0_15, sfg_10, \
                         sfh1_15, pdg_103, pdg_104, pfg_160, pfg_163, \
                         pfg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_1 * pdg_103[k]
                   + f_4 * pc_x[k] * pfg_163[k];

        t_224[k] = f_1 * pdg_104[k]
                   + f_4 * pc_x[k] * pfg_164[k];

        t_225[k] = pa_y[k] * sfh0_15[k]
                   + f_11 * sfg_10[k]
                   - f_9 * pc_y[k] * sfh1_15[k];

        t_226[k] = f_4 * pc_z[k] * pfg_160[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_y, pc_z, sfg_14, pff0_106, pff0_107, \
                         pff1_106, pff1_107, pfg_161, pfg_162, \
                         pfg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_5 * pff0_106[k]
                   - f_6 * pff1_106[k]
                   + f_4 * pc_z[k] * pfg_161[k];

        t_228[k] = f_7 * pff0_107[k]
                   - f_8 * pff1_107[k]
                   + f_4 * pc_z[k] * pfg_162[k];

        t_229[k] = f_0 * sfg_14[k]
                   + f_4 * pc_y[k] * pfg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pa_y, pc_x, pc_y, sfh0_20, sfh1_20, pdg_105, \
                         pdg_106, pff0_110, pff0_111, pff1_110, pff1_111, pfg_165, \
                         pfg_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pa_y[k] * sfh0_20[k]
                   - f_9 * pc_y[k] * sfh1_20[k];

        t_231[k] = f_10 * pdg_105[k]
                   + f_2 * pff0_110[k]
                   - f_3 * pff1_110[k]
                   + f_4 * pc_x[k] * pfg_165[k];

        t_232[k] = f_10 * pdg_106[k]
                   + f_12 * pff0_111[k]
                   - f_13 * pff1_111[k]
                   + f_4 * pc_x[k] * pfg_166[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pc_x, pc_z, pdg_108, pdg_110, pff0_113, \
                         pff0_115, pff1_113, pff1_115, pfg_165, pfg_166, pfg_168, \
                         pfg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_4 * pc_z[k] * pfg_165[k];

        t_234[k] = f_10 * pdg_108[k]
                   + f_7 * pff0_113[k]
                   - f_8 * pff1_113[k]
                   + f_4 * pc_x[k] * pfg_168[k];

        t_235[k] = f_4 * pc_z[k] * pfg_166[k];

        t_236[k] = f_10 * pdg_110[k]
                   + f_7 * pff0_115[k]
                   - f_8 * pff1_115[k]
                   + f_4 * pc_x[k] * pfg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_x, pc_z, pdg_111, pdg_113, pff0_116, \
                         pff0_118, pff1_116, pff1_118, pfg_168, pfg_171, \
                         pfg_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_10 * pdg_111[k]
                   + f_5 * pff0_116[k]
                   - f_6 * pff1_116[k]
                   + f_4 * pc_x[k] * pfg_171[k];

        t_238[k] = f_4 * pc_z[k] * pfg_168[k];

        t_239[k] = f_10 * pdg_113[k]
                   + f_5 * pff0_118[k]
                   - f_6 * pff1_118[k]
                   + f_4 * pc_x[k] * pfg_173[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, pdg_114, pdg_115, pdg_116, pdg_117, \
                         pff0_119, pff1_119, pfg_174, pfg_175, pfg_176, \
                         pfg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * pdg_114[k]
                   + f_5 * pff0_119[k]
                   - f_6 * pff1_119[k]
                   + f_4 * pc_x[k] * pfg_174[k];

        t_241[k] = f_10 * pdg_115[k]
                   + f_4 * pc_x[k] * pfg_175[k];

        t_242[k] = f_10 * pdg_116[k]
                   + f_4 * pc_x[k] * pfg_176[k];

        t_243[k] = f_10 * pdg_117[k]
                   + f_4 * pc_x[k] * pfg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pc_x, pc_z, pph0_99, pph1_99, \
                         pdh0_162, pdg_118, pdg_119, pdh1_162, pfg_175, pfg_178, \
                         pfg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_10 * pdg_118[k]
                   + f_4 * pc_x[k] * pfg_178[k];

        t_245[k] = f_10 * pdg_119[k]
                   + f_4 * pc_x[k] * pfg_179[k];

        t_246[k] = f_14 * pph0_99[k]
                   - f_15 * pph1_99[k]
                   + pb_x[k] * pdh0_162[k]
                   - f_9 * pc_x[k] * pdh1_162[k];

        t_247[k] = f_4 * pc_z[k] * pfg_175[k];
    }
}

static auto
compute_prim_pfh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfh0, const size_t sfg,
                                                          const size_t sfh1, const size_t pdh0,
                                                          const size_t pdg, const size_t pdh1,
                                                          const size_t pff0, const size_t pff1,
                                                          const size_t pfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.0 / q;

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
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_42 = buffer.data(sfh0 + 42);
    const auto *sfh0_47 = buffer.data(sfh0 + 47);
    const auto *sfh0_50 = buffer.data(sfh0 + 50);
    const auto *sfh0_51 = buffer.data(sfh0 + 51);
    const auto *sfh0_59 = buffer.data(sfh0 + 59);
    const auto *sfh0_60 = buffer.data(sfh0 + 60);
    const auto *sfh0_62 = buffer.data(sfh0 + 62);
    const auto *sfh0_105 = buffer.data(sfh0 + 105);
    const auto *sfh0_106 = buffer.data(sfh0 + 106);
    const auto *sfh0_108 = buffer.data(sfh0 + 108);
    const auto *sfh0_110 = buffer.data(sfh0 + 110);
    const auto *sfh0_111 = buffer.data(sfh0 + 111);
    const auto *sfh0_113 = buffer.data(sfh0 + 113);
    const auto *sfh0_114 = buffer.data(sfh0 + 114);
    const auto *sfh0_125 = buffer.data(sfh0 + 125);

    const auto *sfg_29 = buffer.data(sfg + 29);
    const auto *sfg_35 = buffer.data(sfg + 35);
    const auto *sfg_42 = buffer.data(sfg + 42);
    const auto *sfg_43 = buffer.data(sfg + 43);
    const auto *sfg_44 = buffer.data(sfg + 44);
    const auto *sfg_75 = buffer.data(sfg + 75);
    const auto *sfg_76 = buffer.data(sfg + 76);
    const auto *sfg_78 = buffer.data(sfg + 78);
    const auto *sfg_80 = buffer.data(sfg + 80);
    const auto *sfg_89 = buffer.data(sfg + 89);
    const auto *sfg_100 = buffer.data(sfg + 100);
    const auto *sfg_104 = buffer.data(sfg + 104);

    const auto *sfh1_42 = buffer.data(sfh1 + 42);
    const auto *sfh1_47 = buffer.data(sfh1 + 47);
    const auto *sfh1_50 = buffer.data(sfh1 + 50);
    const auto *sfh1_51 = buffer.data(sfh1 + 51);
    const auto *sfh1_59 = buffer.data(sfh1 + 59);
    const auto *sfh1_60 = buffer.data(sfh1 + 60);
    const auto *sfh1_62 = buffer.data(sfh1 + 62);
    const auto *sfh1_105 = buffer.data(sfh1 + 105);
    const auto *sfh1_106 = buffer.data(sfh1 + 106);
    const auto *sfh1_108 = buffer.data(sfh1 + 108);
    const auto *sfh1_110 = buffer.data(sfh1 + 110);
    const auto *sfh1_111 = buffer.data(sfh1 + 111);
    const auto *sfh1_113 = buffer.data(sfh1 + 113);
    const auto *sfh1_114 = buffer.data(sfh1 + 114);
    const auto *sfh1_125 = buffer.data(sfh1 + 125);

    const auto *pdh0_127 = buffer.data(pdh0 + 127);
    const auto *pdh0_129 = buffer.data(pdh0 + 129);
    const auto *pdh0_132 = buffer.data(pdh0 + 132);
    const auto *pdh0_141 = buffer.data(pdh0 + 141);
    const auto *pdh0_148 = buffer.data(pdh0 + 148);
    const auto *pdh0_150 = buffer.data(pdh0 + 150);
    const auto *pdh0_153 = buffer.data(pdh0 + 153);
    const auto *pdh0_189 = buffer.data(pdh0 + 189);
    const auto *pdh0_190 = buffer.data(pdh0 + 190);
    const auto *pdh0_192 = buffer.data(pdh0 + 192);
    const auto *pdh0_194 = buffer.data(pdh0 + 194);
    const auto *pdh0_195 = buffer.data(pdh0 + 195);
    const auto *pdh0_197 = buffer.data(pdh0 + 197);
    const auto *pdh0_198 = buffer.data(pdh0 + 198);
    const auto *pdh0_204 = buffer.data(pdh0 + 204);
    const auto *pdh0_206 = buffer.data(pdh0 + 206);
    const auto *pdh0_207 = buffer.data(pdh0 + 207);
    const auto *pdh0_208 = buffer.data(pdh0 + 208);
    const auto *pdh0_209 = buffer.data(pdh0 + 209);
    const auto *pdh0_218 = buffer.data(pdh0 + 218);
    const auto *pdh0_225 = buffer.data(pdh0 + 225);
    const auto *pdh0_227 = buffer.data(pdh0 + 227);
    const auto *pdh0_228 = buffer.data(pdh0 + 228);
    const auto *pdh0_229 = buffer.data(pdh0 + 229);
    const auto *pdh0_230 = buffer.data(pdh0 + 230);
    const auto *pdh0_246 = buffer.data(pdh0 + 246);
    const auto *pdh0_248 = buffer.data(pdh0 + 248);
    const auto *pdh0_249 = buffer.data(pdh0 + 249);

    const auto *pdg_90 = buffer.data(pdg + 90);
    const auto *pdg_91 = buffer.data(pdg + 91);
    const auto *pdg_93 = buffer.data(pdg + 93);
    const auto *pdg_100 = buffer.data(pdg + 100);
    const auto *pdg_104 = buffer.data(pdg + 104);
    const auto *pdg_105 = buffer.data(pdg + 105);
    const auto *pdg_106 = buffer.data(pdg + 106);
    const auto *pdg_108 = buffer.data(pdg + 108);
    const auto *pdg_115 = buffer.data(pdg + 115);
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
    const auto *pdg_155 = buffer.data(pdg + 155);
    const auto *pdg_158 = buffer.data(pdg + 158);
    const auto *pdg_159 = buffer.data(pdg + 159);
    const auto *pdg_160 = buffer.data(pdg + 160);
    const auto *pdg_161 = buffer.data(pdg + 161);
    const auto *pdg_162 = buffer.data(pdg + 162);
    const auto *pdg_163 = buffer.data(pdg + 163);
    const auto *pdg_164 = buffer.data(pdg + 164);
    const auto *pdg_175 = buffer.data(pdg + 175);
    const auto *pdg_176 = buffer.data(pdg + 176);
    const auto *pdg_177 = buffer.data(pdg + 177);
    const auto *pdg_178 = buffer.data(pdg + 178);
    const auto *pdg_179 = buffer.data(pdg + 179);

    const auto *pdh1_127 = buffer.data(pdh1 + 127);
    const auto *pdh1_129 = buffer.data(pdh1 + 129);
    const auto *pdh1_132 = buffer.data(pdh1 + 132);
    const auto *pdh1_141 = buffer.data(pdh1 + 141);
    const auto *pdh1_148 = buffer.data(pdh1 + 148);
    const auto *pdh1_150 = buffer.data(pdh1 + 150);
    const auto *pdh1_153 = buffer.data(pdh1 + 153);
    const auto *pdh1_189 = buffer.data(pdh1 + 189);
    const auto *pdh1_190 = buffer.data(pdh1 + 190);
    const auto *pdh1_192 = buffer.data(pdh1 + 192);
    const auto *pdh1_194 = buffer.data(pdh1 + 194);
    const auto *pdh1_195 = buffer.data(pdh1 + 195);
    const auto *pdh1_197 = buffer.data(pdh1 + 197);
    const auto *pdh1_198 = buffer.data(pdh1 + 198);
    const auto *pdh1_204 = buffer.data(pdh1 + 204);
    const auto *pdh1_206 = buffer.data(pdh1 + 206);
    const auto *pdh1_207 = buffer.data(pdh1 + 207);
    const auto *pdh1_208 = buffer.data(pdh1 + 208);
    const auto *pdh1_209 = buffer.data(pdh1 + 209);
    const auto *pdh1_218 = buffer.data(pdh1 + 218);
    const auto *pdh1_225 = buffer.data(pdh1 + 225);
    const auto *pdh1_227 = buffer.data(pdh1 + 227);
    const auto *pdh1_228 = buffer.data(pdh1 + 228);
    const auto *pdh1_229 = buffer.data(pdh1 + 229);
    const auto *pdh1_230 = buffer.data(pdh1 + 230);
    const auto *pdh1_246 = buffer.data(pdh1 + 246);
    const auto *pdh1_248 = buffer.data(pdh1 + 248);
    const auto *pdh1_249 = buffer.data(pdh1 + 249);

    const auto *pff0_116 = buffer.data(pff0 + 116);
    const auto *pff0_117 = buffer.data(pff0 + 117);
    const auto *pff0_119 = buffer.data(pff0 + 119);
    const auto *pff0_140 = buffer.data(pff0 + 140);
    const auto *pff0_145 = buffer.data(pff0 + 145);
    const auto *pff0_149 = buffer.data(pff0 + 149);
    const auto *pff0_160 = buffer.data(pff0 + 160);
    const auto *pff0_161 = buffer.data(pff0 + 161);
    const auto *pff0_163 = buffer.data(pff0 + 163);
    const auto *pff0_165 = buffer.data(pff0 + 165);
    const auto *pff0_166 = buffer.data(pff0 + 166);
    const auto *pff0_167 = buffer.data(pff0 + 167);
    const auto *pff0_168 = buffer.data(pff0 + 168);
    const auto *pff0_169 = buffer.data(pff0 + 169);
    const auto *pff0_175 = buffer.data(pff0 + 175);
    const auto *pff0_178 = buffer.data(pff0 + 178);
    const auto *pff0_179 = buffer.data(pff0 + 179);

    const auto *pff1_116 = buffer.data(pff1 + 116);
    const auto *pff1_117 = buffer.data(pff1 + 117);
    const auto *pff1_119 = buffer.data(pff1 + 119);
    const auto *pff1_140 = buffer.data(pff1 + 140);
    const auto *pff1_145 = buffer.data(pff1 + 145);
    const auto *pff1_149 = buffer.data(pff1 + 149);
    const auto *pff1_160 = buffer.data(pff1 + 160);
    const auto *pff1_161 = buffer.data(pff1 + 161);
    const auto *pff1_163 = buffer.data(pff1 + 163);
    const auto *pff1_165 = buffer.data(pff1 + 165);
    const auto *pff1_166 = buffer.data(pff1 + 166);
    const auto *pff1_167 = buffer.data(pff1 + 167);
    const auto *pff1_168 = buffer.data(pff1 + 168);
    const auto *pff1_169 = buffer.data(pff1 + 169);
    const auto *pff1_175 = buffer.data(pff1 + 175);
    const auto *pff1_178 = buffer.data(pff1 + 178);
    const auto *pff1_179 = buffer.data(pff1 + 179);

    const auto *pfg_176 = buffer.data(pfg + 176);
    const auto *pfg_177 = buffer.data(pfg + 177);
    const auto *pfg_179 = buffer.data(pfg + 179);
    const auto *pfg_180 = buffer.data(pfg + 180);
    const auto *pfg_181 = buffer.data(pfg + 181);
    const auto *pfg_183 = buffer.data(pfg + 183);
    const auto *pfg_190 = buffer.data(pfg + 190);
    const auto *pfg_191 = buffer.data(pfg + 191);
    const auto *pfg_192 = buffer.data(pfg + 192);
    const auto *pfg_193 = buffer.data(pfg + 193);
    const auto *pfg_194 = buffer.data(pfg + 194);
    const auto *pfg_195 = buffer.data(pfg + 195);
    const auto *pfg_196 = buffer.data(pfg + 196);
    const auto *pfg_198 = buffer.data(pfg + 198);
    const auto *pfg_205 = buffer.data(pfg + 205);
    const auto *pfg_206 = buffer.data(pfg + 206);
    const auto *pfg_207 = buffer.data(pfg + 207);
    const auto *pfg_208 = buffer.data(pfg + 208);
    const auto *pfg_209 = buffer.data(pfg + 209);
    const auto *pfg_210 = buffer.data(pfg + 210);
    const auto *pfg_211 = buffer.data(pfg + 211);
    const auto *pfg_213 = buffer.data(pfg + 213);
    const auto *pfg_215 = buffer.data(pfg + 215);
    const auto *pfg_219 = buffer.data(pfg + 219);
    const auto *pfg_220 = buffer.data(pfg + 220);
    const auto *pfg_221 = buffer.data(pfg + 221);
    const auto *pfg_222 = buffer.data(pfg + 222);
    const auto *pfg_223 = buffer.data(pfg + 223);
    const auto *pfg_224 = buffer.data(pfg + 224);
    const auto *pfg_225 = buffer.data(pfg + 225);
    const auto *pfg_226 = buffer.data(pfg + 226);
    const auto *pfg_228 = buffer.data(pfg + 228);
    const auto *pfg_235 = buffer.data(pfg + 235);
    const auto *pfg_236 = buffer.data(pfg + 236);
    const auto *pfg_237 = buffer.data(pfg + 237);
    const auto *pfg_238 = buffer.data(pfg + 238);
    const auto *pfg_239 = buffer.data(pfg + 239);
    const auto *pfg_240 = buffer.data(pfg + 240);
    const auto *pfg_241 = buffer.data(pfg + 241);
    const auto *pfg_243 = buffer.data(pfg + 243);
    const auto *pfg_245 = buffer.data(pfg + 245);
    const auto *pfg_246 = buffer.data(pfg + 246);
    const auto *pfg_248 = buffer.data(pfg + 248);
    const auto *pfg_249 = buffer.data(pfg + 249);
    const auto *pfg_250 = buffer.data(pfg + 250);
    const auto *pfg_251 = buffer.data(pfg + 251);
    const auto *pfg_252 = buffer.data(pfg + 252);
    const auto *pfg_253 = buffer.data(pfg + 253);
    const auto *pfg_254 = buffer.data(pfg + 254);
    const auto *pfg_255 = buffer.data(pfg + 255);
    const auto *pfg_256 = buffer.data(pfg + 256);
    const auto *pfg_258 = buffer.data(pfg + 258);
    const auto *pfg_260 = buffer.data(pfg + 260);
    const auto *pfg_263 = buffer.data(pfg + 263);
    const auto *pfg_264 = buffer.data(pfg + 264);
    const auto *pfg_265 = buffer.data(pfg + 265);
    const auto *pfg_266 = buffer.data(pfg + 266);
    const auto *pfg_267 = buffer.data(pfg + 267);
    const auto *pfg_268 = buffer.data(pfg + 268);
    const auto *pfg_269 = buffer.data(pfg + 269);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, pc_z, sfg_29, pdg_104, pff0_116, pff0_117, \
                         pff1_116, pff1_117, pfg_176, pfg_177, \
                         pfg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_5 * pff0_116[k]
                   - f_6 * pff1_116[k]
                   + f_4 * pc_z[k] * pfg_176[k];

        t_249[k] = f_7 * pff0_117[k]
                   - f_8 * pff1_117[k]
                   + f_4 * pc_z[k] * pfg_177[k];

        t_250[k] = f_0 * sfg_29[k]
                   + f_0 * pdg_104[k]
                   + f_4 * pc_y[k] * pfg_179[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, pa_y, pb_z, pc_y, pc_z, sfh0_42, sfh1_42, \
                         pdh0_127, pdh1_127, pff0_119, pff1_119, \
                         pfg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_2 * pff0_119[k]
                   - f_3 * pff1_119[k]
                   + f_4 * pc_z[k] * pfg_179[k];

        t_252[k] = pa_y[k] * sfh0_42[k]
                   - f_9 * pc_y[k] * sfh1_42[k];

        t_253[k] = pb_z[k] * pdh0_127[k]
                   - f_9 * pc_z[k] * pdh1_127[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pa_y, pb_z, pc_y, pc_z, sfh0_47, sfh1_47, \
                         pdh0_129, pdg_90, pdg_91, pdh1_129, pfg_180, \
                         pfg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_0 * pdg_90[k]
                   + f_4 * pc_z[k] * pfg_180[k];

        t_255[k] = pb_z[k] * pdh0_129[k]
                   - f_9 * pc_z[k] * pdh1_129[k];

        t_256[k] = f_0 * pdg_91[k]
                   + f_4 * pc_z[k] * pfg_181[k];

        t_257[k] = pa_y[k] * sfh0_47[k]
                   - f_9 * pc_y[k] * sfh1_47[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_y, pb_z, pc_y, pc_z, sfh0_50, sfg_35, \
                         sfh1_50, pdh0_132, pdg_93, pdh1_132, pfg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pb_z[k] * pdh0_132[k]
                   - f_9 * pc_z[k] * pdh1_132[k];

        t_259[k] = f_0 * pdg_93[k]
                   + f_4 * pc_z[k] * pfg_183[k];

        t_260[k] = pa_y[k] * sfh0_50[k]
                   + f_0 * sfg_35[k]
                   - f_9 * pc_y[k] * sfh1_50[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pc_x, pc_y, sfh0_51, sfh1_51, \
                         pdg_130, pdg_131, pdg_132, pfg_190, pfg_191, \
                         pfg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_y[k] * sfh0_51[k]
                   - f_9 * pc_y[k] * sfh1_51[k];

        t_262[k] = f_10 * pdg_130[k]
                   + f_4 * pc_x[k] * pfg_190[k];

        t_263[k] = f_10 * pdg_131[k]
                   + f_4 * pc_x[k] * pfg_191[k];

        t_264[k] = f_10 * pdg_132[k]
                   + f_4 * pc_x[k] * pfg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_z, pc_x, pc_z, pdh0_141, pdg_100, \
                         pdg_133, pdg_134, pdh1_141, pfg_190, pfg_193, \
                         pfg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_10 * pdg_133[k]
                   + f_4 * pc_x[k] * pfg_193[k];

        t_266[k] = f_10 * pdg_134[k]
                   + f_4 * pc_x[k] * pfg_194[k];

        t_267[k] = pb_z[k] * pdh0_141[k]
                   - f_9 * pc_z[k] * pdh1_141[k];

        t_268[k] = f_0 * pdg_100[k]
                   + f_4 * pc_z[k] * pfg_190[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_y, pc_y, sfh0_59, sfh0_60, sfh0_62, \
                         sfg_42, sfg_43, sfg_44, sfh1_59, sfh1_60, sfh1_62, \
                         pfg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_y[k] * sfh0_59[k]
                   + f_1 * sfg_42[k]
                   - f_9 * pc_y[k] * sfh1_59[k];

        t_270[k] = pa_y[k] * sfh0_60[k]
                   + f_10 * sfg_43[k]
                   - f_9 * pc_y[k] * sfh1_60[k];

        t_271[k] = f_0 * sfg_44[k]
                   + f_4 * pc_y[k] * pfg_194[k];

        t_272[k] = pa_y[k] * sfh0_62[k]
                   - f_9 * pc_y[k] * sfh1_62[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pb_x, pc_x, pc_z, pdh0_189, pdh0_190, pdg_135, \
                         pdg_136, pdh1_189, pdh1_190, pfg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = pb_x[k] * pdh0_189[k]
                   + f_11 * pdg_135[k]
                   - f_9 * pc_x[k] * pdh1_189[k];

        t_274[k] = pb_x[k] * pdh0_190[k]
                   + f_16 * pdg_136[k]
                   - f_9 * pc_x[k] * pdh1_190[k];

        t_275[k] = f_4 * pc_z[k] * pfg_195[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pb_x, pc_x, pc_z, pdh0_192, pdh0_194, pdg_138, \
                         pdg_140, pdh1_192, pdh1_194, pfg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_x[k] * pdh0_192[k]
                   + f_1 * pdg_138[k]
                   - f_9 * pc_x[k] * pdh1_192[k];

        t_277[k] = f_4 * pc_z[k] * pfg_196[k];

        t_278[k] = pb_x[k] * pdh0_194[k]
                   + f_1 * pdg_140[k]
                   - f_9 * pc_x[k] * pdh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pb_x, pc_x, pc_z, pdh0_195, pdh0_197, pdg_141, \
                         pdg_143, pdh1_195, pdh1_197, pfg_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_x[k] * pdh0_195[k]
                   + f_10 * pdg_141[k]
                   - f_9 * pc_x[k] * pdh1_195[k];

        t_280[k] = f_4 * pc_z[k] * pfg_198[k];

        t_281[k] = pb_x[k] * pdh0_197[k]
                   + f_10 * pdg_143[k]
                   - f_9 * pc_x[k] * pdh1_197[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_x, pc_x, pdh0_198, pdg_144, pdg_145, \
                         pdg_146, pdg_147, pdh1_198, pfg_205, pfg_206, \
                         pfg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pb_x[k] * pdh0_198[k]
                   + f_10 * pdg_144[k]
                   - f_9 * pc_x[k] * pdh1_198[k];

        t_283[k] = f_0 * pdg_145[k]
                   + f_4 * pc_x[k] * pfg_205[k];

        t_284[k] = f_0 * pdg_146[k]
                   + f_4 * pc_x[k] * pfg_206[k];

        t_285[k] = f_0 * pdg_147[k]
                   + f_4 * pc_x[k] * pfg_207[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_x, pc_x, pc_z, pdh0_204, pdg_148, \
                         pdg_149, pdh1_204, pfg_205, pfg_208, pfg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * pdg_148[k]
                   + f_4 * pc_x[k] * pfg_208[k];

        t_287[k] = f_0 * pdg_149[k]
                   + f_4 * pc_x[k] * pfg_209[k];

        t_288[k] = pb_x[k] * pdh0_204[k]
                   - f_9 * pc_x[k] * pdh1_204[k];

        t_289[k] = f_4 * pc_z[k] * pfg_205[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_x, pc_x, pdh0_206, pdh0_207, pdh0_208, \
                         pdh0_209, pdh1_206, pdh1_207, pdh1_208, \
                         pdh1_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * pdh0_206[k]
                   - f_9 * pc_x[k] * pdh1_206[k];

        t_291[k] = pb_x[k] * pdh0_207[k]
                   - f_9 * pc_x[k] * pdh1_207[k];

        t_292[k] = pb_x[k] * pdh0_208[k]
                   - f_9 * pc_x[k] * pdh1_208[k];

        t_293[k] = pb_x[k] * pdh0_209[k]
                   - f_9 * pc_x[k] * pdh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_z, pc_x, pc_z, pdh0_148, pdh0_150, \
                         pdg_105, pdg_150, pdh1_148, pdh1_150, pff0_140, pff1_140, \
                         pfg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_0 * pdg_150[k]
                   + f_2 * pff0_140[k]
                   - f_3 * pff1_140[k]
                   + f_4 * pc_x[k] * pfg_210[k];

        t_295[k] = pb_z[k] * pdh0_148[k]
                   - f_9 * pc_z[k] * pdh1_148[k];

        t_296[k] = f_0 * pdg_105[k]
                   + f_4 * pc_z[k] * pfg_210[k];

        t_297[k] = pb_z[k] * pdh0_150[k]
                   - f_9 * pc_z[k] * pdh1_150[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_z, pc_x, pc_z, pdh0_153, pdg_106, pdg_155, \
                         pdh1_153, pff0_145, pff1_145, pfg_211, \
                         pfg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * pdg_106[k]
                   + f_4 * pc_z[k] * pfg_211[k];

        t_299[k] = f_0 * pdg_155[k]
                   + f_7 * pff0_145[k]
                   - f_8 * pff1_145[k]
                   + f_4 * pc_x[k] * pfg_215[k];

        t_300[k] = pb_z[k] * pdh0_153[k]
                   - f_9 * pc_z[k] * pdh1_153[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, pb_x, pc_x, pc_z, pdh0_218, pdg_108, pdg_158, \
                         pdg_159, pdh1_218, pff0_149, pff1_149, pfg_213, \
                         pfg_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_0 * pdg_108[k]
                   + f_4 * pc_z[k] * pfg_213[k];

        t_302[k] = pb_x[k] * pdh0_218[k]
                   + f_10 * pdg_158[k]
                   - f_9 * pc_x[k] * pdh1_218[k];

        t_303[k] = f_0 * pdg_159[k]
                   + f_5 * pff0_149[k]
                   - f_6 * pff1_149[k]
                   + f_4 * pc_x[k] * pfg_219[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pc_x, pdg_160, pdg_161, pdg_162, \
                         pdg_163, pdg_164, pfg_220, pfg_221, pfg_222, pfg_223, \
                         pfg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_0 * pdg_160[k]
                   + f_4 * pc_x[k] * pfg_220[k];

        t_305[k] = f_0 * pdg_161[k]
                   + f_4 * pc_x[k] * pfg_221[k];

        t_306[k] = f_0 * pdg_162[k]
                   + f_4 * pc_x[k] * pfg_222[k];

        t_307[k] = f_0 * pdg_163[k]
                   + f_4 * pc_x[k] * pfg_223[k];

        t_308[k] = f_0 * pdg_164[k]
                   + f_4 * pc_x[k] * pfg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pc_x, pc_z, pdh0_225, pdh0_227, \
                         pdh0_228, pdg_115, pdh1_225, pdh1_227, pdh1_228, \
                         pfg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pb_x[k] * pdh0_225[k]
                   - f_9 * pc_x[k] * pdh1_225[k];

        t_310[k] = f_0 * pdg_115[k]
                   + f_4 * pc_z[k] * pfg_220[k];

        t_311[k] = pb_x[k] * pdh0_227[k]
                   - f_9 * pc_x[k] * pdh1_227[k];

        t_312[k] = pb_x[k] * pdh0_228[k]
                   - f_9 * pc_x[k] * pdh1_228[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_y, pb_x, pc_x, pc_y, sfh0_105, sfh1_105, \
                         pdh0_229, pdh0_230, pdh1_229, pdh1_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_x[k] * pdh0_229[k]
                   - f_9 * pc_x[k] * pdh1_229[k];

        t_314[k] = pb_x[k] * pdh0_230[k]
                   - f_9 * pc_x[k] * pdh1_230[k];

        t_315[k] = pa_y[k] * sfh0_105[k]
                   - f_9 * pc_y[k] * sfh1_105[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_y, pc_y, pc_z, sfh0_106, sfh0_108, sfg_75, \
                         sfg_76, sfh1_106, sfh1_108, pdg_120, pfg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_y[k] * sfh0_106[k]
                   + f_0 * sfg_75[k]
                   - f_9 * pc_y[k] * sfh1_106[k];

        t_317[k] = f_10 * pdg_120[k]
                   + f_4 * pc_z[k] * pfg_225[k];

        t_318[k] = pa_y[k] * sfh0_108[k]
                   + f_10 * sfg_76[k]
                   - f_9 * pc_y[k] * sfh1_108[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pa_y, pc_y, pc_z, sfh0_110, sfh0_111, \
                         sfg_78, sfh1_110, sfh1_111, pdg_121, pdg_123, pfg_226, \
                         pfg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_10 * pdg_121[k]
                   + f_4 * pc_z[k] * pfg_226[k];

        t_320[k] = pa_y[k] * sfh0_110[k]
                   - f_9 * pc_y[k] * sfh1_110[k];

        t_321[k] = pa_y[k] * sfh0_111[k]
                   + f_1 * sfg_78[k]
                   - f_9 * pc_y[k] * sfh1_111[k];

        t_322[k] = f_10 * pdg_123[k]
                   + f_4 * pc_z[k] * pfg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pa_y, pc_x, pc_y, sfh0_113, sfh0_114, \
                         sfg_80, sfh1_113, sfh1_114, pdg_175, pdg_176, pfg_235, \
                         pfg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = pa_y[k] * sfh0_113[k]
                   + f_0 * sfg_80[k]
                   - f_9 * pc_y[k] * sfh1_113[k];

        t_324[k] = pa_y[k] * sfh0_114[k]
                   - f_9 * pc_y[k] * sfh1_114[k];

        t_325[k] = f_0 * pdg_175[k]
                   + f_4 * pc_x[k] * pfg_235[k];

        t_326[k] = f_0 * pdg_176[k]
                   + f_4 * pc_x[k] * pfg_236[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pb_x, pc_x, pdh0_246, pdg_177, pdg_178, \
                         pdg_179, pdh1_246, pfg_237, pfg_238, pfg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_0 * pdg_177[k]
                   + f_4 * pc_x[k] * pfg_237[k];

        t_328[k] = f_0 * pdg_178[k]
                   + f_4 * pc_x[k] * pfg_238[k];

        t_329[k] = f_0 * pdg_179[k]
                   + f_4 * pc_x[k] * pfg_239[k];

        t_330[k] = pb_x[k] * pdh0_246[k]
                   - f_9 * pc_x[k] * pdh1_246[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_x, pc_x, pc_y, pc_z, sfg_89, pdh0_248, \
                         pdh0_249, pdg_130, pdh1_248, pdh1_249, pfg_235, \
                         pfg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_10 * pdg_130[k]
                   + f_4 * pc_z[k] * pfg_235[k];

        t_332[k] = pb_x[k] * pdh0_248[k]
                   - f_9 * pc_x[k] * pdh1_248[k];

        t_333[k] = pb_x[k] * pdh0_249[k]
                   - f_9 * pc_x[k] * pdh1_249[k];

        t_334[k] = f_0 * sfg_89[k]
                   + f_4 * pc_y[k] * pfg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_y, pc_x, pc_y, pc_z, sfh0_125, \
                         sfh1_125, pff0_160, pff0_161, pff1_160, pff1_161, pfg_240, \
                         pfg_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pa_y[k] * sfh0_125[k]
                   - f_9 * pc_y[k] * sfh1_125[k];

        t_336[k] = f_2 * pff0_160[k]
                   - f_3 * pff1_160[k]
                   + f_4 * pc_x[k] * pfg_240[k];

        t_337[k] = f_12 * pff0_161[k]
                   - f_13 * pff1_161[k]
                   + f_4 * pc_x[k] * pfg_241[k];

        t_338[k] = f_4 * pc_z[k] * pfg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pc_x, pc_z, pff0_163, pff0_165, pff0_166, \
                         pff1_163, pff1_165, pff1_166, pfg_241, pfg_243, pfg_245, \
                         pfg_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_7 * pff0_163[k]
                   - f_8 * pff1_163[k]
                   + f_4 * pc_x[k] * pfg_243[k];

        t_340[k] = f_4 * pc_z[k] * pfg_241[k];

        t_341[k] = f_7 * pff0_165[k]
                   - f_8 * pff1_165[k]
                   + f_4 * pc_x[k] * pfg_245[k];

        t_342[k] = f_5 * pff0_166[k]
                   - f_6 * pff1_166[k]
                   + f_4 * pc_x[k] * pfg_246[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, pc_x, pc_z, pff0_168, pff0_169, \
                         pff1_168, pff1_169, pfg_243, pfg_248, pfg_249, pfg_250, \
                         pfg_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_4 * pc_z[k] * pfg_243[k];

        t_344[k] = f_5 * pff0_168[k]
                   - f_6 * pff1_168[k]
                   + f_4 * pc_x[k] * pfg_248[k];

        t_345[k] = f_5 * pff0_169[k]
                   - f_6 * pff1_169[k]
                   + f_4 * pc_x[k] * pfg_249[k];

        t_346[k] = f_4 * pc_x[k] * pfg_250[k];

        t_347[k] = f_4 * pc_x[k] * pfg_251[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, pc_x, pc_y, pc_z, sfg_100, \
                         pdg_145, pff0_166, pff1_166, pfg_250, pfg_252, pfg_253, \
                         pfg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_4 * pc_x[k] * pfg_252[k];

        t_349[k] = f_4 * pc_x[k] * pfg_253[k];

        t_350[k] = f_4 * pc_x[k] * pfg_254[k];

        t_351[k] = f_0 * sfg_100[k]
                   + f_1 * pdg_145[k]
                   + f_2 * pff0_166[k]
                   - f_3 * pff1_166[k]
                   + f_4 * pc_y[k] * pfg_250[k];

        t_352[k] = f_4 * pc_z[k] * pfg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, pc_z, sfg_104, pdg_149, pff0_166, \
                         pff0_167, pff1_166, pff1_167, pfg_251, pfg_252, \
                         pfg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_5 * pff0_166[k]
                   - f_6 * pff1_166[k]
                   + f_4 * pc_z[k] * pfg_251[k];

        t_354[k] = f_7 * pff0_167[k]
                   - f_8 * pff1_167[k]
                   + f_4 * pc_z[k] * pfg_252[k];

        t_355[k] = f_0 * sfg_104[k]
                   + f_1 * pdg_149[k]
                   + f_4 * pc_y[k] * pfg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pb_z, pc_z, pdh0_189, pdh0_190, pdg_135, \
                         pdh1_189, pdh1_190, pff0_169, pff1_169, pfg_254, \
                         pfg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_2 * pff0_169[k]
                   - f_3 * pff1_169[k]
                   + f_4 * pc_z[k] * pfg_254[k];

        t_357[k] = pb_z[k] * pdh0_189[k]
                   - f_9 * pc_z[k] * pdh1_189[k];

        t_358[k] = pb_z[k] * pdh0_190[k]
                   - f_9 * pc_z[k] * pdh1_190[k];

        t_359[k] = f_0 * pdg_135[k]
                   + f_4 * pc_z[k] * pfg_255[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pb_z, pc_x, pc_z, pdh0_192, pdh0_195, \
                         pdg_136, pdh1_192, pdh1_195, pff0_175, pff1_175, pfg_256, \
                         pfg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pb_z[k] * pdh0_192[k]
                   - f_9 * pc_z[k] * pdh1_192[k];

        t_361[k] = f_0 * pdg_136[k]
                   + f_4 * pc_z[k] * pfg_256[k];

        t_362[k] = f_7 * pff0_175[k]
                   - f_8 * pff1_175[k]
                   + f_4 * pc_x[k] * pfg_260[k];

        t_363[k] = pb_z[k] * pdh0_195[k]
                   - f_9 * pc_z[k] * pdh1_195[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pc_x, pc_z, pdg_138, pff0_178, pff0_179, \
                         pff1_178, pff1_179, pfg_258, pfg_263, pfg_264, \
                         pfg_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_0 * pdg_138[k]
                   + f_4 * pc_z[k] * pfg_258[k];

        t_365[k] = f_5 * pff0_178[k]
                   - f_6 * pff1_178[k]
                   + f_4 * pc_x[k] * pfg_263[k];

        t_366[k] = f_5 * pff0_179[k]
                   - f_6 * pff1_179[k]
                   + f_4 * pc_x[k] * pfg_264[k];

        t_367[k] = f_4 * pc_x[k] * pfg_265[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, pb_z, pc_x, pc_z, pdh0_204, \
                         pdh1_204, pfg_266, pfg_267, pfg_268, pfg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_4 * pc_x[k] * pfg_266[k];

        t_369[k] = f_4 * pc_x[k] * pfg_267[k];

        t_370[k] = f_4 * pc_x[k] * pfg_268[k];

        t_371[k] = f_4 * pc_x[k] * pfg_269[k];

        t_372[k] = pb_z[k] * pdh0_204[k]
                   - f_9 * pc_z[k] * pdh1_204[k];
    }
}

static auto
compute_prim_pfh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfh0, const size_t sfg,
                                                          const size_t sfh1, const size_t pph0,
                                                          const size_t pph1, const size_t pdh0,
                                                          const size_t pdg, const size_t pdh1,
                                                          const size_t pff0, const size_t pff1,
                                                          const size_t pfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_0 = buffer.data(sfh0 + 0);
    const auto *sfh0_2 = buffer.data(sfh0 + 2);
    const auto *sfh0_3 = buffer.data(sfh0 + 3);
    const auto *sfh0_5 = buffer.data(sfh0 + 5);
    const auto *sfh0_6 = buffer.data(sfh0 + 6);
    const auto *sfh0_7 = buffer.data(sfh0 + 7);
    const auto *sfh0_9 = buffer.data(sfh0 + 9);
    const auto *sfh0_15 = buffer.data(sfh0 + 15);
    const auto *sfh0_20 = buffer.data(sfh0 + 20);
    const auto *sfh0_21 = buffer.data(sfh0 + 21);
    const auto *sfh0_24 = buffer.data(sfh0 + 24);
    const auto *sfh0_27 = buffer.data(sfh0 + 27);
    const auto *sfh0_28 = buffer.data(sfh0 + 28);
    const auto *sfh0_36 = buffer.data(sfh0 + 36);
    const auto *sfh0_37 = buffer.data(sfh0 + 37);
    const auto *sfh0_38 = buffer.data(sfh0 + 38);
    const auto *sfh0_39 = buffer.data(sfh0 + 39);
    const auto *sfh0_63 = buffer.data(sfh0 + 63);
    const auto *sfh0_65 = buffer.data(sfh0 + 65);
    const auto *sfh0_66 = buffer.data(sfh0 + 66);
    const auto *sfh0_68 = buffer.data(sfh0 + 68);
    const auto *sfh0_69 = buffer.data(sfh0 + 69);
    const auto *sfh0_70 = buffer.data(sfh0 + 70);
    const auto *sfh0_189 = buffer.data(sfh0 + 189);
    const auto *sfh0_194 = buffer.data(sfh0 + 194);
    const auto *sfh0_198 = buffer.data(sfh0 + 198);
    const auto *sfh0_204 = buffer.data(sfh0 + 204);
    const auto *sfh0_206 = buffer.data(sfh0 + 206);
    const auto *sfh0_207 = buffer.data(sfh0 + 207);
    const auto *sfh0_209 = buffer.data(sfh0 + 209);

    const auto *sfg_0 = buffer.data(sfg + 0);
    const auto *sfg_2 = buffer.data(sfg + 2);
    const auto *sfg_3 = buffer.data(sfg + 3);
    const auto *sfg_5 = buffer.data(sfg + 5);
    const auto *sfg_14 = buffer.data(sfg + 14);
    const auto *sfg_18 = buffer.data(sfg + 18);
    const auto *sfg_25 = buffer.data(sfg + 25);
    const auto *sfg_26 = buffer.data(sfg + 26);
    const auto *sfg_27 = buffer.data(sfg + 27);
    const auto *sfg_45 = buffer.data(sfg + 45);
    const auto *sfg_47 = buffer.data(sfg + 47);
    const auto *sfg_48 = buffer.data(sfg + 48);
    const auto *sfg_119 = buffer.data(sfg + 119);
    const auto *sfg_130 = buffer.data(sfg + 130);
    const auto *sfg_134 = buffer.data(sfg + 134);
    const auto *sfg_145 = buffer.data(sfg + 145);
    const auto *sfg_147 = buffer.data(sfg + 147);
    const auto *sfg_148 = buffer.data(sfg + 148);
    const auto *sfg_149 = buffer.data(sfg + 149);

    const auto *sfh1_0 = buffer.data(sfh1 + 0);
    const auto *sfh1_2 = buffer.data(sfh1 + 2);
    const auto *sfh1_3 = buffer.data(sfh1 + 3);
    const auto *sfh1_5 = buffer.data(sfh1 + 5);
    const auto *sfh1_6 = buffer.data(sfh1 + 6);
    const auto *sfh1_7 = buffer.data(sfh1 + 7);
    const auto *sfh1_9 = buffer.data(sfh1 + 9);
    const auto *sfh1_15 = buffer.data(sfh1 + 15);
    const auto *sfh1_20 = buffer.data(sfh1 + 20);
    const auto *sfh1_21 = buffer.data(sfh1 + 21);
    const auto *sfh1_24 = buffer.data(sfh1 + 24);
    const auto *sfh1_27 = buffer.data(sfh1 + 27);
    const auto *sfh1_28 = buffer.data(sfh1 + 28);
    const auto *sfh1_36 = buffer.data(sfh1 + 36);
    const auto *sfh1_37 = buffer.data(sfh1 + 37);
    const auto *sfh1_38 = buffer.data(sfh1 + 38);
    const auto *sfh1_39 = buffer.data(sfh1 + 39);
    const auto *sfh1_63 = buffer.data(sfh1 + 63);
    const auto *sfh1_65 = buffer.data(sfh1 + 65);
    const auto *sfh1_66 = buffer.data(sfh1 + 66);
    const auto *sfh1_68 = buffer.data(sfh1 + 68);
    const auto *sfh1_69 = buffer.data(sfh1 + 69);
    const auto *sfh1_70 = buffer.data(sfh1 + 70);
    const auto *sfh1_189 = buffer.data(sfh1 + 189);
    const auto *sfh1_194 = buffer.data(sfh1 + 194);
    const auto *sfh1_198 = buffer.data(sfh1 + 198);
    const auto *sfh1_204 = buffer.data(sfh1 + 204);
    const auto *sfh1_206 = buffer.data(sfh1 + 206);
    const auto *sfh1_207 = buffer.data(sfh1 + 207);
    const auto *sfh1_209 = buffer.data(sfh1 + 209);

    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *pdh0_206 = buffer.data(pdh0 + 206);
    const auto *pdh0_207 = buffer.data(pdh0 + 207);
    const auto *pdh0_254 = buffer.data(pdh0 + 254);
    const auto *pdh0_257 = buffer.data(pdh0 + 257);
    const auto *pdh0_261 = buffer.data(pdh0 + 261);
    const auto *pdh0_272 = buffer.data(pdh0 + 272);
    const auto *pdh0_314 = buffer.data(pdh0 + 314);

    const auto *pdg_145 = buffer.data(pdg + 145);
    const auto *pdg_146 = buffer.data(pdg + 146);
    const auto *pdg_147 = buffer.data(pdg + 147);
    const auto *pdg_149 = buffer.data(pdg + 149);
    const auto *pdg_150 = buffer.data(pdg + 150);
    const auto *pdg_151 = buffer.data(pdg + 151);
    const auto *pdg_153 = buffer.data(pdg + 153);
    const auto *pdg_160 = buffer.data(pdg + 160);
    const auto *pdg_161 = buffer.data(pdg + 161);
    const auto *pdg_162 = buffer.data(pdg + 162);
    const auto *pdg_164 = buffer.data(pdg + 164);
    const auto *pdg_165 = buffer.data(pdg + 165);
    const auto *pdg_166 = buffer.data(pdg + 166);
    const auto *pdg_168 = buffer.data(pdg + 168);
    const auto *pdg_175 = buffer.data(pdg + 175);
    const auto *pdg_179 = buffer.data(pdg + 179);
    const auto *pdg_180 = buffer.data(pdg + 180);
    const auto *pdg_182 = buffer.data(pdg + 182);
    const auto *pdg_185 = buffer.data(pdg + 185);
    const auto *pdg_190 = buffer.data(pdg + 190);
    const auto *pdg_191 = buffer.data(pdg + 191);
    const auto *pdg_192 = buffer.data(pdg + 192);
    const auto *pdg_193 = buffer.data(pdg + 193);
    const auto *pdg_194 = buffer.data(pdg + 194);
    const auto *pdg_195 = buffer.data(pdg + 195);
    const auto *pdg_197 = buffer.data(pdg + 197);
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
    const auto *pdg_217 = buffer.data(pdg + 217);
    const auto *pdg_219 = buffer.data(pdg + 219);
    const auto *pdg_220 = buffer.data(pdg + 220);
    const auto *pdg_221 = buffer.data(pdg + 221);
    const auto *pdg_222 = buffer.data(pdg + 222);
    const auto *pdg_223 = buffer.data(pdg + 223);
    const auto *pdg_224 = buffer.data(pdg + 224);

    const auto *pdh1_206 = buffer.data(pdh1 + 206);
    const auto *pdh1_207 = buffer.data(pdh1 + 207);
    const auto *pdh1_254 = buffer.data(pdh1 + 254);
    const auto *pdh1_257 = buffer.data(pdh1 + 257);
    const auto *pdh1_261 = buffer.data(pdh1 + 261);
    const auto *pdh1_272 = buffer.data(pdh1 + 272);
    const auto *pdh1_314 = buffer.data(pdh1 + 314);

    const auto *pff0_179 = buffer.data(pff0 + 179);
    const auto *pff0_180 = buffer.data(pff0 + 180);
    const auto *pff0_181 = buffer.data(pff0 + 181);
    const auto *pff0_183 = buffer.data(pff0 + 183);
    const auto *pff0_185 = buffer.data(pff0 + 185);
    const auto *pff0_186 = buffer.data(pff0 + 186);
    const auto *pff0_187 = buffer.data(pff0 + 187);
    const auto *pff0_188 = buffer.data(pff0 + 188);
    const auto *pff0_189 = buffer.data(pff0 + 189);
    const auto *pff0_191 = buffer.data(pff0 + 191);
    const auto *pff0_193 = buffer.data(pff0 + 193);
    const auto *pff0_196 = buffer.data(pff0 + 196);
    const auto *pff0_198 = buffer.data(pff0 + 198);
    const auto *pff0_207 = buffer.data(pff0 + 207);
    const auto *pff0_208 = buffer.data(pff0 + 208);
    const auto *pff0_209 = buffer.data(pff0 + 209);
    const auto *pff0_220 = buffer.data(pff0 + 220);
    const auto *pff0_222 = buffer.data(pff0 + 222);
    const auto *pff0_223 = buffer.data(pff0 + 223);
    const auto *pff0_225 = buffer.data(pff0 + 225);
    const auto *pff0_226 = buffer.data(pff0 + 226);
    const auto *pff0_227 = buffer.data(pff0 + 227);
    const auto *pff0_228 = buffer.data(pff0 + 228);
    const auto *pff0_229 = buffer.data(pff0 + 229);

    const auto *pff1_179 = buffer.data(pff1 + 179);
    const auto *pff1_180 = buffer.data(pff1 + 180);
    const auto *pff1_181 = buffer.data(pff1 + 181);
    const auto *pff1_183 = buffer.data(pff1 + 183);
    const auto *pff1_185 = buffer.data(pff1 + 185);
    const auto *pff1_186 = buffer.data(pff1 + 186);
    const auto *pff1_187 = buffer.data(pff1 + 187);
    const auto *pff1_188 = buffer.data(pff1 + 188);
    const auto *pff1_189 = buffer.data(pff1 + 189);
    const auto *pff1_191 = buffer.data(pff1 + 191);
    const auto *pff1_193 = buffer.data(pff1 + 193);
    const auto *pff1_196 = buffer.data(pff1 + 196);
    const auto *pff1_198 = buffer.data(pff1 + 198);
    const auto *pff1_207 = buffer.data(pff1 + 207);
    const auto *pff1_208 = buffer.data(pff1 + 208);
    const auto *pff1_209 = buffer.data(pff1 + 209);
    const auto *pff1_220 = buffer.data(pff1 + 220);
    const auto *pff1_222 = buffer.data(pff1 + 222);
    const auto *pff1_223 = buffer.data(pff1 + 223);
    const auto *pff1_225 = buffer.data(pff1 + 225);
    const auto *pff1_226 = buffer.data(pff1 + 226);
    const auto *pff1_227 = buffer.data(pff1 + 227);
    const auto *pff1_228 = buffer.data(pff1 + 228);
    const auto *pff1_229 = buffer.data(pff1 + 229);

    const auto *pfg_265 = buffer.data(pfg + 265);
    const auto *pfg_269 = buffer.data(pfg + 269);
    const auto *pfg_270 = buffer.data(pfg + 270);
    const auto *pfg_271 = buffer.data(pfg + 271);
    const auto *pfg_273 = buffer.data(pfg + 273);
    const auto *pfg_275 = buffer.data(pfg + 275);
    const auto *pfg_276 = buffer.data(pfg + 276);
    const auto *pfg_278 = buffer.data(pfg + 278);
    const auto *pfg_279 = buffer.data(pfg + 279);
    const auto *pfg_280 = buffer.data(pfg + 280);
    const auto *pfg_281 = buffer.data(pfg + 281);
    const auto *pfg_282 = buffer.data(pfg + 282);
    const auto *pfg_283 = buffer.data(pfg + 283);
    const auto *pfg_284 = buffer.data(pfg + 284);
    const auto *pfg_285 = buffer.data(pfg + 285);
    const auto *pfg_286 = buffer.data(pfg + 286);
    const auto *pfg_288 = buffer.data(pfg + 288);
    const auto *pfg_291 = buffer.data(pfg + 291);
    const auto *pfg_293 = buffer.data(pfg + 293);
    const auto *pfg_295 = buffer.data(pfg + 295);
    const auto *pfg_296 = buffer.data(pfg + 296);
    const auto *pfg_297 = buffer.data(pfg + 297);
    const auto *pfg_298 = buffer.data(pfg + 298);
    const auto *pfg_299 = buffer.data(pfg + 299);
    const auto *pfg_300 = buffer.data(pfg + 300);
    const auto *pfg_302 = buffer.data(pfg + 302);
    const auto *pfg_305 = buffer.data(pfg + 305);
    const auto *pfg_310 = buffer.data(pfg + 310);
    const auto *pfg_311 = buffer.data(pfg + 311);
    const auto *pfg_312 = buffer.data(pfg + 312);
    const auto *pfg_313 = buffer.data(pfg + 313);
    const auto *pfg_314 = buffer.data(pfg + 314);
    const auto *pfg_315 = buffer.data(pfg + 315);
    const auto *pfg_317 = buffer.data(pfg + 317);
    const auto *pfg_320 = buffer.data(pfg + 320);
    const auto *pfg_325 = buffer.data(pfg + 325);
    const auto *pfg_326 = buffer.data(pfg + 326);
    const auto *pfg_327 = buffer.data(pfg + 327);
    const auto *pfg_328 = buffer.data(pfg + 328);
    const auto *pfg_329 = buffer.data(pfg + 329);
    const auto *pfg_330 = buffer.data(pfg + 330);
    const auto *pfg_332 = buffer.data(pfg + 332);
    const auto *pfg_333 = buffer.data(pfg + 333);
    const auto *pfg_335 = buffer.data(pfg + 335);
    const auto *pfg_336 = buffer.data(pfg + 336);
    const auto *pfg_337 = buffer.data(pfg + 337);
    const auto *pfg_339 = buffer.data(pfg + 339);
    const auto *pfg_340 = buffer.data(pfg + 340);
    const auto *pfg_341 = buffer.data(pfg + 341);
    const auto *pfg_342 = buffer.data(pfg + 342);
    const auto *pfg_343 = buffer.data(pfg + 343);
    const auto *pfg_344 = buffer.data(pfg + 344);
    const auto *pfg_345 = buffer.data(pfg + 345);
    const auto *pfg_347 = buffer.data(pfg + 347);

#pragma omp simd aligned(t_373, t_374, t_375, pb_z, pc_z, pdh0_206, pdh0_207, pdg_145, \
                         pdg_146, pdg_147, pdh1_206, pdh1_207, \
                         pfg_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * pdg_145[k]
                   + f_4 * pc_z[k] * pfg_265[k];

        t_374[k] = pb_z[k] * pdh0_206[k]
                   + f_10 * pdg_146[k]
                   - f_9 * pc_z[k] * pdh1_206[k];

        t_375[k] = pb_z[k] * pdh0_207[k]
                   + f_1 * pdg_147[k]
                   - f_9 * pc_z[k] * pdh1_207[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, sfg_119, pdg_149, pdg_164, \
                         pff0_179, pff0_180, pff1_179, pff1_180, pfg_269, \
                         pfg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_0 * sfg_119[k]
                   + f_10 * pdg_164[k]
                   + f_4 * pc_y[k] * pfg_269[k];

        t_377[k] = f_0 * pdg_149[k]
                   + f_2 * pff0_179[k]
                   - f_3 * pff1_179[k]
                   + f_4 * pc_z[k] * pfg_269[k];

        t_378[k] = f_2 * pff0_180[k]
                   - f_3 * pff1_180[k]
                   + f_4 * pc_x[k] * pfg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_z, pdg_150, pdg_151, pff0_181, \
                         pff0_183, pff1_181, pff1_183, pfg_270, pfg_271, \
                         pfg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_12 * pff0_181[k]
                   - f_13 * pff1_181[k]
                   + f_4 * pc_x[k] * pfg_271[k];

        t_380[k] = f_10 * pdg_150[k]
                   + f_4 * pc_z[k] * pfg_270[k];

        t_381[k] = f_7 * pff0_183[k]
                   - f_8 * pff1_183[k]
                   + f_4 * pc_x[k] * pfg_273[k];

        t_382[k] = f_10 * pdg_151[k]
                   + f_4 * pc_z[k] * pfg_271[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, pdg_153, pff0_185, pff0_186, \
                         pff1_185, pff1_186, pfg_273, pfg_275, \
                         pfg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_7 * pff0_185[k]
                   - f_8 * pff1_185[k]
                   + f_4 * pc_x[k] * pfg_275[k];

        t_384[k] = f_5 * pff0_186[k]
                   - f_6 * pff1_186[k]
                   + f_4 * pc_x[k] * pfg_276[k];

        t_385[k] = f_10 * pdg_153[k]
                   + f_4 * pc_z[k] * pfg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, pc_x, pff0_188, pff0_189, \
                         pff1_188, pff1_189, pfg_278, pfg_279, pfg_280, pfg_281, \
                         pfg_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_5 * pff0_188[k]
                   - f_6 * pff1_188[k]
                   + f_4 * pc_x[k] * pfg_278[k];

        t_387[k] = f_5 * pff0_189[k]
                   - f_6 * pff1_189[k]
                   + f_4 * pc_x[k] * pfg_279[k];

        t_388[k] = f_4 * pc_x[k] * pfg_280[k];

        t_389[k] = f_4 * pc_x[k] * pfg_281[k];

        t_390[k] = f_4 * pc_x[k] * pfg_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, sfg_130, pdg_160, \
                         pdg_175, pff0_186, pff1_186, pfg_280, pfg_283, \
                         pfg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_4 * pc_x[k] * pfg_283[k];

        t_392[k] = f_4 * pc_x[k] * pfg_284[k];

        t_393[k] = f_0 * sfg_130[k]
                   + f_0 * pdg_175[k]
                   + f_2 * pff0_186[k]
                   - f_3 * pff1_186[k]
                   + f_4 * pc_y[k] * pfg_280[k];

        t_394[k] = f_10 * pdg_160[k]
                   + f_4 * pc_z[k] * pfg_280[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_y, pc_z, sfg_134, pdg_161, pdg_162, pdg_179, \
                         pff0_186, pff0_187, pff1_186, pff1_187, pfg_281, pfg_282, \
                         pfg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_10 * pdg_161[k]
                   + f_5 * pff0_186[k]
                   - f_6 * pff1_186[k]
                   + f_4 * pc_z[k] * pfg_281[k];

        t_396[k] = f_10 * pdg_162[k]
                   + f_7 * pff0_187[k]
                   - f_8 * pff1_187[k]
                   + f_4 * pc_z[k] * pfg_282[k];

        t_397[k] = f_0 * sfg_134[k]
                   + f_0 * pdg_179[k]
                   + f_4 * pc_y[k] * pfg_284[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_y, pc_x, pc_y, pc_z, sfh0_189, sfh1_189, \
                         pdg_164, pff0_189, pff0_191, pff1_189, pff1_191, pfg_284, \
                         pfg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_10 * pdg_164[k]
                   + f_2 * pff0_189[k]
                   - f_3 * pff1_189[k]
                   + f_4 * pc_z[k] * pfg_284[k];

        t_399[k] = pa_y[k] * sfh0_189[k]
                   - f_9 * pc_y[k] * sfh1_189[k];

        t_400[k] = f_12 * pff0_191[k]
                   - f_13 * pff1_191[k]
                   + f_4 * pc_x[k] * pfg_286[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_z, pdg_165, pdg_166, pff0_193, \
                         pff1_193, pfg_285, pfg_286, pfg_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_1 * pdg_165[k]
                   + f_4 * pc_z[k] * pfg_285[k];

        t_402[k] = f_7 * pff0_193[k]
                   - f_8 * pff1_193[k]
                   + f_4 * pc_x[k] * pfg_288[k];

        t_403[k] = f_1 * pdg_166[k]
                   + f_4 * pc_z[k] * pfg_286[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_y, pc_x, pc_y, pc_z, sfh0_194, sfh1_194, \
                         pdg_168, pff0_196, pff1_196, pfg_288, \
                         pfg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * sfh0_194[k]
                   - f_9 * pc_y[k] * sfh1_194[k];

        t_405[k] = f_5 * pff0_196[k]
                   - f_6 * pff1_196[k]
                   + f_4 * pc_x[k] * pfg_291[k];

        t_406[k] = f_1 * pdg_168[k]
                   + f_4 * pc_z[k] * pfg_288[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pa_y, pc_x, pc_y, sfh0_198, \
                         sfh1_198, pff0_198, pff1_198, pfg_293, pfg_295, pfg_296, \
                         pfg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_5 * pff0_198[k]
                   - f_6 * pff1_198[k]
                   + f_4 * pc_x[k] * pfg_293[k];

        t_408[k] = pa_y[k] * sfh0_198[k]
                   - f_9 * pc_y[k] * sfh1_198[k];

        t_409[k] = f_4 * pc_x[k] * pfg_295[k];

        t_410[k] = f_4 * pc_x[k] * pfg_296[k];

        t_411[k] = f_4 * pc_x[k] * pfg_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_y, pc_x, pc_y, pc_z, sfh0_204, \
                         sfg_145, sfh1_204, pdg_175, pfg_295, pfg_298, \
                         pfg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_4 * pc_x[k] * pfg_298[k];

        t_413[k] = f_4 * pc_x[k] * pfg_299[k];

        t_414[k] = pa_y[k] * sfh0_204[k]
                   + f_11 * sfg_145[k]
                   - f_9 * pc_y[k] * sfh1_204[k];

        t_415[k] = f_1 * pdg_175[k]
                   + f_4 * pc_z[k] * pfg_295[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pa_y, pc_y, sfh0_206, sfh0_207, sfh0_209, \
                         sfg_147, sfg_148, sfg_149, sfh1_206, sfh1_207, sfh1_209, \
                         pfg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = pa_y[k] * sfh0_206[k]
                   + f_1 * sfg_147[k]
                   - f_9 * pc_y[k] * sfh1_206[k];

        t_417[k] = pa_y[k] * sfh0_207[k]
                   + f_10 * sfg_148[k]
                   - f_9 * pc_y[k] * sfh1_207[k];

        t_418[k] = f_0 * sfg_149[k]
                   + f_4 * pc_y[k] * pfg_299[k];

        t_419[k] = pa_y[k] * sfh0_209[k]
                   - f_9 * pc_y[k] * sfh1_209[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pa_z, pc_y, pc_z, sfh0_0, sfh0_2, sfh0_3, \
                         sfg_0, sfh1_0, sfh1_2, sfh1_3, pfg_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pa_z[k] * sfh0_0[k]
                   - f_9 * pc_z[k] * sfh1_0[k];

        t_421[k] = f_4 * pc_y[k] * pfg_300[k];

        t_422[k] = pa_z[k] * sfh0_2[k]
                   + f_0 * sfg_0[k]
                   - f_9 * pc_z[k] * sfh1_2[k];

        t_423[k] = pa_z[k] * sfh0_3[k]
                   - f_9 * pc_z[k] * sfh1_3[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pa_z, pc_y, pc_z, sfh0_5, sfh0_6, sfh0_7, \
                         sfg_2, sfg_3, sfh1_5, sfh1_6, sfh1_7, \
                         pfg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_4 * pc_y[k] * pfg_302[k];

        t_425[k] = pa_z[k] * sfh0_5[k]
                   + f_10 * sfg_2[k]
                   - f_9 * pc_z[k] * sfh1_5[k];

        t_426[k] = pa_z[k] * sfh0_6[k]
                   - f_9 * pc_z[k] * sfh1_6[k];

        t_427[k] = pa_z[k] * sfh0_7[k]
                   + f_0 * sfg_3[k]
                   - f_9 * pc_z[k] * sfh1_7[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pa_z, pc_x, pc_y, pc_z, sfh0_9, sfg_5, \
                         sfh1_9, pdg_190, pdg_191, pfg_305, pfg_310, \
                         pfg_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_4 * pc_y[k] * pfg_305[k];

        t_429[k] = pa_z[k] * sfh0_9[k]
                   + f_1 * sfg_5[k]
                   - f_9 * pc_z[k] * sfh1_9[k];

        t_430[k] = f_1 * pdg_190[k]
                   + f_4 * pc_x[k] * pfg_310[k];

        t_431[k] = f_1 * pdg_191[k]
                   + f_4 * pc_x[k] * pfg_311[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pa_z, pc_x, pc_z, sfh0_15, sfh1_15, \
                         pdg_192, pdg_193, pdg_194, pfg_312, pfg_313, \
                         pfg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_1 * pdg_192[k]
                   + f_4 * pc_x[k] * pfg_312[k];

        t_433[k] = f_1 * pdg_193[k]
                   + f_4 * pc_x[k] * pfg_313[k];

        t_434[k] = f_1 * pdg_194[k]
                   + f_4 * pc_x[k] * pfg_314[k];

        t_435[k] = pa_z[k] * sfh0_15[k]
                   - f_9 * pc_z[k] * sfh1_15[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pc_y, pff0_207, pff0_208, pff0_209, \
                         pff1_207, pff1_208, pff1_209, pfg_311, pfg_312, pfg_313, \
                         pfg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_12 * pff0_207[k]
                   - f_13 * pff1_207[k]
                   + f_4 * pc_y[k] * pfg_311[k];

        t_437[k] = f_7 * pff0_208[k]
                   - f_8 * pff1_208[k]
                   + f_4 * pc_y[k] * pfg_312[k];

        t_438[k] = f_5 * pff0_209[k]
                   - f_6 * pff1_209[k]
                   + f_4 * pc_y[k] * pfg_313[k];

        t_439[k] = f_4 * pc_y[k] * pfg_314[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, pa_z, pc_y, pc_z, sfh0_20, sfh0_21, sfg_14, \
                         sfh1_20, sfh1_21, pdg_180, pfg_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = pa_z[k] * sfh0_20[k]
                   + f_11 * sfg_14[k]
                   - f_9 * pc_z[k] * sfh1_20[k];

        t_441[k] = pa_z[k] * sfh0_21[k]
                   - f_9 * pc_z[k] * sfh1_21[k];

        t_442[k] = f_0 * pdg_180[k]
                   + f_4 * pc_y[k] * pfg_315[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pa_z, pb_y, pc_y, pc_z, sfh0_24, sfh1_24, \
                         pdh0_254, pdh0_257, pdg_182, pdh1_254, pdh1_257, \
                         pfg_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pb_y[k] * pdh0_254[k]
                   - f_9 * pc_y[k] * pdh1_254[k];

        t_444[k] = pa_z[k] * sfh0_24[k]
                   - f_9 * pc_z[k] * sfh1_24[k];

        t_445[k] = f_0 * pdg_182[k]
                   + f_4 * pc_y[k] * pfg_317[k];

        t_446[k] = pb_y[k] * pdh0_257[k]
                   - f_9 * pc_y[k] * pdh1_257[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pa_z, pc_y, pc_z, sfh0_27, sfh0_28, sfg_18, \
                         sfh1_27, sfh1_28, pdg_185, pfg_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = pa_z[k] * sfh0_27[k]
                   - f_9 * pc_z[k] * sfh1_27[k];

        t_448[k] = pa_z[k] * sfh0_28[k]
                   + f_0 * sfg_18[k]
                   - f_9 * pc_z[k] * sfh1_28[k];

        t_449[k] = f_0 * pdg_185[k]
                   + f_4 * pc_y[k] * pfg_320[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pb_y, pc_x, pc_y, pdh0_261, pdg_205, \
                         pdg_206, pdg_207, pdh1_261, pfg_325, pfg_326, \
                         pfg_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = pb_y[k] * pdh0_261[k]
                   - f_9 * pc_y[k] * pdh1_261[k];

        t_451[k] = f_10 * pdg_205[k]
                   + f_4 * pc_x[k] * pfg_325[k];

        t_452[k] = f_10 * pdg_206[k]
                   + f_4 * pc_x[k] * pfg_326[k];

        t_453[k] = f_10 * pdg_207[k]
                   + f_4 * pc_x[k] * pfg_327[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, pa_z, pc_x, pc_z, sfh0_36, sfh0_37, \
                         sfg_25, sfh1_36, sfh1_37, pdg_208, pdg_209, pfg_328, \
                         pfg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_10 * pdg_208[k]
                   + f_4 * pc_x[k] * pfg_328[k];

        t_455[k] = f_10 * pdg_209[k]
                   + f_4 * pc_x[k] * pfg_329[k];

        t_456[k] = pa_z[k] * sfh0_36[k]
                   - f_9 * pc_z[k] * sfh1_36[k];

        t_457[k] = pa_z[k] * sfh0_37[k]
                   + f_0 * sfg_25[k]
                   - f_9 * pc_z[k] * sfh1_37[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, pa_z, pc_y, pc_z, sfh0_38, sfh0_39, sfg_26, \
                         sfg_27, sfh1_38, sfh1_39, pdg_194, pfg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pa_z[k] * sfh0_38[k]
                   + f_10 * sfg_26[k]
                   - f_9 * pc_z[k] * sfh1_38[k];

        t_459[k] = pa_z[k] * sfh0_39[k]
                   + f_1 * sfg_27[k]
                   - f_9 * pc_z[k] * sfh1_39[k];

        t_460[k] = f_0 * pdg_194[k]
                   + f_4 * pc_y[k] * pfg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, pb_y, pc_x, pc_y, pdh0_272, pdg_210, pdh1_272, \
                         pff0_220, pff1_220, pfg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = pb_y[k] * pdh0_272[k]
                   - f_9 * pc_y[k] * pdh1_272[k];

        t_462[k] = f_10 * pdg_210[k]
                   + f_2 * pff0_220[k]
                   - f_3 * pff1_220[k]
                   + f_4 * pc_x[k] * pfg_330[k];

        t_463[k] = f_4 * pc_y[k] * pfg_330[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pc_x, pc_y, pdg_212, pdg_213, pff0_222, \
                         pff0_223, pff1_222, pff1_223, pfg_332, \
                         pfg_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_10 * pdg_212[k]
                   + f_12 * pff0_222[k]
                   - f_13 * pff1_222[k]
                   + f_4 * pc_x[k] * pfg_332[k];

        t_465[k] = f_10 * pdg_213[k]
                   + f_7 * pff0_223[k]
                   - f_8 * pff1_223[k]
                   + f_4 * pc_x[k] * pfg_333[k];

        t_466[k] = f_4 * pc_y[k] * pfg_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pc_x, pdg_215, pdg_216, pdg_217, pff0_225, \
                         pff0_226, pff0_227, pff1_225, pff1_226, pff1_227, pfg_335, pfg_336, \
                         pfg_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_10 * pdg_215[k]
                   + f_7 * pff0_225[k]
                   - f_8 * pff1_225[k]
                   + f_4 * pc_x[k] * pfg_335[k];

        t_468[k] = f_10 * pdg_216[k]
                   + f_5 * pff0_226[k]
                   - f_6 * pff1_226[k]
                   + f_4 * pc_x[k] * pfg_336[k];

        t_469[k] = f_10 * pdg_217[k]
                   + f_5 * pff0_227[k]
                   - f_6 * pff1_227[k]
                   + f_4 * pc_x[k] * pfg_337[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, pdg_219, pdg_220, pdg_221, \
                         pff0_229, pff1_229, pfg_335, pfg_339, pfg_340, \
                         pfg_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_4 * pc_y[k] * pfg_335[k];

        t_471[k] = f_10 * pdg_219[k]
                   + f_5 * pff0_229[k]
                   - f_6 * pff1_229[k]
                   + f_4 * pc_x[k] * pfg_339[k];

        t_472[k] = f_10 * pdg_220[k]
                   + f_4 * pc_x[k] * pfg_340[k];

        t_473[k] = f_10 * pdg_221[k]
                   + f_4 * pc_x[k] * pfg_341[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_x, pc_y, pdg_222, pdg_223, pdg_224, \
                         pff0_226, pff1_226, pfg_340, pfg_342, pfg_343, \
                         pfg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_10 * pdg_222[k]
                   + f_4 * pc_x[k] * pfg_342[k];

        t_475[k] = f_10 * pdg_223[k]
                   + f_4 * pc_x[k] * pfg_343[k];

        t_476[k] = f_10 * pdg_224[k]
                   + f_4 * pc_x[k] * pfg_344[k];

        t_477[k] = f_2 * pff0_226[k]
                   - f_3 * pff1_226[k]
                   + f_4 * pc_y[k] * pfg_340[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pc_y, pff0_227, pff0_228, pff0_229, \
                         pff1_227, pff1_228, pff1_229, pfg_341, pfg_342, pfg_343, \
                         pfg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_12 * pff0_227[k]
                   - f_13 * pff1_227[k]
                   + f_4 * pc_y[k] * pfg_341[k];

        t_479[k] = f_7 * pff0_228[k]
                   - f_8 * pff1_228[k]
                   + f_4 * pc_y[k] * pfg_342[k];

        t_480[k] = f_5 * pff0_229[k]
                   - f_6 * pff1_229[k]
                   + f_4 * pc_y[k] * pfg_343[k];

        t_481[k] = f_4 * pc_y[k] * pfg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pa_z, pb_x, pc_x, pc_y, pc_z, sfh0_63, sfh1_63, \
                         pph0_188, pph1_188, pdh0_314, pdg_195, pdh1_314, \
                         pfg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_14 * pph0_188[k]
                   - f_15 * pph1_188[k]
                   + pb_x[k] * pdh0_314[k]
                   - f_9 * pc_x[k] * pdh1_314[k];

        t_483[k] = pa_z[k] * sfh0_63[k]
                   - f_9 * pc_z[k] * sfh1_63[k];

        t_484[k] = f_10 * pdg_195[k]
                   + f_4 * pc_y[k] * pfg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pa_z, pc_y, pc_z, sfh0_65, sfh0_66, sfg_45, \
                         sfh1_65, sfh1_66, pdg_197, pfg_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_z[k] * sfh0_65[k]
                   + f_0 * sfg_45[k]
                   - f_9 * pc_z[k] * sfh1_65[k];

        t_486[k] = pa_z[k] * sfh0_66[k]
                   - f_9 * pc_z[k] * sfh1_66[k];

        t_487[k] = f_10 * pdg_197[k]
                   + f_4 * pc_y[k] * pfg_347[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pa_z, pc_z, sfh0_68, sfh0_69, sfh0_70, sfg_47, \
                         sfg_48, sfh1_68, sfh1_69, sfh1_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = pa_z[k] * sfh0_68[k]
                   + f_10 * sfg_47[k]
                   - f_9 * pc_z[k] * sfh1_68[k];

        t_489[k] = pa_z[k] * sfh0_69[k]
                   - f_9 * pc_z[k] * sfh1_69[k];

        t_490[k] = pa_z[k] * sfh0_70[k]
                   + f_0 * sfg_48[k]
                   - f_9 * pc_z[k] * sfh1_70[k];
    }
}

static auto
compute_prim_pfh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfh0, const size_t sfg,
                                                          const size_t sfh1, const size_t pph0,
                                                          const size_t pph1, const size_t pdh0,
                                                          const size_t pdg, const size_t pdh1,
                                                          const size_t pff0, const size_t pff1,
                                                          const size_t pfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);
    const auto f_16 = 2.0 / q;

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
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_72 = buffer.data(sfh0 + 72);
    const auto *sfh0_78 = buffer.data(sfh0 + 78);
    const auto *sfh0_126 = buffer.data(sfh0 + 126);
    const auto *sfh0_129 = buffer.data(sfh0 + 129);
    const auto *sfh0_132 = buffer.data(sfh0 + 132);
    const auto *sfh0_141 = buffer.data(sfh0 + 141);
    const auto *sfh0_142 = buffer.data(sfh0 + 142);
    const auto *sfh0_143 = buffer.data(sfh0 + 143);
    const auto *sfh0_144 = buffer.data(sfh0 + 144);
    const auto *sfh0_146 = buffer.data(sfh0 + 146);

    const auto *sfg_50 = buffer.data(sfg + 50);
    const auto *sfg_100 = buffer.data(sfg + 100);
    const auto *sfg_101 = buffer.data(sfg + 101);
    const auto *sfg_102 = buffer.data(sfg + 102);
    const auto *sfg_104 = buffer.data(sfg + 104);

    const auto *sfh1_72 = buffer.data(sfh1 + 72);
    const auto *sfh1_78 = buffer.data(sfh1 + 78);
    const auto *sfh1_126 = buffer.data(sfh1 + 126);
    const auto *sfh1_129 = buffer.data(sfh1 + 129);
    const auto *sfh1_132 = buffer.data(sfh1 + 132);
    const auto *sfh1_141 = buffer.data(sfh1 + 141);
    const auto *sfh1_142 = buffer.data(sfh1 + 142);
    const auto *sfh1_143 = buffer.data(sfh1 + 143);
    const auto *sfh1_144 = buffer.data(sfh1 + 144);
    const auto *sfh1_146 = buffer.data(sfh1 + 146);

    const auto *pph0_188 = buffer.data(pph0 + 188);

    const auto *pph1_188 = buffer.data(pph1 + 188);

    const auto *pdh0_294 = buffer.data(pdh0 + 294);
    const auto *pdh0_296 = buffer.data(pdh0 + 296);
    const auto *pdh0_299 = buffer.data(pdh0 + 299);
    const auto *pdh0_303 = buffer.data(pdh0 + 303);
    const auto *pdh0_331 = buffer.data(pdh0 + 331);
    const auto *pdh0_332 = buffer.data(pdh0 + 332);
    const auto *pdh0_333 = buffer.data(pdh0 + 333);
    const auto *pdh0_335 = buffer.data(pdh0 + 335);
    const auto *pdh0_343 = buffer.data(pdh0 + 343);
    const auto *pdh0_351 = buffer.data(pdh0 + 351);
    const auto *pdh0_352 = buffer.data(pdh0 + 352);
    const auto *pdh0_353 = buffer.data(pdh0 + 353);
    const auto *pdh0_354 = buffer.data(pdh0 + 354);
    const auto *pdh0_356 = buffer.data(pdh0 + 356);
    const auto *pdh0_357 = buffer.data(pdh0 + 357);
    const auto *pdh0_359 = buffer.data(pdh0 + 359);
    const auto *pdh0_360 = buffer.data(pdh0 + 360);
    const auto *pdh0_362 = buffer.data(pdh0 + 362);
    const auto *pdh0_363 = buffer.data(pdh0 + 363);
    const auto *pdh0_364 = buffer.data(pdh0 + 364);
    const auto *pdh0_366 = buffer.data(pdh0 + 366);
    const auto *pdh0_372 = buffer.data(pdh0 + 372);
    const auto *pdh0_373 = buffer.data(pdh0 + 373);
    const auto *pdh0_374 = buffer.data(pdh0 + 374);
    const auto *pdh0_375 = buffer.data(pdh0 + 375);
    const auto *pdh0_377 = buffer.data(pdh0 + 377);

    const auto *pdg_200 = buffer.data(pdg + 200);
    const auto *pdg_209 = buffer.data(pdg + 209);
    const auto *pdg_210 = buffer.data(pdg + 210);
    const auto *pdg_212 = buffer.data(pdg + 212);
    const auto *pdg_215 = buffer.data(pdg + 215);
    const auto *pdg_224 = buffer.data(pdg + 224);
    const auto *pdg_225 = buffer.data(pdg + 225);
    const auto *pdg_227 = buffer.data(pdg + 227);
    const auto *pdg_230 = buffer.data(pdg + 230);
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

    const auto *pdh1_294 = buffer.data(pdh1 + 294);
    const auto *pdh1_296 = buffer.data(pdh1 + 296);
    const auto *pdh1_299 = buffer.data(pdh1 + 299);
    const auto *pdh1_303 = buffer.data(pdh1 + 303);
    const auto *pdh1_331 = buffer.data(pdh1 + 331);
    const auto *pdh1_332 = buffer.data(pdh1 + 332);
    const auto *pdh1_333 = buffer.data(pdh1 + 333);
    const auto *pdh1_335 = buffer.data(pdh1 + 335);
    const auto *pdh1_343 = buffer.data(pdh1 + 343);
    const auto *pdh1_351 = buffer.data(pdh1 + 351);
    const auto *pdh1_352 = buffer.data(pdh1 + 352);
    const auto *pdh1_353 = buffer.data(pdh1 + 353);
    const auto *pdh1_354 = buffer.data(pdh1 + 354);
    const auto *pdh1_356 = buffer.data(pdh1 + 356);
    const auto *pdh1_357 = buffer.data(pdh1 + 357);
    const auto *pdh1_359 = buffer.data(pdh1 + 359);
    const auto *pdh1_360 = buffer.data(pdh1 + 360);
    const auto *pdh1_362 = buffer.data(pdh1 + 362);
    const auto *pdh1_363 = buffer.data(pdh1 + 363);
    const auto *pdh1_364 = buffer.data(pdh1 + 364);
    const auto *pdh1_366 = buffer.data(pdh1 + 366);
    const auto *pdh1_372 = buffer.data(pdh1 + 372);
    const auto *pdh1_373 = buffer.data(pdh1 + 373);
    const auto *pdh1_374 = buffer.data(pdh1 + 374);
    const auto *pdh1_375 = buffer.data(pdh1 + 375);
    const auto *pdh1_377 = buffer.data(pdh1 + 377);

    const auto *pff0_243 = buffer.data(pff0 + 243);
    const auto *pff0_246 = buffer.data(pff0 + 246);
    const auto *pff0_262 = buffer.data(pff0 + 262);
    const auto *pff0_265 = buffer.data(pff0 + 265);
    const auto *pff0_267 = buffer.data(pff0 + 267);
    const auto *pff0_269 = buffer.data(pff0 + 269);
    const auto *pff0_270 = buffer.data(pff0 + 270);
    const auto *pff0_272 = buffer.data(pff0 + 272);
    const auto *pff0_273 = buffer.data(pff0 + 273);
    const auto *pff0_275 = buffer.data(pff0 + 275);
    const auto *pff0_276 = buffer.data(pff0 + 276);
    const auto *pff0_277 = buffer.data(pff0 + 277);
    const auto *pff0_278 = buffer.data(pff0 + 278);
    const auto *pff0_279 = buffer.data(pff0 + 279);
    const auto *pff0_283 = buffer.data(pff0 + 283);
    const auto *pff0_286 = buffer.data(pff0 + 286);
    const auto *pff0_287 = buffer.data(pff0 + 287);
    const auto *pff0_290 = buffer.data(pff0 + 290);
    const auto *pff0_292 = buffer.data(pff0 + 292);
    const auto *pff0_293 = buffer.data(pff0 + 293);

    const auto *pff1_243 = buffer.data(pff1 + 243);
    const auto *pff1_246 = buffer.data(pff1 + 246);
    const auto *pff1_262 = buffer.data(pff1 + 262);
    const auto *pff1_265 = buffer.data(pff1 + 265);
    const auto *pff1_267 = buffer.data(pff1 + 267);
    const auto *pff1_269 = buffer.data(pff1 + 269);
    const auto *pff1_270 = buffer.data(pff1 + 270);
    const auto *pff1_272 = buffer.data(pff1 + 272);
    const auto *pff1_273 = buffer.data(pff1 + 273);
    const auto *pff1_275 = buffer.data(pff1 + 275);
    const auto *pff1_276 = buffer.data(pff1 + 276);
    const auto *pff1_277 = buffer.data(pff1 + 277);
    const auto *pff1_278 = buffer.data(pff1 + 278);
    const auto *pff1_279 = buffer.data(pff1 + 279);
    const auto *pff1_283 = buffer.data(pff1 + 283);
    const auto *pff1_286 = buffer.data(pff1 + 286);
    const auto *pff1_287 = buffer.data(pff1 + 287);
    const auto *pff1_290 = buffer.data(pff1 + 290);
    const auto *pff1_292 = buffer.data(pff1 + 292);
    const auto *pff1_293 = buffer.data(pff1 + 293);

    const auto *pfg_350 = buffer.data(pfg + 350);
    const auto *pfg_355 = buffer.data(pfg + 355);
    const auto *pfg_356 = buffer.data(pfg + 356);
    const auto *pfg_357 = buffer.data(pfg + 357);
    const auto *pfg_358 = buffer.data(pfg + 358);
    const auto *pfg_359 = buffer.data(pfg + 359);
    const auto *pfg_360 = buffer.data(pfg + 360);
    const auto *pfg_362 = buffer.data(pfg + 362);
    const auto *pfg_363 = buffer.data(pfg + 363);
    const auto *pfg_365 = buffer.data(pfg + 365);
    const auto *pfg_366 = buffer.data(pfg + 366);
    const auto *pfg_370 = buffer.data(pfg + 370);
    const auto *pfg_371 = buffer.data(pfg + 371);
    const auto *pfg_372 = buffer.data(pfg + 372);
    const auto *pfg_373 = buffer.data(pfg + 373);
    const auto *pfg_374 = buffer.data(pfg + 374);
    const auto *pfg_375 = buffer.data(pfg + 375);
    const auto *pfg_377 = buffer.data(pfg + 377);
    const auto *pfg_380 = buffer.data(pfg + 380);
    const auto *pfg_385 = buffer.data(pfg + 385);
    const auto *pfg_386 = buffer.data(pfg + 386);
    const auto *pfg_387 = buffer.data(pfg + 387);
    const auto *pfg_388 = buffer.data(pfg + 388);
    const auto *pfg_389 = buffer.data(pfg + 389);
    const auto *pfg_390 = buffer.data(pfg + 390);
    const auto *pfg_392 = buffer.data(pfg + 392);
    const auto *pfg_395 = buffer.data(pfg + 395);
    const auto *pfg_397 = buffer.data(pfg + 397);
    const auto *pfg_399 = buffer.data(pfg + 399);
    const auto *pfg_400 = buffer.data(pfg + 400);
    const auto *pfg_401 = buffer.data(pfg + 401);
    const auto *pfg_402 = buffer.data(pfg + 402);
    const auto *pfg_403 = buffer.data(pfg + 403);
    const auto *pfg_404 = buffer.data(pfg + 404);
    const auto *pfg_405 = buffer.data(pfg + 405);
    const auto *pfg_407 = buffer.data(pfg + 407);
    const auto *pfg_408 = buffer.data(pfg + 408);
    const auto *pfg_410 = buffer.data(pfg + 410);
    const auto *pfg_411 = buffer.data(pfg + 411);
    const auto *pfg_412 = buffer.data(pfg + 412);
    const auto *pfg_414 = buffer.data(pfg + 414);
    const auto *pfg_415 = buffer.data(pfg + 415);
    const auto *pfg_416 = buffer.data(pfg + 416);
    const auto *pfg_417 = buffer.data(pfg + 417);
    const auto *pfg_418 = buffer.data(pfg + 418);
    const auto *pfg_419 = buffer.data(pfg + 419);
    const auto *pfg_420 = buffer.data(pfg + 420);
    const auto *pfg_422 = buffer.data(pfg + 422);
    const auto *pfg_423 = buffer.data(pfg + 423);
    const auto *pfg_425 = buffer.data(pfg + 425);
    const auto *pfg_426 = buffer.data(pfg + 426);
    const auto *pfg_427 = buffer.data(pfg + 427);
    const auto *pfg_430 = buffer.data(pfg + 430);
    const auto *pfg_431 = buffer.data(pfg + 431);
    const auto *pfg_432 = buffer.data(pfg + 432);
    const auto *pfg_433 = buffer.data(pfg + 433);
    const auto *pfg_434 = buffer.data(pfg + 434);
    const auto *pfg_435 = buffer.data(pfg + 435);
    const auto *pfg_437 = buffer.data(pfg + 437);
    const auto *pfg_438 = buffer.data(pfg + 438);

#pragma omp simd aligned(t_491, t_492, t_493, pa_z, pc_x, pc_y, pc_z, sfh0_72, sfg_50, \
                         sfh1_72, pdg_200, pdg_235, pfg_350, pfg_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * pdg_200[k]
                   + f_4 * pc_y[k] * pfg_350[k];

        t_492[k] = pa_z[k] * sfh0_72[k]
                   + f_1 * sfg_50[k]
                   - f_9 * pc_z[k] * sfh1_72[k];

        t_493[k] = f_0 * pdg_235[k]
                   + f_4 * pc_x[k] * pfg_355[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pc_x, pdg_236, pdg_237, pdg_238, pdg_239, \
                         pfg_356, pfg_357, pfg_358, pfg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_0 * pdg_236[k]
                   + f_4 * pc_x[k] * pfg_356[k];

        t_495[k] = f_0 * pdg_237[k]
                   + f_4 * pc_x[k] * pfg_357[k];

        t_496[k] = f_0 * pdg_238[k]
                   + f_4 * pc_x[k] * pfg_358[k];

        t_497[k] = f_0 * pdg_239[k]
                   + f_4 * pc_x[k] * pfg_359[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pa_z, pb_x, pc_x, pc_z, sfh0_78, sfh1_78, \
                         pdh0_331, pdh0_332, pdh0_333, pdh1_331, pdh1_332, \
                         pdh1_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = pa_z[k] * sfh0_78[k]
                   - f_9 * pc_z[k] * sfh1_78[k];

        t_499[k] = pb_x[k] * pdh0_331[k]
                   - f_9 * pc_x[k] * pdh1_331[k];

        t_500[k] = pb_x[k] * pdh0_332[k]
                   - f_9 * pc_x[k] * pdh1_332[k];

        t_501[k] = pb_x[k] * pdh0_333[k]
                   - f_9 * pc_x[k] * pdh1_333[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, pb_x, pb_y, pc_x, pc_y, pdh0_294, \
                         pdh0_335, pdg_209, pdg_210, pdh1_294, pdh1_335, pfg_359, \
                         pfg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_10 * pdg_209[k]
                   + f_4 * pc_y[k] * pfg_359[k];

        t_503[k] = pb_x[k] * pdh0_335[k]
                   - f_9 * pc_x[k] * pdh1_335[k];

        t_504[k] = pb_y[k] * pdh0_294[k]
                   - f_9 * pc_y[k] * pdh1_294[k];

        t_505[k] = f_0 * pdg_210[k]
                   + f_4 * pc_y[k] * pfg_360[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_y, pc_x, pc_y, pdh0_296, pdg_212, pdg_243, \
                         pdh1_296, pff0_243, pff1_243, pfg_362, \
                         pfg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = pb_y[k] * pdh0_296[k]
                   - f_9 * pc_y[k] * pdh1_296[k];

        t_507[k] = f_0 * pdg_243[k]
                   + f_7 * pff0_243[k]
                   - f_8 * pff1_243[k]
                   + f_4 * pc_x[k] * pfg_363[k];

        t_508[k] = f_0 * pdg_212[k]
                   + f_4 * pc_y[k] * pfg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_x, pb_y, pc_x, pc_y, pdh0_299, pdh0_343, \
                         pdg_246, pdg_247, pdh1_299, pdh1_343, pff0_246, pff1_246, \
                         pfg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pb_y[k] * pdh0_299[k]
                   - f_9 * pc_y[k] * pdh1_299[k];

        t_510[k] = f_0 * pdg_246[k]
                   + f_5 * pff0_246[k]
                   - f_6 * pff1_246[k]
                   + f_4 * pc_x[k] * pfg_366[k];

        t_511[k] = pb_x[k] * pdh0_343[k]
                   + f_10 * pdg_247[k]
                   - f_9 * pc_x[k] * pdh1_343[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pb_y, pc_x, pc_y, pdh0_303, pdg_215, \
                         pdg_250, pdg_251, pdh1_303, pfg_365, pfg_370, \
                         pfg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_0 * pdg_215[k]
                   + f_4 * pc_y[k] * pfg_365[k];

        t_513[k] = pb_y[k] * pdh0_303[k]
                   - f_9 * pc_y[k] * pdh1_303[k];

        t_514[k] = f_0 * pdg_250[k]
                   + f_4 * pc_x[k] * pfg_370[k];

        t_515[k] = f_0 * pdg_251[k]
                   + f_4 * pc_x[k] * pfg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pb_x, pc_x, pdh0_351, pdg_252, pdg_253, \
                         pdg_254, pdh1_351, pfg_372, pfg_373, pfg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_0 * pdg_252[k]
                   + f_4 * pc_x[k] * pfg_372[k];

        t_517[k] = f_0 * pdg_253[k]
                   + f_4 * pc_x[k] * pfg_373[k];

        t_518[k] = f_0 * pdg_254[k]
                   + f_4 * pc_x[k] * pfg_374[k];

        t_519[k] = pb_x[k] * pdh0_351[k]
                   - f_9 * pc_x[k] * pdh1_351[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pb_x, pc_x, pc_y, pdh0_352, pdh0_353, \
                         pdh0_354, pdg_224, pdh1_352, pdh1_353, pdh1_354, \
                         pfg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = pb_x[k] * pdh0_352[k]
                   - f_9 * pc_x[k] * pdh1_352[k];

        t_521[k] = pb_x[k] * pdh0_353[k]
                   - f_9 * pc_x[k] * pdh1_353[k];

        t_522[k] = pb_x[k] * pdh0_354[k]
                   - f_9 * pc_x[k] * pdh1_354[k];

        t_523[k] = f_0 * pdg_224[k]
                   + f_4 * pc_y[k] * pfg_374[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, pb_x, pc_x, pc_y, pdh0_356, pdh0_357, \
                         pdh0_359, pdg_255, pdg_257, pdh1_356, pdh1_357, pdh1_359, \
                         pfg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = pb_x[k] * pdh0_356[k]
                   - f_9 * pc_x[k] * pdh1_356[k];

        t_525[k] = pb_x[k] * pdh0_357[k]
                   + f_11 * pdg_255[k]
                   - f_9 * pc_x[k] * pdh1_357[k];

        t_526[k] = f_4 * pc_y[k] * pfg_375[k];

        t_527[k] = pb_x[k] * pdh0_359[k]
                   + f_16 * pdg_257[k]
                   - f_9 * pc_x[k] * pdh1_359[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pb_x, pc_x, pc_y, pdh0_360, pdh0_362, pdg_258, \
                         pdg_260, pdh1_360, pdh1_362, pfg_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = pb_x[k] * pdh0_360[k]
                   + f_1 * pdg_258[k]
                   - f_9 * pc_x[k] * pdh1_360[k];

        t_529[k] = f_4 * pc_y[k] * pfg_377[k];

        t_530[k] = pb_x[k] * pdh0_362[k]
                   + f_1 * pdg_260[k]
                   - f_9 * pc_x[k] * pdh1_362[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pb_x, pc_x, pc_y, pdh0_363, pdh0_364, pdg_261, \
                         pdg_262, pdh1_363, pdh1_364, pfg_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = pb_x[k] * pdh0_363[k]
                   + f_10 * pdg_261[k]
                   - f_9 * pc_x[k] * pdh1_363[k];

        t_532[k] = pb_x[k] * pdh0_364[k]
                   + f_10 * pdg_262[k]
                   - f_9 * pc_x[k] * pdh1_364[k];

        t_533[k] = f_4 * pc_y[k] * pfg_380[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pb_x, pc_x, pdh0_366, pdg_264, pdg_265, \
                         pdg_266, pdg_267, pdh1_366, pfg_385, pfg_386, \
                         pfg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = pb_x[k] * pdh0_366[k]
                   + f_10 * pdg_264[k]
                   - f_9 * pc_x[k] * pdh1_366[k];

        t_535[k] = f_0 * pdg_265[k]
                   + f_4 * pc_x[k] * pfg_385[k];

        t_536[k] = f_0 * pdg_266[k]
                   + f_4 * pc_x[k] * pfg_386[k];

        t_537[k] = f_0 * pdg_267[k]
                   + f_4 * pc_x[k] * pfg_387[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pb_x, pc_x, pdh0_372, pdh0_373, pdg_268, \
                         pdg_269, pdh1_372, pdh1_373, pfg_388, \
                         pfg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_0 * pdg_268[k]
                   + f_4 * pc_x[k] * pfg_388[k];

        t_539[k] = f_0 * pdg_269[k]
                   + f_4 * pc_x[k] * pfg_389[k];

        t_540[k] = pb_x[k] * pdh0_372[k]
                   - f_9 * pc_x[k] * pdh1_372[k];

        t_541[k] = pb_x[k] * pdh0_373[k]
                   - f_9 * pc_x[k] * pdh1_373[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pb_x, pc_x, pc_y, pdh0_374, pdh0_375, \
                         pdh0_377, pdh1_374, pdh1_375, pdh1_377, \
                         pfg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = pb_x[k] * pdh0_374[k]
                   - f_9 * pc_x[k] * pdh1_374[k];

        t_543[k] = pb_x[k] * pdh0_375[k]
                   - f_9 * pc_x[k] * pdh1_375[k];

        t_544[k] = f_4 * pc_y[k] * pfg_389[k];

        t_545[k] = pb_x[k] * pdh0_377[k]
                   - f_9 * pc_x[k] * pdh1_377[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pa_z, pc_x, pc_y, pc_z, sfh0_126, sfh1_126, \
                         pdg_225, pff0_262, pff1_262, pfg_390, \
                         pfg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_z[k] * sfh0_126[k]
                   - f_9 * pc_z[k] * sfh1_126[k];

        t_547[k] = f_1 * pdg_225[k]
                   + f_4 * pc_y[k] * pfg_390[k];

        t_548[k] = f_12 * pff0_262[k]
                   - f_13 * pff1_262[k]
                   + f_4 * pc_x[k] * pfg_392[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pa_z, pc_x, pc_y, pc_z, sfh0_129, sfh1_129, \
                         pdg_227, pff0_265, pff1_265, pfg_392, \
                         pfg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = pa_z[k] * sfh0_129[k]
                   - f_9 * pc_z[k] * sfh1_129[k];

        t_550[k] = f_1 * pdg_227[k]
                   + f_4 * pc_y[k] * pfg_392[k];

        t_551[k] = f_7 * pff0_265[k]
                   - f_8 * pff1_265[k]
                   + f_4 * pc_x[k] * pfg_395[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_z, pc_x, pc_y, pc_z, sfh0_132, sfh1_132, \
                         pdg_230, pff0_267, pff1_267, pfg_395, \
                         pfg_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_z[k] * sfh0_132[k]
                   - f_9 * pc_z[k] * sfh1_132[k];

        t_553[k] = f_5 * pff0_267[k]
                   - f_6 * pff1_267[k]
                   + f_4 * pc_x[k] * pfg_397[k];

        t_554[k] = f_1 * pdg_230[k]
                   + f_4 * pc_y[k] * pfg_395[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, t_560, pc_x, pff0_269, pff1_269, \
                         pfg_399, pfg_400, pfg_401, pfg_402, pfg_403, \
                         pfg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_5 * pff0_269[k]
                   - f_6 * pff1_269[k]
                   + f_4 * pc_x[k] * pfg_399[k];

        t_556[k] = f_4 * pc_x[k] * pfg_400[k];

        t_557[k] = f_4 * pc_x[k] * pfg_401[k];

        t_558[k] = f_4 * pc_x[k] * pfg_402[k];

        t_559[k] = f_4 * pc_x[k] * pfg_403[k];

        t_560[k] = f_4 * pc_x[k] * pfg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pa_z, pc_z, sfh0_141, sfh0_142, sfh0_143, \
                         sfg_100, sfg_101, sfh1_141, sfh1_142, \
                         sfh1_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = pa_z[k] * sfh0_141[k]
                   - f_9 * pc_z[k] * sfh1_141[k];

        t_562[k] = pa_z[k] * sfh0_142[k]
                   + f_0 * sfg_100[k]
                   - f_9 * pc_z[k] * sfh1_142[k];

        t_563[k] = pa_z[k] * sfh0_143[k]
                   + f_10 * sfg_101[k]
                   - f_9 * pc_z[k] * sfh1_143[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_z, pc_y, pc_z, sfh0_144, sfh0_146, sfg_102, \
                         sfg_104, sfh1_144, sfh1_146, pdg_239, \
                         pfg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = pa_z[k] * sfh0_144[k]
                   + f_1 * sfg_102[k]
                   - f_9 * pc_z[k] * sfh1_144[k];

        t_565[k] = f_1 * pdg_239[k]
                   + f_4 * pc_y[k] * pfg_404[k];

        t_566[k] = pa_z[k] * sfh0_146[k]
                   + f_11 * sfg_104[k]
                   - f_9 * pc_z[k] * sfh1_146[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pc_x, pc_y, pdg_240, pff0_270, pff0_272, \
                         pff0_273, pff1_270, pff1_272, pff1_273, pfg_405, pfg_407, \
                         pfg_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_2 * pff0_270[k]
                   - f_3 * pff1_270[k]
                   + f_4 * pc_x[k] * pfg_405[k];

        t_568[k] = f_10 * pdg_240[k]
                   + f_4 * pc_y[k] * pfg_405[k];

        t_569[k] = f_12 * pff0_272[k]
                   - f_13 * pff1_272[k]
                   + f_4 * pc_x[k] * pfg_407[k];

        t_570[k] = f_7 * pff0_273[k]
                   - f_8 * pff1_273[k]
                   + f_4 * pc_x[k] * pfg_408[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_x, pc_y, pdg_242, pff0_275, pff0_276, \
                         pff1_275, pff1_276, pfg_407, pfg_410, \
                         pfg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_10 * pdg_242[k]
                   + f_4 * pc_y[k] * pfg_407[k];

        t_572[k] = f_7 * pff0_275[k]
                   - f_8 * pff1_275[k]
                   + f_4 * pc_x[k] * pfg_410[k];

        t_573[k] = f_5 * pff0_276[k]
                   - f_6 * pff1_276[k]
                   + f_4 * pc_x[k] * pfg_411[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pc_x, pc_y, pdg_245, pff0_277, pff0_279, \
                         pff1_277, pff1_279, pfg_410, pfg_412, pfg_414, \
                         pfg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_5 * pff0_277[k]
                   - f_6 * pff1_277[k]
                   + f_4 * pc_x[k] * pfg_412[k];

        t_575[k] = f_10 * pdg_245[k]
                   + f_4 * pc_y[k] * pfg_410[k];

        t_576[k] = f_5 * pff0_279[k]
                   - f_6 * pff1_279[k]
                   + f_4 * pc_x[k] * pfg_414[k];

        t_577[k] = f_4 * pc_x[k] * pfg_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, t_582, pc_x, pc_y, pdg_250, pff0_276, \
                         pff1_276, pfg_415, pfg_416, pfg_417, pfg_418, \
                         pfg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_4 * pc_x[k] * pfg_416[k];

        t_579[k] = f_4 * pc_x[k] * pfg_417[k];

        t_580[k] = f_4 * pc_x[k] * pfg_418[k];

        t_581[k] = f_4 * pc_x[k] * pfg_419[k];

        t_582[k] = f_10 * pdg_250[k]
                   + f_2 * pff0_276[k]
                   - f_3 * pff1_276[k]
                   + f_4 * pc_y[k] * pfg_415[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pc_y, pdg_251, pdg_252, pdg_253, pff0_277, \
                         pff0_278, pff0_279, pff1_277, pff1_278, pff1_279, pfg_416, pfg_417, \
                         pfg_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_10 * pdg_251[k]
                   + f_12 * pff0_277[k]
                   - f_13 * pff1_277[k]
                   + f_4 * pc_y[k] * pfg_416[k];

        t_584[k] = f_10 * pdg_252[k]
                   + f_7 * pff0_278[k]
                   - f_8 * pff1_278[k]
                   + f_4 * pc_y[k] * pfg_417[k];

        t_585[k] = f_10 * pdg_253[k]
                   + f_5 * pff0_279[k]
                   - f_6 * pff1_279[k]
                   + f_4 * pc_y[k] * pfg_418[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pb_y, pc_y, pph0_188, pph1_188, pdh0_356, \
                         pdh0_357, pdg_254, pdg_255, pdh1_356, pdh1_357, pfg_419, \
                         pfg_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_10 * pdg_254[k]
                   + f_4 * pc_y[k] * pfg_419[k];

        t_587[k] = f_14 * pph0_188[k]
                   - f_15 * pph1_188[k]
                   + pb_y[k] * pdh0_356[k]
                   - f_9 * pc_y[k] * pdh1_356[k];

        t_588[k] = pb_y[k] * pdh0_357[k]
                   - f_9 * pc_y[k] * pdh1_357[k];

        t_589[k] = f_0 * pdg_255[k]
                   + f_4 * pc_y[k] * pfg_420[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pb_y, pc_x, pc_y, pdh0_359, pdh0_362, \
                         pdg_257, pdh1_359, pdh1_362, pff0_283, pff1_283, pfg_422, \
                         pfg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pb_y[k] * pdh0_359[k]
                   - f_9 * pc_y[k] * pdh1_359[k];

        t_591[k] = f_7 * pff0_283[k]
                   - f_8 * pff1_283[k]
                   + f_4 * pc_x[k] * pfg_423[k];

        t_592[k] = f_0 * pdg_257[k]
                   + f_4 * pc_y[k] * pfg_422[k];

        t_593[k] = pb_y[k] * pdh0_362[k]
                   - f_9 * pc_y[k] * pdh1_362[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pc_x, pc_y, pdg_260, pff0_286, pff0_287, \
                         pff1_286, pff1_287, pfg_425, pfg_426, \
                         pfg_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_5 * pff0_286[k]
                   - f_6 * pff1_286[k]
                   + f_4 * pc_x[k] * pfg_426[k];

        t_595[k] = f_5 * pff0_287[k]
                   - f_6 * pff1_287[k]
                   + f_4 * pc_x[k] * pfg_427[k];

        t_596[k] = f_0 * pdg_260[k]
                   + f_4 * pc_y[k] * pfg_425[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, t_600, t_601, t_602, pb_y, pc_x, pc_y, pdh0_366, \
                         pdh1_366, pfg_430, pfg_431, pfg_432, pfg_433, \
                         pfg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = pb_y[k] * pdh0_366[k]
                   - f_9 * pc_y[k] * pdh1_366[k];

        t_598[k] = f_4 * pc_x[k] * pfg_430[k];

        t_599[k] = f_4 * pc_x[k] * pfg_431[k];

        t_600[k] = f_4 * pc_x[k] * pfg_432[k];

        t_601[k] = f_4 * pc_x[k] * pfg_433[k];

        t_602[k] = f_4 * pc_x[k] * pfg_434[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, pb_y, pc_y, pdh0_372, pdh0_373, pdh0_374, \
                         pdg_265, pdg_266, pdg_267, pdh1_372, pdh1_373, \
                         pdh1_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pb_y[k] * pdh0_372[k]
                   + f_11 * pdg_265[k]
                   - f_9 * pc_y[k] * pdh1_372[k];

        t_604[k] = pb_y[k] * pdh0_373[k]
                   + f_16 * pdg_266[k]
                   - f_9 * pc_y[k] * pdh1_373[k];

        t_605[k] = pb_y[k] * pdh0_374[k]
                   + f_1 * pdg_267[k]
                   - f_9 * pc_y[k] * pdh1_374[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, pb_y, pc_y, pdh0_375, pdh0_377, pdg_268, \
                         pdg_269, pdh1_375, pdh1_377, pfg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = pb_y[k] * pdh0_375[k]
                   + f_10 * pdg_268[k]
                   - f_9 * pc_y[k] * pdh1_375[k];

        t_607[k] = f_0 * pdg_269[k]
                   + f_4 * pc_y[k] * pfg_434[k];

        t_608[k] = pb_y[k] * pdh0_377[k]
                   - f_9 * pc_y[k] * pdh1_377[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, t_613, pc_x, pc_y, pff0_290, pff0_292, \
                         pff0_293, pff1_290, pff1_292, pff1_293, pfg_435, pfg_437, \
                         pfg_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_2 * pff0_290[k]
                   - f_3 * pff1_290[k]
                   + f_4 * pc_x[k] * pfg_435[k];

        t_610[k] = f_4 * pc_y[k] * pfg_435[k];

        t_611[k] = f_12 * pff0_292[k]
                   - f_13 * pff1_292[k]
                   + f_4 * pc_x[k] * pfg_437[k];

        t_612[k] = f_7 * pff0_293[k]
                   - f_8 * pff1_293[k]
                   + f_4 * pc_x[k] * pfg_438[k];

        t_613[k] = f_4 * pc_y[k] * pfg_437[k];
    }
}

static auto
compute_prim_pfh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t sfg, const size_t pdg,
                                                          const size_t pff0, const size_t pff1,
                                                          const size_t pfg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_12 = 1.5 / gamma;
    const auto f_13 = 1.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfg_149 = buffer.data(sfg + 149);

    const auto *pdg_269 = buffer.data(pdg + 269);

    const auto *pff0_295 = buffer.data(pff0 + 295);
    const auto *pff0_296 = buffer.data(pff0 + 296);
    const auto *pff0_297 = buffer.data(pff0 + 297);
    const auto *pff0_298 = buffer.data(pff0 + 298);
    const auto *pff0_299 = buffer.data(pff0 + 299);

    const auto *pff1_295 = buffer.data(pff1 + 295);
    const auto *pff1_296 = buffer.data(pff1 + 296);
    const auto *pff1_297 = buffer.data(pff1 + 297);
    const auto *pff1_298 = buffer.data(pff1 + 298);
    const auto *pff1_299 = buffer.data(pff1 + 299);

    const auto *pfg_440 = buffer.data(pfg + 440);
    const auto *pfg_441 = buffer.data(pfg + 441);
    const auto *pfg_442 = buffer.data(pfg + 442);
    const auto *pfg_444 = buffer.data(pfg + 444);
    const auto *pfg_445 = buffer.data(pfg + 445);
    const auto *pfg_446 = buffer.data(pfg + 446);
    const auto *pfg_447 = buffer.data(pfg + 447);
    const auto *pfg_448 = buffer.data(pfg + 448);
    const auto *pfg_449 = buffer.data(pfg + 449);

#pragma omp simd aligned(t_614, t_615, t_616, t_617, pc_x, pc_y, pff0_295, pff0_296, pff0_297, \
                         pff1_295, pff1_296, pff1_297, pfg_440, pfg_441, \
                         pfg_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_7 * pff0_295[k]
                   - f_8 * pff1_295[k]
                   + f_4 * pc_x[k] * pfg_440[k];

        t_615[k] = f_5 * pff0_296[k]
                   - f_6 * pff1_296[k]
                   + f_4 * pc_x[k] * pfg_441[k];

        t_616[k] = f_5 * pff0_297[k]
                   - f_6 * pff1_297[k]
                   + f_4 * pc_x[k] * pfg_442[k];

        t_617[k] = f_4 * pc_y[k] * pfg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, t_622, t_623, pc_x, pff0_299, pff1_299, \
                         pfg_444, pfg_445, pfg_446, pfg_447, pfg_448, \
                         pfg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_5 * pff0_299[k]
                   - f_6 * pff1_299[k]
                   + f_4 * pc_x[k] * pfg_444[k];

        t_619[k] = f_4 * pc_x[k] * pfg_445[k];

        t_620[k] = f_4 * pc_x[k] * pfg_446[k];

        t_621[k] = f_4 * pc_x[k] * pfg_447[k];

        t_622[k] = f_4 * pc_x[k] * pfg_448[k];

        t_623[k] = f_4 * pc_x[k] * pfg_449[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pc_y, pff0_296, pff0_297, pff0_298, pff1_296, \
                         pff1_297, pff1_298, pfg_445, pfg_446, \
                         pfg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_2 * pff0_296[k]
                   - f_3 * pff1_296[k]
                   + f_4 * pc_y[k] * pfg_445[k];

        t_625[k] = f_12 * pff0_297[k]
                   - f_13 * pff1_297[k]
                   + f_4 * pc_y[k] * pfg_446[k];

        t_626[k] = f_7 * pff0_298[k]
                   - f_8 * pff1_298[k]
                   + f_4 * pc_y[k] * pfg_447[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_y, pc_z, sfg_149, pdg_269, pff0_299, \
                         pff1_299, pfg_448, pfg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_5 * pff0_299[k]
                   - f_6 * pff1_299[k]
                   + f_4 * pc_y[k] * pfg_448[k];

        t_628[k] = f_4 * pc_y[k] * pfg_449[k];

        t_629[k] = f_0 * sfg_149[k]
                   + f_1 * pdg_269[k]
                   + f_2 * pff0_299[k]
                   - f_3 * pff1_299[k]
                   + f_4 * pc_z[k] * pfg_449[k];
    }
}

auto
compute_prim_pfh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sfh0,
                                                   const size_t sfg, const size_t sfh1,
                                                   const size_t pph0, const size_t pph1,
                                                   const size_t pdh0, const size_t pdg,
                                                   const size_t pdh1, const size_t pff0,
                                                   const size_t pff1, const size_t pfg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pfh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfg, pdh0,
                                                              pdg, pdh1, pff0, pff1, pfg, ncols,
                                                              gamma, p, q);

    compute_prim_pfh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sfh0,
                                                              sfg, sfh1, pph0, pph1, pdh0, pdg,
                                                              pdh1, pff0, pff1, pfg, ncols,
                                                              gamma, p, q);

    compute_prim_pfh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sfh0,
                                                              sfg, sfh1, pdh0, pdg, pdh1, pff0,
                                                              pff1, pfg, ncols, gamma, p, q);

    compute_prim_pfh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sfh0,
                                                              sfg, sfh1, pph0, pph1, pdh0, pdg,
                                                              pdh1, pff0, pff1, pfg, ncols,
                                                              gamma, p, q);

    compute_prim_pfh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, sfh0,
                                                              sfg, sfh1, pph0, pph1, pdh0, pdg,
                                                              pdh1, pff0, pff1, pfg, ncols,
                                                              gamma, p, q);

    compute_prim_pfh_three_center_electron_repulsion_0_piece5(buffer, target, pc, sfg, pdg,
                                                              pff0, pff1, pfg, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
