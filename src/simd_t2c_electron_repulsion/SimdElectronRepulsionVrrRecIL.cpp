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


#include "SimdElectronRepulsionVrrRecIL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_il_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_0 = buffer.data(gl0 + 0);

    const auto *gl1_0 = buffer.data(gl1 + 0);

    const auto *hk_0 = buffer.data(hk + 0);
    const auto *hk_1 = buffer.data(hk + 1);
    const auto *hk_2 = buffer.data(hk + 2);
    const auto *hk_3 = buffer.data(hk + 3);
    const auto *hk_5 = buffer.data(hk + 5);
    const auto *hk_6 = buffer.data(hk + 6);
    const auto *hk_7 = buffer.data(hk + 7);
    const auto *hk_8 = buffer.data(hk + 8);
    const auto *hk_9 = buffer.data(hk + 9);
    const auto *hk_10 = buffer.data(hk + 10);
    const auto *hk_11 = buffer.data(hk + 11);
    const auto *hk_12 = buffer.data(hk + 12);
    const auto *hk_13 = buffer.data(hk + 13);
    const auto *hk_14 = buffer.data(hk + 14);
    const auto *hk_15 = buffer.data(hk + 15);
    const auto *hk_16 = buffer.data(hk + 16);
    const auto *hk_17 = buffer.data(hk + 17);
    const auto *hk_18 = buffer.data(hk + 18);
    const auto *hk_19 = buffer.data(hk + 19);
    const auto *hk_20 = buffer.data(hk + 20);
    const auto *hk_28 = buffer.data(hk + 28);
    const auto *hk_29 = buffer.data(hk + 29);
    const auto *hk_30 = buffer.data(hk + 30);
    const auto *hk_31 = buffer.data(hk + 31);
    const auto *hk_32 = buffer.data(hk + 32);
    const auto *hk_33 = buffer.data(hk + 33);
    const auto *hk_34 = buffer.data(hk + 34);
    const auto *hk_35 = buffer.data(hk + 35);
    const auto *hk_36 = buffer.data(hk + 36);
    const auto *hk_41 = buffer.data(hk + 41);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_66 = buffer.data(hk + 66);
    const auto *hk_67 = buffer.data(hk + 67);
    const auto *hk_68 = buffer.data(hk + 68);
    const auto *hk_69 = buffer.data(hk + 69);
    const auto *hk_70 = buffer.data(hk + 70);
    const auto *hk_101 = buffer.data(hk + 101);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_114 = buffer.data(hk + 114);

    const auto *hl_0 = buffer.data(hl + 0);
    const auto *hl_3 = buffer.data(hl + 3);
    const auto *hl_5 = buffer.data(hl + 5);
    const auto *hl_6 = buffer.data(hl + 6);
    const auto *hl_9 = buffer.data(hl + 9);
    const auto *hl_10 = buffer.data(hl + 10);
    const auto *hl_12 = buffer.data(hl + 12);
    const auto *hl_14 = buffer.data(hl + 14);
    const auto *hl_15 = buffer.data(hl + 15);
    const auto *hl_17 = buffer.data(hl + 17);
    const auto *hl_18 = buffer.data(hl + 18);
    const auto *hl_20 = buffer.data(hl + 20);
    const auto *hl_21 = buffer.data(hl + 21);
    const auto *hl_23 = buffer.data(hl + 23);
    const auto *hl_24 = buffer.data(hl + 24);
    const auto *hl_25 = buffer.data(hl + 25);
    const auto *hl_27 = buffer.data(hl + 27);
    const auto *hl_28 = buffer.data(hl + 28);
    const auto *hl_35 = buffer.data(hl + 35);
    const auto *hl_36 = buffer.data(hl + 36);
    const auto *hl_38 = buffer.data(hl + 38);
    const auto *hl_39 = buffer.data(hl + 39);
    const auto *hl_40 = buffer.data(hl + 40);
    const auto *hl_41 = buffer.data(hl + 41);
    const auto *hl_42 = buffer.data(hl + 42);
    const auto *hl_44 = buffer.data(hl + 44);
    const auto *hl_45 = buffer.data(hl + 45);

    const auto *ii0_0 = buffer.data(ii0 + 0);
    const auto *ii0_1 = buffer.data(ii0 + 1);
    const auto *ii0_2 = buffer.data(ii0 + 2);
    const auto *ii0_3 = buffer.data(ii0 + 3);
    const auto *ii0_5 = buffer.data(ii0 + 5);
    const auto *ii0_6 = buffer.data(ii0 + 6);
    const auto *ii0_8 = buffer.data(ii0 + 8);
    const auto *ii0_9 = buffer.data(ii0 + 9);
    const auto *ii0_10 = buffer.data(ii0 + 10);
    const auto *ii0_12 = buffer.data(ii0 + 12);
    const auto *ii0_13 = buffer.data(ii0 + 13);
    const auto *ii0_14 = buffer.data(ii0 + 14);
    const auto *ii0_21 = buffer.data(ii0 + 21);
    const auto *ii0_23 = buffer.data(ii0 + 23);
    const auto *ii0_24 = buffer.data(ii0 + 24);
    const auto *ii0_25 = buffer.data(ii0 + 25);
    const auto *ii0_26 = buffer.data(ii0 + 26);
    const auto *ii0_27 = buffer.data(ii0 + 27);
    const auto *ii0_84 = buffer.data(ii0 + 84);
    const auto *ii0_86 = buffer.data(ii0 + 86);
    const auto *ii0_87 = buffer.data(ii0 + 87);
    const auto *ii0_90 = buffer.data(ii0 + 90);

    const auto *ii1_0 = buffer.data(ii1 + 0);
    const auto *ii1_1 = buffer.data(ii1 + 1);
    const auto *ii1_2 = buffer.data(ii1 + 2);
    const auto *ii1_3 = buffer.data(ii1 + 3);
    const auto *ii1_5 = buffer.data(ii1 + 5);
    const auto *ii1_6 = buffer.data(ii1 + 6);
    const auto *ii1_8 = buffer.data(ii1 + 8);
    const auto *ii1_9 = buffer.data(ii1 + 9);
    const auto *ii1_10 = buffer.data(ii1 + 10);
    const auto *ii1_12 = buffer.data(ii1 + 12);
    const auto *ii1_13 = buffer.data(ii1 + 13);
    const auto *ii1_14 = buffer.data(ii1 + 14);
    const auto *ii1_21 = buffer.data(ii1 + 21);
    const auto *ii1_23 = buffer.data(ii1 + 23);
    const auto *ii1_24 = buffer.data(ii1 + 24);
    const auto *ii1_25 = buffer.data(ii1 + 25);
    const auto *ii1_26 = buffer.data(ii1 + 26);
    const auto *ii1_27 = buffer.data(ii1 + 27);
    const auto *ii1_84 = buffer.data(ii1 + 84);
    const auto *ii1_86 = buffer.data(ii1 + 86);
    const auto *ii1_87 = buffer.data(ii1 + 87);
    const auto *ii1_90 = buffer.data(ii1 + 90);

    const auto *ik_0 = buffer.data(ik + 0);
    const auto *ik_1 = buffer.data(ik + 1);
    const auto *ik_2 = buffer.data(ik + 2);
    const auto *ik_3 = buffer.data(ik + 3);
    const auto *ik_5 = buffer.data(ik + 5);
    const auto *ik_6 = buffer.data(ik + 6);
    const auto *ik_8 = buffer.data(ik + 8);
    const auto *ik_9 = buffer.data(ik + 9);
    const auto *ik_10 = buffer.data(ik + 10);
    const auto *ik_12 = buffer.data(ik + 12);
    const auto *ik_13 = buffer.data(ik + 13);
    const auto *ik_14 = buffer.data(ik + 14);
    const auto *ik_15 = buffer.data(ik + 15);
    const auto *ik_17 = buffer.data(ik + 17);
    const auto *ik_18 = buffer.data(ik + 18);
    const auto *ik_19 = buffer.data(ik + 19);
    const auto *ik_20 = buffer.data(ik + 20);
    const auto *ik_21 = buffer.data(ik + 21);
    const auto *ik_27 = buffer.data(ik + 27);
    const auto *ik_28 = buffer.data(ik + 28);
    const auto *ik_30 = buffer.data(ik + 30);
    const auto *ik_31 = buffer.data(ik + 31);
    const auto *ik_32 = buffer.data(ik + 32);
    const auto *ik_33 = buffer.data(ik + 33);
    const auto *ik_34 = buffer.data(ik + 34);
    const auto *ik_35 = buffer.data(ik + 35);
    const auto *ik_36 = buffer.data(ik + 36);
    const auto *ik_37 = buffer.data(ik + 37);
    const auto *ik_39 = buffer.data(ik + 39);
    const auto *ik_41 = buffer.data(ik + 41);
    const auto *ik_42 = buffer.data(ik + 42);
    const auto *ik_45 = buffer.data(ik + 45);
    const auto *ik_46 = buffer.data(ik + 46);
    const auto *ik_50 = buffer.data(ik + 50);
    const auto *ik_51 = buffer.data(ik + 51);
    const auto *ik_56 = buffer.data(ik + 56);
    const auto *ik_57 = buffer.data(ik + 57);
    const auto *ik_64 = buffer.data(ik + 64);
    const auto *ik_66 = buffer.data(ik + 66);
    const auto *ik_67 = buffer.data(ik + 67);
    const auto *ik_68 = buffer.data(ik + 68);
    const auto *ik_69 = buffer.data(ik + 69);
    const auto *ik_70 = buffer.data(ik + 70);
    const auto *ik_71 = buffer.data(ik + 71);
    const auto *ik_72 = buffer.data(ik + 72);
    const auto *ik_74 = buffer.data(ik + 74);
    const auto *ik_75 = buffer.data(ik + 75);
    const auto *ik_77 = buffer.data(ik + 77);
    const auto *ik_78 = buffer.data(ik + 78);
    const auto *ik_81 = buffer.data(ik + 81);
    const auto *ik_82 = buffer.data(ik + 82);
    const auto *ik_86 = buffer.data(ik + 86);
    const auto *ik_87 = buffer.data(ik + 87);
    const auto *ik_92 = buffer.data(ik + 92);
    const auto *ik_99 = buffer.data(ik + 99);
    const auto *ik_100 = buffer.data(ik + 100);
    const auto *ik_101 = buffer.data(ik + 101);
    const auto *ik_102 = buffer.data(ik + 102);
    const auto *ik_103 = buffer.data(ik + 103);
    const auto *ik_104 = buffer.data(ik + 104);
    const auto *ik_105 = buffer.data(ik + 105);
    const auto *ik_107 = buffer.data(ik + 107);
    const auto *ik_108 = buffer.data(ik + 108);
    const auto *ik_109 = buffer.data(ik + 109);
    const auto *ik_110 = buffer.data(ik + 110);
    const auto *ik_111 = buffer.data(ik + 111);
    const auto *ik_113 = buffer.data(ik + 113);
    const auto *ik_114 = buffer.data(ik + 114);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hk_0, ii0_0, ii1_0, \
                         ik_0, ik_1, ik_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hk_0[k]
                 + f_1 * ii0_0[k]
                 - f_2 * ii1_0[k]
                 + pb_x[k] * ik_0[k];

        t_1[k] = pb_y[k] * ik_0[k];

        t_2[k] = pb_z[k] * ik_0[k];

        t_3[k] = f_3 * ii0_0[k]
                 - f_4 * ii1_0[k]
                 + pb_y[k] * ik_1[k];

        t_4[k] = pb_y[k] * ik_2[k];

        t_5[k] = f_3 * ii0_0[k]
                 - f_4 * ii1_0[k]
                 + pb_z[k] * ik_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, ii0_1, ii0_2, ii0_3, ii1_1, \
                         ii1_2, ii1_3, ik_3, ik_5, ik_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * ii0_1[k]
                 - f_6 * ii1_1[k]
                 + pb_y[k] * ik_3[k];

        t_7[k] = pb_z[k] * ik_3[k];

        t_8[k] = pb_y[k] * ik_5[k];

        t_9[k] = f_5 * ii0_2[k]
                 - f_6 * ii1_2[k]
                 + pb_z[k] * ik_5[k];

        t_10[k] = f_7 * ii0_3[k]
                  - f_8 * ii1_3[k]
                  + pb_y[k] * ik_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, ii0_5, ii0_6, ii1_5, \
                         ii1_6, ik_6, ik_8, ik_9, ik_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * ik_6[k];

        t_12[k] = f_3 * ii0_5[k]
                  - f_4 * ii1_5[k]
                  + pb_y[k] * ik_8[k];

        t_13[k] = pb_y[k] * ik_9[k];

        t_14[k] = f_7 * ii0_5[k]
                  - f_8 * ii1_5[k]
                  + pb_z[k] * ik_9[k];

        t_15[k] = f_9 * ii0_6[k]
                  - f_10 * ii1_6[k]
                  + pb_y[k] * ik_10[k];

        t_16[k] = pb_z[k] * ik_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, ii0_8, ii0_9, ii1_8, ii1_9, \
                         ik_12, ik_13, ik_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * ii0_8[k]
                  - f_6 * ii1_8[k]
                  + pb_y[k] * ik_12[k];

        t_18[k] = f_3 * ii0_9[k]
                  - f_4 * ii1_9[k]
                  + pb_y[k] * ik_13[k];

        t_19[k] = pb_y[k] * ik_14[k];

        t_20[k] = f_9 * ii0_9[k]
                  - f_10 * ii1_9[k]
                  + pb_z[k] * ik_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, ii0_10, ii0_12, ii0_13, ii1_10, \
                         ii1_12, ii1_13, ik_15, ik_17, ik_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * ii0_10[k]
                  - f_12 * ii1_10[k]
                  + pb_y[k] * ik_15[k];

        t_22[k] = pb_z[k] * ik_15[k];

        t_23[k] = f_7 * ii0_12[k]
                  - f_8 * ii1_12[k]
                  + pb_y[k] * ik_17[k];

        t_24[k] = f_5 * ii0_13[k]
                  - f_6 * ii1_13[k]
                  + pb_y[k] * ik_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, hk_28, ii0_14, \
                         ii1_14, ik_19, ik_20, ik_21, ik_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * ii0_14[k]
                  - f_4 * ii1_14[k]
                  + pb_y[k] * ik_19[k];

        t_26[k] = pb_y[k] * ik_20[k];

        t_27[k] = f_11 * ii0_14[k]
                  - f_12 * ii1_14[k]
                  + pb_z[k] * ik_20[k];

        t_28[k] = f_0 * hk_28[k]
                  + pb_x[k] * ik_28[k];

        t_29[k] = pb_z[k] * ik_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, hk_30, hk_31, hk_32, hk_33, \
                         ik_27, ik_30, ik_31, ik_32, ik_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * hk_30[k]
                  + pb_x[k] * ik_30[k];

        t_31[k] = f_0 * hk_31[k]
                  + pb_x[k] * ik_31[k];

        t_32[k] = f_0 * hk_32[k]
                  + pb_x[k] * ik_32[k];

        t_33[k] = f_0 * hk_33[k]
                  + pb_x[k] * ik_33[k];

        t_34[k] = pb_y[k] * ik_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, hk_35, ii0_21, ii0_23, \
                         ii1_21, ii1_23, ik_28, ik_30, ik_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * hk_35[k]
                  + pb_x[k] * ik_35[k];

        t_36[k] = f_1 * ii0_21[k]
                  - f_2 * ii1_21[k]
                  + pb_y[k] * ik_28[k];

        t_37[k] = pb_z[k] * ik_28[k];

        t_38[k] = f_11 * ii0_23[k]
                  - f_12 * ii1_23[k]
                  + pb_y[k] * ik_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, ii0_24, ii0_25, ii0_26, ii1_24, ii1_25, \
                         ii1_26, ik_31, ik_32, ik_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * ii0_24[k]
                  - f_10 * ii1_24[k]
                  + pb_y[k] * ik_31[k];

        t_40[k] = f_7 * ii0_25[k]
                  - f_8 * ii1_25[k]
                  + pb_y[k] * ik_32[k];

        t_41[k] = f_5 * ii0_26[k]
                  - f_6 * ii1_26[k]
                  + pb_y[k] * ik_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, hk_0, hl_0, \
                         ii0_27, ii1_27, ik_34, ik_35, ik_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * ii0_27[k]
                  - f_4 * ii1_27[k]
                  + pb_y[k] * ik_34[k];

        t_43[k] = pb_y[k] * ik_35[k];

        t_44[k] = f_1 * ii0_27[k]
                  - f_2 * ii1_27[k]
                  + pb_z[k] * ik_35[k];

        t_45[k] = pa_y[k] * hl_0[k];

        t_46[k] = f_13 * hk_0[k]
                  + pb_y[k] * ik_36[k];

        t_47[k] = pb_z[k] * ik_36[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_z, hk_1, hk_3, hl_3, hl_5, \
                         hl_6, ik_37, ik_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * hk_1[k]
                  + pa_y[k] * hl_3[k];

        t_49[k] = pb_z[k] * ik_37[k];

        t_50[k] = pa_y[k] * hl_5[k];

        t_51[k] = f_15 * hk_3[k]
                  + pa_y[k] * hl_6[k];

        t_52[k] = pb_z[k] * ik_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_y, pb_z, hk_5, hk_6, hk_8, \
                         hl_9, hl_10, hl_12, ik_41, ik_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * hk_5[k]
                  + pb_y[k] * ik_41[k];

        t_54[k] = pa_y[k] * hl_9[k];

        t_55[k] = f_16 * hk_6[k]
                  + pa_y[k] * hl_10[k];

        t_56[k] = pb_z[k] * ik_42[k];

        t_57[k] = f_14 * hk_8[k]
                  + pa_y[k] * hl_12[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, hk_9, hk_10, hk_12, \
                         hl_14, hl_15, hl_17, ik_45, ik_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_13 * hk_9[k]
                  + pb_y[k] * ik_45[k];

        t_59[k] = pa_y[k] * hl_14[k];

        t_60[k] = f_17 * hk_10[k]
                  + pa_y[k] * hl_15[k];

        t_61[k] = pb_z[k] * ik_46[k];

        t_62[k] = f_15 * hk_12[k]
                  + pa_y[k] * hl_17[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, hk_13, hk_14, hk_15, \
                         hl_18, hl_20, hl_21, ik_50, ik_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * hk_13[k]
                  + pa_y[k] * hl_18[k];

        t_64[k] = f_13 * hk_14[k]
                  + pb_y[k] * ik_50[k];

        t_65[k] = pa_y[k] * hl_20[k];

        t_66[k] = f_0 * hk_15[k]
                  + pa_y[k] * hl_21[k];

        t_67[k] = pb_z[k] * ik_51[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pb_y, hk_17, hk_18, hk_19, hk_20, \
                         hl_23, hl_24, hl_25, hl_27, ik_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_16 * hk_17[k]
                  + pa_y[k] * hl_23[k];

        t_69[k] = f_15 * hk_18[k]
                  + pa_y[k] * hl_24[k];

        t_70[k] = f_14 * hk_19[k]
                  + pa_y[k] * hl_25[k];

        t_71[k] = f_13 * hk_20[k]
                  + pb_y[k] * ik_56[k];

        t_72[k] = pa_y[k] * hl_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, hk_64, hk_66, hk_67, hk_68, \
                         ik_57, ik_64, ik_66, ik_67, ik_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_17 * hk_64[k]
                  + pb_x[k] * ik_64[k];

        t_74[k] = pb_z[k] * ik_57[k];

        t_75[k] = f_17 * hk_66[k]
                  + pb_x[k] * ik_66[k];

        t_76[k] = f_17 * hk_67[k]
                  + pb_x[k] * ik_67[k];

        t_77[k] = f_17 * hk_68[k]
                  + pb_x[k] * ik_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, hk_28, hk_69, hk_70, \
                         hl_35, hl_36, ik_64, ik_69, ik_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_17 * hk_69[k]
                  + pb_x[k] * ik_69[k];

        t_79[k] = f_17 * hk_70[k]
                  + pb_x[k] * ik_70[k];

        t_80[k] = pa_y[k] * hl_35[k];

        t_81[k] = f_18 * hk_28[k]
                  + pa_y[k] * hl_36[k];

        t_82[k] = pb_z[k] * ik_64[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, hk_30, hk_31, hk_32, hk_33, \
                         hk_34, hl_38, hl_39, hl_40, hl_41, hl_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * hk_30[k]
                  + pa_y[k] * hl_38[k];

        t_84[k] = f_17 * hk_31[k]
                  + pa_y[k] * hl_39[k];

        t_85[k] = f_16 * hk_32[k]
                  + pa_y[k] * hl_40[k];

        t_86[k] = f_15 * hk_33[k]
                  + pa_y[k] * hl_41[k];

        t_87[k] = f_14 * hk_34[k]
                  + pa_y[k] * hl_42[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, hk_0, hk_35, \
                         hl_0, hl_44, ik_71, ik_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_13 * hk_35[k]
                  + pb_y[k] * ik_71[k];

        t_89[k] = pa_y[k] * hl_44[k];

        t_90[k] = pa_z[k] * hl_0[k];

        t_91[k] = pb_y[k] * ik_72[k];

        t_92[k] = f_13 * hk_0[k]
                  + pb_z[k] * ik_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, hk_2, hk_3, hl_3, \
                         hl_5, hl_6, ik_74, ik_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * hl_3[k];

        t_94[k] = pb_y[k] * ik_74[k];

        t_95[k] = f_14 * hk_2[k]
                  + pa_z[k] * hl_5[k];

        t_96[k] = pa_z[k] * hl_6[k];

        t_97[k] = f_13 * hk_3[k]
                  + pb_z[k] * ik_75[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, hk_5, hk_6, hk_7, \
                         hl_9, hl_10, hl_12, ik_77, ik_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * ik_77[k];

        t_99[k] = f_15 * hk_5[k]
                  + pa_z[k] * hl_9[k];

        t_100[k] = pa_z[k] * hl_10[k];

        t_101[k] = f_13 * hk_6[k]
                   + pb_z[k] * ik_78[k];

        t_102[k] = f_14 * hk_7[k]
                   + pa_z[k] * hl_12[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, hk_9, hk_10, \
                         hk_11, hl_14, hl_15, hl_17, ik_81, ik_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * ik_81[k];

        t_104[k] = f_16 * hk_9[k]
                   + pa_z[k] * hl_14[k];

        t_105[k] = pa_z[k] * hl_15[k];

        t_106[k] = f_13 * hk_10[k]
                   + pb_z[k] * ik_82[k];

        t_107[k] = f_14 * hk_11[k]
                   + pa_z[k] * hl_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, hk_12, hk_14, \
                         hk_15, hl_18, hl_20, hl_21, ik_86, ik_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * hk_12[k]
                   + pa_z[k] * hl_18[k];

        t_109[k] = pb_y[k] * ik_86[k];

        t_110[k] = f_17 * hk_14[k]
                   + pa_z[k] * hl_20[k];

        t_111[k] = pa_z[k] * hl_21[k];

        t_112[k] = f_13 * hk_15[k]
                   + pb_z[k] * ik_87[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, hk_16, hk_17, hk_18, \
                         hk_20, hl_23, hl_24, hl_25, hl_27, ik_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * hk_16[k]
                   + pa_z[k] * hl_23[k];

        t_114[k] = f_15 * hk_17[k]
                   + pa_z[k] * hl_24[k];

        t_115[k] = f_16 * hk_18[k]
                   + pa_z[k] * hl_25[k];

        t_116[k] = pb_y[k] * ik_92[k];

        t_117[k] = f_0 * hk_20[k]
                   + pa_z[k] * hl_27[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, hk_101, hk_102, \
                         hk_103, hk_104, hl_28, ik_101, ik_102, ik_103, \
                         ik_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * hl_28[k];

        t_119[k] = f_17 * hk_101[k]
                   + pb_x[k] * ik_101[k];

        t_120[k] = f_17 * hk_102[k]
                   + pb_x[k] * ik_102[k];

        t_121[k] = f_17 * hk_103[k]
                   + pb_x[k] * ik_103[k];

        t_122[k] = f_17 * hk_104[k]
                   + pb_x[k] * ik_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, pb_y, hk_105, hk_107, hl_36, \
                         ik_99, ik_105, ik_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_17 * hk_105[k]
                   + pb_x[k] * ik_105[k];

        t_124[k] = pb_y[k] * ik_99[k];

        t_125[k] = f_17 * hk_107[k]
                   + pb_x[k] * ik_107[k];

        t_126[k] = pa_z[k] * hl_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_z, pb_z, hk_28, hk_29, hk_30, hk_31, \
                         hl_38, hl_39, hl_40, ik_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * hk_28[k]
                   + pb_z[k] * ik_100[k];

        t_128[k] = f_14 * hk_29[k]
                   + pa_z[k] * hl_38[k];

        t_129[k] = f_15 * hk_30[k]
                   + pa_z[k] * hl_39[k];

        t_130[k] = f_16 * hk_31[k]
                   + pa_z[k] * hl_40[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_z, pb_y, hk_32, hk_33, hk_35, hl_41, \
                         hl_42, hl_44, ik_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_17 * hk_32[k]
                   + pa_z[k] * hl_41[k];

        t_132[k] = f_0 * hk_33[k]
                   + pa_z[k] * hl_42[k];

        t_133[k] = pb_y[k] * ik_107[k];

        t_134[k] = f_18 * hk_35[k]
                   + pa_z[k] * hl_44[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pb_y, pb_z, gl0_0, gl1_0, hk_36, hl_45, \
                         ik_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_19 * gl0_0[k]
                   - f_20 * gl1_0[k]
                   + pa_y[k] * hl_45[k];

        t_136[k] = f_14 * hk_36[k]
                   + pb_y[k] * ik_108[k];

        t_137[k] = pb_z[k] * ik_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_x, pb_z, hk_111, ii0_84, ii0_87, ii1_84, \
                         ii1_87, ik_109, ik_110, ik_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * hk_111[k]
                   + f_11 * ii0_87[k]
                   - f_12 * ii1_87[k]
                   + pb_x[k] * ik_111[k];

        t_139[k] = pb_z[k] * ik_109[k];

        t_140[k] = f_3 * ii0_84[k]
                   - f_4 * ii1_84[k]
                   + pb_z[k] * ik_110[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_x, pb_y, pb_z, hk_41, hk_114, ii0_86, \
                         ii0_90, ii1_86, ii1_90, ik_111, ik_113, \
                         ik_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_16 * hk_114[k]
                   + f_9 * ii0_90[k]
                   - f_10 * ii1_90[k]
                   + pb_x[k] * ik_114[k];

        t_142[k] = pb_z[k] * ik_111[k];

        t_143[k] = f_14 * hk_41[k]
                   + pb_y[k] * ik_113[k];

        t_144[k] = f_5 * ii0_86[k]
                   - f_6 * ii1_86[k]
                   + pb_z[k] * ik_113[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);
    const auto f_23 = 1.0 / alpha;
    const auto f_24 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_0 = buffer.data(gl0 + 0);
    const auto *gl0_45 = buffer.data(gl0 + 45);
    const auto *gl0_171 = buffer.data(gl0 + 171);
    const auto *gl0_269 = buffer.data(gl0 + 269);

    const auto *gl1_0 = buffer.data(gl1 + 0);
    const auto *gl1_45 = buffer.data(gl1 + 45);
    const auto *gl1_171 = buffer.data(gl1 + 171);
    const auto *gl1_269 = buffer.data(gl1 + 269);

    const auto *hk_39 = buffer.data(hk + 39);
    const auto *hk_42 = buffer.data(hk + 42);
    const auto *hk_45 = buffer.data(hk + 45);
    const auto *hk_46 = buffer.data(hk + 46);
    const auto *hk_50 = buffer.data(hk + 50);
    const auto *hk_51 = buffer.data(hk + 51);
    const auto *hk_56 = buffer.data(hk + 56);
    const auto *hk_64 = buffer.data(hk + 64);
    const auto *hk_71 = buffer.data(hk + 71);
    const auto *hk_72 = buffer.data(hk + 72);
    const auto *hk_74 = buffer.data(hk + 74);
    const auto *hk_75 = buffer.data(hk + 75);
    const auto *hk_77 = buffer.data(hk + 77);
    const auto *hk_78 = buffer.data(hk + 78);
    const auto *hk_80 = buffer.data(hk + 80);
    const auto *hk_81 = buffer.data(hk + 81);
    const auto *hk_82 = buffer.data(hk + 82);
    const auto *hk_84 = buffer.data(hk + 84);
    const auto *hk_85 = buffer.data(hk + 85);
    const auto *hk_86 = buffer.data(hk + 86);
    const auto *hk_87 = buffer.data(hk + 87);
    const auto *hk_89 = buffer.data(hk + 89);
    const auto *hk_90 = buffer.data(hk + 90);
    const auto *hk_91 = buffer.data(hk + 91);
    const auto *hk_92 = buffer.data(hk + 92);
    const auto *hk_100 = buffer.data(hk + 100);
    const auto *hk_102 = buffer.data(hk + 102);
    const auto *hk_103 = buffer.data(hk + 103);
    const auto *hk_104 = buffer.data(hk + 104);
    const auto *hk_105 = buffer.data(hk + 105);
    const auto *hk_106 = buffer.data(hk + 106);
    const auto *hk_107 = buffer.data(hk + 107);
    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_129 = buffer.data(hk + 129);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_142 = buffer.data(hk + 142);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_173 = buffer.data(hk + 173);
    const auto *hk_174 = buffer.data(hk + 174);
    const auto *hk_175 = buffer.data(hk + 175);
    const auto *hk_176 = buffer.data(hk + 176);
    const auto *hk_177 = buffer.data(hk + 177);
    const auto *hk_178 = buffer.data(hk + 178);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_207 = buffer.data(hk + 207);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_209 = buffer.data(hk + 209);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_219 = buffer.data(hk + 219);

    const auto *hl_46 = buffer.data(hl + 46);
    const auto *hl_48 = buffer.data(hl + 48);
    const auto *hl_51 = buffer.data(hl + 51);
    const auto *hl_55 = buffer.data(hl + 55);
    const auto *hl_60 = buffer.data(hl + 60);
    const auto *hl_66 = buffer.data(hl + 66);
    const auto *hl_73 = buffer.data(hl + 73);
    const auto *hl_81 = buffer.data(hl + 81);
    const auto *hl_90 = buffer.data(hl + 90);
    const auto *hl_92 = buffer.data(hl + 92);
    const auto *hl_95 = buffer.data(hl + 95);
    const auto *hl_99 = buffer.data(hl + 99);
    const auto *hl_102 = buffer.data(hl + 102);
    const auto *hl_104 = buffer.data(hl + 104);
    const auto *hl_107 = buffer.data(hl + 107);
    const auto *hl_108 = buffer.data(hl + 108);
    const auto *hl_110 = buffer.data(hl + 110);
    const auto *hl_113 = buffer.data(hl + 113);
    const auto *hl_114 = buffer.data(hl + 114);
    const auto *hl_115 = buffer.data(hl + 115);
    const auto *hl_117 = buffer.data(hl + 117);
    const auto *hl_125 = buffer.data(hl + 125);
    const auto *hl_128 = buffer.data(hl + 128);
    const auto *hl_129 = buffer.data(hl + 129);
    const auto *hl_130 = buffer.data(hl + 130);
    const auto *hl_131 = buffer.data(hl + 131);
    const auto *hl_132 = buffer.data(hl + 132);
    const auto *hl_134 = buffer.data(hl + 134);
    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_269 = buffer.data(hl + 269);

    const auto *ii0_87 = buffer.data(ii0 + 87);
    const auto *ii0_89 = buffer.data(ii0 + 89);
    const auto *ii0_90 = buffer.data(ii0 + 90);
    const auto *ii0_91 = buffer.data(ii0 + 91);
    const auto *ii0_93 = buffer.data(ii0 + 93);
    const auto *ii0_94 = buffer.data(ii0 + 94);
    const auto *ii0_95 = buffer.data(ii0 + 95);
    const auto *ii0_96 = buffer.data(ii0 + 96);
    const auto *ii0_98 = buffer.data(ii0 + 98);
    const auto *ii0_99 = buffer.data(ii0 + 99);
    const auto *ii0_105 = buffer.data(ii0 + 105);
    const auto *ii0_106 = buffer.data(ii0 + 106);
    const auto *ii0_107 = buffer.data(ii0 + 107);
    const auto *ii0_108 = buffer.data(ii0 + 108);
    const auto *ii0_109 = buffer.data(ii0 + 109);
    const auto *ii0_111 = buffer.data(ii0 + 111);
    const auto *ii0_140 = buffer.data(ii0 + 140);
    const auto *ii0_141 = buffer.data(ii0 + 141);
    const auto *ii0_143 = buffer.data(ii0 + 143);
    const auto *ii0_145 = buffer.data(ii0 + 145);
    const auto *ii0_146 = buffer.data(ii0 + 146);
    const auto *ii0_148 = buffer.data(ii0 + 148);
    const auto *ii0_149 = buffer.data(ii0 + 149);
    const auto *ii0_150 = buffer.data(ii0 + 150);
    const auto *ii0_152 = buffer.data(ii0 + 152);
    const auto *ii0_153 = buffer.data(ii0 + 153);
    const auto *ii0_154 = buffer.data(ii0 + 154);
    const auto *ii0_160 = buffer.data(ii0 + 160);
    const auto *ii0_161 = buffer.data(ii0 + 161);
    const auto *ii0_163 = buffer.data(ii0 + 163);
    const auto *ii0_164 = buffer.data(ii0 + 164);
    const auto *ii0_165 = buffer.data(ii0 + 165);
    const auto *ii0_166 = buffer.data(ii0 + 166);
    const auto *ii0_167 = buffer.data(ii0 + 167);
    const auto *ii0_168 = buffer.data(ii0 + 168);
    const auto *ii0_171 = buffer.data(ii0 + 171);

    const auto *ii1_87 = buffer.data(ii1 + 87);
    const auto *ii1_89 = buffer.data(ii1 + 89);
    const auto *ii1_90 = buffer.data(ii1 + 90);
    const auto *ii1_91 = buffer.data(ii1 + 91);
    const auto *ii1_93 = buffer.data(ii1 + 93);
    const auto *ii1_94 = buffer.data(ii1 + 94);
    const auto *ii1_95 = buffer.data(ii1 + 95);
    const auto *ii1_96 = buffer.data(ii1 + 96);
    const auto *ii1_98 = buffer.data(ii1 + 98);
    const auto *ii1_99 = buffer.data(ii1 + 99);
    const auto *ii1_105 = buffer.data(ii1 + 105);
    const auto *ii1_106 = buffer.data(ii1 + 106);
    const auto *ii1_107 = buffer.data(ii1 + 107);
    const auto *ii1_108 = buffer.data(ii1 + 108);
    const auto *ii1_109 = buffer.data(ii1 + 109);
    const auto *ii1_111 = buffer.data(ii1 + 111);
    const auto *ii1_140 = buffer.data(ii1 + 140);
    const auto *ii1_141 = buffer.data(ii1 + 141);
    const auto *ii1_143 = buffer.data(ii1 + 143);
    const auto *ii1_145 = buffer.data(ii1 + 145);
    const auto *ii1_146 = buffer.data(ii1 + 146);
    const auto *ii1_148 = buffer.data(ii1 + 148);
    const auto *ii1_149 = buffer.data(ii1 + 149);
    const auto *ii1_150 = buffer.data(ii1 + 150);
    const auto *ii1_152 = buffer.data(ii1 + 152);
    const auto *ii1_153 = buffer.data(ii1 + 153);
    const auto *ii1_154 = buffer.data(ii1 + 154);
    const auto *ii1_160 = buffer.data(ii1 + 160);
    const auto *ii1_161 = buffer.data(ii1 + 161);
    const auto *ii1_163 = buffer.data(ii1 + 163);
    const auto *ii1_164 = buffer.data(ii1 + 164);
    const auto *ii1_165 = buffer.data(ii1 + 165);
    const auto *ii1_166 = buffer.data(ii1 + 166);
    const auto *ii1_167 = buffer.data(ii1 + 167);
    const auto *ii1_168 = buffer.data(ii1 + 168);
    const auto *ii1_171 = buffer.data(ii1 + 171);

    const auto *ik_114 = buffer.data(ik + 114);
    const auto *ik_115 = buffer.data(ik + 115);
    const auto *ik_117 = buffer.data(ik + 117);
    const auto *ik_118 = buffer.data(ik + 118);
    const auto *ik_119 = buffer.data(ik + 119);
    const auto *ik_120 = buffer.data(ik + 120);
    const auto *ik_122 = buffer.data(ik + 122);
    const auto *ik_123 = buffer.data(ik + 123);
    const auto *ik_124 = buffer.data(ik + 124);
    const auto *ik_125 = buffer.data(ik + 125);
    const auto *ik_126 = buffer.data(ik + 126);
    const auto *ik_128 = buffer.data(ik + 128);
    const auto *ik_129 = buffer.data(ik + 129);
    const auto *ik_136 = buffer.data(ik + 136);
    const auto *ik_137 = buffer.data(ik + 137);
    const auto *ik_138 = buffer.data(ik + 138);
    const auto *ik_139 = buffer.data(ik + 139);
    const auto *ik_140 = buffer.data(ik + 140);
    const auto *ik_141 = buffer.data(ik + 141);
    const auto *ik_142 = buffer.data(ik + 142);
    const auto *ik_143 = buffer.data(ik + 143);
    const auto *ik_146 = buffer.data(ik + 146);
    const auto *ik_147 = buffer.data(ik + 147);
    const auto *ik_149 = buffer.data(ik + 149);
    const auto *ik_150 = buffer.data(ik + 150);
    const auto *ik_153 = buffer.data(ik + 153);
    const auto *ik_154 = buffer.data(ik + 154);
    const auto *ik_158 = buffer.data(ik + 158);
    const auto *ik_159 = buffer.data(ik + 159);
    const auto *ik_164 = buffer.data(ik + 164);
    const auto *ik_172 = buffer.data(ik + 172);
    const auto *ik_173 = buffer.data(ik + 173);
    const auto *ik_174 = buffer.data(ik + 174);
    const auto *ik_175 = buffer.data(ik + 175);
    const auto *ik_176 = buffer.data(ik + 176);
    const auto *ik_177 = buffer.data(ik + 177);
    const auto *ik_178 = buffer.data(ik + 178);
    const auto *ik_179 = buffer.data(ik + 179);
    const auto *ik_180 = buffer.data(ik + 180);
    const auto *ik_181 = buffer.data(ik + 181);
    const auto *ik_182 = buffer.data(ik + 182);
    const auto *ik_183 = buffer.data(ik + 183);
    const auto *ik_185 = buffer.data(ik + 185);
    const auto *ik_186 = buffer.data(ik + 186);
    const auto *ik_188 = buffer.data(ik + 188);
    const auto *ik_189 = buffer.data(ik + 189);
    const auto *ik_190 = buffer.data(ik + 190);
    const auto *ik_192 = buffer.data(ik + 192);
    const auto *ik_193 = buffer.data(ik + 193);
    const auto *ik_194 = buffer.data(ik + 194);
    const auto *ik_195 = buffer.data(ik + 195);
    const auto *ik_197 = buffer.data(ik + 197);
    const auto *ik_198 = buffer.data(ik + 198);
    const auto *ik_199 = buffer.data(ik + 199);
    const auto *ik_200 = buffer.data(ik + 200);
    const auto *ik_207 = buffer.data(ik + 207);
    const auto *ik_208 = buffer.data(ik + 208);
    const auto *ik_209 = buffer.data(ik + 209);
    const auto *ik_210 = buffer.data(ik + 210);
    const auto *ik_211 = buffer.data(ik + 211);
    const auto *ik_212 = buffer.data(ik + 212);
    const auto *ik_213 = buffer.data(ik + 213);
    const auto *ik_214 = buffer.data(ik + 214);
    const auto *ik_215 = buffer.data(ik + 215);
    const auto *ik_216 = buffer.data(ik + 216);
    const auto *ik_217 = buffer.data(ik + 217);
    const auto *ik_218 = buffer.data(ik + 218);
    const auto *ik_219 = buffer.data(ik + 219);

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_z, hk_118, ii0_87, ii0_94, ii1_87, \
                         ii1_94, ik_114, ik_115, ik_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_16 * hk_118[k]
                   + f_7 * ii0_94[k]
                   - f_8 * ii1_94[k]
                   + pb_x[k] * ik_118[k];

        t_146[k] = pb_z[k] * ik_114[k];

        t_147[k] = f_3 * ii0_87[k]
                   - f_4 * ii1_87[k]
                   + pb_z[k] * ik_115[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pb_y, pb_z, hk_45, hk_123, ii0_89, \
                         ii0_99, ii1_89, ii1_99, ik_117, ik_118, \
                         ik_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * hk_45[k]
                   + pb_y[k] * ik_117[k];

        t_149[k] = f_7 * ii0_89[k]
                   - f_8 * ii1_89[k]
                   + pb_z[k] * ik_117[k];

        t_150[k] = f_16 * hk_123[k]
                   + f_5 * ii0_99[k]
                   - f_6 * ii1_99[k]
                   + pb_x[k] * ik_123[k];

        t_151[k] = pb_z[k] * ik_118[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, hk_50, ii0_90, ii0_91, \
                         ii0_93, ii1_90, ii1_91, ii1_93, ik_119, ik_120, \
                         ik_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * ii0_90[k]
                   - f_4 * ii1_90[k]
                   + pb_z[k] * ik_119[k];

        t_153[k] = f_5 * ii0_91[k]
                   - f_6 * ii1_91[k]
                   + pb_z[k] * ik_120[k];

        t_154[k] = f_14 * hk_50[k]
                   + pb_y[k] * ik_122[k];

        t_155[k] = f_9 * ii0_93[k]
                   - f_10 * ii1_93[k]
                   + pb_z[k] * ik_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_z, hk_129, ii0_94, ii0_105, ii1_94, \
                         ii1_105, ik_123, ik_124, ik_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * hk_129[k]
                   + f_3 * ii0_105[k]
                   - f_4 * ii1_105[k]
                   + pb_x[k] * ik_129[k];

        t_157[k] = pb_z[k] * ik_123[k];

        t_158[k] = f_3 * ii0_94[k]
                   - f_4 * ii1_94[k]
                   + pb_z[k] * ik_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_y, pb_z, hk_56, ii0_95, ii0_96, \
                         ii0_98, ii1_95, ii1_96, ii1_98, ik_125, ik_126, \
                         ik_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * ii0_95[k]
                   - f_6 * ii1_95[k]
                   + pb_z[k] * ik_125[k];

        t_160[k] = f_7 * ii0_96[k]
                   - f_8 * ii1_96[k]
                   + pb_z[k] * ik_126[k];

        t_161[k] = f_14 * hk_56[k]
                   + pb_y[k] * ik_128[k];

        t_162[k] = f_11 * ii0_98[k]
                   - f_12 * ii1_98[k]
                   + pb_z[k] * ik_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_x, pb_z, hk_136, hk_138, \
                         hk_139, hk_140, ik_129, ik_136, ik_138, ik_139, \
                         ik_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_16 * hk_136[k]
                   + pb_x[k] * ik_136[k];

        t_164[k] = pb_z[k] * ik_129[k];

        t_165[k] = f_16 * hk_138[k]
                   + pb_x[k] * ik_138[k];

        t_166[k] = f_16 * hk_139[k]
                   + pb_x[k] * ik_139[k];

        t_167[k] = f_16 * hk_140[k]
                   + pb_x[k] * ik_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pb_x, gl0_171, gl1_171, hk_141, \
                         hk_142, hk_143, hl_171, ik_141, ik_142, \
                         ik_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_16 * hk_141[k]
                   + pb_x[k] * ik_141[k];

        t_169[k] = f_16 * hk_142[k]
                   + pb_x[k] * ik_142[k];

        t_170[k] = f_16 * hk_143[k]
                   + pb_x[k] * ik_143[k];

        t_171[k] = f_21 * gl0_171[k]
                   - f_22 * gl1_171[k]
                   + pa_x[k] * hl_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, ii0_105, ii0_106, ii0_107, ii1_105, \
                         ii1_106, ii1_107, ik_136, ik_137, ik_138, \
                         ik_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * ik_136[k];

        t_173[k] = f_3 * ii0_105[k]
                   - f_4 * ii1_105[k]
                   + pb_z[k] * ik_137[k];

        t_174[k] = f_5 * ii0_106[k]
                   - f_6 * ii1_106[k]
                   + pb_z[k] * ik_138[k];

        t_175[k] = f_7 * ii0_107[k]
                   - f_8 * ii1_107[k]
                   + pb_z[k] * ik_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pb_z, hk_71, ii0_108, ii0_109, \
                         ii0_111, ii1_108, ii1_109, ii1_111, ik_140, ik_141, \
                         ik_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * ii0_108[k]
                   - f_10 * ii1_108[k]
                   + pb_z[k] * ik_140[k];

        t_177[k] = f_11 * ii0_109[k]
                   - f_12 * ii1_109[k]
                   + pb_z[k] * ik_141[k];

        t_178[k] = f_14 * hk_71[k]
                   + pb_y[k] * ik_143[k];

        t_179[k] = f_1 * ii0_111[k]
                   - f_2 * ii1_111[k]
                   + pb_z[k] * ik_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pa_y, pa_z, pb_y, hk_74, \
                         hl_46, hl_48, hl_90, hl_92, hl_95, ik_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * hl_90[k];

        t_181[k] = pa_z[k] * hl_46[k];

        t_182[k] = pa_y[k] * hl_92[k];

        t_183[k] = pa_z[k] * hl_48[k];

        t_184[k] = f_13 * hk_74[k]
                   + pb_y[k] * ik_146[k];

        t_185[k] = pa_y[k] * hl_95[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, hk_39, \
                         hk_77, hl_51, hl_55, hl_99, ik_147, ik_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * hl_51[k];

        t_187[k] = f_13 * hk_39[k]
                   + pb_z[k] * ik_147[k];

        t_188[k] = f_13 * hk_77[k]
                   + pb_y[k] * ik_149[k];

        t_189[k] = pa_y[k] * hl_99[k];

        t_190[k] = pa_z[k] * hl_55[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, hk_42, hk_80, hk_81, \
                         hl_102, hl_104, ik_150, ik_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * hk_42[k]
                   + pb_z[k] * ik_150[k];

        t_192[k] = f_14 * hk_80[k]
                   + pa_y[k] * hl_102[k];

        t_193[k] = f_13 * hk_81[k]
                   + pb_y[k] * ik_153[k];

        t_194[k] = pa_y[k] * hl_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pa_z, pb_z, hk_46, hk_84, hk_85, \
                         hl_60, hl_107, hl_108, ik_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * hl_60[k];

        t_196[k] = f_13 * hk_46[k]
                   + pb_z[k] * ik_154[k];

        t_197[k] = f_15 * hk_84[k]
                   + pa_y[k] * hl_107[k];

        t_198[k] = f_14 * hk_85[k]
                   + pa_y[k] * hl_108[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, hk_51, hk_86, \
                         hl_66, hl_110, ik_158, ik_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * hk_86[k]
                   + pb_y[k] * ik_158[k];

        t_200[k] = pa_y[k] * hl_110[k];

        t_201[k] = pa_z[k] * hl_66[k];

        t_202[k] = f_13 * hk_51[k]
                   + pb_z[k] * ik_159[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pa_y, pb_y, hk_89, hk_90, hk_91, \
                         hk_92, hl_113, hl_114, hl_115, hl_117, \
                         ik_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_16 * hk_89[k]
                   + pa_y[k] * hl_113[k];

        t_204[k] = f_15 * hk_90[k]
                   + pa_y[k] * hl_114[k];

        t_205[k] = f_14 * hk_91[k]
                   + pa_y[k] * hl_115[k];

        t_206[k] = f_13 * hk_92[k]
                   + pb_y[k] * ik_164[k];

        t_207[k] = pa_y[k] * hl_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pa_z, pb_x, hk_173, hk_174, \
                         hk_175, hk_176, hl_73, ik_173, ik_174, ik_175, \
                         ik_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_z[k] * hl_73[k];

        t_209[k] = f_16 * hk_173[k]
                   + pb_x[k] * ik_173[k];

        t_210[k] = f_16 * hk_174[k]
                   + pb_x[k] * ik_174[k];

        t_211[k] = f_16 * hk_175[k]
                   + pb_x[k] * ik_175[k];

        t_212[k] = f_16 * hk_176[k]
                   + pb_x[k] * ik_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pa_z, pb_x, hk_177, hk_178, hl_81, \
                         hl_125, ik_177, ik_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_16 * hk_177[k]
                   + pb_x[k] * ik_177[k];

        t_214[k] = f_16 * hk_178[k]
                   + pb_x[k] * ik_178[k];

        t_215[k] = pa_y[k] * hl_125[k];

        t_216[k] = pa_z[k] * hl_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_y, pb_z, hk_64, hk_102, hk_103, \
                         hk_104, hl_128, hl_129, hl_130, ik_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_13 * hk_64[k]
                   + pb_z[k] * ik_172[k];

        t_218[k] = f_0 * hk_102[k]
                   + pa_y[k] * hl_128[k];

        t_219[k] = f_17 * hk_103[k]
                   + pa_y[k] * hl_129[k];

        t_220[k] = f_16 * hk_104[k]
                   + pa_y[k] * hl_130[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pb_y, hk_105, hk_106, hk_107, \
                         hl_131, hl_132, hl_134, ik_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_15 * hk_105[k]
                   + pa_y[k] * hl_131[k];

        t_222[k] = f_14 * hk_106[k]
                   + pa_y[k] * hl_132[k];

        t_223[k] = f_13 * hk_107[k]
                   + pb_y[k] * ik_179[k];

        t_224[k] = pa_y[k] * hl_134[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_z, pb_y, pb_z, gl0_0, gl1_0, hk_72, \
                         hl_90, ii0_140, ii1_140, ik_180, ik_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_19 * gl0_0[k]
                   - f_20 * gl1_0[k]
                   + pa_z[k] * hl_90[k];

        t_226[k] = pb_y[k] * ik_180[k];

        t_227[k] = f_14 * hk_72[k]
                   + pb_z[k] * ik_180[k];

        t_228[k] = f_3 * ii0_140[k]
                   - f_4 * ii1_140[k]
                   + pb_y[k] * ik_181[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, hk_75, hk_185, ii0_141, \
                         ii0_145, ii1_141, ii1_145, ik_182, ik_183, \
                         ik_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * ik_182[k];

        t_230[k] = f_16 * hk_185[k]
                   + f_11 * ii0_145[k]
                   - f_12 * ii1_145[k]
                   + pb_x[k] * ik_185[k];

        t_231[k] = f_5 * ii0_141[k]
                   - f_6 * ii1_141[k]
                   + pb_y[k] * ik_183[k];

        t_232[k] = f_14 * hk_75[k]
                   + pb_z[k] * ik_183[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pb_y, pb_z, hk_78, hk_189, ii0_143, \
                         ii0_149, ii1_143, ii1_149, ik_185, ik_186, \
                         ik_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * ik_185[k];

        t_234[k] = f_16 * hk_189[k]
                   + f_9 * ii0_149[k]
                   - f_10 * ii1_149[k]
                   + pb_x[k] * ik_189[k];

        t_235[k] = f_7 * ii0_143[k]
                   - f_8 * ii1_143[k]
                   + pb_y[k] * ik_186[k];

        t_236[k] = f_14 * hk_78[k]
                   + pb_z[k] * ik_186[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_x, pb_y, hk_194, ii0_145, ii0_154, ii1_145, \
                         ii1_154, ik_188, ik_189, ik_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_3 * ii0_145[k]
                   - f_4 * ii1_145[k]
                   + pb_y[k] * ik_188[k];

        t_238[k] = pb_y[k] * ik_189[k];

        t_239[k] = f_16 * hk_194[k]
                   + f_7 * ii0_154[k]
                   - f_8 * ii1_154[k]
                   + pb_x[k] * ik_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_y, pb_z, hk_82, ii0_146, ii0_148, \
                         ii0_149, ii1_146, ii1_148, ii1_149, ik_190, ik_192, \
                         ik_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_9 * ii0_146[k]
                   - f_10 * ii1_146[k]
                   + pb_y[k] * ik_190[k];

        t_241[k] = f_14 * hk_82[k]
                   + pb_z[k] * ik_190[k];

        t_242[k] = f_5 * ii0_148[k]
                   - f_6 * ii1_148[k]
                   + pb_y[k] * ik_192[k];

        t_243[k] = f_3 * ii0_149[k]
                   - f_4 * ii1_149[k]
                   + pb_y[k] * ik_193[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pb_y, pb_z, hk_87, hk_200, ii0_150, \
                         ii0_160, ii1_150, ii1_160, ik_194, ik_195, \
                         ik_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_y[k] * ik_194[k];

        t_245[k] = f_16 * hk_200[k]
                   + f_5 * ii0_160[k]
                   - f_6 * ii1_160[k]
                   + pb_x[k] * ik_200[k];

        t_246[k] = f_11 * ii0_150[k]
                   - f_12 * ii1_150[k]
                   + pb_y[k] * ik_195[k];

        t_247[k] = f_14 * hk_87[k]
                   + pb_z[k] * ik_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_y, ii0_152, ii0_153, ii0_154, ii1_152, \
                         ii1_153, ii1_154, ik_197, ik_198, ik_199, \
                         ik_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * ii0_152[k]
                   - f_8 * ii1_152[k]
                   + pb_y[k] * ik_197[k];

        t_249[k] = f_5 * ii0_153[k]
                   - f_6 * ii1_153[k]
                   + pb_y[k] * ik_198[k];

        t_250[k] = f_3 * ii0_154[k]
                   - f_4 * ii1_154[k]
                   + pb_y[k] * ik_199[k];

        t_251[k] = pb_y[k] * ik_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, hk_207, hk_208, hk_209, hk_210, \
                         ii0_167, ii1_167, ik_207, ik_208, ik_209, \
                         ik_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_16 * hk_207[k]
                   + f_3 * ii0_167[k]
                   - f_4 * ii1_167[k]
                   + pb_x[k] * ik_207[k];

        t_253[k] = f_16 * hk_208[k]
                   + pb_x[k] * ik_208[k];

        t_254[k] = f_16 * hk_209[k]
                   + pb_x[k] * ik_209[k];

        t_255[k] = f_16 * hk_210[k]
                   + pb_x[k] * ik_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, hk_211, hk_212, \
                         hk_213, hk_215, ik_207, ik_211, ik_212, ik_213, \
                         ik_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_16 * hk_211[k]
                   + pb_x[k] * ik_211[k];

        t_257[k] = f_16 * hk_212[k]
                   + pb_x[k] * ik_212[k];

        t_258[k] = f_16 * hk_213[k]
                   + pb_x[k] * ik_213[k];

        t_259[k] = pb_y[k] * ik_207[k];

        t_260[k] = f_16 * hk_215[k]
                   + pb_x[k] * ik_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_y, pb_z, hk_100, ii0_161, ii0_163, \
                         ii0_164, ii1_161, ii1_163, ii1_164, ik_208, ik_210, \
                         ik_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * ii0_161[k]
                   - f_2 * ii1_161[k]
                   + pb_y[k] * ik_208[k];

        t_262[k] = f_14 * hk_100[k]
                   + pb_z[k] * ik_208[k];

        t_263[k] = f_11 * ii0_163[k]
                   - f_12 * ii1_163[k]
                   + pb_y[k] * ik_210[k];

        t_264[k] = f_9 * ii0_164[k]
                   - f_10 * ii1_164[k]
                   + pb_y[k] * ik_211[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_y, ii0_165, ii0_166, ii0_167, ii1_165, \
                         ii1_166, ii1_167, ik_212, ik_213, ik_214, \
                         ik_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_7 * ii0_165[k]
                   - f_8 * ii1_165[k]
                   + pb_y[k] * ik_212[k];

        t_266[k] = f_5 * ii0_166[k]
                   - f_6 * ii1_166[k]
                   + pb_y[k] * ik_213[k];

        t_267[k] = f_3 * ii0_167[k]
                   - f_4 * ii1_167[k]
                   + pb_y[k] * ik_214[k];

        t_268[k] = pb_y[k] * ik_215[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pa_y, pb_y, pb_z, gl0_45, gl0_269, \
                         gl1_45, gl1_269, hk_108, hl_135, hl_269, \
                         ik_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_21 * gl0_269[k]
                   - f_22 * gl1_269[k]
                   + pa_x[k] * hl_269[k];

        t_270[k] = f_23 * gl0_45[k]
                   - f_24 * gl1_45[k]
                   + pa_y[k] * hl_135[k];

        t_271[k] = f_15 * hk_108[k]
                   + pb_y[k] * ik_216[k];

        t_272[k] = pb_z[k] * ik_216[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pb_x, pb_z, hk_219, ii0_168, ii0_171, ii1_168, \
                         ii1_171, ik_217, ik_218, ik_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_15 * hk_219[k]
                   + f_11 * ii0_171[k]
                   - f_12 * ii1_171[k]
                   + pb_x[k] * ik_219[k];

        t_274[k] = pb_z[k] * ik_217[k];

        t_275[k] = f_3 * ii0_168[k]
                   - f_4 * ii1_168[k]
                   + pb_z[k] * ik_218[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_23 = 1.0 / alpha;
    const auto f_24 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_90 = buffer.data(gl0 + 90);
    const auto *gl0_306 = buffer.data(gl0 + 306);

    const auto *gl1_90 = buffer.data(gl1 + 90);
    const auto *gl1_306 = buffer.data(gl1 + 306);

    const auto *hk_108 = buffer.data(hk + 108);
    const auto *hk_110 = buffer.data(hk + 110);
    const auto *hk_111 = buffer.data(hk + 111);
    const auto *hk_113 = buffer.data(hk + 113);
    const auto *hk_114 = buffer.data(hk + 114);
    const auto *hk_115 = buffer.data(hk + 115);
    const auto *hk_117 = buffer.data(hk + 117);
    const auto *hk_118 = buffer.data(hk + 118);
    const auto *hk_119 = buffer.data(hk + 119);
    const auto *hk_120 = buffer.data(hk + 120);
    const auto *hk_122 = buffer.data(hk + 122);
    const auto *hk_123 = buffer.data(hk + 123);
    const auto *hk_124 = buffer.data(hk + 124);
    const auto *hk_125 = buffer.data(hk + 125);
    const auto *hk_126 = buffer.data(hk + 126);
    const auto *hk_128 = buffer.data(hk + 128);
    const auto *hk_136 = buffer.data(hk + 136);
    const auto *hk_137 = buffer.data(hk + 137);
    const auto *hk_138 = buffer.data(hk + 138);
    const auto *hk_139 = buffer.data(hk + 139);
    const auto *hk_140 = buffer.data(hk + 140);
    const auto *hk_141 = buffer.data(hk + 141);
    const auto *hk_143 = buffer.data(hk + 143);
    const auto *hk_146 = buffer.data(hk + 146);
    const auto *hk_147 = buffer.data(hk + 147);
    const auto *hk_149 = buffer.data(hk + 149);
    const auto *hk_150 = buffer.data(hk + 150);
    const auto *hk_153 = buffer.data(hk + 153);
    const auto *hk_154 = buffer.data(hk + 154);
    const auto *hk_158 = buffer.data(hk + 158);
    const auto *hk_159 = buffer.data(hk + 159);
    const auto *hk_164 = buffer.data(hk + 164);
    const auto *hk_172 = buffer.data(hk + 172);
    const auto *hk_179 = buffer.data(hk + 179);
    const auto *hk_180 = buffer.data(hk + 180);
    const auto *hk_181 = buffer.data(hk + 181);
    const auto *hk_182 = buffer.data(hk + 182);
    const auto *hk_183 = buffer.data(hk + 183);
    const auto *hk_185 = buffer.data(hk + 185);
    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_188 = buffer.data(hk + 188);
    const auto *hk_189 = buffer.data(hk + 189);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_192 = buffer.data(hk + 192);
    const auto *hk_193 = buffer.data(hk + 193);
    const auto *hk_194 = buffer.data(hk + 194);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_197 = buffer.data(hk + 197);
    const auto *hk_198 = buffer.data(hk + 198);
    const auto *hk_199 = buffer.data(hk + 199);
    const auto *hk_200 = buffer.data(hk + 200);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_210 = buffer.data(hk + 210);
    const auto *hk_211 = buffer.data(hk + 211);
    const auto *hk_212 = buffer.data(hk + 212);
    const auto *hk_213 = buffer.data(hk + 213);
    const auto *hk_214 = buffer.data(hk + 214);
    const auto *hk_215 = buffer.data(hk + 215);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_237 = buffer.data(hk + 237);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_250 = buffer.data(hk + 250);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_281 = buffer.data(hk + 281);
    const auto *hk_282 = buffer.data(hk + 282);
    const auto *hk_283 = buffer.data(hk + 283);
    const auto *hk_284 = buffer.data(hk + 284);
    const auto *hk_285 = buffer.data(hk + 285);
    const auto *hk_286 = buffer.data(hk + 286);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_317 = buffer.data(hk + 317);
    const auto *hk_318 = buffer.data(hk + 318);
    const auto *hk_319 = buffer.data(hk + 319);
    const auto *hk_320 = buffer.data(hk + 320);
    const auto *hk_321 = buffer.data(hk + 321);
    const auto *hk_322 = buffer.data(hk + 322);
    const auto *hk_329 = buffer.data(hk + 329);

    const auto *hl_135 = buffer.data(hl + 135);
    const auto *hl_136 = buffer.data(hl + 136);
    const auto *hl_138 = buffer.data(hl + 138);
    const auto *hl_140 = buffer.data(hl + 140);
    const auto *hl_141 = buffer.data(hl + 141);
    const auto *hl_144 = buffer.data(hl + 144);
    const auto *hl_145 = buffer.data(hl + 145);
    const auto *hl_147 = buffer.data(hl + 147);
    const auto *hl_149 = buffer.data(hl + 149);
    const auto *hl_150 = buffer.data(hl + 150);
    const auto *hl_152 = buffer.data(hl + 152);
    const auto *hl_153 = buffer.data(hl + 153);
    const auto *hl_155 = buffer.data(hl + 155);
    const auto *hl_156 = buffer.data(hl + 156);
    const auto *hl_158 = buffer.data(hl + 158);
    const auto *hl_159 = buffer.data(hl + 159);
    const auto *hl_160 = buffer.data(hl + 160);
    const auto *hl_162 = buffer.data(hl + 162);
    const auto *hl_163 = buffer.data(hl + 163);
    const auto *hl_171 = buffer.data(hl + 171);
    const auto *hl_173 = buffer.data(hl + 173);
    const auto *hl_174 = buffer.data(hl + 174);
    const auto *hl_175 = buffer.data(hl + 175);
    const auto *hl_176 = buffer.data(hl + 176);
    const auto *hl_177 = buffer.data(hl + 177);
    const auto *hl_179 = buffer.data(hl + 179);
    const auto *hl_225 = buffer.data(hl + 225);
    const auto *hl_227 = buffer.data(hl + 227);
    const auto *hl_228 = buffer.data(hl + 228);
    const auto *hl_230 = buffer.data(hl + 230);
    const auto *hl_231 = buffer.data(hl + 231);
    const auto *hl_234 = buffer.data(hl + 234);
    const auto *hl_235 = buffer.data(hl + 235);
    const auto *hl_237 = buffer.data(hl + 237);
    const auto *hl_239 = buffer.data(hl + 239);
    const auto *hl_240 = buffer.data(hl + 240);
    const auto *hl_242 = buffer.data(hl + 242);
    const auto *hl_243 = buffer.data(hl + 243);
    const auto *hl_245 = buffer.data(hl + 245);
    const auto *hl_246 = buffer.data(hl + 246);
    const auto *hl_248 = buffer.data(hl + 248);
    const auto *hl_249 = buffer.data(hl + 249);
    const auto *hl_250 = buffer.data(hl + 250);
    const auto *hl_252 = buffer.data(hl + 252);
    const auto *hl_260 = buffer.data(hl + 260);
    const auto *hl_261 = buffer.data(hl + 261);
    const auto *hl_263 = buffer.data(hl + 263);
    const auto *hl_264 = buffer.data(hl + 264);
    const auto *hl_265 = buffer.data(hl + 265);
    const auto *hl_266 = buffer.data(hl + 266);
    const auto *hl_267 = buffer.data(hl + 267);
    const auto *hl_269 = buffer.data(hl + 269);
    const auto *hl_306 = buffer.data(hl + 306);

    const auto *ii0_170 = buffer.data(ii0 + 170);
    const auto *ii0_171 = buffer.data(ii0 + 171);
    const auto *ii0_173 = buffer.data(ii0 + 173);
    const auto *ii0_174 = buffer.data(ii0 + 174);
    const auto *ii0_175 = buffer.data(ii0 + 175);
    const auto *ii0_177 = buffer.data(ii0 + 177);
    const auto *ii0_178 = buffer.data(ii0 + 178);
    const auto *ii0_179 = buffer.data(ii0 + 179);
    const auto *ii0_180 = buffer.data(ii0 + 180);
    const auto *ii0_182 = buffer.data(ii0 + 182);
    const auto *ii0_183 = buffer.data(ii0 + 183);
    const auto *ii0_189 = buffer.data(ii0 + 189);
    const auto *ii0_190 = buffer.data(ii0 + 190);
    const auto *ii0_191 = buffer.data(ii0 + 191);
    const auto *ii0_192 = buffer.data(ii0 + 192);
    const auto *ii0_193 = buffer.data(ii0 + 193);
    const auto *ii0_195 = buffer.data(ii0 + 195);
    const auto *ii0_252 = buffer.data(ii0 + 252);
    const auto *ii0_253 = buffer.data(ii0 + 253);
    const auto *ii0_257 = buffer.data(ii0 + 257);

    const auto *ii1_170 = buffer.data(ii1 + 170);
    const auto *ii1_171 = buffer.data(ii1 + 171);
    const auto *ii1_173 = buffer.data(ii1 + 173);
    const auto *ii1_174 = buffer.data(ii1 + 174);
    const auto *ii1_175 = buffer.data(ii1 + 175);
    const auto *ii1_177 = buffer.data(ii1 + 177);
    const auto *ii1_178 = buffer.data(ii1 + 178);
    const auto *ii1_179 = buffer.data(ii1 + 179);
    const auto *ii1_180 = buffer.data(ii1 + 180);
    const auto *ii1_182 = buffer.data(ii1 + 182);
    const auto *ii1_183 = buffer.data(ii1 + 183);
    const auto *ii1_189 = buffer.data(ii1 + 189);
    const auto *ii1_190 = buffer.data(ii1 + 190);
    const auto *ii1_191 = buffer.data(ii1 + 191);
    const auto *ii1_192 = buffer.data(ii1 + 192);
    const auto *ii1_193 = buffer.data(ii1 + 193);
    const auto *ii1_195 = buffer.data(ii1 + 195);
    const auto *ii1_252 = buffer.data(ii1 + 252);
    const auto *ii1_253 = buffer.data(ii1 + 253);
    const auto *ii1_257 = buffer.data(ii1 + 257);

    const auto *ik_219 = buffer.data(ik + 219);
    const auto *ik_221 = buffer.data(ik + 221);
    const auto *ik_222 = buffer.data(ik + 222);
    const auto *ik_223 = buffer.data(ik + 223);
    const auto *ik_225 = buffer.data(ik + 225);
    const auto *ik_226 = buffer.data(ik + 226);
    const auto *ik_227 = buffer.data(ik + 227);
    const auto *ik_228 = buffer.data(ik + 228);
    const auto *ik_230 = buffer.data(ik + 230);
    const auto *ik_231 = buffer.data(ik + 231);
    const auto *ik_232 = buffer.data(ik + 232);
    const auto *ik_233 = buffer.data(ik + 233);
    const auto *ik_234 = buffer.data(ik + 234);
    const auto *ik_236 = buffer.data(ik + 236);
    const auto *ik_237 = buffer.data(ik + 237);
    const auto *ik_244 = buffer.data(ik + 244);
    const auto *ik_245 = buffer.data(ik + 245);
    const auto *ik_246 = buffer.data(ik + 246);
    const auto *ik_247 = buffer.data(ik + 247);
    const auto *ik_248 = buffer.data(ik + 248);
    const auto *ik_249 = buffer.data(ik + 249);
    const auto *ik_250 = buffer.data(ik + 250);
    const auto *ik_251 = buffer.data(ik + 251);
    const auto *ik_252 = buffer.data(ik + 252);
    const auto *ik_254 = buffer.data(ik + 254);
    const auto *ik_255 = buffer.data(ik + 255);
    const auto *ik_257 = buffer.data(ik + 257);
    const auto *ik_258 = buffer.data(ik + 258);
    const auto *ik_261 = buffer.data(ik + 261);
    const auto *ik_262 = buffer.data(ik + 262);
    const auto *ik_266 = buffer.data(ik + 266);
    const auto *ik_267 = buffer.data(ik + 267);
    const auto *ik_272 = buffer.data(ik + 272);
    const auto *ik_280 = buffer.data(ik + 280);
    const auto *ik_281 = buffer.data(ik + 281);
    const auto *ik_282 = buffer.data(ik + 282);
    const auto *ik_283 = buffer.data(ik + 283);
    const auto *ik_284 = buffer.data(ik + 284);
    const auto *ik_285 = buffer.data(ik + 285);
    const auto *ik_286 = buffer.data(ik + 286);
    const auto *ik_287 = buffer.data(ik + 287);
    const auto *ik_288 = buffer.data(ik + 288);
    const auto *ik_290 = buffer.data(ik + 290);
    const auto *ik_291 = buffer.data(ik + 291);
    const auto *ik_293 = buffer.data(ik + 293);
    const auto *ik_294 = buffer.data(ik + 294);
    const auto *ik_297 = buffer.data(ik + 297);
    const auto *ik_298 = buffer.data(ik + 298);
    const auto *ik_302 = buffer.data(ik + 302);
    const auto *ik_303 = buffer.data(ik + 303);
    const auto *ik_308 = buffer.data(ik + 308);
    const auto *ik_316 = buffer.data(ik + 316);
    const auto *ik_317 = buffer.data(ik + 317);
    const auto *ik_318 = buffer.data(ik + 318);
    const auto *ik_319 = buffer.data(ik + 319);
    const auto *ik_320 = buffer.data(ik + 320);
    const auto *ik_321 = buffer.data(ik + 321);
    const auto *ik_322 = buffer.data(ik + 322);
    const auto *ik_323 = buffer.data(ik + 323);
    const auto *ik_324 = buffer.data(ik + 324);
    const auto *ik_325 = buffer.data(ik + 325);
    const auto *ik_326 = buffer.data(ik + 326);
    const auto *ik_327 = buffer.data(ik + 327);
    const auto *ik_329 = buffer.data(ik + 329);

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, pb_y, pb_z, hk_113, hk_222, \
                         ii0_170, ii0_174, ii1_170, ii1_174, ik_219, ik_221, \
                         ik_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_15 * hk_222[k]
                   + f_9 * ii0_174[k]
                   - f_10 * ii1_174[k]
                   + pb_x[k] * ik_222[k];

        t_277[k] = pb_z[k] * ik_219[k];

        t_278[k] = f_15 * hk_113[k]
                   + pb_y[k] * ik_221[k];

        t_279[k] = f_5 * ii0_170[k]
                   - f_6 * ii1_170[k]
                   + pb_z[k] * ik_221[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_x, pb_z, hk_226, ii0_171, ii0_178, ii1_171, \
                         ii1_178, ik_222, ik_223, ik_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_15 * hk_226[k]
                   + f_7 * ii0_178[k]
                   - f_8 * ii1_178[k]
                   + pb_x[k] * ik_226[k];

        t_281[k] = pb_z[k] * ik_222[k];

        t_282[k] = f_3 * ii0_171[k]
                   - f_4 * ii1_171[k]
                   + pb_z[k] * ik_223[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_x, pb_y, pb_z, hk_117, hk_231, \
                         ii0_173, ii0_183, ii1_173, ii1_183, ik_225, ik_226, \
                         ik_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * hk_117[k]
                   + pb_y[k] * ik_225[k];

        t_284[k] = f_7 * ii0_173[k]
                   - f_8 * ii1_173[k]
                   + pb_z[k] * ik_225[k];

        t_285[k] = f_15 * hk_231[k]
                   + f_5 * ii0_183[k]
                   - f_6 * ii1_183[k]
                   + pb_x[k] * ik_231[k];

        t_286[k] = pb_z[k] * ik_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_y, pb_z, hk_122, ii0_174, ii0_175, \
                         ii0_177, ii1_174, ii1_175, ii1_177, ik_227, ik_228, \
                         ik_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * ii0_174[k]
                   - f_4 * ii1_174[k]
                   + pb_z[k] * ik_227[k];

        t_288[k] = f_5 * ii0_175[k]
                   - f_6 * ii1_175[k]
                   + pb_z[k] * ik_228[k];

        t_289[k] = f_15 * hk_122[k]
                   + pb_y[k] * ik_230[k];

        t_290[k] = f_9 * ii0_177[k]
                   - f_10 * ii1_177[k]
                   + pb_z[k] * ik_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_x, pb_z, hk_237, ii0_178, ii0_189, ii1_178, \
                         ii1_189, ik_231, ik_232, ik_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_15 * hk_237[k]
                   + f_3 * ii0_189[k]
                   - f_4 * ii1_189[k]
                   + pb_x[k] * ik_237[k];

        t_292[k] = pb_z[k] * ik_231[k];

        t_293[k] = f_3 * ii0_178[k]
                   - f_4 * ii1_178[k]
                   + pb_z[k] * ik_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pb_z, hk_128, ii0_179, ii0_180, \
                         ii0_182, ii1_179, ii1_180, ii1_182, ik_233, ik_234, \
                         ik_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_5 * ii0_179[k]
                   - f_6 * ii1_179[k]
                   + pb_z[k] * ik_233[k];

        t_295[k] = f_7 * ii0_180[k]
                   - f_8 * ii1_180[k]
                   + pb_z[k] * ik_234[k];

        t_296[k] = f_15 * hk_128[k]
                   + pb_y[k] * ik_236[k];

        t_297[k] = f_11 * ii0_182[k]
                   - f_12 * ii1_182[k]
                   + pb_z[k] * ik_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pb_x, pb_z, hk_244, hk_246, \
                         hk_247, hk_248, ik_237, ik_244, ik_246, ik_247, \
                         ik_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_15 * hk_244[k]
                   + pb_x[k] * ik_244[k];

        t_299[k] = pb_z[k] * ik_237[k];

        t_300[k] = f_15 * hk_246[k]
                   + pb_x[k] * ik_246[k];

        t_301[k] = f_15 * hk_247[k]
                   + pb_x[k] * ik_247[k];

        t_302[k] = f_15 * hk_248[k]
                   + pb_x[k] * ik_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pb_x, gl0_306, gl1_306, hk_249, \
                         hk_250, hk_251, hl_306, ik_249, ik_250, \
                         ik_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_15 * hk_249[k]
                   + pb_x[k] * ik_249[k];

        t_304[k] = f_15 * hk_250[k]
                   + pb_x[k] * ik_250[k];

        t_305[k] = f_15 * hk_251[k]
                   + pb_x[k] * ik_251[k];

        t_306[k] = f_23 * gl0_306[k]
                   - f_24 * gl1_306[k]
                   + pa_x[k] * hl_306[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pb_z, ii0_189, ii0_190, ii0_191, ii1_189, \
                         ii1_190, ii1_191, ik_244, ik_245, ik_246, \
                         ik_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = pb_z[k] * ik_244[k];

        t_308[k] = f_3 * ii0_189[k]
                   - f_4 * ii1_189[k]
                   + pb_z[k] * ik_245[k];

        t_309[k] = f_5 * ii0_190[k]
                   - f_6 * ii1_190[k]
                   + pb_z[k] * ik_246[k];

        t_310[k] = f_7 * ii0_191[k]
                   - f_8 * ii1_191[k]
                   + pb_z[k] * ik_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_y, pb_z, hk_143, ii0_192, ii0_193, \
                         ii0_195, ii1_192, ii1_193, ii1_195, ik_248, ik_249, \
                         ik_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_9 * ii0_192[k]
                   - f_10 * ii1_192[k]
                   + pb_z[k] * ik_248[k];

        t_312[k] = f_11 * ii0_193[k]
                   - f_12 * ii1_193[k]
                   + pb_z[k] * ik_249[k];

        t_313[k] = f_15 * hk_143[k]
                   + pb_y[k] * ik_251[k];

        t_314[k] = f_1 * ii0_195[k]
                   - f_2 * ii1_195[k]
                   + pb_z[k] * ik_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, pa_z, pb_y, pb_z, hk_108, hk_146, \
                         hl_135, hl_136, hl_138, ik_252, ik_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * hl_135[k];

        t_316[k] = pa_z[k] * hl_136[k];

        t_317[k] = f_13 * hk_108[k]
                   + pb_z[k] * ik_252[k];

        t_318[k] = pa_z[k] * hl_138[k];

        t_319[k] = f_14 * hk_146[k]
                   + pb_y[k] * ik_254[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_z, pb_y, pb_z, hk_110, hk_111, hk_149, \
                         hl_140, hl_141, ik_255, ik_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_14 * hk_110[k]
                   + pa_z[k] * hl_140[k];

        t_321[k] = pa_z[k] * hl_141[k];

        t_322[k] = f_13 * hk_111[k]
                   + pb_z[k] * ik_255[k];

        t_323[k] = f_14 * hk_149[k]
                   + pb_y[k] * ik_257[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_z, hk_113, hk_114, hk_115, \
                         hl_144, hl_145, hl_147, ik_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_15 * hk_113[k]
                   + pa_z[k] * hl_144[k];

        t_325[k] = pa_z[k] * hl_145[k];

        t_326[k] = f_13 * hk_114[k]
                   + pb_z[k] * ik_258[k];

        t_327[k] = f_14 * hk_115[k]
                   + pa_z[k] * hl_147[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_z, pb_y, pb_z, hk_117, hk_118, hk_153, \
                         hl_149, hl_150, ik_261, ik_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_14 * hk_153[k]
                   + pb_y[k] * ik_261[k];

        t_329[k] = f_16 * hk_117[k]
                   + pa_z[k] * hl_149[k];

        t_330[k] = pa_z[k] * hl_150[k];

        t_331[k] = f_13 * hk_118[k]
                   + pb_z[k] * ik_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pa_z, pb_y, hk_119, hk_120, \
                         hk_122, hk_158, hl_152, hl_153, hl_155, hl_156, \
                         ik_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_14 * hk_119[k]
                   + pa_z[k] * hl_152[k];

        t_333[k] = f_15 * hk_120[k]
                   + pa_z[k] * hl_153[k];

        t_334[k] = f_14 * hk_158[k]
                   + pb_y[k] * ik_266[k];

        t_335[k] = f_17 * hk_122[k]
                   + pa_z[k] * hl_155[k];

        t_336[k] = pa_z[k] * hl_156[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pb_z, hk_123, hk_124, hk_125, \
                         hk_126, hl_158, hl_159, hl_160, ik_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_13 * hk_123[k]
                   + pb_z[k] * ik_267[k];

        t_338[k] = f_14 * hk_124[k]
                   + pa_z[k] * hl_158[k];

        t_339[k] = f_15 * hk_125[k]
                   + pa_z[k] * hl_159[k];

        t_340[k] = f_16 * hk_126[k]
                   + pa_z[k] * hl_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_z, pb_x, pb_y, hk_128, hk_164, hk_281, \
                         hl_162, hl_163, ik_272, ik_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_14 * hk_164[k]
                   + pb_y[k] * ik_272[k];

        t_342[k] = f_0 * hk_128[k]
                   + pa_z[k] * hl_162[k];

        t_343[k] = pa_z[k] * hl_163[k];

        t_344[k] = f_15 * hk_281[k]
                   + pb_x[k] * ik_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pb_x, hk_282, hk_283, hk_284, \
                         hk_285, hk_286, ik_282, ik_283, ik_284, ik_285, \
                         ik_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_15 * hk_282[k]
                   + pb_x[k] * ik_282[k];

        t_346[k] = f_15 * hk_283[k]
                   + pb_x[k] * ik_283[k];

        t_347[k] = f_15 * hk_284[k]
                   + pb_x[k] * ik_284[k];

        t_348[k] = f_15 * hk_285[k]
                   + pb_x[k] * ik_285[k];

        t_349[k] = f_15 * hk_286[k]
                   + pb_x[k] * ik_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_z, pb_x, pb_z, hk_136, hk_137, hk_287, \
                         hl_171, hl_173, ik_280, ik_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_15 * hk_287[k]
                   + pb_x[k] * ik_287[k];

        t_351[k] = pa_z[k] * hl_171[k];

        t_352[k] = f_13 * hk_136[k]
                   + pb_z[k] * ik_280[k];

        t_353[k] = f_14 * hk_137[k]
                   + pa_z[k] * hl_173[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_z, hk_138, hk_139, hk_140, hk_141, \
                         hl_174, hl_175, hl_176, hl_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_15 * hk_138[k]
                   + pa_z[k] * hl_174[k];

        t_355[k] = f_16 * hk_139[k]
                   + pa_z[k] * hl_175[k];

        t_356[k] = f_17 * hk_140[k]
                   + pa_z[k] * hl_176[k];

        t_357[k] = f_0 * hk_141[k]
                   + pa_z[k] * hl_177[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, pa_y, pa_z, pb_y, hk_143, hk_179, \
                         hk_180, hl_179, hl_225, hl_227, ik_287, \
                         ik_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_14 * hk_179[k]
                   + pb_y[k] * ik_287[k];

        t_359[k] = f_18 * hk_143[k]
                   + pa_z[k] * hl_179[k];

        t_360[k] = pa_y[k] * hl_225[k];

        t_361[k] = f_13 * hk_180[k]
                   + pb_y[k] * ik_288[k];

        t_362[k] = pa_y[k] * hl_227[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pb_y, hk_181, hk_182, hk_183, \
                         hl_228, hl_230, hl_231, ik_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * hk_181[k]
                   + pa_y[k] * hl_228[k];

        t_364[k] = f_13 * hk_182[k]
                   + pb_y[k] * ik_290[k];

        t_365[k] = pa_y[k] * hl_230[k];

        t_366[k] = f_15 * hk_183[k]
                   + pa_y[k] * hl_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pb_y, pb_z, hk_147, hk_185, hk_186, \
                         hl_234, hl_235, ik_291, ik_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_14 * hk_147[k]
                   + pb_z[k] * ik_291[k];

        t_368[k] = f_13 * hk_185[k]
                   + pb_y[k] * ik_293[k];

        t_369[k] = pa_y[k] * hl_234[k];

        t_370[k] = f_16 * hk_186[k]
                   + pa_y[k] * hl_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pb_y, pb_z, hk_150, hk_188, hk_189, \
                         hl_237, hl_239, ik_294, ik_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_14 * hk_150[k]
                   + pb_z[k] * ik_294[k];

        t_372[k] = f_14 * hk_188[k]
                   + pa_y[k] * hl_237[k];

        t_373[k] = f_13 * hk_189[k]
                   + pb_y[k] * ik_297[k];

        t_374[k] = pa_y[k] * hl_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pa_y, pb_z, hk_154, hk_190, hk_192, \
                         hk_193, hl_240, hl_242, hl_243, ik_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_17 * hk_190[k]
                   + pa_y[k] * hl_240[k];

        t_376[k] = f_14 * hk_154[k]
                   + pb_z[k] * ik_298[k];

        t_377[k] = f_15 * hk_192[k]
                   + pa_y[k] * hl_242[k];

        t_378[k] = f_14 * hk_193[k]
                   + pa_y[k] * hl_243[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pa_y, pb_y, pb_z, hk_159, hk_194, hk_195, \
                         hl_245, hl_246, ik_302, ik_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_13 * hk_194[k]
                   + pb_y[k] * ik_302[k];

        t_380[k] = pa_y[k] * hl_245[k];

        t_381[k] = f_0 * hk_195[k]
                   + pa_y[k] * hl_246[k];

        t_382[k] = f_14 * hk_159[k]
                   + pb_z[k] * ik_303[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pa_y, pb_y, hk_197, hk_198, \
                         hk_199, hk_200, hl_248, hl_249, hl_250, hl_252, \
                         ik_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_16 * hk_197[k]
                   + pa_y[k] * hl_248[k];

        t_384[k] = f_15 * hk_198[k]
                   + pa_y[k] * hl_249[k];

        t_385[k] = f_14 * hk_199[k]
                   + pa_y[k] * hl_250[k];

        t_386[k] = f_13 * hk_200[k]
                   + pb_y[k] * ik_308[k];

        t_387[k] = pa_y[k] * hl_252[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, pb_x, hk_316, hk_317, hk_318, \
                         hk_319, hk_320, ik_316, ik_317, ik_318, ik_319, \
                         ik_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_15 * hk_316[k]
                   + pb_x[k] * ik_316[k];

        t_389[k] = f_15 * hk_317[k]
                   + pb_x[k] * ik_317[k];

        t_390[k] = f_15 * hk_318[k]
                   + pb_x[k] * ik_318[k];

        t_391[k] = f_15 * hk_319[k]
                   + pb_x[k] * ik_319[k];

        t_392[k] = f_15 * hk_320[k]
                   + pb_x[k] * ik_320[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_y, pb_x, hk_208, hk_321, hk_322, \
                         hl_260, hl_261, ik_321, ik_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_15 * hk_321[k]
                   + pb_x[k] * ik_321[k];

        t_394[k] = f_15 * hk_322[k]
                   + pb_x[k] * ik_322[k];

        t_395[k] = pa_y[k] * hl_260[k];

        t_396[k] = f_18 * hk_208[k]
                   + pa_y[k] * hl_261[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pb_z, hk_172, hk_210, hk_211, \
                         hk_212, hl_263, hl_264, hl_265, ik_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_14 * hk_172[k]
                   + pb_z[k] * ik_316[k];

        t_398[k] = f_0 * hk_210[k]
                   + pa_y[k] * hl_263[k];

        t_399[k] = f_17 * hk_211[k]
                   + pa_y[k] * hl_264[k];

        t_400[k] = f_16 * hk_212[k]
                   + pa_y[k] * hl_265[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pb_y, hk_213, hk_214, hk_215, \
                         hl_266, hl_267, hl_269, ik_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * hk_213[k]
                   + pa_y[k] * hl_266[k];

        t_402[k] = f_14 * hk_214[k]
                   + pa_y[k] * hl_267[k];

        t_403[k] = f_13 * hk_215[k]
                   + pb_y[k] * ik_323[k];

        t_404[k] = pa_y[k] * hl_269[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_z, pb_y, pb_z, gl0_90, gl1_90, hk_180, \
                         hl_225, ii0_252, ii1_252, ik_324, ik_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_23 * gl0_90[k]
                   - f_24 * gl1_90[k]
                   + pa_z[k] * hl_225[k];

        t_406[k] = pb_y[k] * ik_324[k];

        t_407[k] = f_15 * hk_180[k]
                   + pb_z[k] * ik_324[k];

        t_408[k] = f_3 * ii0_252[k]
                   - f_4 * ii1_252[k]
                   + pb_y[k] * ik_325[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, pb_x, pb_y, pb_z, hk_183, hk_329, \
                         ii0_253, ii0_257, ii1_253, ii1_257, ik_326, ik_327, \
                         ik_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * ik_326[k];

        t_410[k] = f_15 * hk_329[k]
                   + f_11 * ii0_257[k]
                   - f_12 * ii1_257[k]
                   + pb_x[k] * ik_329[k];

        t_411[k] = f_5 * ii0_253[k]
                   - f_6 * ii1_253[k]
                   + pb_y[k] * ik_327[k];

        t_412[k] = f_15 * hk_183[k]
                   + pb_z[k] * ik_327[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);
    const auto f_23 = 1.0 / alpha;
    const auto f_24 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_135 = buffer.data(gl0 + 135);
    const auto *gl0_225 = buffer.data(gl0 + 225);
    const auto *gl0_449 = buffer.data(gl0 + 449);
    const auto *gl0_486 = buffer.data(gl0 + 486);

    const auto *gl1_135 = buffer.data(gl1 + 135);
    const auto *gl1_225 = buffer.data(gl1 + 225);
    const auto *gl1_449 = buffer.data(gl1 + 449);
    const auto *gl1_486 = buffer.data(gl1 + 486);

    const auto *hk_186 = buffer.data(hk + 186);
    const auto *hk_190 = buffer.data(hk + 190);
    const auto *hk_195 = buffer.data(hk + 195);
    const auto *hk_208 = buffer.data(hk + 208);
    const auto *hk_216 = buffer.data(hk + 216);
    const auto *hk_218 = buffer.data(hk + 218);
    const auto *hk_219 = buffer.data(hk + 219);
    const auto *hk_221 = buffer.data(hk + 221);
    const auto *hk_222 = buffer.data(hk + 222);
    const auto *hk_223 = buffer.data(hk + 223);
    const auto *hk_225 = buffer.data(hk + 225);
    const auto *hk_226 = buffer.data(hk + 226);
    const auto *hk_227 = buffer.data(hk + 227);
    const auto *hk_228 = buffer.data(hk + 228);
    const auto *hk_230 = buffer.data(hk + 230);
    const auto *hk_231 = buffer.data(hk + 231);
    const auto *hk_232 = buffer.data(hk + 232);
    const auto *hk_233 = buffer.data(hk + 233);
    const auto *hk_234 = buffer.data(hk + 234);
    const auto *hk_236 = buffer.data(hk + 236);
    const auto *hk_244 = buffer.data(hk + 244);
    const auto *hk_245 = buffer.data(hk + 245);
    const auto *hk_246 = buffer.data(hk + 246);
    const auto *hk_247 = buffer.data(hk + 247);
    const auto *hk_248 = buffer.data(hk + 248);
    const auto *hk_249 = buffer.data(hk + 249);
    const auto *hk_251 = buffer.data(hk + 251);
    const auto *hk_254 = buffer.data(hk + 254);
    const auto *hk_257 = buffer.data(hk + 257);
    const auto *hk_261 = buffer.data(hk + 261);
    const auto *hk_266 = buffer.data(hk + 266);
    const auto *hk_272 = buffer.data(hk + 272);
    const auto *hk_287 = buffer.data(hk + 287);
    const auto *hk_288 = buffer.data(hk + 288);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_351 = buffer.data(hk + 351);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_353 = buffer.data(hk + 353);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_381 = buffer.data(hk + 381);
    const auto *hk_388 = buffer.data(hk + 388);
    const auto *hk_390 = buffer.data(hk + 390);
    const auto *hk_391 = buffer.data(hk + 391);
    const auto *hk_392 = buffer.data(hk + 392);
    const auto *hk_393 = buffer.data(hk + 393);
    const auto *hk_394 = buffer.data(hk + 394);
    const auto *hk_395 = buffer.data(hk + 395);
    const auto *hk_425 = buffer.data(hk + 425);
    const auto *hk_426 = buffer.data(hk + 426);
    const auto *hk_427 = buffer.data(hk + 427);
    const auto *hk_428 = buffer.data(hk + 428);
    const auto *hk_429 = buffer.data(hk + 429);
    const auto *hk_430 = buffer.data(hk + 430);
    const auto *hk_431 = buffer.data(hk + 431);

    const auto *hl_270 = buffer.data(hl + 270);
    const auto *hl_271 = buffer.data(hl + 271);
    const auto *hl_273 = buffer.data(hl + 273);
    const auto *hl_275 = buffer.data(hl + 275);
    const auto *hl_276 = buffer.data(hl + 276);
    const auto *hl_279 = buffer.data(hl + 279);
    const auto *hl_280 = buffer.data(hl + 280);
    const auto *hl_282 = buffer.data(hl + 282);
    const auto *hl_284 = buffer.data(hl + 284);
    const auto *hl_285 = buffer.data(hl + 285);
    const auto *hl_287 = buffer.data(hl + 287);
    const auto *hl_288 = buffer.data(hl + 288);
    const auto *hl_290 = buffer.data(hl + 290);
    const auto *hl_291 = buffer.data(hl + 291);
    const auto *hl_293 = buffer.data(hl + 293);
    const auto *hl_294 = buffer.data(hl + 294);
    const auto *hl_295 = buffer.data(hl + 295);
    const auto *hl_297 = buffer.data(hl + 297);
    const auto *hl_298 = buffer.data(hl + 298);
    const auto *hl_306 = buffer.data(hl + 306);
    const auto *hl_308 = buffer.data(hl + 308);
    const auto *hl_309 = buffer.data(hl + 309);
    const auto *hl_310 = buffer.data(hl + 310);
    const auto *hl_311 = buffer.data(hl + 311);
    const auto *hl_312 = buffer.data(hl + 312);
    const auto *hl_314 = buffer.data(hl + 314);
    const auto *hl_360 = buffer.data(hl + 360);
    const auto *hl_449 = buffer.data(hl + 449);
    const auto *hl_486 = buffer.data(hl + 486);

    const auto *ii0_255 = buffer.data(ii0 + 255);
    const auto *ii0_257 = buffer.data(ii0 + 257);
    const auto *ii0_258 = buffer.data(ii0 + 258);
    const auto *ii0_260 = buffer.data(ii0 + 260);
    const auto *ii0_261 = buffer.data(ii0 + 261);
    const auto *ii0_262 = buffer.data(ii0 + 262);
    const auto *ii0_264 = buffer.data(ii0 + 264);
    const auto *ii0_265 = buffer.data(ii0 + 265);
    const auto *ii0_266 = buffer.data(ii0 + 266);
    const auto *ii0_272 = buffer.data(ii0 + 272);
    const auto *ii0_273 = buffer.data(ii0 + 273);
    const auto *ii0_275 = buffer.data(ii0 + 275);
    const auto *ii0_276 = buffer.data(ii0 + 276);
    const auto *ii0_277 = buffer.data(ii0 + 277);
    const auto *ii0_278 = buffer.data(ii0 + 278);
    const auto *ii0_279 = buffer.data(ii0 + 279);
    const auto *ii0_280 = buffer.data(ii0 + 280);
    const auto *ii0_282 = buffer.data(ii0 + 282);
    const auto *ii0_283 = buffer.data(ii0 + 283);
    const auto *ii0_285 = buffer.data(ii0 + 285);
    const auto *ii0_286 = buffer.data(ii0 + 286);
    const auto *ii0_287 = buffer.data(ii0 + 287);
    const auto *ii0_289 = buffer.data(ii0 + 289);
    const auto *ii0_290 = buffer.data(ii0 + 290);
    const auto *ii0_291 = buffer.data(ii0 + 291);
    const auto *ii0_292 = buffer.data(ii0 + 292);
    const auto *ii0_294 = buffer.data(ii0 + 294);
    const auto *ii0_295 = buffer.data(ii0 + 295);
    const auto *ii0_301 = buffer.data(ii0 + 301);
    const auto *ii0_302 = buffer.data(ii0 + 302);
    const auto *ii0_303 = buffer.data(ii0 + 303);
    const auto *ii0_304 = buffer.data(ii0 + 304);
    const auto *ii0_305 = buffer.data(ii0 + 305);
    const auto *ii0_307 = buffer.data(ii0 + 307);

    const auto *ii1_255 = buffer.data(ii1 + 255);
    const auto *ii1_257 = buffer.data(ii1 + 257);
    const auto *ii1_258 = buffer.data(ii1 + 258);
    const auto *ii1_260 = buffer.data(ii1 + 260);
    const auto *ii1_261 = buffer.data(ii1 + 261);
    const auto *ii1_262 = buffer.data(ii1 + 262);
    const auto *ii1_264 = buffer.data(ii1 + 264);
    const auto *ii1_265 = buffer.data(ii1 + 265);
    const auto *ii1_266 = buffer.data(ii1 + 266);
    const auto *ii1_272 = buffer.data(ii1 + 272);
    const auto *ii1_273 = buffer.data(ii1 + 273);
    const auto *ii1_275 = buffer.data(ii1 + 275);
    const auto *ii1_276 = buffer.data(ii1 + 276);
    const auto *ii1_277 = buffer.data(ii1 + 277);
    const auto *ii1_278 = buffer.data(ii1 + 278);
    const auto *ii1_279 = buffer.data(ii1 + 279);
    const auto *ii1_280 = buffer.data(ii1 + 280);
    const auto *ii1_282 = buffer.data(ii1 + 282);
    const auto *ii1_283 = buffer.data(ii1 + 283);
    const auto *ii1_285 = buffer.data(ii1 + 285);
    const auto *ii1_286 = buffer.data(ii1 + 286);
    const auto *ii1_287 = buffer.data(ii1 + 287);
    const auto *ii1_289 = buffer.data(ii1 + 289);
    const auto *ii1_290 = buffer.data(ii1 + 290);
    const auto *ii1_291 = buffer.data(ii1 + 291);
    const auto *ii1_292 = buffer.data(ii1 + 292);
    const auto *ii1_294 = buffer.data(ii1 + 294);
    const auto *ii1_295 = buffer.data(ii1 + 295);
    const auto *ii1_301 = buffer.data(ii1 + 301);
    const auto *ii1_302 = buffer.data(ii1 + 302);
    const auto *ii1_303 = buffer.data(ii1 + 303);
    const auto *ii1_304 = buffer.data(ii1 + 304);
    const auto *ii1_305 = buffer.data(ii1 + 305);
    const auto *ii1_307 = buffer.data(ii1 + 307);

    const auto *ik_329 = buffer.data(ik + 329);
    const auto *ik_330 = buffer.data(ik + 330);
    const auto *ik_332 = buffer.data(ik + 332);
    const auto *ik_333 = buffer.data(ik + 333);
    const auto *ik_334 = buffer.data(ik + 334);
    const auto *ik_336 = buffer.data(ik + 336);
    const auto *ik_337 = buffer.data(ik + 337);
    const auto *ik_338 = buffer.data(ik + 338);
    const auto *ik_339 = buffer.data(ik + 339);
    const auto *ik_341 = buffer.data(ik + 341);
    const auto *ik_342 = buffer.data(ik + 342);
    const auto *ik_343 = buffer.data(ik + 343);
    const auto *ik_344 = buffer.data(ik + 344);
    const auto *ik_351 = buffer.data(ik + 351);
    const auto *ik_352 = buffer.data(ik + 352);
    const auto *ik_353 = buffer.data(ik + 353);
    const auto *ik_354 = buffer.data(ik + 354);
    const auto *ik_355 = buffer.data(ik + 355);
    const auto *ik_356 = buffer.data(ik + 356);
    const auto *ik_357 = buffer.data(ik + 357);
    const auto *ik_358 = buffer.data(ik + 358);
    const auto *ik_359 = buffer.data(ik + 359);
    const auto *ik_360 = buffer.data(ik + 360);
    const auto *ik_361 = buffer.data(ik + 361);
    const auto *ik_362 = buffer.data(ik + 362);
    const auto *ik_363 = buffer.data(ik + 363);
    const auto *ik_365 = buffer.data(ik + 365);
    const auto *ik_366 = buffer.data(ik + 366);
    const auto *ik_367 = buffer.data(ik + 367);
    const auto *ik_369 = buffer.data(ik + 369);
    const auto *ik_370 = buffer.data(ik + 370);
    const auto *ik_371 = buffer.data(ik + 371);
    const auto *ik_372 = buffer.data(ik + 372);
    const auto *ik_374 = buffer.data(ik + 374);
    const auto *ik_375 = buffer.data(ik + 375);
    const auto *ik_376 = buffer.data(ik + 376);
    const auto *ik_377 = buffer.data(ik + 377);
    const auto *ik_378 = buffer.data(ik + 378);
    const auto *ik_380 = buffer.data(ik + 380);
    const auto *ik_381 = buffer.data(ik + 381);
    const auto *ik_388 = buffer.data(ik + 388);
    const auto *ik_389 = buffer.data(ik + 389);
    const auto *ik_390 = buffer.data(ik + 390);
    const auto *ik_391 = buffer.data(ik + 391);
    const auto *ik_392 = buffer.data(ik + 392);
    const auto *ik_393 = buffer.data(ik + 393);
    const auto *ik_394 = buffer.data(ik + 394);
    const auto *ik_395 = buffer.data(ik + 395);
    const auto *ik_396 = buffer.data(ik + 396);
    const auto *ik_398 = buffer.data(ik + 398);
    const auto *ik_399 = buffer.data(ik + 399);
    const auto *ik_401 = buffer.data(ik + 401);
    const auto *ik_402 = buffer.data(ik + 402);
    const auto *ik_405 = buffer.data(ik + 405);
    const auto *ik_406 = buffer.data(ik + 406);
    const auto *ik_410 = buffer.data(ik + 410);
    const auto *ik_411 = buffer.data(ik + 411);
    const auto *ik_416 = buffer.data(ik + 416);
    const auto *ik_424 = buffer.data(ik + 424);
    const auto *ik_425 = buffer.data(ik + 425);
    const auto *ik_426 = buffer.data(ik + 426);
    const auto *ik_427 = buffer.data(ik + 427);
    const auto *ik_428 = buffer.data(ik + 428);
    const auto *ik_429 = buffer.data(ik + 429);
    const auto *ik_430 = buffer.data(ik + 430);
    const auto *ik_431 = buffer.data(ik + 431);
    const auto *ik_432 = buffer.data(ik + 432);

#pragma omp simd aligned(t_413, t_414, t_415, t_416, pb_x, pb_y, pb_z, hk_186, hk_333, \
                         ii0_255, ii0_261, ii1_255, ii1_261, ik_329, ik_330, \
                         ik_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_y[k] * ik_329[k];

        t_414[k] = f_15 * hk_333[k]
                   + f_9 * ii0_261[k]
                   - f_10 * ii1_261[k]
                   + pb_x[k] * ik_333[k];

        t_415[k] = f_7 * ii0_255[k]
                   - f_8 * ii1_255[k]
                   + pb_y[k] * ik_330[k];

        t_416[k] = f_15 * hk_186[k]
                   + pb_z[k] * ik_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_x, pb_y, hk_338, ii0_257, ii0_266, ii1_257, \
                         ii1_266, ik_332, ik_333, ik_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_3 * ii0_257[k]
                   - f_4 * ii1_257[k]
                   + pb_y[k] * ik_332[k];

        t_418[k] = pb_y[k] * ik_333[k];

        t_419[k] = f_15 * hk_338[k]
                   + f_7 * ii0_266[k]
                   - f_8 * ii1_266[k]
                   + pb_x[k] * ik_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_y, pb_z, hk_190, ii0_258, ii0_260, \
                         ii0_261, ii1_258, ii1_260, ii1_261, ik_334, ik_336, \
                         ik_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_9 * ii0_258[k]
                   - f_10 * ii1_258[k]
                   + pb_y[k] * ik_334[k];

        t_421[k] = f_15 * hk_190[k]
                   + pb_z[k] * ik_334[k];

        t_422[k] = f_5 * ii0_260[k]
                   - f_6 * ii1_260[k]
                   + pb_y[k] * ik_336[k];

        t_423[k] = f_3 * ii0_261[k]
                   - f_4 * ii1_261[k]
                   + pb_y[k] * ik_337[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pb_x, pb_y, pb_z, hk_195, hk_344, \
                         ii0_262, ii0_272, ii1_262, ii1_272, ik_338, ik_339, \
                         ik_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_y[k] * ik_338[k];

        t_425[k] = f_15 * hk_344[k]
                   + f_5 * ii0_272[k]
                   - f_6 * ii1_272[k]
                   + pb_x[k] * ik_344[k];

        t_426[k] = f_11 * ii0_262[k]
                   - f_12 * ii1_262[k]
                   + pb_y[k] * ik_339[k];

        t_427[k] = f_15 * hk_195[k]
                   + pb_z[k] * ik_339[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pb_y, ii0_264, ii0_265, ii0_266, ii1_264, \
                         ii1_265, ii1_266, ik_341, ik_342, ik_343, \
                         ik_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_7 * ii0_264[k]
                   - f_8 * ii1_264[k]
                   + pb_y[k] * ik_341[k];

        t_429[k] = f_5 * ii0_265[k]
                   - f_6 * ii1_265[k]
                   + pb_y[k] * ik_342[k];

        t_430[k] = f_3 * ii0_266[k]
                   - f_4 * ii1_266[k]
                   + pb_y[k] * ik_343[k];

        t_431[k] = pb_y[k] * ik_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_x, hk_351, hk_352, hk_353, hk_354, \
                         ii0_279, ii1_279, ik_351, ik_352, ik_353, \
                         ik_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_15 * hk_351[k]
                   + f_3 * ii0_279[k]
                   - f_4 * ii1_279[k]
                   + pb_x[k] * ik_351[k];

        t_433[k] = f_15 * hk_352[k]
                   + pb_x[k] * ik_352[k];

        t_434[k] = f_15 * hk_353[k]
                   + pb_x[k] * ik_353[k];

        t_435[k] = f_15 * hk_354[k]
                   + pb_x[k] * ik_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pb_x, pb_y, hk_355, hk_356, \
                         hk_357, hk_359, ik_351, ik_355, ik_356, ik_357, \
                         ik_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_15 * hk_355[k]
                   + pb_x[k] * ik_355[k];

        t_437[k] = f_15 * hk_356[k]
                   + pb_x[k] * ik_356[k];

        t_438[k] = f_15 * hk_357[k]
                   + pb_x[k] * ik_357[k];

        t_439[k] = pb_y[k] * ik_351[k];

        t_440[k] = f_15 * hk_359[k]
                   + pb_x[k] * ik_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_y, pb_z, hk_208, ii0_273, ii0_275, \
                         ii0_276, ii1_273, ii1_275, ii1_276, ik_352, ik_354, \
                         ik_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * ii0_273[k]
                   - f_2 * ii1_273[k]
                   + pb_y[k] * ik_352[k];

        t_442[k] = f_15 * hk_208[k]
                   + pb_z[k] * ik_352[k];

        t_443[k] = f_11 * ii0_275[k]
                   - f_12 * ii1_275[k]
                   + pb_y[k] * ik_354[k];

        t_444[k] = f_9 * ii0_276[k]
                   - f_10 * ii1_276[k]
                   + pb_y[k] * ik_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pb_y, ii0_277, ii0_278, ii0_279, ii1_277, \
                         ii1_278, ii1_279, ik_356, ik_357, ik_358, \
                         ik_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_7 * ii0_277[k]
                   - f_8 * ii1_277[k]
                   + pb_y[k] * ik_356[k];

        t_446[k] = f_5 * ii0_278[k]
                   - f_6 * ii1_278[k]
                   + pb_y[k] * ik_357[k];

        t_447[k] = f_3 * ii0_279[k]
                   - f_4 * ii1_279[k]
                   + pb_y[k] * ik_358[k];

        t_448[k] = pb_y[k] * ik_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_x, pa_y, pb_y, pb_z, gl0_135, gl0_449, \
                         gl1_135, gl1_449, hk_216, hl_270, hl_449, \
                         ik_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_23 * gl0_449[k]
                   - f_24 * gl1_449[k]
                   + pa_x[k] * hl_449[k];

        t_450[k] = f_21 * gl0_135[k]
                   - f_22 * gl1_135[k]
                   + pa_y[k] * hl_270[k];

        t_451[k] = f_16 * hk_216[k]
                   + pb_y[k] * ik_360[k];

        t_452[k] = pb_z[k] * ik_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pb_x, pb_z, hk_363, ii0_280, ii0_283, ii1_280, \
                         ii1_283, ik_361, ik_362, ik_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_14 * hk_363[k]
                   + f_11 * ii0_283[k]
                   - f_12 * ii1_283[k]
                   + pb_x[k] * ik_363[k];

        t_454[k] = pb_z[k] * ik_361[k];

        t_455[k] = f_3 * ii0_280[k]
                   - f_4 * ii1_280[k]
                   + pb_z[k] * ik_362[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_x, pb_y, pb_z, hk_221, hk_366, \
                         ii0_282, ii0_286, ii1_282, ii1_286, ik_363, ik_365, \
                         ik_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * hk_366[k]
                   + f_9 * ii0_286[k]
                   - f_10 * ii1_286[k]
                   + pb_x[k] * ik_366[k];

        t_457[k] = pb_z[k] * ik_363[k];

        t_458[k] = f_16 * hk_221[k]
                   + pb_y[k] * ik_365[k];

        t_459[k] = f_5 * ii0_282[k]
                   - f_6 * ii1_282[k]
                   + pb_z[k] * ik_365[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pb_x, pb_z, hk_370, ii0_283, ii0_290, ii1_283, \
                         ii1_290, ik_366, ik_367, ik_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * hk_370[k]
                   + f_7 * ii0_290[k]
                   - f_8 * ii1_290[k]
                   + pb_x[k] * ik_370[k];

        t_461[k] = pb_z[k] * ik_366[k];

        t_462[k] = f_3 * ii0_283[k]
                   - f_4 * ii1_283[k]
                   + pb_z[k] * ik_367[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pb_x, pb_y, pb_z, hk_225, hk_375, \
                         ii0_285, ii0_295, ii1_285, ii1_295, ik_369, ik_370, \
                         ik_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_16 * hk_225[k]
                   + pb_y[k] * ik_369[k];

        t_464[k] = f_7 * ii0_285[k]
                   - f_8 * ii1_285[k]
                   + pb_z[k] * ik_369[k];

        t_465[k] = f_14 * hk_375[k]
                   + f_5 * ii0_295[k]
                   - f_6 * ii1_295[k]
                   + pb_x[k] * ik_375[k];

        t_466[k] = pb_z[k] * ik_370[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, pb_y, pb_z, hk_230, ii0_286, ii0_287, \
                         ii0_289, ii1_286, ii1_287, ii1_289, ik_371, ik_372, \
                         ik_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_3 * ii0_286[k]
                   - f_4 * ii1_286[k]
                   + pb_z[k] * ik_371[k];

        t_468[k] = f_5 * ii0_287[k]
                   - f_6 * ii1_287[k]
                   + pb_z[k] * ik_372[k];

        t_469[k] = f_16 * hk_230[k]
                   + pb_y[k] * ik_374[k];

        t_470[k] = f_9 * ii0_289[k]
                   - f_10 * ii1_289[k]
                   + pb_z[k] * ik_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pb_x, pb_z, hk_381, ii0_290, ii0_301, ii1_290, \
                         ii1_301, ik_375, ik_376, ik_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_14 * hk_381[k]
                   + f_3 * ii0_301[k]
                   - f_4 * ii1_301[k]
                   + pb_x[k] * ik_381[k];

        t_472[k] = pb_z[k] * ik_375[k];

        t_473[k] = f_3 * ii0_290[k]
                   - f_4 * ii1_290[k]
                   + pb_z[k] * ik_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pb_z, hk_236, ii0_291, ii0_292, \
                         ii0_294, ii1_291, ii1_292, ii1_294, ik_377, ik_378, \
                         ik_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_5 * ii0_291[k]
                   - f_6 * ii1_291[k]
                   + pb_z[k] * ik_377[k];

        t_475[k] = f_7 * ii0_292[k]
                   - f_8 * ii1_292[k]
                   + pb_z[k] * ik_378[k];

        t_476[k] = f_16 * hk_236[k]
                   + pb_y[k] * ik_380[k];

        t_477[k] = f_11 * ii0_294[k]
                   - f_12 * ii1_294[k]
                   + pb_z[k] * ik_380[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pb_x, pb_z, hk_388, hk_390, \
                         hk_391, hk_392, ik_381, ik_388, ik_390, ik_391, \
                         ik_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_14 * hk_388[k]
                   + pb_x[k] * ik_388[k];

        t_479[k] = pb_z[k] * ik_381[k];

        t_480[k] = f_14 * hk_390[k]
                   + pb_x[k] * ik_390[k];

        t_481[k] = f_14 * hk_391[k]
                   + pb_x[k] * ik_391[k];

        t_482[k] = f_14 * hk_392[k]
                   + pb_x[k] * ik_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_x, pb_x, gl0_486, gl1_486, hk_393, \
                         hk_394, hk_395, hl_486, ik_393, ik_394, \
                         ik_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_14 * hk_393[k]
                   + pb_x[k] * ik_393[k];

        t_484[k] = f_14 * hk_394[k]
                   + pb_x[k] * ik_394[k];

        t_485[k] = f_14 * hk_395[k]
                   + pb_x[k] * ik_395[k];

        t_486[k] = f_19 * gl0_486[k]
                   - f_20 * gl1_486[k]
                   + pa_x[k] * hl_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pb_z, ii0_301, ii0_302, ii0_303, ii1_301, \
                         ii1_302, ii1_303, ik_388, ik_389, ik_390, \
                         ik_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = pb_z[k] * ik_388[k];

        t_488[k] = f_3 * ii0_301[k]
                   - f_4 * ii1_301[k]
                   + pb_z[k] * ik_389[k];

        t_489[k] = f_5 * ii0_302[k]
                   - f_6 * ii1_302[k]
                   + pb_z[k] * ik_390[k];

        t_490[k] = f_7 * ii0_303[k]
                   - f_8 * ii1_303[k]
                   + pb_z[k] * ik_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pb_y, pb_z, hk_251, ii0_304, ii0_305, \
                         ii0_307, ii1_304, ii1_305, ii1_307, ik_392, ik_393, \
                         ik_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_9 * ii0_304[k]
                   - f_10 * ii1_304[k]
                   + pb_z[k] * ik_392[k];

        t_492[k] = f_11 * ii0_305[k]
                   - f_12 * ii1_305[k]
                   + pb_z[k] * ik_393[k];

        t_493[k] = f_16 * hk_251[k]
                   + pb_y[k] * ik_395[k];

        t_494[k] = f_1 * ii0_307[k]
                   - f_2 * ii1_307[k]
                   + pb_z[k] * ik_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, pa_z, pb_y, pb_z, hk_216, hk_254, \
                         hl_270, hl_271, hl_273, ik_396, ik_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * hl_270[k];

        t_496[k] = pa_z[k] * hl_271[k];

        t_497[k] = f_13 * hk_216[k]
                   + pb_z[k] * ik_396[k];

        t_498[k] = pa_z[k] * hl_273[k];

        t_499[k] = f_15 * hk_254[k]
                   + pb_y[k] * ik_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_z, pb_y, pb_z, hk_218, hk_219, hk_257, \
                         hl_275, hl_276, ik_399, ik_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_14 * hk_218[k]
                   + pa_z[k] * hl_275[k];

        t_501[k] = pa_z[k] * hl_276[k];

        t_502[k] = f_13 * hk_219[k]
                   + pb_z[k] * ik_399[k];

        t_503[k] = f_15 * hk_257[k]
                   + pb_y[k] * ik_401[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_z, pb_z, hk_221, hk_222, hk_223, \
                         hl_279, hl_280, hl_282, ik_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_15 * hk_221[k]
                   + pa_z[k] * hl_279[k];

        t_505[k] = pa_z[k] * hl_280[k];

        t_506[k] = f_13 * hk_222[k]
                   + pb_z[k] * ik_402[k];

        t_507[k] = f_14 * hk_223[k]
                   + pa_z[k] * hl_282[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_z, pb_y, pb_z, hk_225, hk_226, hk_261, \
                         hl_284, hl_285, ik_405, ik_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_15 * hk_261[k]
                   + pb_y[k] * ik_405[k];

        t_509[k] = f_16 * hk_225[k]
                   + pa_z[k] * hl_284[k];

        t_510[k] = pa_z[k] * hl_285[k];

        t_511[k] = f_13 * hk_226[k]
                   + pb_z[k] * ik_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pa_z, pb_y, hk_227, hk_228, \
                         hk_230, hk_266, hl_287, hl_288, hl_290, hl_291, \
                         ik_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_14 * hk_227[k]
                   + pa_z[k] * hl_287[k];

        t_513[k] = f_15 * hk_228[k]
                   + pa_z[k] * hl_288[k];

        t_514[k] = f_15 * hk_266[k]
                   + pb_y[k] * ik_410[k];

        t_515[k] = f_17 * hk_230[k]
                   + pa_z[k] * hl_290[k];

        t_516[k] = pa_z[k] * hl_291[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_z, pb_z, hk_231, hk_232, hk_233, \
                         hk_234, hl_293, hl_294, hl_295, ik_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_13 * hk_231[k]
                   + pb_z[k] * ik_411[k];

        t_518[k] = f_14 * hk_232[k]
                   + pa_z[k] * hl_293[k];

        t_519[k] = f_15 * hk_233[k]
                   + pa_z[k] * hl_294[k];

        t_520[k] = f_16 * hk_234[k]
                   + pa_z[k] * hl_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pa_z, pb_x, pb_y, hk_236, hk_272, hk_425, \
                         hl_297, hl_298, ik_416, ik_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_15 * hk_272[k]
                   + pb_y[k] * ik_416[k];

        t_522[k] = f_0 * hk_236[k]
                   + pa_z[k] * hl_297[k];

        t_523[k] = pa_z[k] * hl_298[k];

        t_524[k] = f_14 * hk_425[k]
                   + pb_x[k] * ik_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pb_x, hk_426, hk_427, hk_428, \
                         hk_429, hk_430, ik_426, ik_427, ik_428, ik_429, \
                         ik_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_14 * hk_426[k]
                   + pb_x[k] * ik_426[k];

        t_526[k] = f_14 * hk_427[k]
                   + pb_x[k] * ik_427[k];

        t_527[k] = f_14 * hk_428[k]
                   + pb_x[k] * ik_428[k];

        t_528[k] = f_14 * hk_429[k]
                   + pb_x[k] * ik_429[k];

        t_529[k] = f_14 * hk_430[k]
                   + pb_x[k] * ik_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pa_z, pb_x, pb_z, hk_244, hk_245, hk_431, \
                         hl_306, hl_308, ik_424, ik_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_14 * hk_431[k]
                   + pb_x[k] * ik_431[k];

        t_531[k] = pa_z[k] * hl_306[k];

        t_532[k] = f_13 * hk_244[k]
                   + pb_z[k] * ik_424[k];

        t_533[k] = f_14 * hk_245[k]
                   + pa_z[k] * hl_308[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pa_z, hk_246, hk_247, hk_248, hk_249, \
                         hl_309, hl_310, hl_311, hl_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_15 * hk_246[k]
                   + pa_z[k] * hl_309[k];

        t_535[k] = f_16 * hk_247[k]
                   + pa_z[k] * hl_310[k];

        t_536[k] = f_17 * hk_248[k]
                   + pa_z[k] * hl_311[k];

        t_537[k] = f_0 * hk_249[k]
                   + pa_z[k] * hl_312[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_y, pa_z, pb_y, gl0_225, gl1_225, \
                         hk_251, hk_287, hk_288, hl_314, hl_360, ik_431, \
                         ik_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_15 * hk_287[k]
                   + pb_y[k] * ik_431[k];

        t_539[k] = f_18 * hk_251[k]
                   + pa_z[k] * hl_314[k];

        t_540[k] = f_19 * gl0_225[k]
                   - f_20 * gl1_225[k]
                   + pa_y[k] * hl_360[k];

        t_541[k] = f_14 * hk_288[k]
                   + pb_y[k] * ik_432[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_138 = buffer.data(gl0 + 138);
    const auto *gl0_141 = buffer.data(gl0 + 141);
    const auto *gl0_145 = buffer.data(gl0 + 145);
    const auto *gl0_150 = buffer.data(gl0 + 150);
    const auto *gl0_156 = buffer.data(gl0 + 156);
    const auto *gl0_225 = buffer.data(gl0 + 225);
    const auto *gl0_230 = buffer.data(gl0 + 230);
    const auto *gl0_234 = buffer.data(gl0 + 234);
    const auto *gl0_239 = buffer.data(gl0 + 239);
    const auto *gl0_245 = buffer.data(gl0 + 245);
    const auto *gl0_252 = buffer.data(gl0 + 252);
    const auto *gl0_576 = buffer.data(gl0 + 576);
    const auto *gl0_578 = buffer.data(gl0 + 578);
    const auto *gl0_579 = buffer.data(gl0 + 579);
    const auto *gl0_580 = buffer.data(gl0 + 580);
    const auto *gl0_581 = buffer.data(gl0 + 581);
    const auto *gl0_582 = buffer.data(gl0 + 582);
    const auto *gl0_584 = buffer.data(gl0 + 584);

    const auto *gl1_138 = buffer.data(gl1 + 138);
    const auto *gl1_141 = buffer.data(gl1 + 141);
    const auto *gl1_145 = buffer.data(gl1 + 145);
    const auto *gl1_150 = buffer.data(gl1 + 150);
    const auto *gl1_156 = buffer.data(gl1 + 156);
    const auto *gl1_225 = buffer.data(gl1 + 225);
    const auto *gl1_230 = buffer.data(gl1 + 230);
    const auto *gl1_234 = buffer.data(gl1 + 234);
    const auto *gl1_239 = buffer.data(gl1 + 239);
    const auto *gl1_245 = buffer.data(gl1 + 245);
    const auto *gl1_252 = buffer.data(gl1 + 252);
    const auto *gl1_576 = buffer.data(gl1 + 576);
    const auto *gl1_578 = buffer.data(gl1 + 578);
    const auto *gl1_579 = buffer.data(gl1 + 579);
    const auto *gl1_580 = buffer.data(gl1 + 580);
    const auto *gl1_581 = buffer.data(gl1 + 581);
    const auto *gl1_582 = buffer.data(gl1 + 582);
    const auto *gl1_584 = buffer.data(gl1 + 584);

    const auto *hk_252 = buffer.data(hk + 252);
    const auto *hk_255 = buffer.data(hk + 255);
    const auto *hk_258 = buffer.data(hk + 258);
    const auto *hk_262 = buffer.data(hk + 262);
    const auto *hk_267 = buffer.data(hk + 267);
    const auto *hk_280 = buffer.data(hk + 280);
    const auto *hk_290 = buffer.data(hk + 290);
    const auto *hk_291 = buffer.data(hk + 291);
    const auto *hk_293 = buffer.data(hk + 293);
    const auto *hk_294 = buffer.data(hk + 294);
    const auto *hk_297 = buffer.data(hk + 297);
    const auto *hk_298 = buffer.data(hk + 298);
    const auto *hk_302 = buffer.data(hk + 302);
    const auto *hk_303 = buffer.data(hk + 303);
    const auto *hk_308 = buffer.data(hk + 308);
    const auto *hk_316 = buffer.data(hk + 316);
    const auto *hk_323 = buffer.data(hk + 323);
    const auto *hk_324 = buffer.data(hk + 324);
    const auto *hk_325 = buffer.data(hk + 325);
    const auto *hk_326 = buffer.data(hk + 326);
    const auto *hk_327 = buffer.data(hk + 327);
    const auto *hk_329 = buffer.data(hk + 329);
    const auto *hk_330 = buffer.data(hk + 330);
    const auto *hk_332 = buffer.data(hk + 332);
    const auto *hk_333 = buffer.data(hk + 333);
    const auto *hk_334 = buffer.data(hk + 334);
    const auto *hk_336 = buffer.data(hk + 336);
    const auto *hk_337 = buffer.data(hk + 337);
    const auto *hk_338 = buffer.data(hk + 338);
    const auto *hk_339 = buffer.data(hk + 339);
    const auto *hk_341 = buffer.data(hk + 341);
    const auto *hk_342 = buffer.data(hk + 342);
    const auto *hk_343 = buffer.data(hk + 343);
    const auto *hk_344 = buffer.data(hk + 344);
    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_354 = buffer.data(hk + 354);
    const auto *hk_355 = buffer.data(hk + 355);
    const auto *hk_356 = buffer.data(hk + 356);
    const auto *hk_357 = buffer.data(hk + 357);
    const auto *hk_358 = buffer.data(hk + 358);
    const auto *hk_359 = buffer.data(hk + 359);
    const auto *hk_444 = buffer.data(hk + 444);
    const auto *hk_449 = buffer.data(hk + 449);
    const auto *hk_450 = buffer.data(hk + 450);
    const auto *hk_455 = buffer.data(hk + 455);
    const auto *hk_456 = buffer.data(hk + 456);
    const auto *hk_457 = buffer.data(hk + 457);
    const auto *hk_460 = buffer.data(hk + 460);
    const auto *hk_461 = buffer.data(hk + 461);
    const auto *hk_462 = buffer.data(hk + 462);
    const auto *hk_463 = buffer.data(hk + 463);
    const auto *hk_464 = buffer.data(hk + 464);
    const auto *hk_465 = buffer.data(hk + 465);
    const auto *hk_466 = buffer.data(hk + 466);
    const auto *hk_467 = buffer.data(hk + 467);
    const auto *hk_496 = buffer.data(hk + 496);
    const auto *hk_497 = buffer.data(hk + 497);
    const auto *hk_498 = buffer.data(hk + 498);
    const auto *hk_499 = buffer.data(hk + 499);
    const auto *hk_500 = buffer.data(hk + 500);
    const auto *hk_501 = buffer.data(hk + 501);
    const auto *hk_502 = buffer.data(hk + 502);
    const auto *hk_509 = buffer.data(hk + 509);
    const auto *hk_513 = buffer.data(hk + 513);
    const auto *hk_518 = buffer.data(hk + 518);
    const auto *hk_524 = buffer.data(hk + 524);
    const auto *hk_531 = buffer.data(hk + 531);
    const auto *hk_532 = buffer.data(hk + 532);
    const auto *hk_533 = buffer.data(hk + 533);
    const auto *hk_534 = buffer.data(hk + 534);
    const auto *hk_535 = buffer.data(hk + 535);
    const auto *hk_536 = buffer.data(hk + 536);
    const auto *hk_537 = buffer.data(hk + 537);
    const auto *hk_539 = buffer.data(hk + 539);

    const auto *hl_318 = buffer.data(hl + 318);
    const auto *hl_321 = buffer.data(hl + 321);
    const auto *hl_325 = buffer.data(hl + 325);
    const auto *hl_330 = buffer.data(hl + 330);
    const auto *hl_336 = buffer.data(hl + 336);
    const auto *hl_365 = buffer.data(hl + 365);
    const auto *hl_369 = buffer.data(hl + 369);
    const auto *hl_374 = buffer.data(hl + 374);
    const auto *hl_380 = buffer.data(hl + 380);
    const auto *hl_387 = buffer.data(hl + 387);
    const auto *hl_405 = buffer.data(hl + 405);
    const auto *hl_407 = buffer.data(hl + 407);
    const auto *hl_408 = buffer.data(hl + 408);
    const auto *hl_410 = buffer.data(hl + 410);
    const auto *hl_411 = buffer.data(hl + 411);
    const auto *hl_414 = buffer.data(hl + 414);
    const auto *hl_415 = buffer.data(hl + 415);
    const auto *hl_417 = buffer.data(hl + 417);
    const auto *hl_419 = buffer.data(hl + 419);
    const auto *hl_420 = buffer.data(hl + 420);
    const auto *hl_422 = buffer.data(hl + 422);
    const auto *hl_423 = buffer.data(hl + 423);
    const auto *hl_425 = buffer.data(hl + 425);
    const auto *hl_426 = buffer.data(hl + 426);
    const auto *hl_428 = buffer.data(hl + 428);
    const auto *hl_429 = buffer.data(hl + 429);
    const auto *hl_430 = buffer.data(hl + 430);
    const auto *hl_432 = buffer.data(hl + 432);
    const auto *hl_440 = buffer.data(hl + 440);
    const auto *hl_441 = buffer.data(hl + 441);
    const auto *hl_443 = buffer.data(hl + 443);
    const auto *hl_444 = buffer.data(hl + 444);
    const auto *hl_445 = buffer.data(hl + 445);
    const auto *hl_446 = buffer.data(hl + 446);
    const auto *hl_447 = buffer.data(hl + 447);
    const auto *hl_449 = buffer.data(hl + 449);
    const auto *hl_576 = buffer.data(hl + 576);
    const auto *hl_578 = buffer.data(hl + 578);
    const auto *hl_579 = buffer.data(hl + 579);
    const auto *hl_580 = buffer.data(hl + 580);
    const auto *hl_581 = buffer.data(hl + 581);
    const auto *hl_582 = buffer.data(hl + 582);
    const auto *hl_584 = buffer.data(hl + 584);

    const auto *ii0_348 = buffer.data(ii0 + 348);
    const auto *ii0_353 = buffer.data(ii0 + 353);
    const auto *ii0_354 = buffer.data(ii0 + 354);
    const auto *ii0_359 = buffer.data(ii0 + 359);
    const auto *ii0_360 = buffer.data(ii0 + 360);
    const auto *ii0_361 = buffer.data(ii0 + 361);
    const auto *ii0_392 = buffer.data(ii0 + 392);
    const auto *ii0_393 = buffer.data(ii0 + 393);
    const auto *ii0_395 = buffer.data(ii0 + 395);
    const auto *ii0_397 = buffer.data(ii0 + 397);
    const auto *ii0_398 = buffer.data(ii0 + 398);
    const auto *ii0_400 = buffer.data(ii0 + 400);
    const auto *ii0_401 = buffer.data(ii0 + 401);
    const auto *ii0_402 = buffer.data(ii0 + 402);
    const auto *ii0_404 = buffer.data(ii0 + 404);
    const auto *ii0_405 = buffer.data(ii0 + 405);
    const auto *ii0_406 = buffer.data(ii0 + 406);
    const auto *ii0_412 = buffer.data(ii0 + 412);
    const auto *ii0_419 = buffer.data(ii0 + 419);

    const auto *ii1_348 = buffer.data(ii1 + 348);
    const auto *ii1_353 = buffer.data(ii1 + 353);
    const auto *ii1_354 = buffer.data(ii1 + 354);
    const auto *ii1_359 = buffer.data(ii1 + 359);
    const auto *ii1_360 = buffer.data(ii1 + 360);
    const auto *ii1_361 = buffer.data(ii1 + 361);
    const auto *ii1_392 = buffer.data(ii1 + 392);
    const auto *ii1_393 = buffer.data(ii1 + 393);
    const auto *ii1_395 = buffer.data(ii1 + 395);
    const auto *ii1_397 = buffer.data(ii1 + 397);
    const auto *ii1_398 = buffer.data(ii1 + 398);
    const auto *ii1_400 = buffer.data(ii1 + 400);
    const auto *ii1_401 = buffer.data(ii1 + 401);
    const auto *ii1_402 = buffer.data(ii1 + 402);
    const auto *ii1_404 = buffer.data(ii1 + 404);
    const auto *ii1_405 = buffer.data(ii1 + 405);
    const auto *ii1_406 = buffer.data(ii1 + 406);
    const auto *ii1_412 = buffer.data(ii1 + 412);
    const auto *ii1_419 = buffer.data(ii1 + 419);

    const auto *ik_432 = buffer.data(ik + 432);
    const auto *ik_434 = buffer.data(ik + 434);
    const auto *ik_435 = buffer.data(ik + 435);
    const auto *ik_437 = buffer.data(ik + 437);
    const auto *ik_438 = buffer.data(ik + 438);
    const auto *ik_441 = buffer.data(ik + 441);
    const auto *ik_442 = buffer.data(ik + 442);
    const auto *ik_444 = buffer.data(ik + 444);
    const auto *ik_446 = buffer.data(ik + 446);
    const auto *ik_447 = buffer.data(ik + 447);
    const auto *ik_449 = buffer.data(ik + 449);
    const auto *ik_450 = buffer.data(ik + 450);
    const auto *ik_452 = buffer.data(ik + 452);
    const auto *ik_455 = buffer.data(ik + 455);
    const auto *ik_456 = buffer.data(ik + 456);
    const auto *ik_457 = buffer.data(ik + 457);
    const auto *ik_460 = buffer.data(ik + 460);
    const auto *ik_461 = buffer.data(ik + 461);
    const auto *ik_462 = buffer.data(ik + 462);
    const auto *ik_463 = buffer.data(ik + 463);
    const auto *ik_464 = buffer.data(ik + 464);
    const auto *ik_465 = buffer.data(ik + 465);
    const auto *ik_466 = buffer.data(ik + 466);
    const auto *ik_467 = buffer.data(ik + 467);
    const auto *ik_468 = buffer.data(ik + 468);
    const auto *ik_470 = buffer.data(ik + 470);
    const auto *ik_471 = buffer.data(ik + 471);
    const auto *ik_473 = buffer.data(ik + 473);
    const auto *ik_474 = buffer.data(ik + 474);
    const auto *ik_477 = buffer.data(ik + 477);
    const auto *ik_478 = buffer.data(ik + 478);
    const auto *ik_482 = buffer.data(ik + 482);
    const auto *ik_483 = buffer.data(ik + 483);
    const auto *ik_488 = buffer.data(ik + 488);
    const auto *ik_496 = buffer.data(ik + 496);
    const auto *ik_497 = buffer.data(ik + 497);
    const auto *ik_498 = buffer.data(ik + 498);
    const auto *ik_499 = buffer.data(ik + 499);
    const auto *ik_500 = buffer.data(ik + 500);
    const auto *ik_501 = buffer.data(ik + 501);
    const auto *ik_502 = buffer.data(ik + 502);
    const auto *ik_503 = buffer.data(ik + 503);
    const auto *ik_504 = buffer.data(ik + 504);
    const auto *ik_505 = buffer.data(ik + 505);
    const auto *ik_506 = buffer.data(ik + 506);
    const auto *ik_507 = buffer.data(ik + 507);
    const auto *ik_509 = buffer.data(ik + 509);
    const auto *ik_510 = buffer.data(ik + 510);
    const auto *ik_512 = buffer.data(ik + 512);
    const auto *ik_513 = buffer.data(ik + 513);
    const auto *ik_514 = buffer.data(ik + 514);
    const auto *ik_516 = buffer.data(ik + 516);
    const auto *ik_517 = buffer.data(ik + 517);
    const auto *ik_518 = buffer.data(ik + 518);
    const auto *ik_519 = buffer.data(ik + 519);
    const auto *ik_521 = buffer.data(ik + 521);
    const auto *ik_522 = buffer.data(ik + 522);
    const auto *ik_523 = buffer.data(ik + 523);
    const auto *ik_524 = buffer.data(ik + 524);
    const auto *ik_531 = buffer.data(ik + 531);
    const auto *ik_532 = buffer.data(ik + 532);
    const auto *ik_533 = buffer.data(ik + 533);
    const auto *ik_534 = buffer.data(ik + 534);
    const auto *ik_535 = buffer.data(ik + 535);
    const auto *ik_536 = buffer.data(ik + 536);
    const auto *ik_537 = buffer.data(ik + 537);
    const auto *ik_539 = buffer.data(ik + 539);

#pragma omp simd aligned(t_542, t_543, t_544, pa_z, pb_y, pb_z, gl0_138, gl1_138, hk_252, \
                         hk_290, hl_318, ik_432, ik_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_14 * hk_252[k]
                   + pb_z[k] * ik_432[k];

        t_543[k] = f_19 * gl0_138[k]
                   - f_20 * gl1_138[k]
                   + pa_z[k] * hl_318[k];

        t_544[k] = f_14 * hk_290[k]
                   + pb_y[k] * ik_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pa_y, pa_z, pb_z, gl0_141, gl0_230, gl1_141, \
                         gl1_230, hk_255, hl_321, hl_365, ik_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_19 * gl0_230[k]
                   - f_20 * gl1_230[k]
                   + pa_y[k] * hl_365[k];

        t_546[k] = f_19 * gl0_141[k]
                   - f_20 * gl1_141[k]
                   + pa_z[k] * hl_321[k];

        t_547[k] = f_14 * hk_255[k]
                   + pb_z[k] * ik_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pa_y, pa_z, pb_y, gl0_145, gl0_234, gl1_145, \
                         gl1_234, hk_293, hl_325, hl_369, ik_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_14 * hk_293[k]
                   + pb_y[k] * ik_437[k];

        t_549[k] = f_19 * gl0_234[k]
                   - f_20 * gl1_234[k]
                   + pa_y[k] * hl_369[k];

        t_550[k] = f_19 * gl0_145[k]
                   - f_20 * gl1_145[k]
                   + pa_z[k] * hl_325[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pb_x, pb_y, pb_z, hk_258, hk_297, hk_444, \
                         ii0_348, ii1_348, ik_438, ik_441, ik_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_14 * hk_258[k]
                   + pb_z[k] * ik_438[k];

        t_552[k] = f_14 * hk_444[k]
                   + f_7 * ii0_348[k]
                   - f_8 * ii1_348[k]
                   + pb_x[k] * ik_444[k];

        t_553[k] = f_14 * hk_297[k]
                   + pb_y[k] * ik_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pa_y, pa_z, pb_z, gl0_150, gl0_239, gl1_150, \
                         gl1_239, hk_262, hl_330, hl_374, ik_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_19 * gl0_239[k]
                   - f_20 * gl1_239[k]
                   + pa_y[k] * hl_374[k];

        t_555[k] = f_19 * gl0_150[k]
                   - f_20 * gl1_150[k]
                   + pa_z[k] * hl_330[k];

        t_556[k] = f_14 * hk_262[k]
                   + pb_z[k] * ik_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pb_x, pb_y, hk_302, hk_449, hk_450, ii0_353, \
                         ii0_354, ii1_353, ii1_354, ik_446, ik_449, \
                         ik_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_14 * hk_449[k]
                   + f_5 * ii0_353[k]
                   - f_6 * ii1_353[k]
                   + pb_x[k] * ik_449[k];

        t_558[k] = f_14 * hk_450[k]
                   + f_5 * ii0_354[k]
                   - f_6 * ii1_354[k]
                   + pb_x[k] * ik_450[k];

        t_559[k] = f_14 * hk_302[k]
                   + pb_y[k] * ik_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pa_y, pa_z, pb_z, gl0_156, gl0_245, gl1_156, \
                         gl1_245, hk_267, hl_336, hl_380, ik_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_19 * gl0_245[k]
                   - f_20 * gl1_245[k]
                   + pa_y[k] * hl_380[k];

        t_561[k] = f_19 * gl0_156[k]
                   - f_20 * gl1_156[k]
                   + pa_z[k] * hl_336[k];

        t_562[k] = f_14 * hk_267[k]
                   + pb_z[k] * ik_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pb_x, hk_455, hk_456, hk_457, ii0_359, ii0_360, \
                         ii0_361, ii1_359, ii1_360, ii1_361, ik_455, ik_456, \
                         ik_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_14 * hk_455[k]
                   + f_3 * ii0_359[k]
                   - f_4 * ii1_359[k]
                   + pb_x[k] * ik_455[k];

        t_564[k] = f_14 * hk_456[k]
                   + f_3 * ii0_360[k]
                   - f_4 * ii1_360[k]
                   + pb_x[k] * ik_456[k];

        t_565[k] = f_14 * hk_457[k]
                   + f_3 * ii0_361[k]
                   - f_4 * ii1_361[k]
                   + pb_x[k] * ik_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_y, pb_x, pb_y, gl0_252, gl1_252, \
                         hk_308, hk_460, hk_461, hl_387, ik_452, ik_460, \
                         ik_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_14 * hk_308[k]
                   + pb_y[k] * ik_452[k];

        t_567[k] = f_19 * gl0_252[k]
                   - f_20 * gl1_252[k]
                   + pa_y[k] * hl_387[k];

        t_568[k] = f_14 * hk_460[k]
                   + pb_x[k] * ik_460[k];

        t_569[k] = f_14 * hk_461[k]
                   + pb_x[k] * ik_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pb_x, hk_462, hk_463, hk_464, \
                         hk_465, hk_466, ik_462, ik_463, ik_464, ik_465, \
                         ik_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_14 * hk_462[k]
                   + pb_x[k] * ik_462[k];

        t_571[k] = f_14 * hk_463[k]
                   + pb_x[k] * ik_463[k];

        t_572[k] = f_14 * hk_464[k]
                   + pb_x[k] * ik_464[k];

        t_573[k] = f_14 * hk_465[k]
                   + pb_x[k] * ik_465[k];

        t_574[k] = f_14 * hk_466[k]
                   + pb_x[k] * ik_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pa_x, pb_x, pb_z, gl0_576, gl1_576, hk_280, \
                         hk_467, hl_576, ik_460, ik_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_14 * hk_467[k]
                   + pb_x[k] * ik_467[k];

        t_576[k] = f_19 * gl0_576[k]
                   - f_20 * gl1_576[k]
                   + pa_x[k] * hl_576[k];

        t_577[k] = f_14 * hk_280[k]
                   + pb_z[k] * ik_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_x, gl0_578, gl0_579, gl0_580, gl1_578, \
                         gl1_579, gl1_580, hl_578, hl_579, hl_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_19 * gl0_578[k]
                   - f_20 * gl1_578[k]
                   + pa_x[k] * hl_578[k];

        t_579[k] = f_19 * gl0_579[k]
                   - f_20 * gl1_579[k]
                   + pa_x[k] * hl_579[k];

        t_580[k] = f_19 * gl0_580[k]
                   - f_20 * gl1_580[k]
                   + pa_x[k] * hl_580[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_x, pb_y, gl0_581, gl0_582, gl1_581, gl1_582, \
                         hk_323, hl_581, hl_582, ik_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_19 * gl0_581[k]
                   - f_20 * gl1_581[k]
                   + pa_x[k] * hl_581[k];

        t_582[k] = f_19 * gl0_582[k]
                   - f_20 * gl1_582[k]
                   + pa_x[k] * hl_582[k];

        t_583[k] = f_14 * hk_323[k]
                   + pb_y[k] * ik_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_x, pa_y, pb_y, gl0_584, gl1_584, \
                         hk_324, hl_405, hl_407, hl_584, ik_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_19 * gl0_584[k]
                   - f_20 * gl1_584[k]
                   + pa_x[k] * hl_584[k];

        t_585[k] = pa_y[k] * hl_405[k];

        t_586[k] = f_13 * hk_324[k]
                   + pb_y[k] * ik_468[k];

        t_587[k] = pa_y[k] * hl_407[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_y, pb_y, hk_325, hk_326, hk_327, \
                         hl_408, hl_410, hl_411, ik_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_14 * hk_325[k]
                   + pa_y[k] * hl_408[k];

        t_589[k] = f_13 * hk_326[k]
                   + pb_y[k] * ik_470[k];

        t_590[k] = pa_y[k] * hl_410[k];

        t_591[k] = f_15 * hk_327[k]
                   + pa_y[k] * hl_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_y, pb_y, pb_z, hk_291, hk_329, hk_330, \
                         hl_414, hl_415, ik_471, ik_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_15 * hk_291[k]
                   + pb_z[k] * ik_471[k];

        t_593[k] = f_13 * hk_329[k]
                   + pb_y[k] * ik_473[k];

        t_594[k] = pa_y[k] * hl_414[k];

        t_595[k] = f_16 * hk_330[k]
                   + pa_y[k] * hl_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_y, pb_y, pb_z, hk_294, hk_332, hk_333, \
                         hl_417, hl_419, ik_474, ik_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * hk_294[k]
                   + pb_z[k] * ik_474[k];

        t_597[k] = f_14 * hk_332[k]
                   + pa_y[k] * hl_417[k];

        t_598[k] = f_13 * hk_333[k]
                   + pb_y[k] * ik_477[k];

        t_599[k] = pa_y[k] * hl_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pa_y, pb_z, hk_298, hk_334, hk_336, \
                         hk_337, hl_420, hl_422, hl_423, ik_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_17 * hk_334[k]
                   + pa_y[k] * hl_420[k];

        t_601[k] = f_15 * hk_298[k]
                   + pb_z[k] * ik_478[k];

        t_602[k] = f_15 * hk_336[k]
                   + pa_y[k] * hl_422[k];

        t_603[k] = f_14 * hk_337[k]
                   + pa_y[k] * hl_423[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_y, pb_y, pb_z, hk_303, hk_338, hk_339, \
                         hl_425, hl_426, ik_482, ik_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_13 * hk_338[k]
                   + pb_y[k] * ik_482[k];

        t_605[k] = pa_y[k] * hl_425[k];

        t_606[k] = f_0 * hk_339[k]
                   + pa_y[k] * hl_426[k];

        t_607[k] = f_15 * hk_303[k]
                   + pb_z[k] * ik_483[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, t_612, pa_y, pb_y, hk_341, hk_342, \
                         hk_343, hk_344, hl_428, hl_429, hl_430, hl_432, \
                         ik_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_16 * hk_341[k]
                   + pa_y[k] * hl_428[k];

        t_609[k] = f_15 * hk_342[k]
                   + pa_y[k] * hl_429[k];

        t_610[k] = f_14 * hk_343[k]
                   + pa_y[k] * hl_430[k];

        t_611[k] = f_13 * hk_344[k]
                   + pb_y[k] * ik_488[k];

        t_612[k] = pa_y[k] * hl_432[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, t_617, pb_x, hk_496, hk_497, hk_498, \
                         hk_499, hk_500, ik_496, ik_497, ik_498, ik_499, \
                         ik_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_14 * hk_496[k]
                   + pb_x[k] * ik_496[k];

        t_614[k] = f_14 * hk_497[k]
                   + pb_x[k] * ik_497[k];

        t_615[k] = f_14 * hk_498[k]
                   + pb_x[k] * ik_498[k];

        t_616[k] = f_14 * hk_499[k]
                   + pb_x[k] * ik_499[k];

        t_617[k] = f_14 * hk_500[k]
                   + pb_x[k] * ik_500[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pa_y, pb_x, hk_352, hk_501, hk_502, \
                         hl_440, hl_441, ik_501, ik_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_14 * hk_501[k]
                   + pb_x[k] * ik_501[k];

        t_619[k] = f_14 * hk_502[k]
                   + pb_x[k] * ik_502[k];

        t_620[k] = pa_y[k] * hl_440[k];

        t_621[k] = f_18 * hk_352[k]
                   + pa_y[k] * hl_441[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_y, pb_z, hk_316, hk_354, hk_355, \
                         hk_356, hl_443, hl_444, hl_445, ik_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_15 * hk_316[k]
                   + pb_z[k] * ik_496[k];

        t_623[k] = f_0 * hk_354[k]
                   + pa_y[k] * hl_443[k];

        t_624[k] = f_17 * hk_355[k]
                   + pa_y[k] * hl_444[k];

        t_625[k] = f_16 * hk_356[k]
                   + pa_y[k] * hl_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pa_y, pb_y, hk_357, hk_358, hk_359, \
                         hl_446, hl_447, hl_449, ik_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * hk_357[k]
                   + pa_y[k] * hl_446[k];

        t_627[k] = f_14 * hk_358[k]
                   + pa_y[k] * hl_447[k];

        t_628[k] = f_13 * hk_359[k]
                   + pb_y[k] * ik_503[k];

        t_629[k] = pa_y[k] * hl_449[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pa_z, pb_y, pb_z, gl0_225, gl1_225, \
                         hk_324, hl_405, ii0_392, ii1_392, ik_504, \
                         ik_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_21 * gl0_225[k]
                   - f_22 * gl1_225[k]
                   + pa_z[k] * hl_405[k];

        t_631[k] = pb_y[k] * ik_504[k];

        t_632[k] = f_16 * hk_324[k]
                   + pb_z[k] * ik_504[k];

        t_633[k] = f_3 * ii0_392[k]
                   - f_4 * ii1_392[k]
                   + pb_y[k] * ik_505[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_x, pb_y, pb_z, hk_327, hk_509, \
                         ii0_393, ii0_397, ii1_393, ii1_397, ik_506, ik_507, \
                         ik_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = pb_y[k] * ik_506[k];

        t_635[k] = f_14 * hk_509[k]
                   + f_11 * ii0_397[k]
                   - f_12 * ii1_397[k]
                   + pb_x[k] * ik_509[k];

        t_636[k] = f_5 * ii0_393[k]
                   - f_6 * ii1_393[k]
                   + pb_y[k] * ik_507[k];

        t_637[k] = f_16 * hk_327[k]
                   + pb_z[k] * ik_507[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pb_x, pb_y, pb_z, hk_330, hk_513, \
                         ii0_395, ii0_401, ii1_395, ii1_401, ik_509, ik_510, \
                         ik_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = pb_y[k] * ik_509[k];

        t_639[k] = f_14 * hk_513[k]
                   + f_9 * ii0_401[k]
                   - f_10 * ii1_401[k]
                   + pb_x[k] * ik_513[k];

        t_640[k] = f_7 * ii0_395[k]
                   - f_8 * ii1_395[k]
                   + pb_y[k] * ik_510[k];

        t_641[k] = f_16 * hk_330[k]
                   + pb_z[k] * ik_510[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pb_x, pb_y, hk_518, ii0_397, ii0_406, ii1_397, \
                         ii1_406, ik_512, ik_513, ik_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_3 * ii0_397[k]
                   - f_4 * ii1_397[k]
                   + pb_y[k] * ik_512[k];

        t_643[k] = pb_y[k] * ik_513[k];

        t_644[k] = f_14 * hk_518[k]
                   + f_7 * ii0_406[k]
                   - f_8 * ii1_406[k]
                   + pb_x[k] * ik_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, pb_y, pb_z, hk_334, ii0_398, ii0_400, \
                         ii0_401, ii1_398, ii1_400, ii1_401, ik_514, ik_516, \
                         ik_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_9 * ii0_398[k]
                   - f_10 * ii1_398[k]
                   + pb_y[k] * ik_514[k];

        t_646[k] = f_16 * hk_334[k]
                   + pb_z[k] * ik_514[k];

        t_647[k] = f_5 * ii0_400[k]
                   - f_6 * ii1_400[k]
                   + pb_y[k] * ik_516[k];

        t_648[k] = f_3 * ii0_401[k]
                   - f_4 * ii1_401[k]
                   + pb_y[k] * ik_517[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pb_x, pb_y, pb_z, hk_339, hk_524, \
                         ii0_402, ii0_412, ii1_402, ii1_412, ik_518, ik_519, \
                         ik_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pb_y[k] * ik_518[k];

        t_650[k] = f_14 * hk_524[k]
                   + f_5 * ii0_412[k]
                   - f_6 * ii1_412[k]
                   + pb_x[k] * ik_524[k];

        t_651[k] = f_11 * ii0_402[k]
                   - f_12 * ii1_402[k]
                   + pb_y[k] * ik_519[k];

        t_652[k] = f_16 * hk_339[k]
                   + pb_z[k] * ik_519[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, pb_y, ii0_404, ii0_405, ii0_406, ii1_404, \
                         ii1_405, ii1_406, ik_521, ik_522, ik_523, \
                         ik_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_7 * ii0_404[k]
                   - f_8 * ii1_404[k]
                   + pb_y[k] * ik_521[k];

        t_654[k] = f_5 * ii0_405[k]
                   - f_6 * ii1_405[k]
                   + pb_y[k] * ik_522[k];

        t_655[k] = f_3 * ii0_406[k]
                   - f_4 * ii1_406[k]
                   + pb_y[k] * ik_523[k];

        t_656[k] = pb_y[k] * ik_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pb_x, hk_531, hk_532, hk_533, hk_534, \
                         ii0_419, ii1_419, ik_531, ik_532, ik_533, \
                         ik_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_14 * hk_531[k]
                   + f_3 * ii0_419[k]
                   - f_4 * ii1_419[k]
                   + pb_x[k] * ik_531[k];

        t_658[k] = f_14 * hk_532[k]
                   + pb_x[k] * ik_532[k];

        t_659[k] = f_14 * hk_533[k]
                   + pb_x[k] * ik_533[k];

        t_660[k] = f_14 * hk_534[k]
                   + pb_x[k] * ik_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pb_x, pb_y, hk_535, hk_536, \
                         hk_537, hk_539, ik_531, ik_535, ik_536, ik_537, \
                         ik_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_14 * hk_535[k]
                   + pb_x[k] * ik_535[k];

        t_662[k] = f_14 * hk_536[k]
                   + pb_x[k] * ik_536[k];

        t_663[k] = f_14 * hk_537[k]
                   + pb_x[k] * ik_537[k];

        t_664[k] = pb_y[k] * ik_531[k];

        t_665[k] = f_14 * hk_539[k]
                   + pb_x[k] * ik_539[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_674 = buffer.data(gl0 + 674);

    const auto *gl1_674 = buffer.data(gl1 + 674);

    const auto *hk_352 = buffer.data(hk + 352);
    const auto *hk_360 = buffer.data(hk + 360);
    const auto *hk_363 = buffer.data(hk + 363);
    const auto *hk_365 = buffer.data(hk + 365);
    const auto *hk_366 = buffer.data(hk + 366);
    const auto *hk_369 = buffer.data(hk + 369);
    const auto *hk_370 = buffer.data(hk + 370);
    const auto *hk_374 = buffer.data(hk + 374);
    const auto *hk_375 = buffer.data(hk + 375);
    const auto *hk_380 = buffer.data(hk + 380);
    const auto *hk_396 = buffer.data(hk + 396);
    const auto *hk_398 = buffer.data(hk + 398);
    const auto *hk_399 = buffer.data(hk + 399);
    const auto *hk_401 = buffer.data(hk + 401);
    const auto *hk_402 = buffer.data(hk + 402);
    const auto *hk_405 = buffer.data(hk + 405);
    const auto *hk_406 = buffer.data(hk + 406);
    const auto *hk_410 = buffer.data(hk + 410);
    const auto *hk_411 = buffer.data(hk + 411);
    const auto *hk_416 = buffer.data(hk + 416);
    const auto *hk_432 = buffer.data(hk + 432);
    const auto *hk_434 = buffer.data(hk + 434);
    const auto *hk_437 = buffer.data(hk + 437);
    const auto *hk_441 = buffer.data(hk + 441);
    const auto *hk_446 = buffer.data(hk + 446);
    const auto *hk_452 = buffer.data(hk + 452);
    const auto *hk_468 = buffer.data(hk + 468);
    const auto *hk_470 = buffer.data(hk + 470);
    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_543 = buffer.data(hk + 543);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_546 = buffer.data(hk + 546);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_550 = buffer.data(hk + 550);
    const auto *hk_552 = buffer.data(hk + 552);
    const auto *hk_554 = buffer.data(hk + 554);
    const auto *hk_555 = buffer.data(hk + 555);
    const auto *hk_557 = buffer.data(hk + 557);
    const auto *hk_558 = buffer.data(hk + 558);
    const auto *hk_560 = buffer.data(hk + 560);
    const auto *hk_561 = buffer.data(hk + 561);
    const auto *hk_563 = buffer.data(hk + 563);
    const auto *hk_564 = buffer.data(hk + 564);
    const auto *hk_565 = buffer.data(hk + 565);
    const auto *hk_567 = buffer.data(hk + 567);
    const auto *hk_568 = buffer.data(hk + 568);
    const auto *hk_570 = buffer.data(hk + 570);
    const auto *hk_571 = buffer.data(hk + 571);
    const auto *hk_572 = buffer.data(hk + 572);
    const auto *hk_573 = buffer.data(hk + 573);
    const auto *hk_574 = buffer.data(hk + 574);
    const auto *hk_575 = buffer.data(hk + 575);
    const auto *hk_581 = buffer.data(hk + 581);
    const auto *hk_585 = buffer.data(hk + 585);
    const auto *hk_588 = buffer.data(hk + 588);
    const auto *hk_590 = buffer.data(hk + 590);
    const auto *hk_593 = buffer.data(hk + 593);
    const auto *hk_594 = buffer.data(hk + 594);
    const auto *hk_596 = buffer.data(hk + 596);
    const auto *hk_599 = buffer.data(hk + 599);
    const auto *hk_600 = buffer.data(hk + 600);
    const auto *hk_601 = buffer.data(hk + 601);
    const auto *hk_603 = buffer.data(hk + 603);
    const auto *hk_605 = buffer.data(hk + 605);
    const auto *hk_606 = buffer.data(hk + 606);
    const auto *hk_607 = buffer.data(hk + 607);
    const auto *hk_608 = buffer.data(hk + 608);
    const auto *hk_609 = buffer.data(hk + 609);
    const auto *hk_610 = buffer.data(hk + 610);
    const auto *hk_611 = buffer.data(hk + 611);
    const auto *hk_612 = buffer.data(hk + 612);
    const auto *hk_615 = buffer.data(hk + 615);
    const auto *hk_617 = buffer.data(hk + 617);
    const auto *hk_618 = buffer.data(hk + 618);
    const auto *hk_621 = buffer.data(hk + 621);
    const auto *hk_622 = buffer.data(hk + 622);
    const auto *hk_624 = buffer.data(hk + 624);
    const auto *hk_626 = buffer.data(hk + 626);
    const auto *hk_627 = buffer.data(hk + 627);
    const auto *hk_629 = buffer.data(hk + 629);
    const auto *hk_630 = buffer.data(hk + 630);
    const auto *hk_632 = buffer.data(hk + 632);
    const auto *hk_633 = buffer.data(hk + 633);
    const auto *hk_635 = buffer.data(hk + 635);
    const auto *hk_636 = buffer.data(hk + 636);
    const auto *hk_637 = buffer.data(hk + 637);
    const auto *hk_639 = buffer.data(hk + 639);
    const auto *hk_640 = buffer.data(hk + 640);
    const auto *hk_641 = buffer.data(hk + 641);
    const auto *hk_642 = buffer.data(hk + 642);
    const auto *hk_643 = buffer.data(hk + 643);
    const auto *hk_644 = buffer.data(hk + 644);
    const auto *hk_645 = buffer.data(hk + 645);
    const auto *hk_646 = buffer.data(hk + 646);
    const auto *hk_647 = buffer.data(hk + 647);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_651 = buffer.data(hk + 651);

    const auto *hl_450 = buffer.data(hl + 450);
    const auto *hl_451 = buffer.data(hl + 451);
    const auto *hl_453 = buffer.data(hl + 453);
    const auto *hl_456 = buffer.data(hl + 456);
    const auto *hl_460 = buffer.data(hl + 460);
    const auto *hl_465 = buffer.data(hl + 465);
    const auto *hl_471 = buffer.data(hl + 471);
    const auto *hl_478 = buffer.data(hl + 478);
    const auto *hl_674 = buffer.data(hl + 674);
    const auto *hl_675 = buffer.data(hl + 675);
    const auto *hl_678 = buffer.data(hl + 678);
    const auto *hl_680 = buffer.data(hl + 680);
    const auto *hl_681 = buffer.data(hl + 681);
    const auto *hl_684 = buffer.data(hl + 684);
    const auto *hl_685 = buffer.data(hl + 685);
    const auto *hl_687 = buffer.data(hl + 687);
    const auto *hl_689 = buffer.data(hl + 689);
    const auto *hl_690 = buffer.data(hl + 690);
    const auto *hl_692 = buffer.data(hl + 692);
    const auto *hl_693 = buffer.data(hl + 693);
    const auto *hl_695 = buffer.data(hl + 695);
    const auto *hl_696 = buffer.data(hl + 696);
    const auto *hl_698 = buffer.data(hl + 698);
    const auto *hl_699 = buffer.data(hl + 699);
    const auto *hl_700 = buffer.data(hl + 700);
    const auto *hl_702 = buffer.data(hl + 702);
    const auto *hl_711 = buffer.data(hl + 711);
    const auto *hl_713 = buffer.data(hl + 713);
    const auto *hl_714 = buffer.data(hl + 714);
    const auto *hl_715 = buffer.data(hl + 715);
    const auto *hl_716 = buffer.data(hl + 716);
    const auto *hl_717 = buffer.data(hl + 717);
    const auto *hl_718 = buffer.data(hl + 718);
    const auto *hl_719 = buffer.data(hl + 719);
    const auto *hl_725 = buffer.data(hl + 725);
    const auto *hl_729 = buffer.data(hl + 729);
    const auto *hl_732 = buffer.data(hl + 732);
    const auto *hl_734 = buffer.data(hl + 734);
    const auto *hl_737 = buffer.data(hl + 737);
    const auto *hl_738 = buffer.data(hl + 738);
    const auto *hl_740 = buffer.data(hl + 740);
    const auto *hl_743 = buffer.data(hl + 743);
    const auto *hl_744 = buffer.data(hl + 744);
    const auto *hl_745 = buffer.data(hl + 745);
    const auto *hl_747 = buffer.data(hl + 747);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_757 = buffer.data(hl + 757);
    const auto *hl_758 = buffer.data(hl + 758);
    const auto *hl_759 = buffer.data(hl + 759);
    const auto *hl_760 = buffer.data(hl + 760);
    const auto *hl_761 = buffer.data(hl + 761);
    const auto *hl_762 = buffer.data(hl + 762);
    const auto *hl_763 = buffer.data(hl + 763);
    const auto *hl_764 = buffer.data(hl + 764);
    const auto *hl_765 = buffer.data(hl + 765);
    const auto *hl_768 = buffer.data(hl + 768);
    const auto *hl_770 = buffer.data(hl + 770);
    const auto *hl_771 = buffer.data(hl + 771);
    const auto *hl_774 = buffer.data(hl + 774);
    const auto *hl_775 = buffer.data(hl + 775);
    const auto *hl_777 = buffer.data(hl + 777);
    const auto *hl_779 = buffer.data(hl + 779);
    const auto *hl_780 = buffer.data(hl + 780);
    const auto *hl_782 = buffer.data(hl + 782);
    const auto *hl_783 = buffer.data(hl + 783);
    const auto *hl_785 = buffer.data(hl + 785);
    const auto *hl_786 = buffer.data(hl + 786);
    const auto *hl_788 = buffer.data(hl + 788);
    const auto *hl_789 = buffer.data(hl + 789);
    const auto *hl_790 = buffer.data(hl + 790);
    const auto *hl_792 = buffer.data(hl + 792);
    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_802 = buffer.data(hl + 802);
    const auto *hl_803 = buffer.data(hl + 803);
    const auto *hl_804 = buffer.data(hl + 804);
    const auto *hl_805 = buffer.data(hl + 805);
    const auto *hl_806 = buffer.data(hl + 806);
    const auto *hl_807 = buffer.data(hl + 807);
    const auto *hl_808 = buffer.data(hl + 808);
    const auto *hl_809 = buffer.data(hl + 809);
    const auto *hl_810 = buffer.data(hl + 810);
    const auto *hl_813 = buffer.data(hl + 813);

    const auto *ii0_413 = buffer.data(ii0 + 413);
    const auto *ii0_415 = buffer.data(ii0 + 415);
    const auto *ii0_416 = buffer.data(ii0 + 416);
    const auto *ii0_417 = buffer.data(ii0 + 417);
    const auto *ii0_418 = buffer.data(ii0 + 418);
    const auto *ii0_419 = buffer.data(ii0 + 419);

    const auto *ii1_413 = buffer.data(ii1 + 413);
    const auto *ii1_415 = buffer.data(ii1 + 415);
    const auto *ii1_416 = buffer.data(ii1 + 416);
    const auto *ii1_417 = buffer.data(ii1 + 417);
    const auto *ii1_418 = buffer.data(ii1 + 418);
    const auto *ii1_419 = buffer.data(ii1 + 419);

    const auto *ik_532 = buffer.data(ik + 532);
    const auto *ik_534 = buffer.data(ik + 534);
    const auto *ik_535 = buffer.data(ik + 535);
    const auto *ik_536 = buffer.data(ik + 536);
    const auto *ik_537 = buffer.data(ik + 537);
    const auto *ik_538 = buffer.data(ik + 538);
    const auto *ik_539 = buffer.data(ik + 539);
    const auto *ik_540 = buffer.data(ik + 540);
    const auto *ik_541 = buffer.data(ik + 541);
    const auto *ik_543 = buffer.data(ik + 543);
    const auto *ik_545 = buffer.data(ik + 545);
    const auto *ik_546 = buffer.data(ik + 546);
    const auto *ik_549 = buffer.data(ik + 549);
    const auto *ik_550 = buffer.data(ik + 550);
    const auto *ik_554 = buffer.data(ik + 554);
    const auto *ik_555 = buffer.data(ik + 555);
    const auto *ik_560 = buffer.data(ik + 560);
    const auto *ik_561 = buffer.data(ik + 561);
    const auto *ik_568 = buffer.data(ik + 568);
    const auto *ik_570 = buffer.data(ik + 570);
    const auto *ik_571 = buffer.data(ik + 571);
    const auto *ik_572 = buffer.data(ik + 572);
    const auto *ik_573 = buffer.data(ik + 573);
    const auto *ik_574 = buffer.data(ik + 574);
    const auto *ik_575 = buffer.data(ik + 575);
    const auto *ik_576 = buffer.data(ik + 576);
    const auto *ik_578 = buffer.data(ik + 578);
    const auto *ik_579 = buffer.data(ik + 579);
    const auto *ik_581 = buffer.data(ik + 581);
    const auto *ik_582 = buffer.data(ik + 582);
    const auto *ik_585 = buffer.data(ik + 585);
    const auto *ik_586 = buffer.data(ik + 586);
    const auto *ik_590 = buffer.data(ik + 590);
    const auto *ik_591 = buffer.data(ik + 591);
    const auto *ik_596 = buffer.data(ik + 596);
    const auto *ik_605 = buffer.data(ik + 605);
    const auto *ik_606 = buffer.data(ik + 606);
    const auto *ik_607 = buffer.data(ik + 607);
    const auto *ik_608 = buffer.data(ik + 608);
    const auto *ik_609 = buffer.data(ik + 609);
    const auto *ik_610 = buffer.data(ik + 610);
    const auto *ik_611 = buffer.data(ik + 611);
    const auto *ik_612 = buffer.data(ik + 612);
    const auto *ik_614 = buffer.data(ik + 614);
    const auto *ik_615 = buffer.data(ik + 615);
    const auto *ik_617 = buffer.data(ik + 617);
    const auto *ik_618 = buffer.data(ik + 618);
    const auto *ik_621 = buffer.data(ik + 621);
    const auto *ik_622 = buffer.data(ik + 622);
    const auto *ik_626 = buffer.data(ik + 626);
    const auto *ik_627 = buffer.data(ik + 627);
    const auto *ik_632 = buffer.data(ik + 632);
    const auto *ik_640 = buffer.data(ik + 640);
    const auto *ik_641 = buffer.data(ik + 641);
    const auto *ik_642 = buffer.data(ik + 642);
    const auto *ik_643 = buffer.data(ik + 643);
    const auto *ik_644 = buffer.data(ik + 644);
    const auto *ik_645 = buffer.data(ik + 645);
    const auto *ik_646 = buffer.data(ik + 646);
    const auto *ik_647 = buffer.data(ik + 647);
    const auto *ik_648 = buffer.data(ik + 648);
    const auto *ik_650 = buffer.data(ik + 650);

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_y, pb_z, hk_352, ii0_413, ii0_415, \
                         ii0_416, ii1_413, ii1_415, ii1_416, ik_532, ik_534, \
                         ik_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * ii0_413[k]
                   - f_2 * ii1_413[k]
                   + pb_y[k] * ik_532[k];

        t_667[k] = f_16 * hk_352[k]
                   + pb_z[k] * ik_532[k];

        t_668[k] = f_11 * ii0_415[k]
                   - f_12 * ii1_415[k]
                   + pb_y[k] * ik_534[k];

        t_669[k] = f_9 * ii0_416[k]
                   - f_10 * ii1_416[k]
                   + pb_y[k] * ik_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_y, ii0_417, ii0_418, ii0_419, ii1_417, \
                         ii1_418, ii1_419, ik_536, ik_537, ik_538, \
                         ik_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_7 * ii0_417[k]
                   - f_8 * ii1_417[k]
                   + pb_y[k] * ik_536[k];

        t_671[k] = f_5 * ii0_418[k]
                   - f_6 * ii1_418[k]
                   + pb_y[k] * ik_537[k];

        t_672[k] = f_3 * ii0_419[k]
                   - f_4 * ii1_419[k]
                   + pb_y[k] * ik_538[k];

        t_673[k] = pb_y[k] * ik_539[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pa_x, pb_y, pb_z, gl0_674, gl1_674, \
                         hk_360, hk_540, hl_674, hl_675, ik_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_19 * gl0_674[k]
                   - f_20 * gl1_674[k]
                   + pa_x[k] * hl_674[k];

        t_675[k] = f_18 * hk_540[k]
                   + pa_x[k] * hl_675[k];

        t_676[k] = f_17 * hk_360[k]
                   + pb_y[k] * ik_540[k];

        t_677[k] = pb_z[k] * ik_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, t_681, t_682, pa_x, pb_z, hk_543, hk_545, \
                         hk_546, hl_678, hl_680, hl_681, ik_541, \
                         ik_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_0 * hk_543[k]
                   + pa_x[k] * hl_678[k];

        t_679[k] = pb_z[k] * ik_541[k];

        t_680[k] = f_0 * hk_545[k]
                   + pa_x[k] * hl_680[k];

        t_681[k] = f_17 * hk_546[k]
                   + pa_x[k] * hl_681[k];

        t_682[k] = pb_z[k] * ik_543[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, pa_x, pb_y, pb_z, hk_365, hk_549, hk_550, \
                         hl_684, hl_685, ik_545, ik_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_17 * hk_365[k]
                   + pb_y[k] * ik_545[k];

        t_684[k] = f_17 * hk_549[k]
                   + pa_x[k] * hl_684[k];

        t_685[k] = f_16 * hk_550[k]
                   + pa_x[k] * hl_685[k];

        t_686[k] = pb_z[k] * ik_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, pa_x, pb_y, hk_369, hk_552, hk_554, \
                         hk_555, hl_687, hl_689, hl_690, ik_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_16 * hk_552[k]
                   + pa_x[k] * hl_687[k];

        t_688[k] = f_17 * hk_369[k]
                   + pb_y[k] * ik_549[k];

        t_689[k] = f_16 * hk_554[k]
                   + pa_x[k] * hl_689[k];

        t_690[k] = f_15 * hk_555[k]
                   + pa_x[k] * hl_690[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pa_x, pb_y, pb_z, hk_374, hk_557, hk_558, \
                         hl_692, hl_693, ik_550, ik_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = pb_z[k] * ik_550[k];

        t_692[k] = f_15 * hk_557[k]
                   + pa_x[k] * hl_692[k];

        t_693[k] = f_15 * hk_558[k]
                   + pa_x[k] * hl_693[k];

        t_694[k] = f_17 * hk_374[k]
                   + pb_y[k] * ik_554[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, t_699, pa_x, pb_z, hk_560, hk_561, \
                         hk_563, hk_564, hl_695, hl_696, hl_698, hl_699, \
                         ik_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_15 * hk_560[k]
                   + pa_x[k] * hl_695[k];

        t_696[k] = f_14 * hk_561[k]
                   + pa_x[k] * hl_696[k];

        t_697[k] = pb_z[k] * ik_555[k];

        t_698[k] = f_14 * hk_563[k]
                   + pa_x[k] * hl_698[k];

        t_699[k] = f_14 * hk_564[k]
                   + pa_x[k] * hl_699[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pa_x, pb_x, pb_y, hk_380, hk_565, hk_567, \
                         hk_568, hl_700, hl_702, ik_560, ik_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_14 * hk_565[k]
                   + pa_x[k] * hl_700[k];

        t_701[k] = f_17 * hk_380[k]
                   + pb_y[k] * ik_560[k];

        t_702[k] = f_14 * hk_567[k]
                   + pa_x[k] * hl_702[k];

        t_703[k] = f_13 * hk_568[k]
                   + pb_x[k] * ik_568[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, pb_x, pb_z, hk_570, hk_571, \
                         hk_572, hk_573, ik_561, ik_570, ik_571, ik_572, \
                         ik_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = pb_z[k] * ik_561[k];

        t_705[k] = f_13 * hk_570[k]
                   + pb_x[k] * ik_570[k];

        t_706[k] = f_13 * hk_571[k]
                   + pb_x[k] * ik_571[k];

        t_707[k] = f_13 * hk_572[k]
                   + pb_x[k] * ik_572[k];

        t_708[k] = f_13 * hk_573[k]
                   + pb_x[k] * ik_573[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, pa_x, pb_x, pb_z, hk_574, hk_575, \
                         hl_711, hl_713, ik_568, ik_574, ik_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_13 * hk_574[k]
                   + pb_x[k] * ik_574[k];

        t_710[k] = f_13 * hk_575[k]
                   + pb_x[k] * ik_575[k];

        t_711[k] = pa_x[k] * hl_711[k];

        t_712[k] = pb_z[k] * ik_568[k];

        t_713[k] = pa_x[k] * hl_713[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, t_718, t_719, t_720, pa_x, pa_z, hl_450, \
                         hl_714, hl_715, hl_716, hl_717, hl_718, \
                         hl_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pa_x[k] * hl_714[k];

        t_715[k] = pa_x[k] * hl_715[k];

        t_716[k] = pa_x[k] * hl_716[k];

        t_717[k] = pa_x[k] * hl_717[k];

        t_718[k] = pa_x[k] * hl_718[k];

        t_719[k] = pa_x[k] * hl_719[k];

        t_720[k] = pa_z[k] * hl_450[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pa_z, pb_y, pb_z, hk_360, hk_398, hl_451, \
                         hl_453, ik_576, ik_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = pa_z[k] * hl_451[k];

        t_722[k] = f_13 * hk_360[k]
                   + pb_z[k] * ik_576[k];

        t_723[k] = pa_z[k] * hl_453[k];

        t_724[k] = f_16 * hk_398[k]
                   + pb_y[k] * ik_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pa_x, pa_z, pb_y, pb_z, hk_363, hk_401, \
                         hk_581, hl_456, hl_725, ik_579, ik_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_0 * hk_581[k]
                   + pa_x[k] * hl_725[k];

        t_726[k] = pa_z[k] * hl_456[k];

        t_727[k] = f_13 * hk_363[k]
                   + pb_z[k] * ik_579[k];

        t_728[k] = f_16 * hk_401[k]
                   + pb_y[k] * ik_581[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, pa_x, pa_z, pb_z, hk_366, hk_585, hk_588, \
                         hl_460, hl_729, hl_732, ik_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_17 * hk_585[k]
                   + pa_x[k] * hl_729[k];

        t_730[k] = pa_z[k] * hl_460[k];

        t_731[k] = f_13 * hk_366[k]
                   + pb_z[k] * ik_582[k];

        t_732[k] = f_16 * hk_588[k]
                   + pa_x[k] * hl_732[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, pa_x, pa_z, pb_y, pb_z, hk_370, hk_405, \
                         hk_590, hl_465, hl_734, ik_585, ik_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_16 * hk_405[k]
                   + pb_y[k] * ik_585[k];

        t_734[k] = f_16 * hk_590[k]
                   + pa_x[k] * hl_734[k];

        t_735[k] = pa_z[k] * hl_465[k];

        t_736[k] = f_13 * hk_370[k]
                   + pb_z[k] * ik_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, pa_x, pb_y, hk_410, hk_593, hk_594, \
                         hk_596, hl_737, hl_738, hl_740, ik_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_15 * hk_593[k]
                   + pa_x[k] * hl_737[k];

        t_738[k] = f_15 * hk_594[k]
                   + pa_x[k] * hl_738[k];

        t_739[k] = f_16 * hk_410[k]
                   + pb_y[k] * ik_590[k];

        t_740[k] = f_15 * hk_596[k]
                   + pa_x[k] * hl_740[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_x, pa_z, pb_z, hk_375, hk_599, hk_600, \
                         hl_471, hl_743, hl_744, ik_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = pa_z[k] * hl_471[k];

        t_742[k] = f_13 * hk_375[k]
                   + pb_z[k] * ik_591[k];

        t_743[k] = f_14 * hk_599[k]
                   + pa_x[k] * hl_743[k];

        t_744[k] = f_14 * hk_600[k]
                   + pa_x[k] * hl_744[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pa_x, pa_z, pb_y, hk_416, hk_601, hk_603, \
                         hl_478, hl_745, hl_747, ik_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_14 * hk_601[k]
                   + pa_x[k] * hl_745[k];

        t_746[k] = f_16 * hk_416[k]
                   + pb_y[k] * ik_596[k];

        t_747[k] = f_14 * hk_603[k]
                   + pa_x[k] * hl_747[k];

        t_748[k] = pa_z[k] * hl_478[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, t_752, t_753, pb_x, hk_605, hk_606, hk_607, \
                         hk_608, hk_609, ik_605, ik_606, ik_607, ik_608, \
                         ik_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_13 * hk_605[k]
                   + pb_x[k] * ik_605[k];

        t_750[k] = f_13 * hk_606[k]
                   + pb_x[k] * ik_606[k];

        t_751[k] = f_13 * hk_607[k]
                   + pb_x[k] * ik_607[k];

        t_752[k] = f_13 * hk_608[k]
                   + pb_x[k] * ik_608[k];

        t_753[k] = f_13 * hk_609[k]
                   + pb_x[k] * ik_609[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, t_758, t_759, pa_x, pb_x, hk_610, hk_611, \
                         hl_756, hl_757, hl_758, hl_759, ik_610, \
                         ik_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_13 * hk_610[k]
                   + pb_x[k] * ik_610[k];

        t_755[k] = f_13 * hk_611[k]
                   + pb_x[k] * ik_611[k];

        t_756[k] = pa_x[k] * hl_756[k];

        t_757[k] = pa_x[k] * hl_757[k];

        t_758[k] = pa_x[k] * hl_758[k];

        t_759[k] = pa_x[k] * hl_759[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, t_763, t_764, t_765, pa_x, hk_612, hl_760, \
                         hl_761, hl_762, hl_763, hl_764, hl_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = pa_x[k] * hl_760[k];

        t_761[k] = pa_x[k] * hl_761[k];

        t_762[k] = pa_x[k] * hl_762[k];

        t_763[k] = pa_x[k] * hl_763[k];

        t_764[k] = pa_x[k] * hl_764[k];

        t_765[k] = f_18 * hk_612[k]
                   + pa_x[k] * hl_765[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, pa_x, pb_y, pb_z, hk_396, hk_432, hk_434, \
                         hk_615, hl_768, ik_612, ik_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_15 * hk_432[k]
                   + pb_y[k] * ik_612[k];

        t_767[k] = f_14 * hk_396[k]
                   + pb_z[k] * ik_612[k];

        t_768[k] = f_0 * hk_615[k]
                   + pa_x[k] * hl_768[k];

        t_769[k] = f_15 * hk_434[k]
                   + pb_y[k] * ik_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pa_x, pb_y, pb_z, hk_399, hk_437, hk_617, \
                         hk_618, hl_770, hl_771, ik_615, ik_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_0 * hk_617[k]
                   + pa_x[k] * hl_770[k];

        t_771[k] = f_17 * hk_618[k]
                   + pa_x[k] * hl_771[k];

        t_772[k] = f_14 * hk_399[k]
                   + pb_z[k] * ik_615[k];

        t_773[k] = f_15 * hk_437[k]
                   + pb_y[k] * ik_617[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pa_x, pb_z, hk_402, hk_621, hk_622, \
                         hk_624, hl_774, hl_775, hl_777, ik_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_17 * hk_621[k]
                   + pa_x[k] * hl_774[k];

        t_775[k] = f_16 * hk_622[k]
                   + pa_x[k] * hl_775[k];

        t_776[k] = f_14 * hk_402[k]
                   + pb_z[k] * ik_618[k];

        t_777[k] = f_16 * hk_624[k]
                   + pa_x[k] * hl_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, pa_x, pb_y, pb_z, hk_406, hk_441, hk_626, \
                         hk_627, hl_779, hl_780, ik_621, ik_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_15 * hk_441[k]
                   + pb_y[k] * ik_621[k];

        t_779[k] = f_16 * hk_626[k]
                   + pa_x[k] * hl_779[k];

        t_780[k] = f_15 * hk_627[k]
                   + pa_x[k] * hl_780[k];

        t_781[k] = f_14 * hk_406[k]
                   + pb_z[k] * ik_622[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, pa_x, pb_y, hk_446, hk_629, hk_630, \
                         hk_632, hl_782, hl_783, hl_785, ik_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_15 * hk_629[k]
                   + pa_x[k] * hl_782[k];

        t_783[k] = f_15 * hk_630[k]
                   + pa_x[k] * hl_783[k];

        t_784[k] = f_15 * hk_446[k]
                   + pb_y[k] * ik_626[k];

        t_785[k] = f_15 * hk_632[k]
                   + pa_x[k] * hl_785[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pa_x, pb_z, hk_411, hk_633, hk_635, \
                         hk_636, hl_786, hl_788, hl_789, ik_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_14 * hk_633[k]
                   + pa_x[k] * hl_786[k];

        t_787[k] = f_14 * hk_411[k]
                   + pb_z[k] * ik_627[k];

        t_788[k] = f_14 * hk_635[k]
                   + pa_x[k] * hl_788[k];

        t_789[k] = f_14 * hk_636[k]
                   + pa_x[k] * hl_789[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_x, pb_x, pb_y, hk_452, hk_637, hk_639, \
                         hk_640, hl_790, hl_792, ik_632, ik_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_14 * hk_637[k]
                   + pa_x[k] * hl_790[k];

        t_791[k] = f_15 * hk_452[k]
                   + pb_y[k] * ik_632[k];

        t_792[k] = f_14 * hk_639[k]
                   + pa_x[k] * hl_792[k];

        t_793[k] = f_13 * hk_640[k]
                   + pb_x[k] * ik_640[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, t_798, pb_x, hk_641, hk_642, hk_643, \
                         hk_644, hk_645, ik_641, ik_642, ik_643, ik_644, \
                         ik_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_13 * hk_641[k]
                   + pb_x[k] * ik_641[k];

        t_795[k] = f_13 * hk_642[k]
                   + pb_x[k] * ik_642[k];

        t_796[k] = f_13 * hk_643[k]
                   + pb_x[k] * ik_643[k];

        t_797[k] = f_13 * hk_644[k]
                   + pb_x[k] * ik_644[k];

        t_798[k] = f_13 * hk_645[k]
                   + pb_x[k] * ik_645[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, t_803, t_804, pa_x, pb_x, hk_646, hk_647, \
                         hl_801, hl_802, hl_803, hl_804, ik_646, \
                         ik_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_13 * hk_646[k]
                   + pb_x[k] * ik_646[k];

        t_800[k] = f_13 * hk_647[k]
                   + pb_x[k] * ik_647[k];

        t_801[k] = pa_x[k] * hl_801[k];

        t_802[k] = pa_x[k] * hl_802[k];

        t_803[k] = pa_x[k] * hl_803[k];

        t_804[k] = pa_x[k] * hl_804[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, t_809, t_810, pa_x, hk_648, hl_805, \
                         hl_806, hl_807, hl_808, hl_809, hl_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = pa_x[k] * hl_805[k];

        t_806[k] = pa_x[k] * hl_806[k];

        t_807[k] = pa_x[k] * hl_807[k];

        t_808[k] = pa_x[k] * hl_808[k];

        t_809[k] = pa_x[k] * hl_809[k];

        t_810[k] = f_18 * hk_648[k]
                   + pa_x[k] * hl_810[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pa_x, pb_y, pb_z, hk_432, hk_468, hk_470, \
                         hk_651, hl_813, ik_648, ik_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_14 * hk_468[k]
                   + pb_y[k] * ik_648[k];

        t_812[k] = f_15 * hk_432[k]
                   + pb_z[k] * ik_648[k];

        t_813[k] = f_0 * hk_651[k]
                   + pa_x[k] * hl_813[k];

        t_814[k] = f_14 * hk_470[k]
                   + pb_y[k] * ik_650[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t hk,
                                            const size_t hl, const size_t ii0, const size_t ii1,
                                            const size_t ik, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hk_435 = buffer.data(hk + 435);
    const auto *hk_438 = buffer.data(hk + 438);
    const auto *hk_442 = buffer.data(hk + 442);
    const auto *hk_447 = buffer.data(hk + 447);
    const auto *hk_471 = buffer.data(hk + 471);
    const auto *hk_473 = buffer.data(hk + 473);
    const auto *hk_474 = buffer.data(hk + 474);
    const auto *hk_477 = buffer.data(hk + 477);
    const auto *hk_478 = buffer.data(hk + 478);
    const auto *hk_482 = buffer.data(hk + 482);
    const auto *hk_483 = buffer.data(hk + 483);
    const auto *hk_488 = buffer.data(hk + 488);
    const auto *hk_504 = buffer.data(hk + 504);
    const auto *hk_506 = buffer.data(hk + 506);
    const auto *hk_507 = buffer.data(hk + 507);
    const auto *hk_509 = buffer.data(hk + 509);
    const auto *hk_510 = buffer.data(hk + 510);
    const auto *hk_513 = buffer.data(hk + 513);
    const auto *hk_514 = buffer.data(hk + 514);
    const auto *hk_518 = buffer.data(hk + 518);
    const auto *hk_519 = buffer.data(hk + 519);
    const auto *hk_524 = buffer.data(hk + 524);
    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_653 = buffer.data(hk + 653);
    const auto *hk_654 = buffer.data(hk + 654);
    const auto *hk_657 = buffer.data(hk + 657);
    const auto *hk_658 = buffer.data(hk + 658);
    const auto *hk_660 = buffer.data(hk + 660);
    const auto *hk_662 = buffer.data(hk + 662);
    const auto *hk_663 = buffer.data(hk + 663);
    const auto *hk_665 = buffer.data(hk + 665);
    const auto *hk_666 = buffer.data(hk + 666);
    const auto *hk_668 = buffer.data(hk + 668);
    const auto *hk_669 = buffer.data(hk + 669);
    const auto *hk_671 = buffer.data(hk + 671);
    const auto *hk_672 = buffer.data(hk + 672);
    const auto *hk_673 = buffer.data(hk + 673);
    const auto *hk_675 = buffer.data(hk + 675);
    const auto *hk_676 = buffer.data(hk + 676);
    const auto *hk_677 = buffer.data(hk + 677);
    const auto *hk_678 = buffer.data(hk + 678);
    const auto *hk_679 = buffer.data(hk + 679);
    const auto *hk_680 = buffer.data(hk + 680);
    const auto *hk_681 = buffer.data(hk + 681);
    const auto *hk_682 = buffer.data(hk + 682);
    const auto *hk_683 = buffer.data(hk + 683);
    const auto *hk_687 = buffer.data(hk + 687);
    const auto *hk_690 = buffer.data(hk + 690);
    const auto *hk_694 = buffer.data(hk + 694);
    const auto *hk_696 = buffer.data(hk + 696);
    const auto *hk_699 = buffer.data(hk + 699);
    const auto *hk_701 = buffer.data(hk + 701);
    const auto *hk_702 = buffer.data(hk + 702);
    const auto *hk_705 = buffer.data(hk + 705);
    const auto *hk_707 = buffer.data(hk + 707);
    const auto *hk_708 = buffer.data(hk + 708);
    const auto *hk_709 = buffer.data(hk + 709);
    const auto *hk_712 = buffer.data(hk + 712);
    const auto *hk_713 = buffer.data(hk + 713);
    const auto *hk_714 = buffer.data(hk + 714);
    const auto *hk_715 = buffer.data(hk + 715);
    const auto *hk_716 = buffer.data(hk + 716);
    const auto *hk_717 = buffer.data(hk + 717);
    const auto *hk_718 = buffer.data(hk + 718);
    const auto *hk_720 = buffer.data(hk + 720);
    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_725 = buffer.data(hk + 725);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_729 = buffer.data(hk + 729);
    const auto *hk_730 = buffer.data(hk + 730);
    const auto *hk_732 = buffer.data(hk + 732);
    const auto *hk_734 = buffer.data(hk + 734);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_737 = buffer.data(hk + 737);
    const auto *hk_738 = buffer.data(hk + 738);
    const auto *hk_740 = buffer.data(hk + 740);
    const auto *hk_741 = buffer.data(hk + 741);
    const auto *hk_743 = buffer.data(hk + 743);
    const auto *hk_744 = buffer.data(hk + 744);
    const auto *hk_745 = buffer.data(hk + 745);
    const auto *hk_747 = buffer.data(hk + 747);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_749 = buffer.data(hk + 749);
    const auto *hk_750 = buffer.data(hk + 750);
    const auto *hk_751 = buffer.data(hk + 751);
    const auto *hk_752 = buffer.data(hk + 752);
    const auto *hk_753 = buffer.data(hk + 753);
    const auto *hk_755 = buffer.data(hk + 755);

    const auto *hl_630 = buffer.data(hl + 630);
    const auto *hl_632 = buffer.data(hl + 632);
    const auto *hl_635 = buffer.data(hl + 635);
    const auto *hl_639 = buffer.data(hl + 639);
    const auto *hl_644 = buffer.data(hl + 644);
    const auto *hl_650 = buffer.data(hl + 650);
    const auto *hl_657 = buffer.data(hl + 657);
    const auto *hl_665 = buffer.data(hl + 665);
    const auto *hl_815 = buffer.data(hl + 815);
    const auto *hl_816 = buffer.data(hl + 816);
    const auto *hl_819 = buffer.data(hl + 819);
    const auto *hl_820 = buffer.data(hl + 820);
    const auto *hl_822 = buffer.data(hl + 822);
    const auto *hl_824 = buffer.data(hl + 824);
    const auto *hl_825 = buffer.data(hl + 825);
    const auto *hl_827 = buffer.data(hl + 827);
    const auto *hl_828 = buffer.data(hl + 828);
    const auto *hl_830 = buffer.data(hl + 830);
    const auto *hl_831 = buffer.data(hl + 831);
    const auto *hl_833 = buffer.data(hl + 833);
    const auto *hl_834 = buffer.data(hl + 834);
    const auto *hl_835 = buffer.data(hl + 835);
    const auto *hl_837 = buffer.data(hl + 837);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_847 = buffer.data(hl + 847);
    const auto *hl_848 = buffer.data(hl + 848);
    const auto *hl_849 = buffer.data(hl + 849);
    const auto *hl_850 = buffer.data(hl + 850);
    const auto *hl_851 = buffer.data(hl + 851);
    const auto *hl_852 = buffer.data(hl + 852);
    const auto *hl_853 = buffer.data(hl + 853);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_858 = buffer.data(hl + 858);
    const auto *hl_861 = buffer.data(hl + 861);
    const auto *hl_865 = buffer.data(hl + 865);
    const auto *hl_867 = buffer.data(hl + 867);
    const auto *hl_870 = buffer.data(hl + 870);
    const auto *hl_872 = buffer.data(hl + 872);
    const auto *hl_873 = buffer.data(hl + 873);
    const auto *hl_876 = buffer.data(hl + 876);
    const auto *hl_878 = buffer.data(hl + 878);
    const auto *hl_879 = buffer.data(hl + 879);
    const auto *hl_880 = buffer.data(hl + 880);
    const auto *hl_891 = buffer.data(hl + 891);
    const auto *hl_892 = buffer.data(hl + 892);
    const auto *hl_893 = buffer.data(hl + 893);
    const auto *hl_894 = buffer.data(hl + 894);
    const auto *hl_895 = buffer.data(hl + 895);
    const auto *hl_896 = buffer.data(hl + 896);
    const auto *hl_897 = buffer.data(hl + 897);
    const auto *hl_898 = buffer.data(hl + 898);
    const auto *hl_899 = buffer.data(hl + 899);
    const auto *hl_900 = buffer.data(hl + 900);
    const auto *hl_903 = buffer.data(hl + 903);
    const auto *hl_905 = buffer.data(hl + 905);
    const auto *hl_906 = buffer.data(hl + 906);
    const auto *hl_909 = buffer.data(hl + 909);
    const auto *hl_910 = buffer.data(hl + 910);
    const auto *hl_912 = buffer.data(hl + 912);
    const auto *hl_914 = buffer.data(hl + 914);
    const auto *hl_915 = buffer.data(hl + 915);
    const auto *hl_917 = buffer.data(hl + 917);
    const auto *hl_918 = buffer.data(hl + 918);
    const auto *hl_920 = buffer.data(hl + 920);
    const auto *hl_921 = buffer.data(hl + 921);
    const auto *hl_923 = buffer.data(hl + 923);
    const auto *hl_924 = buffer.data(hl + 924);
    const auto *hl_925 = buffer.data(hl + 925);
    const auto *hl_927 = buffer.data(hl + 927);
    const auto *hl_936 = buffer.data(hl + 936);
    const auto *hl_937 = buffer.data(hl + 937);
    const auto *hl_938 = buffer.data(hl + 938);
    const auto *hl_939 = buffer.data(hl + 939);
    const auto *hl_940 = buffer.data(hl + 940);
    const auto *hl_941 = buffer.data(hl + 941);
    const auto *hl_942 = buffer.data(hl + 942);
    const auto *hl_944 = buffer.data(hl + 944);

    const auto *ii0_588 = buffer.data(ii0 + 588);
    const auto *ii0_591 = buffer.data(ii0 + 591);
    const auto *ii0_593 = buffer.data(ii0 + 593);
    const auto *ii0_594 = buffer.data(ii0 + 594);
    const auto *ii0_597 = buffer.data(ii0 + 597);
    const auto *ii0_598 = buffer.data(ii0 + 598);
    const auto *ii0_600 = buffer.data(ii0 + 600);
    const auto *ii0_602 = buffer.data(ii0 + 602);
    const auto *ii0_603 = buffer.data(ii0 + 603);

    const auto *ii1_588 = buffer.data(ii1 + 588);
    const auto *ii1_591 = buffer.data(ii1 + 591);
    const auto *ii1_593 = buffer.data(ii1 + 593);
    const auto *ii1_594 = buffer.data(ii1 + 594);
    const auto *ii1_597 = buffer.data(ii1 + 597);
    const auto *ii1_598 = buffer.data(ii1 + 598);
    const auto *ii1_600 = buffer.data(ii1 + 600);
    const auto *ii1_602 = buffer.data(ii1 + 602);
    const auto *ii1_603 = buffer.data(ii1 + 603);

    const auto *ik_651 = buffer.data(ik + 651);
    const auto *ik_653 = buffer.data(ik + 653);
    const auto *ik_654 = buffer.data(ik + 654);
    const auto *ik_657 = buffer.data(ik + 657);
    const auto *ik_658 = buffer.data(ik + 658);
    const auto *ik_662 = buffer.data(ik + 662);
    const auto *ik_663 = buffer.data(ik + 663);
    const auto *ik_668 = buffer.data(ik + 668);
    const auto *ik_676 = buffer.data(ik + 676);
    const auto *ik_677 = buffer.data(ik + 677);
    const auto *ik_678 = buffer.data(ik + 678);
    const auto *ik_679 = buffer.data(ik + 679);
    const auto *ik_680 = buffer.data(ik + 680);
    const auto *ik_681 = buffer.data(ik + 681);
    const auto *ik_682 = buffer.data(ik + 682);
    const auto *ik_683 = buffer.data(ik + 683);
    const auto *ik_684 = buffer.data(ik + 684);
    const auto *ik_686 = buffer.data(ik + 686);
    const auto *ik_687 = buffer.data(ik + 687);
    const auto *ik_689 = buffer.data(ik + 689);
    const auto *ik_690 = buffer.data(ik + 690);
    const auto *ik_693 = buffer.data(ik + 693);
    const auto *ik_694 = buffer.data(ik + 694);
    const auto *ik_698 = buffer.data(ik + 698);
    const auto *ik_699 = buffer.data(ik + 699);
    const auto *ik_704 = buffer.data(ik + 704);
    const auto *ik_712 = buffer.data(ik + 712);
    const auto *ik_713 = buffer.data(ik + 713);
    const auto *ik_714 = buffer.data(ik + 714);
    const auto *ik_715 = buffer.data(ik + 715);
    const auto *ik_716 = buffer.data(ik + 716);
    const auto *ik_717 = buffer.data(ik + 717);
    const auto *ik_718 = buffer.data(ik + 718);
    const auto *ik_720 = buffer.data(ik + 720);
    const auto *ik_722 = buffer.data(ik + 722);
    const auto *ik_723 = buffer.data(ik + 723);
    const auto *ik_725 = buffer.data(ik + 725);
    const auto *ik_726 = buffer.data(ik + 726);
    const auto *ik_729 = buffer.data(ik + 729);
    const auto *ik_730 = buffer.data(ik + 730);
    const auto *ik_734 = buffer.data(ik + 734);
    const auto *ik_735 = buffer.data(ik + 735);
    const auto *ik_740 = buffer.data(ik + 740);
    const auto *ik_747 = buffer.data(ik + 747);
    const auto *ik_748 = buffer.data(ik + 748);
    const auto *ik_749 = buffer.data(ik + 749);
    const auto *ik_750 = buffer.data(ik + 750);
    const auto *ik_751 = buffer.data(ik + 751);
    const auto *ik_752 = buffer.data(ik + 752);
    const auto *ik_753 = buffer.data(ik + 753);
    const auto *ik_755 = buffer.data(ik + 755);
    const auto *ik_756 = buffer.data(ik + 756);
    const auto *ik_757 = buffer.data(ik + 757);
    const auto *ik_759 = buffer.data(ik + 759);
    const auto *ik_761 = buffer.data(ik + 761);
    const auto *ik_762 = buffer.data(ik + 762);
    const auto *ik_765 = buffer.data(ik + 765);
    const auto *ik_766 = buffer.data(ik + 766);
    const auto *ik_768 = buffer.data(ik + 768);
    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pa_x, pb_y, pb_z, hk_435, hk_473, hk_653, \
                         hk_654, hl_815, hl_816, ik_651, ik_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_0 * hk_653[k]
                   + pa_x[k] * hl_815[k];

        t_816[k] = f_17 * hk_654[k]
                   + pa_x[k] * hl_816[k];

        t_817[k] = f_15 * hk_435[k]
                   + pb_z[k] * ik_651[k];

        t_818[k] = f_14 * hk_473[k]
                   + pb_y[k] * ik_653[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, pa_x, pb_z, hk_438, hk_657, hk_658, \
                         hk_660, hl_819, hl_820, hl_822, ik_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_17 * hk_657[k]
                   + pa_x[k] * hl_819[k];

        t_820[k] = f_16 * hk_658[k]
                   + pa_x[k] * hl_820[k];

        t_821[k] = f_15 * hk_438[k]
                   + pb_z[k] * ik_654[k];

        t_822[k] = f_16 * hk_660[k]
                   + pa_x[k] * hl_822[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, t_826, pa_x, pb_y, pb_z, hk_442, hk_477, hk_662, \
                         hk_663, hl_824, hl_825, ik_657, ik_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_14 * hk_477[k]
                   + pb_y[k] * ik_657[k];

        t_824[k] = f_16 * hk_662[k]
                   + pa_x[k] * hl_824[k];

        t_825[k] = f_15 * hk_663[k]
                   + pa_x[k] * hl_825[k];

        t_826[k] = f_15 * hk_442[k]
                   + pb_z[k] * ik_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, pa_x, pb_y, hk_482, hk_665, hk_666, \
                         hk_668, hl_827, hl_828, hl_830, ik_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_15 * hk_665[k]
                   + pa_x[k] * hl_827[k];

        t_828[k] = f_15 * hk_666[k]
                   + pa_x[k] * hl_828[k];

        t_829[k] = f_14 * hk_482[k]
                   + pb_y[k] * ik_662[k];

        t_830[k] = f_15 * hk_668[k]
                   + pa_x[k] * hl_830[k];
    }

#pragma omp simd aligned(t_831, t_832, t_833, t_834, pa_x, pb_z, hk_447, hk_669, hk_671, \
                         hk_672, hl_831, hl_833, hl_834, ik_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_831[k] = f_14 * hk_669[k]
                   + pa_x[k] * hl_831[k];

        t_832[k] = f_15 * hk_447[k]
                   + pb_z[k] * ik_663[k];

        t_833[k] = f_14 * hk_671[k]
                   + pa_x[k] * hl_833[k];

        t_834[k] = f_14 * hk_672[k]
                   + pa_x[k] * hl_834[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, t_838, pa_x, pb_x, pb_y, hk_488, hk_673, hk_675, \
                         hk_676, hl_835, hl_837, ik_668, ik_676 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_14 * hk_673[k]
                   + pa_x[k] * hl_835[k];

        t_836[k] = f_14 * hk_488[k]
                   + pb_y[k] * ik_668[k];

        t_837[k] = f_14 * hk_675[k]
                   + pa_x[k] * hl_837[k];

        t_838[k] = f_13 * hk_676[k]
                   + pb_x[k] * ik_676[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, t_842, t_843, pb_x, hk_677, hk_678, hk_679, \
                         hk_680, hk_681, ik_677, ik_678, ik_679, ik_680, \
                         ik_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_13 * hk_677[k]
                   + pb_x[k] * ik_677[k];

        t_840[k] = f_13 * hk_678[k]
                   + pb_x[k] * ik_678[k];

        t_841[k] = f_13 * hk_679[k]
                   + pb_x[k] * ik_679[k];

        t_842[k] = f_13 * hk_680[k]
                   + pb_x[k] * ik_680[k];

        t_843[k] = f_13 * hk_681[k]
                   + pb_x[k] * ik_681[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, t_847, t_848, t_849, pa_x, pb_x, hk_682, hk_683, \
                         hl_846, hl_847, hl_848, hl_849, ik_682, \
                         ik_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = f_13 * hk_682[k]
                   + pb_x[k] * ik_682[k];

        t_845[k] = f_13 * hk_683[k]
                   + pb_x[k] * ik_683[k];

        t_846[k] = pa_x[k] * hl_846[k];

        t_847[k] = pa_x[k] * hl_847[k];

        t_848[k] = pa_x[k] * hl_848[k];

        t_849[k] = pa_x[k] * hl_849[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, t_853, t_854, t_855, pa_x, pa_y, hl_630, hl_850, \
                         hl_851, hl_852, hl_853, hl_854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = pa_x[k] * hl_850[k];

        t_851[k] = pa_x[k] * hl_851[k];

        t_852[k] = pa_x[k] * hl_852[k];

        t_853[k] = pa_x[k] * hl_853[k];

        t_854[k] = pa_x[k] * hl_854[k];

        t_855[k] = pa_y[k] * hl_630[k];
    }

#pragma omp simd aligned(t_856, t_857, t_858, t_859, t_860, pa_x, pa_y, pb_y, hk_504, hk_506, \
                         hk_687, hl_632, hl_635, hl_858, ik_684, \
                         ik_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_13 * hk_504[k]
                   + pb_y[k] * ik_684[k];

        t_857[k] = pa_y[k] * hl_632[k];

        t_858[k] = f_0 * hk_687[k]
                   + pa_x[k] * hl_858[k];

        t_859[k] = f_13 * hk_506[k]
                   + pb_y[k] * ik_686[k];

        t_860[k] = pa_y[k] * hl_635[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, t_864, pa_x, pa_y, pb_y, pb_z, hk_471, hk_509, \
                         hk_690, hl_639, hl_861, ik_687, ik_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_17 * hk_690[k]
                   + pa_x[k] * hl_861[k];

        t_862[k] = f_16 * hk_471[k]
                   + pb_z[k] * ik_687[k];

        t_863[k] = f_13 * hk_509[k]
                   + pb_y[k] * ik_689[k];

        t_864[k] = pa_y[k] * hl_639[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, t_868, pa_x, pb_y, pb_z, hk_474, hk_513, hk_694, \
                         hk_696, hl_865, hl_867, ik_690, ik_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_16 * hk_694[k]
                   + pa_x[k] * hl_865[k];

        t_866[k] = f_16 * hk_474[k]
                   + pb_z[k] * ik_690[k];

        t_867[k] = f_16 * hk_696[k]
                   + pa_x[k] * hl_867[k];

        t_868[k] = f_13 * hk_513[k]
                   + pb_y[k] * ik_693[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, pa_x, pa_y, pb_z, hk_478, hk_699, hk_701, \
                         hl_644, hl_870, hl_872, ik_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = pa_y[k] * hl_644[k];

        t_870[k] = f_15 * hk_699[k]
                   + pa_x[k] * hl_870[k];

        t_871[k] = f_16 * hk_478[k]
                   + pb_z[k] * ik_694[k];

        t_872[k] = f_15 * hk_701[k]
                   + pa_x[k] * hl_872[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pa_x, pa_y, pb_y, hk_518, hk_702, hk_705, \
                         hl_650, hl_873, hl_876, ik_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = f_15 * hk_702[k]
                   + pa_x[k] * hl_873[k];

        t_874[k] = f_13 * hk_518[k]
                   + pb_y[k] * ik_698[k];

        t_875[k] = pa_y[k] * hl_650[k];

        t_876[k] = f_14 * hk_705[k]
                   + pa_x[k] * hl_876[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, t_880, pa_x, pb_z, hk_483, hk_707, hk_708, \
                         hk_709, hl_878, hl_879, hl_880, ik_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_16 * hk_483[k]
                   + pb_z[k] * ik_699[k];

        t_878[k] = f_14 * hk_707[k]
                   + pa_x[k] * hl_878[k];

        t_879[k] = f_14 * hk_708[k]
                   + pa_x[k] * hl_879[k];

        t_880[k] = f_14 * hk_709[k]
                   + pa_x[k] * hl_880[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, t_884, pa_y, pb_x, pb_y, hk_524, hk_712, hk_713, \
                         hl_657, ik_704, ik_712, ik_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_13 * hk_524[k]
                   + pb_y[k] * ik_704[k];

        t_882[k] = pa_y[k] * hl_657[k];

        t_883[k] = f_13 * hk_712[k]
                   + pb_x[k] * ik_712[k];

        t_884[k] = f_13 * hk_713[k]
                   + pb_x[k] * ik_713[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, pb_x, hk_714, hk_715, hk_716, \
                         hk_717, hk_718, ik_714, ik_715, ik_716, ik_717, \
                         ik_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = f_13 * hk_714[k]
                   + pb_x[k] * ik_714[k];

        t_886[k] = f_13 * hk_715[k]
                   + pb_x[k] * ik_715[k];

        t_887[k] = f_13 * hk_716[k]
                   + pb_x[k] * ik_716[k];

        t_888[k] = f_13 * hk_717[k]
                   + pb_x[k] * ik_717[k];

        t_889[k] = f_13 * hk_718[k]
                   + pb_x[k] * ik_718[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, t_895, t_896, pa_x, pa_y, hl_665, \
                         hl_891, hl_892, hl_893, hl_894, hl_895, \
                         hl_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = pa_y[k] * hl_665[k];

        t_891[k] = pa_x[k] * hl_891[k];

        t_892[k] = pa_x[k] * hl_892[k];

        t_893[k] = pa_x[k] * hl_893[k];

        t_894[k] = pa_x[k] * hl_894[k];

        t_895[k] = pa_x[k] * hl_895[k];

        t_896[k] = pa_x[k] * hl_896[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, t_900, t_901, t_902, pa_x, pb_y, pb_z, hk_504, \
                         hk_720, hl_897, hl_898, hl_899, hl_900, \
                         ik_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = pa_x[k] * hl_897[k];

        t_898[k] = pa_x[k] * hl_898[k];

        t_899[k] = pa_x[k] * hl_899[k];

        t_900[k] = f_18 * hk_720[k]
                   + pa_x[k] * hl_900[k];

        t_901[k] = pb_y[k] * ik_720[k];

        t_902[k] = f_17 * hk_504[k]
                   + pb_z[k] * ik_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, pa_x, pb_y, hk_723, hk_725, hk_726, \
                         hl_903, hl_905, hl_906, ik_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_0 * hk_723[k]
                   + pa_x[k] * hl_903[k];

        t_904[k] = pb_y[k] * ik_722[k];

        t_905[k] = f_0 * hk_725[k]
                   + pa_x[k] * hl_905[k];

        t_906[k] = f_17 * hk_726[k]
                   + pa_x[k] * hl_906[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, t_910, pa_x, pb_y, pb_z, hk_507, hk_729, hk_730, \
                         hl_909, hl_910, ik_723, ik_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_17 * hk_507[k]
                   + pb_z[k] * ik_723[k];

        t_908[k] = pb_y[k] * ik_725[k];

        t_909[k] = f_17 * hk_729[k]
                   + pa_x[k] * hl_909[k];

        t_910[k] = f_16 * hk_730[k]
                   + pa_x[k] * hl_910[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, pa_x, pb_y, pb_z, hk_510, hk_732, hk_734, \
                         hl_912, hl_914, ik_726, ik_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_17 * hk_510[k]
                   + pb_z[k] * ik_726[k];

        t_912[k] = f_16 * hk_732[k]
                   + pa_x[k] * hl_912[k];

        t_913[k] = pb_y[k] * ik_729[k];

        t_914[k] = f_16 * hk_734[k]
                   + pa_x[k] * hl_914[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, t_918, pa_x, pb_z, hk_514, hk_735, hk_737, \
                         hk_738, hl_915, hl_917, hl_918, ik_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_15 * hk_735[k]
                   + pa_x[k] * hl_915[k];

        t_916[k] = f_17 * hk_514[k]
                   + pb_z[k] * ik_730[k];

        t_917[k] = f_15 * hk_737[k]
                   + pa_x[k] * hl_917[k];

        t_918[k] = f_15 * hk_738[k]
                   + pa_x[k] * hl_918[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, pa_x, pb_y, pb_z, hk_519, hk_740, hk_741, \
                         hl_920, hl_921, ik_734, ik_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = pb_y[k] * ik_734[k];

        t_920[k] = f_15 * hk_740[k]
                   + pa_x[k] * hl_920[k];

        t_921[k] = f_14 * hk_741[k]
                   + pa_x[k] * hl_921[k];

        t_922[k] = f_17 * hk_519[k]
                   + pb_z[k] * ik_735[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, t_927, pa_x, pb_y, hk_743, hk_744, \
                         hk_745, hk_747, hl_923, hl_924, hl_925, hl_927, \
                         ik_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = f_14 * hk_743[k]
                   + pa_x[k] * hl_923[k];

        t_924[k] = f_14 * hk_744[k]
                   + pa_x[k] * hl_924[k];

        t_925[k] = f_14 * hk_745[k]
                   + pa_x[k] * hl_925[k];

        t_926[k] = pb_y[k] * ik_740[k];

        t_927[k] = f_14 * hk_747[k]
                   + pa_x[k] * hl_927[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, t_932, pb_x, hk_748, hk_749, hk_750, \
                         hk_751, hk_752, ik_748, ik_749, ik_750, ik_751, \
                         ik_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_13 * hk_748[k]
                   + pb_x[k] * ik_748[k];

        t_929[k] = f_13 * hk_749[k]
                   + pb_x[k] * ik_749[k];

        t_930[k] = f_13 * hk_750[k]
                   + pb_x[k] * ik_750[k];

        t_931[k] = f_13 * hk_751[k]
                   + pb_x[k] * ik_751[k];

        t_932[k] = f_13 * hk_752[k]
                   + pb_x[k] * ik_752[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, pa_x, pb_x, pb_y, hk_753, hk_755, \
                         hl_936, hl_937, ik_747, ik_753, ik_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_13 * hk_753[k]
                   + pb_x[k] * ik_753[k];

        t_934[k] = pb_y[k] * ik_747[k];

        t_935[k] = f_13 * hk_755[k]
                   + pb_x[k] * ik_755[k];

        t_936[k] = pa_x[k] * hl_936[k];

        t_937[k] = pa_x[k] * hl_937[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, t_942, t_943, t_944, pa_x, pb_y, hl_938, \
                         hl_939, hl_940, hl_941, hl_942, hl_944, \
                         ik_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = pa_x[k] * hl_938[k];

        t_939[k] = pa_x[k] * hl_939[k];

        t_940[k] = pa_x[k] * hl_940[k];

        t_941[k] = pa_x[k] * hl_941[k];

        t_942[k] = pa_x[k] * hl_942[k];

        t_943[k] = pb_y[k] * ik_755[k];

        t_944[k] = pa_x[k] * hl_944[k];
    }

#pragma omp simd aligned(t_945, t_946, t_947, t_948, t_949, pb_x, pb_y, pb_z, hk_540, ii0_588, \
                         ii0_591, ii1_588, ii1_591, ik_756, ik_757, \
                         ik_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_945[k] = f_1 * ii0_588[k]
                   - f_2 * ii1_588[k]
                   + pb_x[k] * ik_756[k];

        t_946[k] = f_0 * hk_540[k]
                   + pb_y[k] * ik_756[k];

        t_947[k] = pb_z[k] * ik_756[k];

        t_948[k] = f_11 * ii0_591[k]
                   - f_12 * ii1_591[k]
                   + pb_x[k] * ik_759[k];

        t_949[k] = pb_z[k] * ik_757[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pb_x, pb_y, pb_z, hk_545, ii0_593, \
                         ii0_594, ii1_593, ii1_594, ik_759, ik_761, \
                         ik_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_11 * ii0_593[k]
                   - f_12 * ii1_593[k]
                   + pb_x[k] * ik_761[k];

        t_951[k] = f_9 * ii0_594[k]
                   - f_10 * ii1_594[k]
                   + pb_x[k] * ik_762[k];

        t_952[k] = pb_z[k] * ik_759[k];

        t_953[k] = f_0 * hk_545[k]
                   + pb_y[k] * ik_761[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, t_957, pb_x, pb_z, ii0_597, ii0_598, ii0_600, \
                         ii1_597, ii1_598, ii1_600, ik_762, ik_765, ik_766, \
                         ik_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_9 * ii0_597[k]
                   - f_10 * ii1_597[k]
                   + pb_x[k] * ik_765[k];

        t_955[k] = f_7 * ii0_598[k]
                   - f_8 * ii1_598[k]
                   + pb_x[k] * ik_766[k];

        t_956[k] = pb_z[k] * ik_762[k];

        t_957[k] = f_7 * ii0_600[k]
                   - f_8 * ii1_600[k]
                   + pb_x[k] * ik_768[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, pb_x, pb_y, pb_z, hk_549, ii0_602, \
                         ii0_603, ii1_602, ii1_603, ik_765, ik_766, ik_770, \
                         ik_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_0 * hk_549[k]
                   + pb_y[k] * ik_765[k];

        t_959[k] = f_7 * ii0_602[k]
                   - f_8 * ii1_602[k]
                   + pb_x[k] * ik_770[k];

        t_960[k] = f_5 * ii0_603[k]
                   - f_6 * ii1_603[k]
                   + pb_x[k] * ik_771[k];

        t_961[k] = pb_z[k] * ik_766[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_486 = buffer.data(gl0 + 486);
    const auto *gl0_584 = buffer.data(gl0 + 584);

    const auto *gl1_486 = buffer.data(gl1 + 486);
    const auto *gl1_584 = buffer.data(gl1 + 584);

    const auto *hk_540 = buffer.data(hk + 540);
    const auto *hk_542 = buffer.data(hk + 542);
    const auto *hk_543 = buffer.data(hk + 543);
    const auto *hk_545 = buffer.data(hk + 545);
    const auto *hk_546 = buffer.data(hk + 546);
    const auto *hk_547 = buffer.data(hk + 547);
    const auto *hk_549 = buffer.data(hk + 549);
    const auto *hk_550 = buffer.data(hk + 550);
    const auto *hk_551 = buffer.data(hk + 551);
    const auto *hk_552 = buffer.data(hk + 552);
    const auto *hk_554 = buffer.data(hk + 554);
    const auto *hk_555 = buffer.data(hk + 555);
    const auto *hk_556 = buffer.data(hk + 556);
    const auto *hk_557 = buffer.data(hk + 557);
    const auto *hk_558 = buffer.data(hk + 558);
    const auto *hk_560 = buffer.data(hk + 560);
    const auto *hk_568 = buffer.data(hk + 568);
    const auto *hk_569 = buffer.data(hk + 569);
    const auto *hk_570 = buffer.data(hk + 570);
    const auto *hk_571 = buffer.data(hk + 571);
    const auto *hk_572 = buffer.data(hk + 572);
    const auto *hk_573 = buffer.data(hk + 573);
    const auto *hk_575 = buffer.data(hk + 575);
    const auto *hk_576 = buffer.data(hk + 576);
    const auto *hk_578 = buffer.data(hk + 578);
    const auto *hk_579 = buffer.data(hk + 579);
    const auto *hk_581 = buffer.data(hk + 581);
    const auto *hk_582 = buffer.data(hk + 582);
    const auto *hk_585 = buffer.data(hk + 585);
    const auto *hk_586 = buffer.data(hk + 586);
    const auto *hk_590 = buffer.data(hk + 590);
    const auto *hk_591 = buffer.data(hk + 591);
    const auto *hk_596 = buffer.data(hk + 596);
    const auto *hk_604 = buffer.data(hk + 604);
    const auto *hk_611 = buffer.data(hk + 611);
    const auto *hk_612 = buffer.data(hk + 612);
    const auto *hk_614 = buffer.data(hk + 614);
    const auto *hk_615 = buffer.data(hk + 615);
    const auto *hk_617 = buffer.data(hk + 617);
    const auto *hk_618 = buffer.data(hk + 618);
    const auto *hk_621 = buffer.data(hk + 621);
    const auto *hk_626 = buffer.data(hk + 626);
    const auto *hk_632 = buffer.data(hk + 632);
    const auto *hk_642 = buffer.data(hk + 642);
    const auto *hk_643 = buffer.data(hk + 643);
    const auto *hk_644 = buffer.data(hk + 644);
    const auto *hk_645 = buffer.data(hk + 645);
    const auto *hk_646 = buffer.data(hk + 646);
    const auto *hk_647 = buffer.data(hk + 647);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_650 = buffer.data(hk + 650);
    const auto *hk_653 = buffer.data(hk + 653);

    const auto *hl_675 = buffer.data(hl + 675);
    const auto *hl_676 = buffer.data(hl + 676);
    const auto *hl_678 = buffer.data(hl + 678);
    const auto *hl_680 = buffer.data(hl + 680);
    const auto *hl_681 = buffer.data(hl + 681);
    const auto *hl_684 = buffer.data(hl + 684);
    const auto *hl_685 = buffer.data(hl + 685);
    const auto *hl_687 = buffer.data(hl + 687);
    const auto *hl_689 = buffer.data(hl + 689);
    const auto *hl_690 = buffer.data(hl + 690);
    const auto *hl_692 = buffer.data(hl + 692);
    const auto *hl_693 = buffer.data(hl + 693);
    const auto *hl_695 = buffer.data(hl + 695);
    const auto *hl_696 = buffer.data(hl + 696);
    const auto *hl_698 = buffer.data(hl + 698);
    const auto *hl_699 = buffer.data(hl + 699);
    const auto *hl_700 = buffer.data(hl + 700);
    const auto *hl_702 = buffer.data(hl + 702);
    const auto *hl_711 = buffer.data(hl + 711);
    const auto *hl_713 = buffer.data(hl + 713);
    const auto *hl_714 = buffer.data(hl + 714);
    const auto *hl_715 = buffer.data(hl + 715);
    const auto *hl_716 = buffer.data(hl + 716);
    const auto *hl_717 = buffer.data(hl + 717);
    const auto *hl_719 = buffer.data(hl + 719);
    const auto *hl_756 = buffer.data(hl + 756);
    const auto *hl_809 = buffer.data(hl + 809);

    const auto *ii0_605 = buffer.data(ii0 + 605);
    const auto *ii0_606 = buffer.data(ii0 + 606);
    const auto *ii0_608 = buffer.data(ii0 + 608);
    const auto *ii0_609 = buffer.data(ii0 + 609);
    const auto *ii0_610 = buffer.data(ii0 + 610);
    const auto *ii0_611 = buffer.data(ii0 + 611);
    const auto *ii0_612 = buffer.data(ii0 + 612);
    const auto *ii0_613 = buffer.data(ii0 + 613);
    const auto *ii0_615 = buffer.data(ii0 + 615);
    const auto *ii0_644 = buffer.data(ii0 + 644);
    const auto *ii0_647 = buffer.data(ii0 + 647);
    const auto *ii0_649 = buffer.data(ii0 + 649);
    const auto *ii0_650 = buffer.data(ii0 + 650);
    const auto *ii0_653 = buffer.data(ii0 + 653);
    const auto *ii0_654 = buffer.data(ii0 + 654);
    const auto *ii0_656 = buffer.data(ii0 + 656);
    const auto *ii0_658 = buffer.data(ii0 + 658);
    const auto *ii0_659 = buffer.data(ii0 + 659);
    const auto *ii0_661 = buffer.data(ii0 + 661);
    const auto *ii0_662 = buffer.data(ii0 + 662);
    const auto *ii0_664 = buffer.data(ii0 + 664);
    const auto *ii0_665 = buffer.data(ii0 + 665);
    const auto *ii0_667 = buffer.data(ii0 + 667);
    const auto *ii0_668 = buffer.data(ii0 + 668);
    const auto *ii0_669 = buffer.data(ii0 + 669);
    const auto *ii0_670 = buffer.data(ii0 + 670);
    const auto *ii0_671 = buffer.data(ii0 + 671);
    const auto *ii0_672 = buffer.data(ii0 + 672);
    const auto *ii0_675 = buffer.data(ii0 + 675);
    const auto *ii0_677 = buffer.data(ii0 + 677);
    const auto *ii0_678 = buffer.data(ii0 + 678);
    const auto *ii0_681 = buffer.data(ii0 + 681);
    const auto *ii0_682 = buffer.data(ii0 + 682);
    const auto *ii0_684 = buffer.data(ii0 + 684);

    const auto *ii1_605 = buffer.data(ii1 + 605);
    const auto *ii1_606 = buffer.data(ii1 + 606);
    const auto *ii1_608 = buffer.data(ii1 + 608);
    const auto *ii1_609 = buffer.data(ii1 + 609);
    const auto *ii1_610 = buffer.data(ii1 + 610);
    const auto *ii1_611 = buffer.data(ii1 + 611);
    const auto *ii1_612 = buffer.data(ii1 + 612);
    const auto *ii1_613 = buffer.data(ii1 + 613);
    const auto *ii1_615 = buffer.data(ii1 + 615);
    const auto *ii1_644 = buffer.data(ii1 + 644);
    const auto *ii1_647 = buffer.data(ii1 + 647);
    const auto *ii1_649 = buffer.data(ii1 + 649);
    const auto *ii1_650 = buffer.data(ii1 + 650);
    const auto *ii1_653 = buffer.data(ii1 + 653);
    const auto *ii1_654 = buffer.data(ii1 + 654);
    const auto *ii1_656 = buffer.data(ii1 + 656);
    const auto *ii1_658 = buffer.data(ii1 + 658);
    const auto *ii1_659 = buffer.data(ii1 + 659);
    const auto *ii1_661 = buffer.data(ii1 + 661);
    const auto *ii1_662 = buffer.data(ii1 + 662);
    const auto *ii1_664 = buffer.data(ii1 + 664);
    const auto *ii1_665 = buffer.data(ii1 + 665);
    const auto *ii1_667 = buffer.data(ii1 + 667);
    const auto *ii1_668 = buffer.data(ii1 + 668);
    const auto *ii1_669 = buffer.data(ii1 + 669);
    const auto *ii1_670 = buffer.data(ii1 + 670);
    const auto *ii1_671 = buffer.data(ii1 + 671);
    const auto *ii1_672 = buffer.data(ii1 + 672);
    const auto *ii1_675 = buffer.data(ii1 + 675);
    const auto *ii1_677 = buffer.data(ii1 + 677);
    const auto *ii1_678 = buffer.data(ii1 + 678);
    const auto *ii1_681 = buffer.data(ii1 + 681);
    const auto *ii1_682 = buffer.data(ii1 + 682);
    const auto *ii1_684 = buffer.data(ii1 + 684);

    const auto *ik_770 = buffer.data(ik + 770);
    const auto *ik_771 = buffer.data(ik + 771);
    const auto *ik_773 = buffer.data(ik + 773);
    const auto *ik_774 = buffer.data(ik + 774);
    const auto *ik_776 = buffer.data(ik + 776);
    const auto *ik_777 = buffer.data(ik + 777);
    const auto *ik_779 = buffer.data(ik + 779);
    const auto *ik_780 = buffer.data(ik + 780);
    const auto *ik_781 = buffer.data(ik + 781);
    const auto *ik_783 = buffer.data(ik + 783);
    const auto *ik_784 = buffer.data(ik + 784);
    const auto *ik_785 = buffer.data(ik + 785);
    const auto *ik_786 = buffer.data(ik + 786);
    const auto *ik_787 = buffer.data(ik + 787);
    const auto *ik_788 = buffer.data(ik + 788);
    const auto *ik_789 = buffer.data(ik + 789);
    const auto *ik_790 = buffer.data(ik + 790);
    const auto *ik_791 = buffer.data(ik + 791);
    const auto *ik_792 = buffer.data(ik + 792);
    const auto *ik_794 = buffer.data(ik + 794);
    const auto *ik_795 = buffer.data(ik + 795);
    const auto *ik_797 = buffer.data(ik + 797);
    const auto *ik_798 = buffer.data(ik + 798);
    const auto *ik_801 = buffer.data(ik + 801);
    const auto *ik_802 = buffer.data(ik + 802);
    const auto *ik_806 = buffer.data(ik + 806);
    const auto *ik_807 = buffer.data(ik + 807);
    const auto *ik_812 = buffer.data(ik + 812);
    const auto *ik_820 = buffer.data(ik + 820);
    const auto *ik_821 = buffer.data(ik + 821);
    const auto *ik_822 = buffer.data(ik + 822);
    const auto *ik_823 = buffer.data(ik + 823);
    const auto *ik_824 = buffer.data(ik + 824);
    const auto *ik_825 = buffer.data(ik + 825);
    const auto *ik_826 = buffer.data(ik + 826);
    const auto *ik_827 = buffer.data(ik + 827);
    const auto *ik_828 = buffer.data(ik + 828);
    const auto *ik_830 = buffer.data(ik + 830);
    const auto *ik_831 = buffer.data(ik + 831);
    const auto *ik_833 = buffer.data(ik + 833);
    const auto *ik_834 = buffer.data(ik + 834);
    const auto *ik_837 = buffer.data(ik + 837);
    const auto *ik_838 = buffer.data(ik + 838);
    const auto *ik_840 = buffer.data(ik + 840);
    const auto *ik_842 = buffer.data(ik + 842);
    const auto *ik_843 = buffer.data(ik + 843);
    const auto *ik_845 = buffer.data(ik + 845);
    const auto *ik_846 = buffer.data(ik + 846);
    const auto *ik_848 = buffer.data(ik + 848);
    const auto *ik_849 = buffer.data(ik + 849);
    const auto *ik_851 = buffer.data(ik + 851);
    const auto *ik_852 = buffer.data(ik + 852);
    const auto *ik_853 = buffer.data(ik + 853);
    const auto *ik_855 = buffer.data(ik + 855);
    const auto *ik_856 = buffer.data(ik + 856);
    const auto *ik_857 = buffer.data(ik + 857);
    const auto *ik_858 = buffer.data(ik + 858);
    const auto *ik_859 = buffer.data(ik + 859);
    const auto *ik_860 = buffer.data(ik + 860);
    const auto *ik_861 = buffer.data(ik + 861);
    const auto *ik_862 = buffer.data(ik + 862);
    const auto *ik_863 = buffer.data(ik + 863);
    const auto *ik_864 = buffer.data(ik + 864);
    const auto *ik_866 = buffer.data(ik + 866);
    const auto *ik_867 = buffer.data(ik + 867);
    const auto *ik_869 = buffer.data(ik + 869);
    const auto *ik_870 = buffer.data(ik + 870);
    const auto *ik_873 = buffer.data(ik + 873);
    const auto *ik_874 = buffer.data(ik + 874);
    const auto *ik_876 = buffer.data(ik + 876);

#pragma omp simd aligned(t_962, t_963, t_964, pb_x, pb_y, hk_554, ii0_605, ii0_606, ii1_605, \
                         ii1_606, ik_770, ik_773, ik_774 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_5 * ii0_605[k]
                   - f_6 * ii1_605[k]
                   + pb_x[k] * ik_773[k];

        t_963[k] = f_5 * ii0_606[k]
                   - f_6 * ii1_606[k]
                   + pb_x[k] * ik_774[k];

        t_964[k] = f_0 * hk_554[k]
                   + pb_y[k] * ik_770[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pb_x, pb_z, ii0_608, ii0_609, ii0_611, \
                         ii1_608, ii1_609, ii1_611, ik_771, ik_776, ik_777, \
                         ik_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_5 * ii0_608[k]
                   - f_6 * ii1_608[k]
                   + pb_x[k] * ik_776[k];

        t_966[k] = f_3 * ii0_609[k]
                   - f_4 * ii1_609[k]
                   + pb_x[k] * ik_777[k];

        t_967[k] = pb_z[k] * ik_771[k];

        t_968[k] = f_3 * ii0_611[k]
                   - f_4 * ii1_611[k]
                   + pb_x[k] * ik_779[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pb_x, pb_y, hk_560, ii0_612, ii0_613, ii1_612, \
                         ii1_613, ik_776, ik_780, ik_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_3 * ii0_612[k]
                   - f_4 * ii1_612[k]
                   + pb_x[k] * ik_780[k];

        t_970[k] = f_3 * ii0_613[k]
                   - f_4 * ii1_613[k]
                   + pb_x[k] * ik_781[k];

        t_971[k] = f_0 * hk_560[k]
                   + pb_y[k] * ik_776[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, t_976, t_977, pb_x, ii0_615, ii1_615, \
                         ik_783, ik_784, ik_785, ik_786, ik_787, \
                         ik_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_3 * ii0_615[k]
                   - f_4 * ii1_615[k]
                   + pb_x[k] * ik_783[k];

        t_973[k] = pb_x[k] * ik_784[k];

        t_974[k] = pb_x[k] * ik_785[k];

        t_975[k] = pb_x[k] * ik_786[k];

        t_976[k] = pb_x[k] * ik_787[k];

        t_977[k] = pb_x[k] * ik_788[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, t_982, pb_x, pb_y, pb_z, hk_568, ii0_609, \
                         ii1_609, ik_784, ik_789, ik_790, ik_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = pb_x[k] * ik_789[k];

        t_979[k] = pb_x[k] * ik_790[k];

        t_980[k] = pb_x[k] * ik_791[k];

        t_981[k] = f_0 * hk_568[k]
                   + f_1 * ii0_609[k]
                   - f_2 * ii1_609[k]
                   + pb_y[k] * ik_784[k];

        t_982[k] = pb_z[k] * ik_784[k];
    }

#pragma omp simd aligned(t_983, t_984, t_985, pb_z, ii0_609, ii0_610, ii0_611, ii1_609, \
                         ii1_610, ii1_611, ik_785, ik_786, ik_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_983[k] = f_3 * ii0_609[k]
                   - f_4 * ii1_609[k]
                   + pb_z[k] * ik_785[k];

        t_984[k] = f_5 * ii0_610[k]
                   - f_6 * ii1_610[k]
                   + pb_z[k] * ik_786[k];

        t_985[k] = f_7 * ii0_611[k]
                   - f_8 * ii1_611[k]
                   + pb_z[k] * ik_787[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, pb_y, pb_z, hk_575, ii0_612, ii0_613, \
                         ii0_615, ii1_612, ii1_613, ii1_615, ik_788, ik_789, \
                         ik_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_9 * ii0_612[k]
                   - f_10 * ii1_612[k]
                   + pb_z[k] * ik_788[k];

        t_987[k] = f_11 * ii0_613[k]
                   - f_12 * ii1_613[k]
                   + pb_z[k] * ik_789[k];

        t_988[k] = f_0 * hk_575[k]
                   + pb_y[k] * ik_791[k];

        t_989[k] = f_1 * ii0_615[k]
                   - f_2 * ii1_615[k]
                   + pb_z[k] * ik_791[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, t_994, pa_z, pb_y, pb_z, hk_540, hk_578, \
                         hl_675, hl_676, hl_678, ik_792, ik_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pa_z[k] * hl_675[k];

        t_991[k] = pa_z[k] * hl_676[k];

        t_992[k] = f_13 * hk_540[k]
                   + pb_z[k] * ik_792[k];

        t_993[k] = pa_z[k] * hl_678[k];

        t_994[k] = f_17 * hk_578[k]
                   + pb_y[k] * ik_794[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pa_z, pb_y, pb_z, hk_542, hk_543, hk_581, \
                         hl_680, hl_681, ik_795, ik_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_14 * hk_542[k]
                   + pa_z[k] * hl_680[k];

        t_996[k] = pa_z[k] * hl_681[k];

        t_997[k] = f_13 * hk_543[k]
                   + pb_z[k] * ik_795[k];

        t_998[k] = f_17 * hk_581[k]
                   + pb_y[k] * ik_797[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, t_1002, pa_z, pb_z, hk_545, hk_546, hk_547, \
                         hl_684, hl_685, hl_687, ik_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_15 * hk_545[k]
                   + pa_z[k] * hl_684[k];

        t_1000[k] = pa_z[k] * hl_685[k];

        t_1001[k] = f_13 * hk_546[k]
                    + pb_z[k] * ik_798[k];

        t_1002[k] = f_14 * hk_547[k]
                    + pa_z[k] * hl_687[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, t_1006, pa_z, pb_y, pb_z, hk_549, hk_550, \
                         hk_585, hl_689, hl_690, ik_801, ik_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_17 * hk_585[k]
                    + pb_y[k] * ik_801[k];

        t_1004[k] = f_16 * hk_549[k]
                    + pa_z[k] * hl_689[k];

        t_1005[k] = pa_z[k] * hl_690[k];

        t_1006[k] = f_13 * hk_550[k]
                    + pb_z[k] * ik_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, t_1010, t_1011, pa_z, pb_y, hk_551, hk_552, \
                         hk_554, hk_590, hl_692, hl_693, hl_695, hl_696, \
                         ik_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = f_14 * hk_551[k]
                    + pa_z[k] * hl_692[k];

        t_1008[k] = f_15 * hk_552[k]
                    + pa_z[k] * hl_693[k];

        t_1009[k] = f_17 * hk_590[k]
                    + pb_y[k] * ik_806[k];

        t_1010[k] = f_17 * hk_554[k]
                    + pa_z[k] * hl_695[k];

        t_1011[k] = pa_z[k] * hl_696[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pa_z, pb_z, hk_555, hk_556, hk_557, \
                         hk_558, hl_698, hl_699, hl_700, ik_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_13 * hk_555[k]
                    + pb_z[k] * ik_807[k];

        t_1013[k] = f_14 * hk_556[k]
                    + pa_z[k] * hl_698[k];

        t_1014[k] = f_15 * hk_557[k]
                    + pa_z[k] * hl_699[k];

        t_1015[k] = f_16 * hk_558[k]
                    + pa_z[k] * hl_700[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, t_1020, pa_z, pb_x, pb_y, hk_560, \
                         hk_596, hl_702, ik_812, ik_820, ik_821, \
                         ik_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_17 * hk_596[k]
                    + pb_y[k] * ik_812[k];

        t_1017[k] = f_0 * hk_560[k]
                    + pa_z[k] * hl_702[k];

        t_1018[k] = pb_x[k] * ik_820[k];

        t_1019[k] = pb_x[k] * ik_821[k];

        t_1020[k] = pb_x[k] * ik_822[k];
    }

#pragma omp simd aligned(t_1021, t_1022, t_1023, t_1024, t_1025, t_1026, pa_z, pb_x, hl_711, \
                         ik_823, ik_824, ik_825, ik_826, ik_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1021[k] = pb_x[k] * ik_823[k];

        t_1022[k] = pb_x[k] * ik_824[k];

        t_1023[k] = pb_x[k] * ik_825[k];

        t_1024[k] = pb_x[k] * ik_826[k];

        t_1025[k] = pb_x[k] * ik_827[k];

        t_1026[k] = pa_z[k] * hl_711[k];
    }

#pragma omp simd aligned(t_1027, t_1028, t_1029, t_1030, pa_z, pb_z, hk_568, hk_569, hk_570, \
                         hk_571, hl_713, hl_714, hl_715, ik_820 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1027[k] = f_13 * hk_568[k]
                    + pb_z[k] * ik_820[k];

        t_1028[k] = f_14 * hk_569[k]
                    + pa_z[k] * hl_713[k];

        t_1029[k] = f_15 * hk_570[k]
                    + pa_z[k] * hl_714[k];

        t_1030[k] = f_16 * hk_571[k]
                    + pa_z[k] * hl_715[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, t_1034, pa_z, pb_y, hk_572, hk_573, hk_575, \
                         hk_611, hl_716, hl_717, hl_719, ik_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_17 * hk_572[k]
                    + pa_z[k] * hl_716[k];

        t_1032[k] = f_0 * hk_573[k]
                    + pa_z[k] * hl_717[k];

        t_1033[k] = f_17 * hk_611[k]
                    + pb_y[k] * ik_827[k];

        t_1034[k] = f_18 * hk_575[k]
                    + pa_z[k] * hl_719[k];
    }

#pragma omp simd aligned(t_1035, t_1036, t_1037, t_1038, pb_x, pb_y, pb_z, hk_576, hk_612, \
                         ii0_644, ii0_647, ii1_644, ii1_647, ik_828, \
                         ik_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1035[k] = f_1 * ii0_644[k]
                    - f_2 * ii1_644[k]
                    + pb_x[k] * ik_828[k];

        t_1036[k] = f_16 * hk_612[k]
                    + pb_y[k] * ik_828[k];

        t_1037[k] = f_14 * hk_576[k]
                    + pb_z[k] * ik_828[k];

        t_1038[k] = f_11 * ii0_647[k]
                    - f_12 * ii1_647[k]
                    + pb_x[k] * ik_831[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, pb_x, pb_y, hk_614, ii0_649, ii0_650, \
                         ii1_649, ii1_650, ik_830, ik_833, ik_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_16 * hk_614[k]
                    + pb_y[k] * ik_830[k];

        t_1040[k] = f_11 * ii0_649[k]
                    - f_12 * ii1_649[k]
                    + pb_x[k] * ik_833[k];

        t_1041[k] = f_9 * ii0_650[k]
                    - f_10 * ii1_650[k]
                    + pb_x[k] * ik_834[k];
    }

#pragma omp simd aligned(t_1042, t_1043, t_1044, pb_x, pb_y, pb_z, hk_579, hk_617, ii0_653, \
                         ii1_653, ik_831, ik_833, ik_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1042[k] = f_14 * hk_579[k]
                    + pb_z[k] * ik_831[k];

        t_1043[k] = f_16 * hk_617[k]
                    + pb_y[k] * ik_833[k];

        t_1044[k] = f_9 * ii0_653[k]
                    - f_10 * ii1_653[k]
                    + pb_x[k] * ik_837[k];
    }

#pragma omp simd aligned(t_1045, t_1046, t_1047, pb_x, pb_z, hk_582, ii0_654, ii0_656, \
                         ii1_654, ii1_656, ik_834, ik_838, ik_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1045[k] = f_7 * ii0_654[k]
                    - f_8 * ii1_654[k]
                    + pb_x[k] * ik_838[k];

        t_1046[k] = f_14 * hk_582[k]
                    + pb_z[k] * ik_834[k];

        t_1047[k] = f_7 * ii0_656[k]
                    - f_8 * ii1_656[k]
                    + pb_x[k] * ik_840[k];
    }

#pragma omp simd aligned(t_1048, t_1049, t_1050, pb_x, pb_y, hk_621, ii0_658, ii0_659, \
                         ii1_658, ii1_659, ik_837, ik_842, ik_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1048[k] = f_16 * hk_621[k]
                    + pb_y[k] * ik_837[k];

        t_1049[k] = f_7 * ii0_658[k]
                    - f_8 * ii1_658[k]
                    + pb_x[k] * ik_842[k];

        t_1050[k] = f_5 * ii0_659[k]
                    - f_6 * ii1_659[k]
                    + pb_x[k] * ik_843[k];
    }

#pragma omp simd aligned(t_1051, t_1052, t_1053, pb_x, pb_z, hk_586, ii0_661, ii0_662, \
                         ii1_661, ii1_662, ik_838, ik_845, ik_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1051[k] = f_14 * hk_586[k]
                    + pb_z[k] * ik_838[k];

        t_1052[k] = f_5 * ii0_661[k]
                    - f_6 * ii1_661[k]
                    + pb_x[k] * ik_845[k];

        t_1053[k] = f_5 * ii0_662[k]
                    - f_6 * ii1_662[k]
                    + pb_x[k] * ik_846[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, pb_x, pb_y, hk_626, ii0_664, ii0_665, \
                         ii1_664, ii1_665, ik_842, ik_848, ik_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_16 * hk_626[k]
                    + pb_y[k] * ik_842[k];

        t_1055[k] = f_5 * ii0_664[k]
                    - f_6 * ii1_664[k]
                    + pb_x[k] * ik_848[k];

        t_1056[k] = f_3 * ii0_665[k]
                    - f_4 * ii1_665[k]
                    + pb_x[k] * ik_849[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, pb_x, pb_z, hk_591, ii0_667, ii0_668, \
                         ii1_667, ii1_668, ik_843, ik_851, ik_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_14 * hk_591[k]
                    + pb_z[k] * ik_843[k];

        t_1058[k] = f_3 * ii0_667[k]
                    - f_4 * ii1_667[k]
                    + pb_x[k] * ik_851[k];

        t_1059[k] = f_3 * ii0_668[k]
                    - f_4 * ii1_668[k]
                    + pb_x[k] * ik_852[k];
    }

#pragma omp simd aligned(t_1060, t_1061, t_1062, t_1063, pb_x, pb_y, hk_632, ii0_669, ii0_671, \
                         ii1_669, ii1_671, ik_848, ik_853, ik_855, \
                         ik_856 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1060[k] = f_3 * ii0_669[k]
                    - f_4 * ii1_669[k]
                    + pb_x[k] * ik_853[k];

        t_1061[k] = f_16 * hk_632[k]
                    + pb_y[k] * ik_848[k];

        t_1062[k] = f_3 * ii0_671[k]
                    - f_4 * ii1_671[k]
                    + pb_x[k] * ik_855[k];

        t_1063[k] = pb_x[k] * ik_856[k];
    }

#pragma omp simd aligned(t_1064, t_1065, t_1066, t_1067, t_1068, t_1069, t_1070, pb_x, ik_857, \
                         ik_858, ik_859, ik_860, ik_861, ik_862, \
                         ik_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1064[k] = pb_x[k] * ik_857[k];

        t_1065[k] = pb_x[k] * ik_858[k];

        t_1066[k] = pb_x[k] * ik_859[k];

        t_1067[k] = pb_x[k] * ik_860[k];

        t_1068[k] = pb_x[k] * ik_861[k];

        t_1069[k] = pb_x[k] * ik_862[k];

        t_1070[k] = pb_x[k] * ik_863[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, pa_z, pb_y, pb_z, gl0_486, gl1_486, hk_604, \
                         hk_642, hl_756, ii0_667, ii1_667, ik_856, \
                         ik_858 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = f_19 * gl0_486[k]
                    - f_20 * gl1_486[k]
                    + pa_z[k] * hl_756[k];

        t_1072[k] = f_14 * hk_604[k]
                    + pb_z[k] * ik_856[k];

        t_1073[k] = f_16 * hk_642[k]
                    + f_11 * ii0_667[k]
                    - f_12 * ii1_667[k]
                    + pb_y[k] * ik_858[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, pb_y, hk_643, hk_644, hk_645, ii0_668, \
                         ii0_669, ii0_670, ii1_668, ii1_669, ii1_670, ik_859, ik_860, \
                         ik_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = f_16 * hk_643[k]
                    + f_9 * ii0_668[k]
                    - f_10 * ii1_668[k]
                    + pb_y[k] * ik_859[k];

        t_1075[k] = f_16 * hk_644[k]
                    + f_7 * ii0_669[k]
                    - f_8 * ii1_669[k]
                    + pb_y[k] * ik_860[k];

        t_1076[k] = f_16 * hk_645[k]
                    + f_5 * ii0_670[k]
                    - f_6 * ii1_670[k]
                    + pb_y[k] * ik_861[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, pa_y, pb_y, gl0_584, gl1_584, hk_646, hk_647, \
                         hl_809, ii0_671, ii1_671, ik_862, ik_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_16 * hk_646[k]
                    + f_3 * ii0_671[k]
                    - f_4 * ii1_671[k]
                    + pb_y[k] * ik_862[k];

        t_1078[k] = f_16 * hk_647[k]
                    + pb_y[k] * ik_863[k];

        t_1079[k] = f_21 * gl0_584[k]
                    - f_22 * gl1_584[k]
                    + pa_y[k] * hl_809[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, t_1083, pb_x, pb_y, pb_z, hk_612, hk_648, \
                         ii0_672, ii0_675, ii1_672, ii1_675, ik_864, \
                         ik_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_1 * ii0_672[k]
                    - f_2 * ii1_672[k]
                    + pb_x[k] * ik_864[k];

        t_1081[k] = f_15 * hk_648[k]
                    + pb_y[k] * ik_864[k];

        t_1082[k] = f_15 * hk_612[k]
                    + pb_z[k] * ik_864[k];

        t_1083[k] = f_11 * ii0_675[k]
                    - f_12 * ii1_675[k]
                    + pb_x[k] * ik_867[k];
    }

#pragma omp simd aligned(t_1084, t_1085, t_1086, pb_x, pb_y, hk_650, ii0_677, ii0_678, \
                         ii1_677, ii1_678, ik_866, ik_869, ik_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1084[k] = f_15 * hk_650[k]
                    + pb_y[k] * ik_866[k];

        t_1085[k] = f_11 * ii0_677[k]
                    - f_12 * ii1_677[k]
                    + pb_x[k] * ik_869[k];

        t_1086[k] = f_9 * ii0_678[k]
                    - f_10 * ii1_678[k]
                    + pb_x[k] * ik_870[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pb_x, pb_y, pb_z, hk_615, hk_653, ii0_681, \
                         ii1_681, ik_867, ik_869, ik_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_15 * hk_615[k]
                    + pb_z[k] * ik_867[k];

        t_1088[k] = f_15 * hk_653[k]
                    + pb_y[k] * ik_869[k];

        t_1089[k] = f_9 * ii0_681[k]
                    - f_10 * ii1_681[k]
                    + pb_x[k] * ik_873[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pb_x, pb_z, hk_618, ii0_682, ii0_684, \
                         ii1_682, ii1_684, ik_870, ik_874, ik_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_7 * ii0_682[k]
                    - f_8 * ii1_682[k]
                    + pb_x[k] * ik_874[k];

        t_1091[k] = f_15 * hk_618[k]
                    + pb_z[k] * ik_870[k];

        t_1092[k] = f_7 * ii0_684[k]
                    - f_8 * ii1_684[k]
                    + pb_x[k] * ik_876[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t gl0,
                                            const size_t gl1, const size_t hk, const size_t hl,
                                            const size_t ii0, const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 4.0 / p;
    const auto f_19 = 0.5 / alpha;
    const auto f_20 = 0.5 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);
    const auto f_23 = 1.0 / alpha;
    const auto f_24 = beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gl0_531 = buffer.data(gl0 + 531);
    const auto *gl0_576 = buffer.data(gl0 + 576);
    const auto *gl0_629 = buffer.data(gl0 + 629);
    const auto *gl0_674 = buffer.data(gl0 + 674);

    const auto *gl1_531 = buffer.data(gl1 + 531);
    const auto *gl1_576 = buffer.data(gl1 + 576);
    const auto *gl1_629 = buffer.data(gl1 + 629);
    const auto *gl1_674 = buffer.data(gl1 + 674);

    const auto *hk_622 = buffer.data(hk + 622);
    const auto *hk_627 = buffer.data(hk + 627);
    const auto *hk_640 = buffer.data(hk + 640);
    const auto *hk_648 = buffer.data(hk + 648);
    const auto *hk_651 = buffer.data(hk + 651);
    const auto *hk_654 = buffer.data(hk + 654);
    const auto *hk_657 = buffer.data(hk + 657);
    const auto *hk_658 = buffer.data(hk + 658);
    const auto *hk_662 = buffer.data(hk + 662);
    const auto *hk_663 = buffer.data(hk + 663);
    const auto *hk_668 = buffer.data(hk + 668);
    const auto *hk_676 = buffer.data(hk + 676);
    const auto *hk_678 = buffer.data(hk + 678);
    const auto *hk_679 = buffer.data(hk + 679);
    const auto *hk_680 = buffer.data(hk + 680);
    const auto *hk_681 = buffer.data(hk + 681);
    const auto *hk_682 = buffer.data(hk + 682);
    const auto *hk_683 = buffer.data(hk + 683);
    const auto *hk_684 = buffer.data(hk + 684);
    const auto *hk_686 = buffer.data(hk + 686);
    const auto *hk_687 = buffer.data(hk + 687);
    const auto *hk_689 = buffer.data(hk + 689);
    const auto *hk_690 = buffer.data(hk + 690);
    const auto *hk_693 = buffer.data(hk + 693);
    const auto *hk_694 = buffer.data(hk + 694);
    const auto *hk_698 = buffer.data(hk + 698);
    const auto *hk_699 = buffer.data(hk + 699);
    const auto *hk_704 = buffer.data(hk + 704);
    const auto *hk_712 = buffer.data(hk + 712);
    const auto *hk_714 = buffer.data(hk + 714);
    const auto *hk_715 = buffer.data(hk + 715);
    const auto *hk_716 = buffer.data(hk + 716);
    const auto *hk_717 = buffer.data(hk + 717);
    const auto *hk_718 = buffer.data(hk + 718);
    const auto *hk_719 = buffer.data(hk + 719);
    const auto *hk_720 = buffer.data(hk + 720);
    const auto *hk_721 = buffer.data(hk + 721);
    const auto *hk_722 = buffer.data(hk + 722);
    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_725 = buffer.data(hk + 725);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_728 = buffer.data(hk + 728);
    const auto *hk_729 = buffer.data(hk + 729);
    const auto *hk_730 = buffer.data(hk + 730);
    const auto *hk_732 = buffer.data(hk + 732);
    const auto *hk_733 = buffer.data(hk + 733);
    const auto *hk_734 = buffer.data(hk + 734);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_737 = buffer.data(hk + 737);
    const auto *hk_738 = buffer.data(hk + 738);
    const auto *hk_739 = buffer.data(hk + 739);
    const auto *hk_740 = buffer.data(hk + 740);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_750 = buffer.data(hk + 750);
    const auto *hk_751 = buffer.data(hk + 751);
    const auto *hk_752 = buffer.data(hk + 752);
    const auto *hk_753 = buffer.data(hk + 753);
    const auto *hk_754 = buffer.data(hk + 754);
    const auto *hk_755 = buffer.data(hk + 755);

    const auto *hl_801 = buffer.data(hl + 801);
    const auto *hl_846 = buffer.data(hl + 846);
    const auto *hl_854 = buffer.data(hl + 854);
    const auto *hl_899 = buffer.data(hl + 899);
    const auto *hl_900 = buffer.data(hl + 900);
    const auto *hl_902 = buffer.data(hl + 902);
    const auto *hl_903 = buffer.data(hl + 903);
    const auto *hl_905 = buffer.data(hl + 905);
    const auto *hl_906 = buffer.data(hl + 906);
    const auto *hl_909 = buffer.data(hl + 909);
    const auto *hl_910 = buffer.data(hl + 910);
    const auto *hl_912 = buffer.data(hl + 912);
    const auto *hl_914 = buffer.data(hl + 914);
    const auto *hl_915 = buffer.data(hl + 915);
    const auto *hl_917 = buffer.data(hl + 917);
    const auto *hl_918 = buffer.data(hl + 918);
    const auto *hl_920 = buffer.data(hl + 920);
    const auto *hl_921 = buffer.data(hl + 921);
    const auto *hl_923 = buffer.data(hl + 923);
    const auto *hl_924 = buffer.data(hl + 924);
    const auto *hl_925 = buffer.data(hl + 925);
    const auto *hl_927 = buffer.data(hl + 927);
    const auto *hl_936 = buffer.data(hl + 936);
    const auto *hl_938 = buffer.data(hl + 938);
    const auto *hl_939 = buffer.data(hl + 939);
    const auto *hl_940 = buffer.data(hl + 940);
    const auto *hl_941 = buffer.data(hl + 941);
    const auto *hl_942 = buffer.data(hl + 942);
    const auto *hl_944 = buffer.data(hl + 944);

    const auto *ii0_686 = buffer.data(ii0 + 686);
    const auto *ii0_687 = buffer.data(ii0 + 687);
    const auto *ii0_689 = buffer.data(ii0 + 689);
    const auto *ii0_690 = buffer.data(ii0 + 690);
    const auto *ii0_692 = buffer.data(ii0 + 692);
    const auto *ii0_693 = buffer.data(ii0 + 693);
    const auto *ii0_695 = buffer.data(ii0 + 695);
    const auto *ii0_696 = buffer.data(ii0 + 696);
    const auto *ii0_697 = buffer.data(ii0 + 697);
    const auto *ii0_698 = buffer.data(ii0 + 698);
    const auto *ii0_699 = buffer.data(ii0 + 699);
    const auto *ii0_700 = buffer.data(ii0 + 700);
    const auto *ii0_703 = buffer.data(ii0 + 703);
    const auto *ii0_705 = buffer.data(ii0 + 705);
    const auto *ii0_706 = buffer.data(ii0 + 706);
    const auto *ii0_709 = buffer.data(ii0 + 709);
    const auto *ii0_710 = buffer.data(ii0 + 710);
    const auto *ii0_712 = buffer.data(ii0 + 712);
    const auto *ii0_714 = buffer.data(ii0 + 714);
    const auto *ii0_715 = buffer.data(ii0 + 715);
    const auto *ii0_717 = buffer.data(ii0 + 717);
    const auto *ii0_718 = buffer.data(ii0 + 718);
    const auto *ii0_720 = buffer.data(ii0 + 720);
    const auto *ii0_721 = buffer.data(ii0 + 721);
    const auto *ii0_723 = buffer.data(ii0 + 723);
    const auto *ii0_724 = buffer.data(ii0 + 724);
    const auto *ii0_725 = buffer.data(ii0 + 725);
    const auto *ii0_726 = buffer.data(ii0 + 726);
    const auto *ii0_727 = buffer.data(ii0 + 727);
    const auto *ii0_756 = buffer.data(ii0 + 756);
    const auto *ii0_759 = buffer.data(ii0 + 759);
    const auto *ii0_761 = buffer.data(ii0 + 761);
    const auto *ii0_762 = buffer.data(ii0 + 762);

    const auto *ii1_686 = buffer.data(ii1 + 686);
    const auto *ii1_687 = buffer.data(ii1 + 687);
    const auto *ii1_689 = buffer.data(ii1 + 689);
    const auto *ii1_690 = buffer.data(ii1 + 690);
    const auto *ii1_692 = buffer.data(ii1 + 692);
    const auto *ii1_693 = buffer.data(ii1 + 693);
    const auto *ii1_695 = buffer.data(ii1 + 695);
    const auto *ii1_696 = buffer.data(ii1 + 696);
    const auto *ii1_697 = buffer.data(ii1 + 697);
    const auto *ii1_698 = buffer.data(ii1 + 698);
    const auto *ii1_699 = buffer.data(ii1 + 699);
    const auto *ii1_700 = buffer.data(ii1 + 700);
    const auto *ii1_703 = buffer.data(ii1 + 703);
    const auto *ii1_705 = buffer.data(ii1 + 705);
    const auto *ii1_706 = buffer.data(ii1 + 706);
    const auto *ii1_709 = buffer.data(ii1 + 709);
    const auto *ii1_710 = buffer.data(ii1 + 710);
    const auto *ii1_712 = buffer.data(ii1 + 712);
    const auto *ii1_714 = buffer.data(ii1 + 714);
    const auto *ii1_715 = buffer.data(ii1 + 715);
    const auto *ii1_717 = buffer.data(ii1 + 717);
    const auto *ii1_718 = buffer.data(ii1 + 718);
    const auto *ii1_720 = buffer.data(ii1 + 720);
    const auto *ii1_721 = buffer.data(ii1 + 721);
    const auto *ii1_723 = buffer.data(ii1 + 723);
    const auto *ii1_724 = buffer.data(ii1 + 724);
    const auto *ii1_725 = buffer.data(ii1 + 725);
    const auto *ii1_726 = buffer.data(ii1 + 726);
    const auto *ii1_727 = buffer.data(ii1 + 727);
    const auto *ii1_756 = buffer.data(ii1 + 756);
    const auto *ii1_759 = buffer.data(ii1 + 759);
    const auto *ii1_761 = buffer.data(ii1 + 761);
    const auto *ii1_762 = buffer.data(ii1 + 762);

    const auto *ik_873 = buffer.data(ik + 873);
    const auto *ik_874 = buffer.data(ik + 874);
    const auto *ik_878 = buffer.data(ik + 878);
    const auto *ik_879 = buffer.data(ik + 879);
    const auto *ik_881 = buffer.data(ik + 881);
    const auto *ik_882 = buffer.data(ik + 882);
    const auto *ik_884 = buffer.data(ik + 884);
    const auto *ik_885 = buffer.data(ik + 885);
    const auto *ik_887 = buffer.data(ik + 887);
    const auto *ik_888 = buffer.data(ik + 888);
    const auto *ik_889 = buffer.data(ik + 889);
    const auto *ik_891 = buffer.data(ik + 891);
    const auto *ik_892 = buffer.data(ik + 892);
    const auto *ik_893 = buffer.data(ik + 893);
    const auto *ik_894 = buffer.data(ik + 894);
    const auto *ik_895 = buffer.data(ik + 895);
    const auto *ik_896 = buffer.data(ik + 896);
    const auto *ik_897 = buffer.data(ik + 897);
    const auto *ik_898 = buffer.data(ik + 898);
    const auto *ik_899 = buffer.data(ik + 899);
    const auto *ik_900 = buffer.data(ik + 900);
    const auto *ik_902 = buffer.data(ik + 902);
    const auto *ik_903 = buffer.data(ik + 903);
    const auto *ik_905 = buffer.data(ik + 905);
    const auto *ik_906 = buffer.data(ik + 906);
    const auto *ik_909 = buffer.data(ik + 909);
    const auto *ik_910 = buffer.data(ik + 910);
    const auto *ik_912 = buffer.data(ik + 912);
    const auto *ik_914 = buffer.data(ik + 914);
    const auto *ik_915 = buffer.data(ik + 915);
    const auto *ik_917 = buffer.data(ik + 917);
    const auto *ik_918 = buffer.data(ik + 918);
    const auto *ik_920 = buffer.data(ik + 920);
    const auto *ik_921 = buffer.data(ik + 921);
    const auto *ik_923 = buffer.data(ik + 923);
    const auto *ik_924 = buffer.data(ik + 924);
    const auto *ik_925 = buffer.data(ik + 925);
    const auto *ik_927 = buffer.data(ik + 927);
    const auto *ik_928 = buffer.data(ik + 928);
    const auto *ik_929 = buffer.data(ik + 929);
    const auto *ik_930 = buffer.data(ik + 930);
    const auto *ik_931 = buffer.data(ik + 931);
    const auto *ik_932 = buffer.data(ik + 932);
    const auto *ik_933 = buffer.data(ik + 933);
    const auto *ik_934 = buffer.data(ik + 934);
    const auto *ik_935 = buffer.data(ik + 935);
    const auto *ik_936 = buffer.data(ik + 936);
    const auto *ik_938 = buffer.data(ik + 938);
    const auto *ik_939 = buffer.data(ik + 939);
    const auto *ik_941 = buffer.data(ik + 941);
    const auto *ik_942 = buffer.data(ik + 942);
    const auto *ik_945 = buffer.data(ik + 945);
    const auto *ik_946 = buffer.data(ik + 946);
    const auto *ik_950 = buffer.data(ik + 950);
    const auto *ik_951 = buffer.data(ik + 951);
    const auto *ik_956 = buffer.data(ik + 956);
    const auto *ik_964 = buffer.data(ik + 964);
    const auto *ik_965 = buffer.data(ik + 965);
    const auto *ik_966 = buffer.data(ik + 966);
    const auto *ik_967 = buffer.data(ik + 967);
    const auto *ik_968 = buffer.data(ik + 968);
    const auto *ik_969 = buffer.data(ik + 969);
    const auto *ik_970 = buffer.data(ik + 970);
    const auto *ik_971 = buffer.data(ik + 971);
    const auto *ik_972 = buffer.data(ik + 972);
    const auto *ik_974 = buffer.data(ik + 974);
    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);

#pragma omp simd aligned(t_1093, t_1094, t_1095, pb_x, pb_y, hk_657, ii0_686, ii0_687, \
                         ii1_686, ii1_687, ik_873, ik_878, ik_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_15 * hk_657[k]
                    + pb_y[k] * ik_873[k];

        t_1094[k] = f_7 * ii0_686[k]
                    - f_8 * ii1_686[k]
                    + pb_x[k] * ik_878[k];

        t_1095[k] = f_5 * ii0_687[k]
                    - f_6 * ii1_687[k]
                    + pb_x[k] * ik_879[k];
    }

#pragma omp simd aligned(t_1096, t_1097, t_1098, pb_x, pb_z, hk_622, ii0_689, ii0_690, \
                         ii1_689, ii1_690, ik_874, ik_881, ik_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1096[k] = f_15 * hk_622[k]
                    + pb_z[k] * ik_874[k];

        t_1097[k] = f_5 * ii0_689[k]
                    - f_6 * ii1_689[k]
                    + pb_x[k] * ik_881[k];

        t_1098[k] = f_5 * ii0_690[k]
                    - f_6 * ii1_690[k]
                    + pb_x[k] * ik_882[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, pb_x, pb_y, hk_662, ii0_692, ii0_693, \
                         ii1_692, ii1_693, ik_878, ik_884, ik_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_15 * hk_662[k]
                    + pb_y[k] * ik_878[k];

        t_1100[k] = f_5 * ii0_692[k]
                    - f_6 * ii1_692[k]
                    + pb_x[k] * ik_884[k];

        t_1101[k] = f_3 * ii0_693[k]
                    - f_4 * ii1_693[k]
                    + pb_x[k] * ik_885[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, pb_x, pb_z, hk_627, ii0_695, ii0_696, \
                         ii1_695, ii1_696, ik_879, ik_887, ik_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = f_15 * hk_627[k]
                    + pb_z[k] * ik_879[k];

        t_1103[k] = f_3 * ii0_695[k]
                    - f_4 * ii1_695[k]
                    + pb_x[k] * ik_887[k];

        t_1104[k] = f_3 * ii0_696[k]
                    - f_4 * ii1_696[k]
                    + pb_x[k] * ik_888[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, pb_x, pb_y, hk_668, ii0_697, ii0_699, \
                         ii1_697, ii1_699, ik_884, ik_889, ik_891, \
                         ik_892 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_3 * ii0_697[k]
                    - f_4 * ii1_697[k]
                    + pb_x[k] * ik_889[k];

        t_1106[k] = f_15 * hk_668[k]
                    + pb_y[k] * ik_884[k];

        t_1107[k] = f_3 * ii0_699[k]
                    - f_4 * ii1_699[k]
                    + pb_x[k] * ik_891[k];

        t_1108[k] = pb_x[k] * ik_892[k];
    }

#pragma omp simd aligned(t_1109, t_1110, t_1111, t_1112, t_1113, t_1114, t_1115, pb_x, ik_893, \
                         ik_894, ik_895, ik_896, ik_897, ik_898, \
                         ik_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1109[k] = pb_x[k] * ik_893[k];

        t_1110[k] = pb_x[k] * ik_894[k];

        t_1111[k] = pb_x[k] * ik_895[k];

        t_1112[k] = pb_x[k] * ik_896[k];

        t_1113[k] = pb_x[k] * ik_897[k];

        t_1114[k] = pb_x[k] * ik_898[k];

        t_1115[k] = pb_x[k] * ik_899[k];
    }

#pragma omp simd aligned(t_1116, t_1117, t_1118, pa_z, pb_y, pb_z, gl0_531, gl1_531, hk_640, \
                         hk_678, hl_801, ii0_695, ii1_695, ik_892, \
                         ik_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1116[k] = f_23 * gl0_531[k]
                    - f_24 * gl1_531[k]
                    + pa_z[k] * hl_801[k];

        t_1117[k] = f_15 * hk_640[k]
                    + pb_z[k] * ik_892[k];

        t_1118[k] = f_15 * hk_678[k]
                    + f_11 * ii0_695[k]
                    - f_12 * ii1_695[k]
                    + pb_y[k] * ik_894[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, pb_y, hk_679, hk_680, hk_681, ii0_696, \
                         ii0_697, ii0_698, ii1_696, ii1_697, ii1_698, ik_895, ik_896, \
                         ik_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = f_15 * hk_679[k]
                    + f_9 * ii0_696[k]
                    - f_10 * ii1_696[k]
                    + pb_y[k] * ik_895[k];

        t_1120[k] = f_15 * hk_680[k]
                    + f_7 * ii0_697[k]
                    - f_8 * ii1_697[k]
                    + pb_y[k] * ik_896[k];

        t_1121[k] = f_15 * hk_681[k]
                    + f_5 * ii0_698[k]
                    - f_6 * ii1_698[k]
                    + pb_y[k] * ik_897[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, pa_y, pb_y, gl0_629, gl1_629, hk_682, hk_683, \
                         hl_854, ii0_699, ii1_699, ik_898, ik_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = f_15 * hk_682[k]
                    + f_3 * ii0_699[k]
                    - f_4 * ii1_699[k]
                    + pb_y[k] * ik_898[k];

        t_1123[k] = f_15 * hk_683[k]
                    + pb_y[k] * ik_899[k];

        t_1124[k] = f_23 * gl0_629[k]
                    - f_24 * gl1_629[k]
                    + pa_y[k] * hl_854[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, t_1128, pb_x, pb_y, pb_z, hk_648, hk_684, \
                         ii0_700, ii0_703, ii1_700, ii1_703, ik_900, \
                         ik_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_1 * ii0_700[k]
                    - f_2 * ii1_700[k]
                    + pb_x[k] * ik_900[k];

        t_1126[k] = f_14 * hk_684[k]
                    + pb_y[k] * ik_900[k];

        t_1127[k] = f_16 * hk_648[k]
                    + pb_z[k] * ik_900[k];

        t_1128[k] = f_11 * ii0_703[k]
                    - f_12 * ii1_703[k]
                    + pb_x[k] * ik_903[k];
    }

#pragma omp simd aligned(t_1129, t_1130, t_1131, pb_x, pb_y, hk_686, ii0_705, ii0_706, \
                         ii1_705, ii1_706, ik_902, ik_905, ik_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1129[k] = f_14 * hk_686[k]
                    + pb_y[k] * ik_902[k];

        t_1130[k] = f_11 * ii0_705[k]
                    - f_12 * ii1_705[k]
                    + pb_x[k] * ik_905[k];

        t_1131[k] = f_9 * ii0_706[k]
                    - f_10 * ii1_706[k]
                    + pb_x[k] * ik_906[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, pb_x, pb_y, pb_z, hk_651, hk_689, ii0_709, \
                         ii1_709, ik_903, ik_905, ik_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_16 * hk_651[k]
                    + pb_z[k] * ik_903[k];

        t_1133[k] = f_14 * hk_689[k]
                    + pb_y[k] * ik_905[k];

        t_1134[k] = f_9 * ii0_709[k]
                    - f_10 * ii1_709[k]
                    + pb_x[k] * ik_909[k];
    }

#pragma omp simd aligned(t_1135, t_1136, t_1137, pb_x, pb_z, hk_654, ii0_710, ii0_712, \
                         ii1_710, ii1_712, ik_906, ik_910, ik_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1135[k] = f_7 * ii0_710[k]
                    - f_8 * ii1_710[k]
                    + pb_x[k] * ik_910[k];

        t_1136[k] = f_16 * hk_654[k]
                    + pb_z[k] * ik_906[k];

        t_1137[k] = f_7 * ii0_712[k]
                    - f_8 * ii1_712[k]
                    + pb_x[k] * ik_912[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, pb_x, pb_y, hk_693, ii0_714, ii0_715, \
                         ii1_714, ii1_715, ik_909, ik_914, ik_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_14 * hk_693[k]
                    + pb_y[k] * ik_909[k];

        t_1139[k] = f_7 * ii0_714[k]
                    - f_8 * ii1_714[k]
                    + pb_x[k] * ik_914[k];

        t_1140[k] = f_5 * ii0_715[k]
                    - f_6 * ii1_715[k]
                    + pb_x[k] * ik_915[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, pb_x, pb_z, hk_658, ii0_717, ii0_718, \
                         ii1_717, ii1_718, ik_910, ik_917, ik_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_16 * hk_658[k]
                    + pb_z[k] * ik_910[k];

        t_1142[k] = f_5 * ii0_717[k]
                    - f_6 * ii1_717[k]
                    + pb_x[k] * ik_917[k];

        t_1143[k] = f_5 * ii0_718[k]
                    - f_6 * ii1_718[k]
                    + pb_x[k] * ik_918[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, pb_x, pb_y, hk_698, ii0_720, ii0_721, \
                         ii1_720, ii1_721, ik_914, ik_920, ik_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_14 * hk_698[k]
                    + pb_y[k] * ik_914[k];

        t_1145[k] = f_5 * ii0_720[k]
                    - f_6 * ii1_720[k]
                    + pb_x[k] * ik_920[k];

        t_1146[k] = f_3 * ii0_721[k]
                    - f_4 * ii1_721[k]
                    + pb_x[k] * ik_921[k];
    }

#pragma omp simd aligned(t_1147, t_1148, t_1149, pb_x, pb_z, hk_663, ii0_723, ii0_724, \
                         ii1_723, ii1_724, ik_915, ik_923, ik_924 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1147[k] = f_16 * hk_663[k]
                    + pb_z[k] * ik_915[k];

        t_1148[k] = f_3 * ii0_723[k]
                    - f_4 * ii1_723[k]
                    + pb_x[k] * ik_923[k];

        t_1149[k] = f_3 * ii0_724[k]
                    - f_4 * ii1_724[k]
                    + pb_x[k] * ik_924[k];
    }

#pragma omp simd aligned(t_1150, t_1151, t_1152, t_1153, pb_x, pb_y, hk_704, ii0_725, ii0_727, \
                         ii1_725, ii1_727, ik_920, ik_925, ik_927, \
                         ik_928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1150[k] = f_3 * ii0_725[k]
                    - f_4 * ii1_725[k]
                    + pb_x[k] * ik_925[k];

        t_1151[k] = f_14 * hk_704[k]
                    + pb_y[k] * ik_920[k];

        t_1152[k] = f_3 * ii0_727[k]
                    - f_4 * ii1_727[k]
                    + pb_x[k] * ik_927[k];

        t_1153[k] = pb_x[k] * ik_928[k];
    }

#pragma omp simd aligned(t_1154, t_1155, t_1156, t_1157, t_1158, t_1159, t_1160, pb_x, ik_929, \
                         ik_930, ik_931, ik_932, ik_933, ik_934, \
                         ik_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1154[k] = pb_x[k] * ik_929[k];

        t_1155[k] = pb_x[k] * ik_930[k];

        t_1156[k] = pb_x[k] * ik_931[k];

        t_1157[k] = pb_x[k] * ik_932[k];

        t_1158[k] = pb_x[k] * ik_933[k];

        t_1159[k] = pb_x[k] * ik_934[k];

        t_1160[k] = pb_x[k] * ik_935[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, pa_z, pb_y, pb_z, gl0_576, gl1_576, hk_676, \
                         hk_714, hl_846, ii0_723, ii1_723, ik_928, \
                         ik_930 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_21 * gl0_576[k]
                    - f_22 * gl1_576[k]
                    + pa_z[k] * hl_846[k];

        t_1162[k] = f_16 * hk_676[k]
                    + pb_z[k] * ik_928[k];

        t_1163[k] = f_14 * hk_714[k]
                    + f_11 * ii0_723[k]
                    - f_12 * ii1_723[k]
                    + pb_y[k] * ik_930[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, pb_y, hk_715, hk_716, hk_717, ii0_724, \
                         ii0_725, ii0_726, ii1_724, ii1_725, ii1_726, ik_931, ik_932, \
                         ik_933 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = f_14 * hk_715[k]
                    + f_9 * ii0_724[k]
                    - f_10 * ii1_724[k]
                    + pb_y[k] * ik_931[k];

        t_1165[k] = f_14 * hk_716[k]
                    + f_7 * ii0_725[k]
                    - f_8 * ii1_725[k]
                    + pb_y[k] * ik_932[k];

        t_1166[k] = f_14 * hk_717[k]
                    + f_5 * ii0_726[k]
                    - f_6 * ii1_726[k]
                    + pb_y[k] * ik_933[k];
    }

#pragma omp simd aligned(t_1167, t_1168, t_1169, t_1170, pa_y, pb_y, gl0_674, gl1_674, hk_718, \
                         hk_719, hl_899, hl_900, ii0_727, ii1_727, ik_934, \
                         ik_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1167[k] = f_14 * hk_718[k]
                    + f_3 * ii0_727[k]
                    - f_4 * ii1_727[k]
                    + pb_y[k] * ik_934[k];

        t_1168[k] = f_14 * hk_719[k]
                    + pb_y[k] * ik_935[k];

        t_1169[k] = f_19 * gl0_674[k]
                    - f_20 * gl1_674[k]
                    + pa_y[k] * hl_899[k];

        t_1170[k] = pa_y[k] * hl_900[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, t_1174, t_1175, pa_y, pb_y, hk_720, hk_721, \
                         hk_722, hl_902, hl_903, hl_905, ik_936, \
                         ik_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_13 * hk_720[k]
                    + pb_y[k] * ik_936[k];

        t_1172[k] = pa_y[k] * hl_902[k];

        t_1173[k] = f_14 * hk_721[k]
                    + pa_y[k] * hl_903[k];

        t_1174[k] = f_13 * hk_722[k]
                    + pb_y[k] * ik_938[k];

        t_1175[k] = pa_y[k] * hl_905[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, t_1179, pa_y, pb_y, pb_z, hk_687, hk_723, \
                         hk_725, hl_906, hl_909, ik_939, ik_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_15 * hk_723[k]
                    + pa_y[k] * hl_906[k];

        t_1177[k] = f_17 * hk_687[k]
                    + pb_z[k] * ik_939[k];

        t_1178[k] = f_13 * hk_725[k]
                    + pb_y[k] * ik_941[k];

        t_1179[k] = pa_y[k] * hl_909[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, t_1183, pa_y, pb_y, pb_z, hk_690, hk_726, \
                         hk_728, hk_729, hl_910, hl_912, ik_942, \
                         ik_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_16 * hk_726[k]
                    + pa_y[k] * hl_910[k];

        t_1181[k] = f_17 * hk_690[k]
                    + pb_z[k] * ik_942[k];

        t_1182[k] = f_14 * hk_728[k]
                    + pa_y[k] * hl_912[k];

        t_1183[k] = f_13 * hk_729[k]
                    + pb_y[k] * ik_945[k];
    }

#pragma omp simd aligned(t_1184, t_1185, t_1186, t_1187, t_1188, pa_y, pb_z, hk_694, hk_730, \
                         hk_732, hk_733, hl_914, hl_915, hl_917, hl_918, \
                         ik_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1184[k] = pa_y[k] * hl_914[k];

        t_1185[k] = f_17 * hk_730[k]
                    + pa_y[k] * hl_915[k];

        t_1186[k] = f_17 * hk_694[k]
                    + pb_z[k] * ik_946[k];

        t_1187[k] = f_15 * hk_732[k]
                    + pa_y[k] * hl_917[k];

        t_1188[k] = f_14 * hk_733[k]
                    + pa_y[k] * hl_918[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_y, pb_y, pb_z, hk_699, hk_734, \
                         hk_735, hl_920, hl_921, ik_950, ik_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_13 * hk_734[k]
                    + pb_y[k] * ik_950[k];

        t_1190[k] = pa_y[k] * hl_920[k];

        t_1191[k] = f_0 * hk_735[k]
                    + pa_y[k] * hl_921[k];

        t_1192[k] = f_17 * hk_699[k]
                    + pb_z[k] * ik_951[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, t_1197, pa_y, pb_y, hk_737, hk_738, \
                         hk_739, hk_740, hl_923, hl_924, hl_925, hl_927, \
                         ik_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_16 * hk_737[k]
                    + pa_y[k] * hl_923[k];

        t_1194[k] = f_15 * hk_738[k]
                    + pa_y[k] * hl_924[k];

        t_1195[k] = f_14 * hk_739[k]
                    + pa_y[k] * hl_925[k];

        t_1196[k] = f_13 * hk_740[k]
                    + pb_y[k] * ik_956[k];

        t_1197[k] = pa_y[k] * hl_927[k];
    }

#pragma omp simd aligned(t_1198, t_1199, t_1200, t_1201, t_1202, t_1203, t_1204, pb_x, ik_964, \
                         ik_965, ik_966, ik_967, ik_968, ik_969, \
                         ik_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1198[k] = pb_x[k] * ik_964[k];

        t_1199[k] = pb_x[k] * ik_965[k];

        t_1200[k] = pb_x[k] * ik_966[k];

        t_1201[k] = pb_x[k] * ik_967[k];

        t_1202[k] = pb_x[k] * ik_968[k];

        t_1203[k] = pb_x[k] * ik_969[k];

        t_1204[k] = pb_x[k] * ik_970[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, pa_y, pb_x, pb_z, hk_712, hk_748, \
                         hk_750, hl_936, hl_938, ik_964, ik_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = pb_x[k] * ik_971[k];

        t_1206[k] = f_18 * hk_748[k]
                    + pa_y[k] * hl_936[k];

        t_1207[k] = f_17 * hk_712[k]
                    + pb_z[k] * ik_964[k];

        t_1208[k] = f_0 * hk_750[k]
                    + pa_y[k] * hl_938[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, t_1212, pa_y, hk_751, hk_752, hk_753, hk_754, \
                         hl_939, hl_940, hl_941, hl_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = f_17 * hk_751[k]
                    + pa_y[k] * hl_939[k];

        t_1210[k] = f_16 * hk_752[k]
                    + pa_y[k] * hl_940[k];

        t_1211[k] = f_15 * hk_753[k]
                    + pa_y[k] * hl_941[k];

        t_1212[k] = f_14 * hk_754[k]
                    + pa_y[k] * hl_942[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, t_1216, t_1217, pa_y, pb_x, pb_y, pb_z, \
                         hk_720, hk_755, hl_944, ii0_756, ii1_756, ik_971, \
                         ik_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = f_13 * hk_755[k]
                    + pb_y[k] * ik_971[k];

        t_1214[k] = pa_y[k] * hl_944[k];

        t_1215[k] = f_1 * ii0_756[k]
                    - f_2 * ii1_756[k]
                    + pb_x[k] * ik_972[k];

        t_1216[k] = pb_y[k] * ik_972[k];

        t_1217[k] = f_0 * hk_720[k]
                    + pb_z[k] * ik_972[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, t_1221, pb_x, pb_y, ii0_759, ii0_761, \
                         ii0_762, ii1_759, ii1_761, ii1_762, ik_974, ik_975, ik_977, \
                         ik_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_11 * ii0_759[k]
                    - f_12 * ii1_759[k]
                    + pb_x[k] * ik_975[k];

        t_1219[k] = pb_y[k] * ik_974[k];

        t_1220[k] = f_11 * ii0_761[k]
                    - f_12 * ii1_761[k]
                    + pb_x[k] * ik_977[k];

        t_1221[k] = f_9 * ii0_762[k]
                    - f_10 * ii1_762[k]
                    + pb_x[k] * ik_978[k];
    }
}

static auto
compute_prim_il_electron_repulsion_0_piece9(CSimdMatrix &buffer, const size_t target,
                                            const size_t pb, const size_t hk, const size_t ii0,
                                            const size_t ii1, const size_t ik,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.5 / beta;
    const auto f_2 = 3.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 2.5 / beta;
    const auto f_12 = 2.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hk_723 = buffer.data(hk + 723);
    const auto *hk_726 = buffer.data(hk + 726);
    const auto *hk_730 = buffer.data(hk + 730);
    const auto *hk_735 = buffer.data(hk + 735);
    const auto *hk_748 = buffer.data(hk + 748);
    const auto *hk_755 = buffer.data(hk + 755);

    const auto *ii0_765 = buffer.data(ii0 + 765);
    const auto *ii0_766 = buffer.data(ii0 + 766);
    const auto *ii0_768 = buffer.data(ii0 + 768);
    const auto *ii0_770 = buffer.data(ii0 + 770);
    const auto *ii0_771 = buffer.data(ii0 + 771);
    const auto *ii0_773 = buffer.data(ii0 + 773);
    const auto *ii0_774 = buffer.data(ii0 + 774);
    const auto *ii0_776 = buffer.data(ii0 + 776);
    const auto *ii0_777 = buffer.data(ii0 + 777);
    const auto *ii0_779 = buffer.data(ii0 + 779);
    const auto *ii0_780 = buffer.data(ii0 + 780);
    const auto *ii0_781 = buffer.data(ii0 + 781);
    const auto *ii0_782 = buffer.data(ii0 + 782);
    const auto *ii0_783 = buffer.data(ii0 + 783);

    const auto *ii1_765 = buffer.data(ii1 + 765);
    const auto *ii1_766 = buffer.data(ii1 + 766);
    const auto *ii1_768 = buffer.data(ii1 + 768);
    const auto *ii1_770 = buffer.data(ii1 + 770);
    const auto *ii1_771 = buffer.data(ii1 + 771);
    const auto *ii1_773 = buffer.data(ii1 + 773);
    const auto *ii1_774 = buffer.data(ii1 + 774);
    const auto *ii1_776 = buffer.data(ii1 + 776);
    const auto *ii1_777 = buffer.data(ii1 + 777);
    const auto *ii1_779 = buffer.data(ii1 + 779);
    const auto *ii1_780 = buffer.data(ii1 + 780);
    const auto *ii1_781 = buffer.data(ii1 + 781);
    const auto *ii1_782 = buffer.data(ii1 + 782);
    const auto *ii1_783 = buffer.data(ii1 + 783);

    const auto *ik_975 = buffer.data(ik + 975);
    const auto *ik_977 = buffer.data(ik + 977);
    const auto *ik_978 = buffer.data(ik + 978);
    const auto *ik_981 = buffer.data(ik + 981);
    const auto *ik_982 = buffer.data(ik + 982);
    const auto *ik_984 = buffer.data(ik + 984);
    const auto *ik_986 = buffer.data(ik + 986);
    const auto *ik_987 = buffer.data(ik + 987);
    const auto *ik_989 = buffer.data(ik + 989);
    const auto *ik_990 = buffer.data(ik + 990);
    const auto *ik_992 = buffer.data(ik + 992);
    const auto *ik_993 = buffer.data(ik + 993);
    const auto *ik_995 = buffer.data(ik + 995);
    const auto *ik_996 = buffer.data(ik + 996);
    const auto *ik_997 = buffer.data(ik + 997);
    const auto *ik_999 = buffer.data(ik + 999);
    const auto *ik_1000 = buffer.data(ik + 1000);
    const auto *ik_1001 = buffer.data(ik + 1001);
    const auto *ik_1002 = buffer.data(ik + 1002);
    const auto *ik_1003 = buffer.data(ik + 1003);
    const auto *ik_1004 = buffer.data(ik + 1004);
    const auto *ik_1005 = buffer.data(ik + 1005);
    const auto *ik_1006 = buffer.data(ik + 1006);
    const auto *ik_1007 = buffer.data(ik + 1007);

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pb_x, pb_y, pb_z, hk_723, ii0_765, \
                         ii0_766, ii1_765, ii1_766, ik_975, ik_977, ik_981, \
                         ik_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_0 * hk_723[k]
                    + pb_z[k] * ik_975[k];

        t_1223[k] = pb_y[k] * ik_977[k];

        t_1224[k] = f_9 * ii0_765[k]
                    - f_10 * ii1_765[k]
                    + pb_x[k] * ik_981[k];

        t_1225[k] = f_7 * ii0_766[k]
                    - f_8 * ii1_766[k]
                    + pb_x[k] * ik_982[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, t_1229, pb_x, pb_y, pb_z, hk_726, ii0_768, \
                         ii0_770, ii1_768, ii1_770, ik_978, ik_981, ik_984, \
                         ik_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_0 * hk_726[k]
                    + pb_z[k] * ik_978[k];

        t_1227[k] = f_7 * ii0_768[k]
                    - f_8 * ii1_768[k]
                    + pb_x[k] * ik_984[k];

        t_1228[k] = pb_y[k] * ik_981[k];

        t_1229[k] = f_7 * ii0_770[k]
                    - f_8 * ii1_770[k]
                    + pb_x[k] * ik_986[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pb_x, pb_z, hk_730, ii0_771, ii0_773, \
                         ii1_771, ii1_773, ik_982, ik_987, ik_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_5 * ii0_771[k]
                    - f_6 * ii1_771[k]
                    + pb_x[k] * ik_987[k];

        t_1231[k] = f_0 * hk_730[k]
                    + pb_z[k] * ik_982[k];

        t_1232[k] = f_5 * ii0_773[k]
                    - f_6 * ii1_773[k]
                    + pb_x[k] * ik_989[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, t_1236, pb_x, pb_y, ii0_774, ii0_776, \
                         ii0_777, ii1_774, ii1_776, ii1_777, ik_986, ik_990, ik_992, \
                         ik_993 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_5 * ii0_774[k]
                    - f_6 * ii1_774[k]
                    + pb_x[k] * ik_990[k];

        t_1234[k] = pb_y[k] * ik_986[k];

        t_1235[k] = f_5 * ii0_776[k]
                    - f_6 * ii1_776[k]
                    + pb_x[k] * ik_992[k];

        t_1236[k] = f_3 * ii0_777[k]
                    - f_4 * ii1_777[k]
                    + pb_x[k] * ik_993[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, pb_x, pb_z, hk_735, ii0_779, ii0_780, \
                         ii1_779, ii1_780, ik_987, ik_995, ik_996 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_0 * hk_735[k]
                    + pb_z[k] * ik_987[k];

        t_1238[k] = f_3 * ii0_779[k]
                    - f_4 * ii1_779[k]
                    + pb_x[k] * ik_995[k];

        t_1239[k] = f_3 * ii0_780[k]
                    - f_4 * ii1_780[k]
                    + pb_x[k] * ik_996[k];
    }

#pragma omp simd aligned(t_1240, t_1241, t_1242, t_1243, t_1244, pb_x, pb_y, ii0_781, ii0_783, \
                         ii1_781, ii1_783, ik_992, ik_997, ik_999, ik_1000, \
                         ik_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1240[k] = f_3 * ii0_781[k]
                    - f_4 * ii1_781[k]
                    + pb_x[k] * ik_997[k];

        t_1241[k] = pb_y[k] * ik_992[k];

        t_1242[k] = f_3 * ii0_783[k]
                    - f_4 * ii1_783[k]
                    + pb_x[k] * ik_999[k];

        t_1243[k] = pb_x[k] * ik_1000[k];

        t_1244[k] = pb_x[k] * ik_1001[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, t_1250, pb_x, ik_1002, \
                         ik_1003, ik_1004, ik_1005, ik_1006, ik_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = pb_x[k] * ik_1002[k];

        t_1246[k] = pb_x[k] * ik_1003[k];

        t_1247[k] = pb_x[k] * ik_1004[k];

        t_1248[k] = pb_x[k] * ik_1005[k];

        t_1249[k] = pb_x[k] * ik_1006[k];

        t_1250[k] = pb_x[k] * ik_1007[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, pb_y, pb_z, hk_748, ii0_777, ii0_779, \
                         ii0_780, ii1_777, ii1_779, ii1_780, ik_1000, ik_1002, \
                         ik_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_1 * ii0_777[k]
                    - f_2 * ii1_777[k]
                    + pb_y[k] * ik_1000[k];

        t_1252[k] = f_0 * hk_748[k]
                    + pb_z[k] * ik_1000[k];

        t_1253[k] = f_11 * ii0_779[k]
                    - f_12 * ii1_779[k]
                    + pb_y[k] * ik_1002[k];

        t_1254[k] = f_9 * ii0_780[k]
                    - f_10 * ii1_780[k]
                    + pb_y[k] * ik_1003[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, t_1258, pb_y, ii0_781, ii0_782, ii0_783, \
                         ii1_781, ii1_782, ii1_783, ik_1004, ik_1005, ik_1006, \
                         ik_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = f_7 * ii0_781[k]
                    - f_8 * ii1_781[k]
                    + pb_y[k] * ik_1004[k];

        t_1256[k] = f_5 * ii0_782[k]
                    - f_6 * ii1_782[k]
                    + pb_y[k] * ik_1005[k];

        t_1257[k] = f_3 * ii0_783[k]
                    - f_4 * ii1_783[k]
                    + pb_y[k] * ik_1006[k];

        t_1258[k] = pb_y[k] * ik_1007[k];
    }

#pragma omp simd aligned(t_1259, pb_z, hk_755, ii0_783, ii1_783, \
                         ik_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_0 * hk_755[k]
                    + f_1 * ii0_783[k]
                    - f_2 * ii1_783[k]
                    + pb_z[k] * ik_1007[k];
    }
}

auto
compute_prim_il_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gl0, const size_t gl1,
                                     const size_t hk, const size_t hl, const size_t ii0,
                                     const size_t ii1, const size_t ik, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_il_electron_repulsion_0_piece0(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece1(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece2(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece3(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece4(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece5(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece6(buffer, target, pa, pb, hk, hl, ii0, ii1, ik,
                                                ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece7(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece8(buffer, target, pa, pb, gl0, gl1, hk, hl, ii0,
                                                ii1, ik, ncols, alpha, beta, p);

    compute_prim_il_electron_repulsion_0_piece9(buffer, target, pb, hk, ii0, ii1, ik, ncols,
                                                alpha, beta, p);
}

}  // namespace simdt2ceri
