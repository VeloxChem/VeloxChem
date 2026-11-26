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


#include "SimdElectronRepulsionVrrRecML.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_ml_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;

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

    const auto *kl0_0 = buffer.data(kl0 + 0);

    const auto *kl1_0 = buffer.data(kl1 + 0);

    const auto *lk_0 = buffer.data(lk + 0);
    const auto *lk_1 = buffer.data(lk + 1);
    const auto *lk_2 = buffer.data(lk + 2);
    const auto *lk_3 = buffer.data(lk + 3);
    const auto *lk_5 = buffer.data(lk + 5);
    const auto *lk_6 = buffer.data(lk + 6);
    const auto *lk_7 = buffer.data(lk + 7);
    const auto *lk_8 = buffer.data(lk + 8);
    const auto *lk_9 = buffer.data(lk + 9);
    const auto *lk_10 = buffer.data(lk + 10);
    const auto *lk_11 = buffer.data(lk + 11);
    const auto *lk_12 = buffer.data(lk + 12);
    const auto *lk_13 = buffer.data(lk + 13);
    const auto *lk_14 = buffer.data(lk + 14);
    const auto *lk_15 = buffer.data(lk + 15);
    const auto *lk_16 = buffer.data(lk + 16);
    const auto *lk_17 = buffer.data(lk + 17);
    const auto *lk_18 = buffer.data(lk + 18);
    const auto *lk_19 = buffer.data(lk + 19);
    const auto *lk_20 = buffer.data(lk + 20);
    const auto *lk_28 = buffer.data(lk + 28);
    const auto *lk_29 = buffer.data(lk + 29);
    const auto *lk_30 = buffer.data(lk + 30);
    const auto *lk_31 = buffer.data(lk + 31);
    const auto *lk_32 = buffer.data(lk + 32);
    const auto *lk_33 = buffer.data(lk + 33);
    const auto *lk_34 = buffer.data(lk + 34);
    const auto *lk_35 = buffer.data(lk + 35);
    const auto *lk_36 = buffer.data(lk + 36);
    const auto *lk_41 = buffer.data(lk + 41);
    const auto *lk_64 = buffer.data(lk + 64);
    const auto *lk_66 = buffer.data(lk + 66);
    const auto *lk_67 = buffer.data(lk + 67);
    const auto *lk_68 = buffer.data(lk + 68);
    const auto *lk_69 = buffer.data(lk + 69);
    const auto *lk_70 = buffer.data(lk + 70);
    const auto *lk_101 = buffer.data(lk + 101);
    const auto *lk_102 = buffer.data(lk + 102);
    const auto *lk_103 = buffer.data(lk + 103);
    const auto *lk_104 = buffer.data(lk + 104);
    const auto *lk_105 = buffer.data(lk + 105);
    const auto *lk_107 = buffer.data(lk + 107);
    const auto *lk_111 = buffer.data(lk + 111);
    const auto *lk_114 = buffer.data(lk + 114);

    const auto *ll_0 = buffer.data(ll + 0);
    const auto *ll_3 = buffer.data(ll + 3);
    const auto *ll_5 = buffer.data(ll + 5);
    const auto *ll_6 = buffer.data(ll + 6);
    const auto *ll_9 = buffer.data(ll + 9);
    const auto *ll_10 = buffer.data(ll + 10);
    const auto *ll_12 = buffer.data(ll + 12);
    const auto *ll_14 = buffer.data(ll + 14);
    const auto *ll_15 = buffer.data(ll + 15);
    const auto *ll_17 = buffer.data(ll + 17);
    const auto *ll_18 = buffer.data(ll + 18);
    const auto *ll_20 = buffer.data(ll + 20);
    const auto *ll_21 = buffer.data(ll + 21);
    const auto *ll_23 = buffer.data(ll + 23);
    const auto *ll_24 = buffer.data(ll + 24);
    const auto *ll_25 = buffer.data(ll + 25);
    const auto *ll_27 = buffer.data(ll + 27);
    const auto *ll_28 = buffer.data(ll + 28);
    const auto *ll_35 = buffer.data(ll + 35);
    const auto *ll_36 = buffer.data(ll + 36);
    const auto *ll_38 = buffer.data(ll + 38);
    const auto *ll_39 = buffer.data(ll + 39);
    const auto *ll_40 = buffer.data(ll + 40);
    const auto *ll_41 = buffer.data(ll + 41);
    const auto *ll_42 = buffer.data(ll + 42);
    const auto *ll_44 = buffer.data(ll + 44);
    const auto *ll_45 = buffer.data(ll + 45);

    const auto *mi0_0 = buffer.data(mi0 + 0);
    const auto *mi0_1 = buffer.data(mi0 + 1);
    const auto *mi0_2 = buffer.data(mi0 + 2);
    const auto *mi0_3 = buffer.data(mi0 + 3);
    const auto *mi0_5 = buffer.data(mi0 + 5);
    const auto *mi0_6 = buffer.data(mi0 + 6);
    const auto *mi0_8 = buffer.data(mi0 + 8);
    const auto *mi0_9 = buffer.data(mi0 + 9);
    const auto *mi0_10 = buffer.data(mi0 + 10);
    const auto *mi0_12 = buffer.data(mi0 + 12);
    const auto *mi0_13 = buffer.data(mi0 + 13);
    const auto *mi0_14 = buffer.data(mi0 + 14);
    const auto *mi0_21 = buffer.data(mi0 + 21);
    const auto *mi0_23 = buffer.data(mi0 + 23);
    const auto *mi0_24 = buffer.data(mi0 + 24);
    const auto *mi0_25 = buffer.data(mi0 + 25);
    const auto *mi0_26 = buffer.data(mi0 + 26);
    const auto *mi0_27 = buffer.data(mi0 + 27);
    const auto *mi0_84 = buffer.data(mi0 + 84);
    const auto *mi0_86 = buffer.data(mi0 + 86);
    const auto *mi0_87 = buffer.data(mi0 + 87);
    const auto *mi0_90 = buffer.data(mi0 + 90);

    const auto *mi1_0 = buffer.data(mi1 + 0);
    const auto *mi1_1 = buffer.data(mi1 + 1);
    const auto *mi1_2 = buffer.data(mi1 + 2);
    const auto *mi1_3 = buffer.data(mi1 + 3);
    const auto *mi1_5 = buffer.data(mi1 + 5);
    const auto *mi1_6 = buffer.data(mi1 + 6);
    const auto *mi1_8 = buffer.data(mi1 + 8);
    const auto *mi1_9 = buffer.data(mi1 + 9);
    const auto *mi1_10 = buffer.data(mi1 + 10);
    const auto *mi1_12 = buffer.data(mi1 + 12);
    const auto *mi1_13 = buffer.data(mi1 + 13);
    const auto *mi1_14 = buffer.data(mi1 + 14);
    const auto *mi1_21 = buffer.data(mi1 + 21);
    const auto *mi1_23 = buffer.data(mi1 + 23);
    const auto *mi1_24 = buffer.data(mi1 + 24);
    const auto *mi1_25 = buffer.data(mi1 + 25);
    const auto *mi1_26 = buffer.data(mi1 + 26);
    const auto *mi1_27 = buffer.data(mi1 + 27);
    const auto *mi1_84 = buffer.data(mi1 + 84);
    const auto *mi1_86 = buffer.data(mi1 + 86);
    const auto *mi1_87 = buffer.data(mi1 + 87);
    const auto *mi1_90 = buffer.data(mi1 + 90);

    const auto *mk_0 = buffer.data(mk + 0);
    const auto *mk_1 = buffer.data(mk + 1);
    const auto *mk_2 = buffer.data(mk + 2);
    const auto *mk_3 = buffer.data(mk + 3);
    const auto *mk_5 = buffer.data(mk + 5);
    const auto *mk_6 = buffer.data(mk + 6);
    const auto *mk_8 = buffer.data(mk + 8);
    const auto *mk_9 = buffer.data(mk + 9);
    const auto *mk_10 = buffer.data(mk + 10);
    const auto *mk_12 = buffer.data(mk + 12);
    const auto *mk_13 = buffer.data(mk + 13);
    const auto *mk_14 = buffer.data(mk + 14);
    const auto *mk_15 = buffer.data(mk + 15);
    const auto *mk_17 = buffer.data(mk + 17);
    const auto *mk_18 = buffer.data(mk + 18);
    const auto *mk_19 = buffer.data(mk + 19);
    const auto *mk_20 = buffer.data(mk + 20);
    const auto *mk_21 = buffer.data(mk + 21);
    const auto *mk_27 = buffer.data(mk + 27);
    const auto *mk_28 = buffer.data(mk + 28);
    const auto *mk_30 = buffer.data(mk + 30);
    const auto *mk_31 = buffer.data(mk + 31);
    const auto *mk_32 = buffer.data(mk + 32);
    const auto *mk_33 = buffer.data(mk + 33);
    const auto *mk_34 = buffer.data(mk + 34);
    const auto *mk_35 = buffer.data(mk + 35);
    const auto *mk_36 = buffer.data(mk + 36);
    const auto *mk_37 = buffer.data(mk + 37);
    const auto *mk_39 = buffer.data(mk + 39);
    const auto *mk_41 = buffer.data(mk + 41);
    const auto *mk_42 = buffer.data(mk + 42);
    const auto *mk_45 = buffer.data(mk + 45);
    const auto *mk_46 = buffer.data(mk + 46);
    const auto *mk_50 = buffer.data(mk + 50);
    const auto *mk_51 = buffer.data(mk + 51);
    const auto *mk_56 = buffer.data(mk + 56);
    const auto *mk_57 = buffer.data(mk + 57);
    const auto *mk_64 = buffer.data(mk + 64);
    const auto *mk_66 = buffer.data(mk + 66);
    const auto *mk_67 = buffer.data(mk + 67);
    const auto *mk_68 = buffer.data(mk + 68);
    const auto *mk_69 = buffer.data(mk + 69);
    const auto *mk_70 = buffer.data(mk + 70);
    const auto *mk_71 = buffer.data(mk + 71);
    const auto *mk_72 = buffer.data(mk + 72);
    const auto *mk_74 = buffer.data(mk + 74);
    const auto *mk_75 = buffer.data(mk + 75);
    const auto *mk_77 = buffer.data(mk + 77);
    const auto *mk_78 = buffer.data(mk + 78);
    const auto *mk_81 = buffer.data(mk + 81);
    const auto *mk_82 = buffer.data(mk + 82);
    const auto *mk_86 = buffer.data(mk + 86);
    const auto *mk_87 = buffer.data(mk + 87);
    const auto *mk_92 = buffer.data(mk + 92);
    const auto *mk_99 = buffer.data(mk + 99);
    const auto *mk_100 = buffer.data(mk + 100);
    const auto *mk_101 = buffer.data(mk + 101);
    const auto *mk_102 = buffer.data(mk + 102);
    const auto *mk_103 = buffer.data(mk + 103);
    const auto *mk_104 = buffer.data(mk + 104);
    const auto *mk_105 = buffer.data(mk + 105);
    const auto *mk_107 = buffer.data(mk + 107);
    const auto *mk_108 = buffer.data(mk + 108);
    const auto *mk_109 = buffer.data(mk + 109);
    const auto *mk_110 = buffer.data(mk + 110);
    const auto *mk_111 = buffer.data(mk + 111);
    const auto *mk_113 = buffer.data(mk + 113);
    const auto *mk_114 = buffer.data(mk + 114);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, lk_0, mi0_0, mi1_0, \
                         mk_0, mk_1, mk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lk_0[k]
                 + f_1 * mi0_0[k]
                 - f_2 * mi1_0[k]
                 + pb_x[k] * mk_0[k];

        t_1[k] = pb_y[k] * mk_0[k];

        t_2[k] = pb_z[k] * mk_0[k];

        t_3[k] = f_3 * mi0_0[k]
                 - f_4 * mi1_0[k]
                 + pb_y[k] * mk_1[k];

        t_4[k] = pb_y[k] * mk_2[k];

        t_5[k] = f_3 * mi0_0[k]
                 - f_4 * mi1_0[k]
                 + pb_z[k] * mk_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, mi0_1, mi0_2, mi0_3, mi1_1, \
                         mi1_2, mi1_3, mk_3, mk_5, mk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * mi0_1[k]
                 - f_6 * mi1_1[k]
                 + pb_y[k] * mk_3[k];

        t_7[k] = pb_z[k] * mk_3[k];

        t_8[k] = pb_y[k] * mk_5[k];

        t_9[k] = f_5 * mi0_2[k]
                 - f_6 * mi1_2[k]
                 + pb_z[k] * mk_5[k];

        t_10[k] = f_7 * mi0_3[k]
                  - f_8 * mi1_3[k]
                  + pb_y[k] * mk_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, mi0_5, mi0_6, mi1_5, \
                         mi1_6, mk_6, mk_8, mk_9, mk_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * mk_6[k];

        t_12[k] = f_3 * mi0_5[k]
                  - f_4 * mi1_5[k]
                  + pb_y[k] * mk_8[k];

        t_13[k] = pb_y[k] * mk_9[k];

        t_14[k] = f_7 * mi0_5[k]
                  - f_8 * mi1_5[k]
                  + pb_z[k] * mk_9[k];

        t_15[k] = f_9 * mi0_6[k]
                  - f_10 * mi1_6[k]
                  + pb_y[k] * mk_10[k];

        t_16[k] = pb_z[k] * mk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, mi0_8, mi0_9, mi1_8, mi1_9, \
                         mk_12, mk_13, mk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * mi0_8[k]
                  - f_6 * mi1_8[k]
                  + pb_y[k] * mk_12[k];

        t_18[k] = f_3 * mi0_9[k]
                  - f_4 * mi1_9[k]
                  + pb_y[k] * mk_13[k];

        t_19[k] = pb_y[k] * mk_14[k];

        t_20[k] = f_9 * mi0_9[k]
                  - f_10 * mi1_9[k]
                  + pb_z[k] * mk_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, mi0_10, mi0_12, mi0_13, mi1_10, \
                         mi1_12, mi1_13, mk_15, mk_17, mk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * mi0_10[k]
                  - f_12 * mi1_10[k]
                  + pb_y[k] * mk_15[k];

        t_22[k] = pb_z[k] * mk_15[k];

        t_23[k] = f_7 * mi0_12[k]
                  - f_8 * mi1_12[k]
                  + pb_y[k] * mk_17[k];

        t_24[k] = f_5 * mi0_13[k]
                  - f_6 * mi1_13[k]
                  + pb_y[k] * mk_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, lk_28, mi0_14, \
                         mi1_14, mk_19, mk_20, mk_21, mk_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * mi0_14[k]
                  - f_4 * mi1_14[k]
                  + pb_y[k] * mk_19[k];

        t_26[k] = pb_y[k] * mk_20[k];

        t_27[k] = f_11 * mi0_14[k]
                  - f_12 * mi1_14[k]
                  + pb_z[k] * mk_20[k];

        t_28[k] = f_0 * lk_28[k]
                  + pb_x[k] * mk_28[k];

        t_29[k] = pb_z[k] * mk_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_x, pb_y, lk_30, lk_31, lk_32, lk_33, \
                         mk_27, mk_30, mk_31, mk_32, mk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * lk_30[k]
                  + pb_x[k] * mk_30[k];

        t_31[k] = f_0 * lk_31[k]
                  + pb_x[k] * mk_31[k];

        t_32[k] = f_0 * lk_32[k]
                  + pb_x[k] * mk_32[k];

        t_33[k] = f_0 * lk_33[k]
                  + pb_x[k] * mk_33[k];

        t_34[k] = pb_y[k] * mk_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, lk_35, mi0_21, mi0_23, \
                         mi1_21, mi1_23, mk_28, mk_30, mk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * lk_35[k]
                  + pb_x[k] * mk_35[k];

        t_36[k] = f_1 * mi0_21[k]
                  - f_2 * mi1_21[k]
                  + pb_y[k] * mk_28[k];

        t_37[k] = pb_z[k] * mk_28[k];

        t_38[k] = f_11 * mi0_23[k]
                  - f_12 * mi1_23[k]
                  + pb_y[k] * mk_30[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, mi0_24, mi0_25, mi0_26, mi1_24, mi1_25, \
                         mi1_26, mk_31, mk_32, mk_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * mi0_24[k]
                  - f_10 * mi1_24[k]
                  + pb_y[k] * mk_31[k];

        t_40[k] = f_7 * mi0_25[k]
                  - f_8 * mi1_25[k]
                  + pb_y[k] * mk_32[k];

        t_41[k] = f_5 * mi0_26[k]
                  - f_6 * mi1_26[k]
                  + pb_y[k] * mk_33[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, t_47, pa_y, pb_y, pb_z, lk_0, ll_0, \
                         mi0_27, mi1_27, mk_34, mk_35, mk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * mi0_27[k]
                  - f_4 * mi1_27[k]
                  + pb_y[k] * mk_34[k];

        t_43[k] = pb_y[k] * mk_35[k];

        t_44[k] = f_1 * mi0_27[k]
                  - f_2 * mi1_27[k]
                  + pb_z[k] * mk_35[k];

        t_45[k] = pa_y[k] * ll_0[k];

        t_46[k] = f_13 * lk_0[k]
                  + pb_y[k] * mk_36[k];

        t_47[k] = pb_z[k] * mk_36[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pb_z, lk_1, lk_3, ll_3, ll_5, \
                         ll_6, mk_37, mk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * lk_1[k]
                  + pa_y[k] * ll_3[k];

        t_49[k] = pb_z[k] * mk_37[k];

        t_50[k] = pa_y[k] * ll_5[k];

        t_51[k] = f_15 * lk_3[k]
                  + pa_y[k] * ll_6[k];

        t_52[k] = pb_z[k] * mk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pb_y, pb_z, lk_5, lk_6, lk_8, \
                         ll_9, ll_10, ll_12, mk_41, mk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_13 * lk_5[k]
                  + pb_y[k] * mk_41[k];

        t_54[k] = pa_y[k] * ll_9[k];

        t_55[k] = f_16 * lk_6[k]
                  + pa_y[k] * ll_10[k];

        t_56[k] = pb_z[k] * mk_42[k];

        t_57[k] = f_14 * lk_8[k]
                  + pa_y[k] * ll_12[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_y, pb_y, pb_z, lk_9, lk_10, lk_12, \
                         ll_14, ll_15, ll_17, mk_45, mk_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_13 * lk_9[k]
                  + pb_y[k] * mk_45[k];

        t_59[k] = pa_y[k] * ll_14[k];

        t_60[k] = f_17 * lk_10[k]
                  + pa_y[k] * ll_15[k];

        t_61[k] = pb_z[k] * mk_46[k];

        t_62[k] = f_15 * lk_12[k]
                  + pa_y[k] * ll_17[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_y, pb_z, lk_13, lk_14, lk_15, \
                         ll_18, ll_20, ll_21, mk_50, mk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_14 * lk_13[k]
                  + pa_y[k] * ll_18[k];

        t_64[k] = f_13 * lk_14[k]
                  + pb_y[k] * mk_50[k];

        t_65[k] = pa_y[k] * ll_20[k];

        t_66[k] = f_18 * lk_15[k]
                  + pa_y[k] * ll_21[k];

        t_67[k] = pb_z[k] * mk_51[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pb_y, lk_17, lk_18, lk_19, lk_20, \
                         ll_23, ll_24, ll_25, ll_27, mk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_16 * lk_17[k]
                  + pa_y[k] * ll_23[k];

        t_69[k] = f_15 * lk_18[k]
                  + pa_y[k] * ll_24[k];

        t_70[k] = f_14 * lk_19[k]
                  + pa_y[k] * ll_25[k];

        t_71[k] = f_13 * lk_20[k]
                  + pb_y[k] * mk_56[k];

        t_72[k] = pa_y[k] * ll_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, lk_64, lk_66, lk_67, lk_68, \
                         mk_57, mk_64, mk_66, mk_67, mk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_19 * lk_64[k]
                  + pb_x[k] * mk_64[k];

        t_74[k] = pb_z[k] * mk_57[k];

        t_75[k] = f_19 * lk_66[k]
                  + pb_x[k] * mk_66[k];

        t_76[k] = f_19 * lk_67[k]
                  + pb_x[k] * mk_67[k];

        t_77[k] = f_19 * lk_68[k]
                  + pb_x[k] * mk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pb_x, pb_z, lk_28, lk_69, lk_70, \
                         ll_35, ll_36, mk_64, mk_69, mk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * lk_69[k]
                  + pb_x[k] * mk_69[k];

        t_79[k] = f_19 * lk_70[k]
                  + pb_x[k] * mk_70[k];

        t_80[k] = pa_y[k] * ll_35[k];

        t_81[k] = f_19 * lk_28[k]
                  + pa_y[k] * ll_36[k];

        t_82[k] = pb_z[k] * mk_64[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, lk_30, lk_31, lk_32, lk_33, \
                         lk_34, ll_38, ll_39, ll_40, ll_41, ll_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_18 * lk_30[k]
                  + pa_y[k] * ll_38[k];

        t_84[k] = f_17 * lk_31[k]
                  + pa_y[k] * ll_39[k];

        t_85[k] = f_16 * lk_32[k]
                  + pa_y[k] * ll_40[k];

        t_86[k] = f_15 * lk_33[k]
                  + pa_y[k] * ll_41[k];

        t_87[k] = f_14 * lk_34[k]
                  + pa_y[k] * ll_42[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, lk_0, lk_35, \
                         ll_0, ll_44, mk_71, mk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_13 * lk_35[k]
                  + pb_y[k] * mk_71[k];

        t_89[k] = pa_y[k] * ll_44[k];

        t_90[k] = pa_z[k] * ll_0[k];

        t_91[k] = pb_y[k] * mk_72[k];

        t_92[k] = f_13 * lk_0[k]
                  + pb_z[k] * mk_72[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, lk_2, lk_3, ll_3, \
                         ll_5, ll_6, mk_74, mk_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_z[k] * ll_3[k];

        t_94[k] = pb_y[k] * mk_74[k];

        t_95[k] = f_14 * lk_2[k]
                  + pa_z[k] * ll_5[k];

        t_96[k] = pa_z[k] * ll_6[k];

        t_97[k] = f_13 * lk_3[k]
                  + pb_z[k] * mk_75[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, lk_5, lk_6, lk_7, \
                         ll_9, ll_10, ll_12, mk_77, mk_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pb_y[k] * mk_77[k];

        t_99[k] = f_15 * lk_5[k]
                  + pa_z[k] * ll_9[k];

        t_100[k] = pa_z[k] * ll_10[k];

        t_101[k] = f_13 * lk_6[k]
                   + pb_z[k] * mk_78[k];

        t_102[k] = f_14 * lk_7[k]
                   + pa_z[k] * ll_12[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, lk_9, lk_10, \
                         lk_11, ll_14, ll_15, ll_17, mk_81, mk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * mk_81[k];

        t_104[k] = f_16 * lk_9[k]
                   + pa_z[k] * ll_14[k];

        t_105[k] = pa_z[k] * ll_15[k];

        t_106[k] = f_13 * lk_10[k]
                   + pb_z[k] * mk_82[k];

        t_107[k] = f_14 * lk_11[k]
                   + pa_z[k] * ll_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, lk_12, lk_14, \
                         lk_15, ll_18, ll_20, ll_21, mk_86, mk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_15 * lk_12[k]
                   + pa_z[k] * ll_18[k];

        t_109[k] = pb_y[k] * mk_86[k];

        t_110[k] = f_17 * lk_14[k]
                   + pa_z[k] * ll_20[k];

        t_111[k] = pa_z[k] * ll_21[k];

        t_112[k] = f_13 * lk_15[k]
                   + pb_z[k] * mk_87[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_y, lk_16, lk_17, lk_18, \
                         lk_20, ll_23, ll_24, ll_25, ll_27, mk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_14 * lk_16[k]
                   + pa_z[k] * ll_23[k];

        t_114[k] = f_15 * lk_17[k]
                   + pa_z[k] * ll_24[k];

        t_115[k] = f_16 * lk_18[k]
                   + pa_z[k] * ll_25[k];

        t_116[k] = pb_y[k] * mk_92[k];

        t_117[k] = f_18 * lk_20[k]
                   + pa_z[k] * ll_27[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pa_z, pb_x, lk_101, lk_102, \
                         lk_103, lk_104, ll_28, mk_101, mk_102, mk_103, \
                         mk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * ll_28[k];

        t_119[k] = f_19 * lk_101[k]
                   + pb_x[k] * mk_101[k];

        t_120[k] = f_19 * lk_102[k]
                   + pb_x[k] * mk_102[k];

        t_121[k] = f_19 * lk_103[k]
                   + pb_x[k] * mk_103[k];

        t_122[k] = f_19 * lk_104[k]
                   + pb_x[k] * mk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_z, pb_x, pb_y, lk_105, lk_107, ll_36, \
                         mk_99, mk_105, mk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_19 * lk_105[k]
                   + pb_x[k] * mk_105[k];

        t_124[k] = pb_y[k] * mk_99[k];

        t_125[k] = f_19 * lk_107[k]
                   + pb_x[k] * mk_107[k];

        t_126[k] = pa_z[k] * ll_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_z, pb_z, lk_28, lk_29, lk_30, lk_31, \
                         ll_38, ll_39, ll_40, mk_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_13 * lk_28[k]
                   + pb_z[k] * mk_100[k];

        t_128[k] = f_14 * lk_29[k]
                   + pa_z[k] * ll_38[k];

        t_129[k] = f_15 * lk_30[k]
                   + pa_z[k] * ll_39[k];

        t_130[k] = f_16 * lk_31[k]
                   + pa_z[k] * ll_40[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_z, pb_y, lk_32, lk_33, lk_35, ll_41, \
                         ll_42, ll_44, mk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_17 * lk_32[k]
                   + pa_z[k] * ll_41[k];

        t_132[k] = f_18 * lk_33[k]
                   + pa_z[k] * ll_42[k];

        t_133[k] = pb_y[k] * mk_107[k];

        t_134[k] = f_19 * lk_35[k]
                   + pa_z[k] * ll_44[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pb_y, pb_z, kl0_0, kl1_0, lk_36, ll_45, \
                         mk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_20 * kl0_0[k]
                   - f_21 * kl1_0[k]
                   + pa_y[k] * ll_45[k];

        t_136[k] = f_14 * lk_36[k]
                   + pb_y[k] * mk_108[k];

        t_137[k] = pb_z[k] * mk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_x, pb_z, lk_111, mi0_84, mi0_87, mi1_84, \
                         mi1_87, mk_109, mk_110, mk_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_22 * lk_111[k]
                   + f_11 * mi0_87[k]
                   - f_12 * mi1_87[k]
                   + pb_x[k] * mk_111[k];

        t_139[k] = pb_z[k] * mk_109[k];

        t_140[k] = f_3 * mi0_84[k]
                   - f_4 * mi1_84[k]
                   + pb_z[k] * mk_110[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_x, pb_y, pb_z, lk_41, lk_114, mi0_86, \
                         mi0_90, mi1_86, mi1_90, mk_111, mk_113, \
                         mk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_22 * lk_114[k]
                   + f_9 * mi0_90[k]
                   - f_10 * mi1_90[k]
                   + pb_x[k] * mk_114[k];

        t_142[k] = pb_z[k] * mk_111[k];

        t_143[k] = f_14 * lk_41[k]
                   + pb_y[k] * mk_113[k];

        t_144[k] = f_5 * mi0_86[k]
                   - f_6 * mi1_86[k]
                   + pb_z[k] * mk_113[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;
    const auto f_23 = 3.0 / alpha;
    const auto f_24 = 3.0 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);

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

    const auto *kl0_0 = buffer.data(kl0 + 0);
    const auto *kl0_45 = buffer.data(kl0 + 45);
    const auto *kl0_171 = buffer.data(kl0 + 171);
    const auto *kl0_269 = buffer.data(kl0 + 269);

    const auto *kl1_0 = buffer.data(kl1 + 0);
    const auto *kl1_45 = buffer.data(kl1 + 45);
    const auto *kl1_171 = buffer.data(kl1 + 171);
    const auto *kl1_269 = buffer.data(kl1 + 269);

    const auto *lk_39 = buffer.data(lk + 39);
    const auto *lk_42 = buffer.data(lk + 42);
    const auto *lk_45 = buffer.data(lk + 45);
    const auto *lk_46 = buffer.data(lk + 46);
    const auto *lk_50 = buffer.data(lk + 50);
    const auto *lk_51 = buffer.data(lk + 51);
    const auto *lk_56 = buffer.data(lk + 56);
    const auto *lk_64 = buffer.data(lk + 64);
    const auto *lk_71 = buffer.data(lk + 71);
    const auto *lk_72 = buffer.data(lk + 72);
    const auto *lk_74 = buffer.data(lk + 74);
    const auto *lk_75 = buffer.data(lk + 75);
    const auto *lk_77 = buffer.data(lk + 77);
    const auto *lk_78 = buffer.data(lk + 78);
    const auto *lk_80 = buffer.data(lk + 80);
    const auto *lk_81 = buffer.data(lk + 81);
    const auto *lk_82 = buffer.data(lk + 82);
    const auto *lk_84 = buffer.data(lk + 84);
    const auto *lk_85 = buffer.data(lk + 85);
    const auto *lk_86 = buffer.data(lk + 86);
    const auto *lk_87 = buffer.data(lk + 87);
    const auto *lk_89 = buffer.data(lk + 89);
    const auto *lk_90 = buffer.data(lk + 90);
    const auto *lk_91 = buffer.data(lk + 91);
    const auto *lk_92 = buffer.data(lk + 92);
    const auto *lk_100 = buffer.data(lk + 100);
    const auto *lk_102 = buffer.data(lk + 102);
    const auto *lk_103 = buffer.data(lk + 103);
    const auto *lk_104 = buffer.data(lk + 104);
    const auto *lk_105 = buffer.data(lk + 105);
    const auto *lk_106 = buffer.data(lk + 106);
    const auto *lk_107 = buffer.data(lk + 107);
    const auto *lk_108 = buffer.data(lk + 108);
    const auto *lk_118 = buffer.data(lk + 118);
    const auto *lk_123 = buffer.data(lk + 123);
    const auto *lk_129 = buffer.data(lk + 129);
    const auto *lk_136 = buffer.data(lk + 136);
    const auto *lk_138 = buffer.data(lk + 138);
    const auto *lk_139 = buffer.data(lk + 139);
    const auto *lk_140 = buffer.data(lk + 140);
    const auto *lk_141 = buffer.data(lk + 141);
    const auto *lk_142 = buffer.data(lk + 142);
    const auto *lk_143 = buffer.data(lk + 143);
    const auto *lk_173 = buffer.data(lk + 173);
    const auto *lk_174 = buffer.data(lk + 174);
    const auto *lk_175 = buffer.data(lk + 175);
    const auto *lk_176 = buffer.data(lk + 176);
    const auto *lk_177 = buffer.data(lk + 177);
    const auto *lk_178 = buffer.data(lk + 178);
    const auto *lk_185 = buffer.data(lk + 185);
    const auto *lk_189 = buffer.data(lk + 189);
    const auto *lk_194 = buffer.data(lk + 194);
    const auto *lk_200 = buffer.data(lk + 200);
    const auto *lk_207 = buffer.data(lk + 207);
    const auto *lk_208 = buffer.data(lk + 208);
    const auto *lk_209 = buffer.data(lk + 209);
    const auto *lk_210 = buffer.data(lk + 210);
    const auto *lk_211 = buffer.data(lk + 211);
    const auto *lk_212 = buffer.data(lk + 212);
    const auto *lk_213 = buffer.data(lk + 213);
    const auto *lk_215 = buffer.data(lk + 215);
    const auto *lk_219 = buffer.data(lk + 219);

    const auto *ll_46 = buffer.data(ll + 46);
    const auto *ll_48 = buffer.data(ll + 48);
    const auto *ll_51 = buffer.data(ll + 51);
    const auto *ll_55 = buffer.data(ll + 55);
    const auto *ll_60 = buffer.data(ll + 60);
    const auto *ll_66 = buffer.data(ll + 66);
    const auto *ll_73 = buffer.data(ll + 73);
    const auto *ll_81 = buffer.data(ll + 81);
    const auto *ll_90 = buffer.data(ll + 90);
    const auto *ll_92 = buffer.data(ll + 92);
    const auto *ll_95 = buffer.data(ll + 95);
    const auto *ll_99 = buffer.data(ll + 99);
    const auto *ll_102 = buffer.data(ll + 102);
    const auto *ll_104 = buffer.data(ll + 104);
    const auto *ll_107 = buffer.data(ll + 107);
    const auto *ll_108 = buffer.data(ll + 108);
    const auto *ll_110 = buffer.data(ll + 110);
    const auto *ll_113 = buffer.data(ll + 113);
    const auto *ll_114 = buffer.data(ll + 114);
    const auto *ll_115 = buffer.data(ll + 115);
    const auto *ll_117 = buffer.data(ll + 117);
    const auto *ll_125 = buffer.data(ll + 125);
    const auto *ll_128 = buffer.data(ll + 128);
    const auto *ll_129 = buffer.data(ll + 129);
    const auto *ll_130 = buffer.data(ll + 130);
    const auto *ll_131 = buffer.data(ll + 131);
    const auto *ll_132 = buffer.data(ll + 132);
    const auto *ll_134 = buffer.data(ll + 134);
    const auto *ll_135 = buffer.data(ll + 135);
    const auto *ll_171 = buffer.data(ll + 171);
    const auto *ll_269 = buffer.data(ll + 269);

    const auto *mi0_87 = buffer.data(mi0 + 87);
    const auto *mi0_89 = buffer.data(mi0 + 89);
    const auto *mi0_90 = buffer.data(mi0 + 90);
    const auto *mi0_91 = buffer.data(mi0 + 91);
    const auto *mi0_93 = buffer.data(mi0 + 93);
    const auto *mi0_94 = buffer.data(mi0 + 94);
    const auto *mi0_95 = buffer.data(mi0 + 95);
    const auto *mi0_96 = buffer.data(mi0 + 96);
    const auto *mi0_98 = buffer.data(mi0 + 98);
    const auto *mi0_99 = buffer.data(mi0 + 99);
    const auto *mi0_105 = buffer.data(mi0 + 105);
    const auto *mi0_106 = buffer.data(mi0 + 106);
    const auto *mi0_107 = buffer.data(mi0 + 107);
    const auto *mi0_108 = buffer.data(mi0 + 108);
    const auto *mi0_109 = buffer.data(mi0 + 109);
    const auto *mi0_111 = buffer.data(mi0 + 111);
    const auto *mi0_140 = buffer.data(mi0 + 140);
    const auto *mi0_141 = buffer.data(mi0 + 141);
    const auto *mi0_143 = buffer.data(mi0 + 143);
    const auto *mi0_145 = buffer.data(mi0 + 145);
    const auto *mi0_146 = buffer.data(mi0 + 146);
    const auto *mi0_148 = buffer.data(mi0 + 148);
    const auto *mi0_149 = buffer.data(mi0 + 149);
    const auto *mi0_150 = buffer.data(mi0 + 150);
    const auto *mi0_152 = buffer.data(mi0 + 152);
    const auto *mi0_153 = buffer.data(mi0 + 153);
    const auto *mi0_154 = buffer.data(mi0 + 154);
    const auto *mi0_160 = buffer.data(mi0 + 160);
    const auto *mi0_161 = buffer.data(mi0 + 161);
    const auto *mi0_163 = buffer.data(mi0 + 163);
    const auto *mi0_164 = buffer.data(mi0 + 164);
    const auto *mi0_165 = buffer.data(mi0 + 165);
    const auto *mi0_166 = buffer.data(mi0 + 166);
    const auto *mi0_167 = buffer.data(mi0 + 167);
    const auto *mi0_168 = buffer.data(mi0 + 168);
    const auto *mi0_171 = buffer.data(mi0 + 171);

    const auto *mi1_87 = buffer.data(mi1 + 87);
    const auto *mi1_89 = buffer.data(mi1 + 89);
    const auto *mi1_90 = buffer.data(mi1 + 90);
    const auto *mi1_91 = buffer.data(mi1 + 91);
    const auto *mi1_93 = buffer.data(mi1 + 93);
    const auto *mi1_94 = buffer.data(mi1 + 94);
    const auto *mi1_95 = buffer.data(mi1 + 95);
    const auto *mi1_96 = buffer.data(mi1 + 96);
    const auto *mi1_98 = buffer.data(mi1 + 98);
    const auto *mi1_99 = buffer.data(mi1 + 99);
    const auto *mi1_105 = buffer.data(mi1 + 105);
    const auto *mi1_106 = buffer.data(mi1 + 106);
    const auto *mi1_107 = buffer.data(mi1 + 107);
    const auto *mi1_108 = buffer.data(mi1 + 108);
    const auto *mi1_109 = buffer.data(mi1 + 109);
    const auto *mi1_111 = buffer.data(mi1 + 111);
    const auto *mi1_140 = buffer.data(mi1 + 140);
    const auto *mi1_141 = buffer.data(mi1 + 141);
    const auto *mi1_143 = buffer.data(mi1 + 143);
    const auto *mi1_145 = buffer.data(mi1 + 145);
    const auto *mi1_146 = buffer.data(mi1 + 146);
    const auto *mi1_148 = buffer.data(mi1 + 148);
    const auto *mi1_149 = buffer.data(mi1 + 149);
    const auto *mi1_150 = buffer.data(mi1 + 150);
    const auto *mi1_152 = buffer.data(mi1 + 152);
    const auto *mi1_153 = buffer.data(mi1 + 153);
    const auto *mi1_154 = buffer.data(mi1 + 154);
    const auto *mi1_160 = buffer.data(mi1 + 160);
    const auto *mi1_161 = buffer.data(mi1 + 161);
    const auto *mi1_163 = buffer.data(mi1 + 163);
    const auto *mi1_164 = buffer.data(mi1 + 164);
    const auto *mi1_165 = buffer.data(mi1 + 165);
    const auto *mi1_166 = buffer.data(mi1 + 166);
    const auto *mi1_167 = buffer.data(mi1 + 167);
    const auto *mi1_168 = buffer.data(mi1 + 168);
    const auto *mi1_171 = buffer.data(mi1 + 171);

    const auto *mk_114 = buffer.data(mk + 114);
    const auto *mk_115 = buffer.data(mk + 115);
    const auto *mk_117 = buffer.data(mk + 117);
    const auto *mk_118 = buffer.data(mk + 118);
    const auto *mk_119 = buffer.data(mk + 119);
    const auto *mk_120 = buffer.data(mk + 120);
    const auto *mk_122 = buffer.data(mk + 122);
    const auto *mk_123 = buffer.data(mk + 123);
    const auto *mk_124 = buffer.data(mk + 124);
    const auto *mk_125 = buffer.data(mk + 125);
    const auto *mk_126 = buffer.data(mk + 126);
    const auto *mk_128 = buffer.data(mk + 128);
    const auto *mk_129 = buffer.data(mk + 129);
    const auto *mk_136 = buffer.data(mk + 136);
    const auto *mk_137 = buffer.data(mk + 137);
    const auto *mk_138 = buffer.data(mk + 138);
    const auto *mk_139 = buffer.data(mk + 139);
    const auto *mk_140 = buffer.data(mk + 140);
    const auto *mk_141 = buffer.data(mk + 141);
    const auto *mk_142 = buffer.data(mk + 142);
    const auto *mk_143 = buffer.data(mk + 143);
    const auto *mk_146 = buffer.data(mk + 146);
    const auto *mk_147 = buffer.data(mk + 147);
    const auto *mk_149 = buffer.data(mk + 149);
    const auto *mk_150 = buffer.data(mk + 150);
    const auto *mk_153 = buffer.data(mk + 153);
    const auto *mk_154 = buffer.data(mk + 154);
    const auto *mk_158 = buffer.data(mk + 158);
    const auto *mk_159 = buffer.data(mk + 159);
    const auto *mk_164 = buffer.data(mk + 164);
    const auto *mk_172 = buffer.data(mk + 172);
    const auto *mk_173 = buffer.data(mk + 173);
    const auto *mk_174 = buffer.data(mk + 174);
    const auto *mk_175 = buffer.data(mk + 175);
    const auto *mk_176 = buffer.data(mk + 176);
    const auto *mk_177 = buffer.data(mk + 177);
    const auto *mk_178 = buffer.data(mk + 178);
    const auto *mk_179 = buffer.data(mk + 179);
    const auto *mk_180 = buffer.data(mk + 180);
    const auto *mk_181 = buffer.data(mk + 181);
    const auto *mk_182 = buffer.data(mk + 182);
    const auto *mk_183 = buffer.data(mk + 183);
    const auto *mk_185 = buffer.data(mk + 185);
    const auto *mk_186 = buffer.data(mk + 186);
    const auto *mk_188 = buffer.data(mk + 188);
    const auto *mk_189 = buffer.data(mk + 189);
    const auto *mk_190 = buffer.data(mk + 190);
    const auto *mk_192 = buffer.data(mk + 192);
    const auto *mk_193 = buffer.data(mk + 193);
    const auto *mk_194 = buffer.data(mk + 194);
    const auto *mk_195 = buffer.data(mk + 195);
    const auto *mk_197 = buffer.data(mk + 197);
    const auto *mk_198 = buffer.data(mk + 198);
    const auto *mk_199 = buffer.data(mk + 199);
    const auto *mk_200 = buffer.data(mk + 200);
    const auto *mk_207 = buffer.data(mk + 207);
    const auto *mk_208 = buffer.data(mk + 208);
    const auto *mk_209 = buffer.data(mk + 209);
    const auto *mk_210 = buffer.data(mk + 210);
    const auto *mk_211 = buffer.data(mk + 211);
    const auto *mk_212 = buffer.data(mk + 212);
    const auto *mk_213 = buffer.data(mk + 213);
    const auto *mk_214 = buffer.data(mk + 214);
    const auto *mk_215 = buffer.data(mk + 215);
    const auto *mk_216 = buffer.data(mk + 216);
    const auto *mk_217 = buffer.data(mk + 217);
    const auto *mk_218 = buffer.data(mk + 218);
    const auto *mk_219 = buffer.data(mk + 219);

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_z, lk_118, mi0_87, mi0_94, mi1_87, \
                         mi1_94, mk_114, mk_115, mk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_22 * lk_118[k]
                   + f_7 * mi0_94[k]
                   - f_8 * mi1_94[k]
                   + pb_x[k] * mk_118[k];

        t_146[k] = pb_z[k] * mk_114[k];

        t_147[k] = f_3 * mi0_87[k]
                   - f_4 * mi1_87[k]
                   + pb_z[k] * mk_115[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pb_y, pb_z, lk_45, lk_123, mi0_89, \
                         mi0_99, mi1_89, mi1_99, mk_117, mk_118, \
                         mk_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_14 * lk_45[k]
                   + pb_y[k] * mk_117[k];

        t_149[k] = f_7 * mi0_89[k]
                   - f_8 * mi1_89[k]
                   + pb_z[k] * mk_117[k];

        t_150[k] = f_22 * lk_123[k]
                   + f_5 * mi0_99[k]
                   - f_6 * mi1_99[k]
                   + pb_x[k] * mk_123[k];

        t_151[k] = pb_z[k] * mk_118[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, lk_50, mi0_90, mi0_91, \
                         mi0_93, mi1_90, mi1_91, mi1_93, mk_119, mk_120, \
                         mk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * mi0_90[k]
                   - f_4 * mi1_90[k]
                   + pb_z[k] * mk_119[k];

        t_153[k] = f_5 * mi0_91[k]
                   - f_6 * mi1_91[k]
                   + pb_z[k] * mk_120[k];

        t_154[k] = f_14 * lk_50[k]
                   + pb_y[k] * mk_122[k];

        t_155[k] = f_9 * mi0_93[k]
                   - f_10 * mi1_93[k]
                   + pb_z[k] * mk_122[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_z, lk_129, mi0_94, mi0_105, mi1_94, \
                         mi1_105, mk_123, mk_124, mk_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_22 * lk_129[k]
                   + f_3 * mi0_105[k]
                   - f_4 * mi1_105[k]
                   + pb_x[k] * mk_129[k];

        t_157[k] = pb_z[k] * mk_123[k];

        t_158[k] = f_3 * mi0_94[k]
                   - f_4 * mi1_94[k]
                   + pb_z[k] * mk_124[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_y, pb_z, lk_56, mi0_95, mi0_96, \
                         mi0_98, mi1_95, mi1_96, mi1_98, mk_125, mk_126, \
                         mk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * mi0_95[k]
                   - f_6 * mi1_95[k]
                   + pb_z[k] * mk_125[k];

        t_160[k] = f_7 * mi0_96[k]
                   - f_8 * mi1_96[k]
                   + pb_z[k] * mk_126[k];

        t_161[k] = f_14 * lk_56[k]
                   + pb_y[k] * mk_128[k];

        t_162[k] = f_11 * mi0_98[k]
                   - f_12 * mi1_98[k]
                   + pb_z[k] * mk_128[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pb_x, pb_z, lk_136, lk_138, \
                         lk_139, lk_140, mk_129, mk_136, mk_138, mk_139, \
                         mk_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_22 * lk_136[k]
                   + pb_x[k] * mk_136[k];

        t_164[k] = pb_z[k] * mk_129[k];

        t_165[k] = f_22 * lk_138[k]
                   + pb_x[k] * mk_138[k];

        t_166[k] = f_22 * lk_139[k]
                   + pb_x[k] * mk_139[k];

        t_167[k] = f_22 * lk_140[k]
                   + pb_x[k] * mk_140[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_x, pb_x, kl0_171, kl1_171, lk_141, \
                         lk_142, lk_143, ll_171, mk_141, mk_142, \
                         mk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_22 * lk_141[k]
                   + pb_x[k] * mk_141[k];

        t_169[k] = f_22 * lk_142[k]
                   + pb_x[k] * mk_142[k];

        t_170[k] = f_22 * lk_143[k]
                   + pb_x[k] * mk_143[k];

        t_171[k] = f_23 * kl0_171[k]
                   - f_24 * kl1_171[k]
                   + pa_x[k] * ll_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, mi0_105, mi0_106, mi0_107, mi1_105, \
                         mi1_106, mi1_107, mk_136, mk_137, mk_138, \
                         mk_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * mk_136[k];

        t_173[k] = f_3 * mi0_105[k]
                   - f_4 * mi1_105[k]
                   + pb_z[k] * mk_137[k];

        t_174[k] = f_5 * mi0_106[k]
                   - f_6 * mi1_106[k]
                   + pb_z[k] * mk_138[k];

        t_175[k] = f_7 * mi0_107[k]
                   - f_8 * mi1_107[k]
                   + pb_z[k] * mk_139[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_y, pb_z, lk_71, mi0_108, mi0_109, \
                         mi0_111, mi1_108, mi1_109, mi1_111, mk_140, mk_141, \
                         mk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_9 * mi0_108[k]
                   - f_10 * mi1_108[k]
                   + pb_z[k] * mk_140[k];

        t_177[k] = f_11 * mi0_109[k]
                   - f_12 * mi1_109[k]
                   + pb_z[k] * mk_141[k];

        t_178[k] = f_14 * lk_71[k]
                   + pb_y[k] * mk_143[k];

        t_179[k] = f_1 * mi0_111[k]
                   - f_2 * mi1_111[k]
                   + pb_z[k] * mk_143[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, t_185, pa_y, pa_z, pb_y, lk_74, \
                         ll_46, ll_48, ll_90, ll_92, ll_95, mk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * ll_90[k];

        t_181[k] = pa_z[k] * ll_46[k];

        t_182[k] = pa_y[k] * ll_92[k];

        t_183[k] = pa_z[k] * ll_48[k];

        t_184[k] = f_13 * lk_74[k]
                   + pb_y[k] * mk_146[k];

        t_185[k] = pa_y[k] * ll_95[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, t_190, pa_y, pa_z, pb_y, pb_z, lk_39, \
                         lk_77, ll_51, ll_55, ll_99, mk_147, mk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_z[k] * ll_51[k];

        t_187[k] = f_13 * lk_39[k]
                   + pb_z[k] * mk_147[k];

        t_188[k] = f_13 * lk_77[k]
                   + pb_y[k] * mk_149[k];

        t_189[k] = pa_y[k] * ll_99[k];

        t_190[k] = pa_z[k] * ll_55[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, lk_42, lk_80, lk_81, \
                         ll_102, ll_104, mk_150, mk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_13 * lk_42[k]
                   + pb_z[k] * mk_150[k];

        t_192[k] = f_14 * lk_80[k]
                   + pa_y[k] * ll_102[k];

        t_193[k] = f_13 * lk_81[k]
                   + pb_y[k] * mk_153[k];

        t_194[k] = pa_y[k] * ll_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_y, pa_z, pb_z, lk_46, lk_84, lk_85, \
                         ll_60, ll_107, ll_108, mk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_z[k] * ll_60[k];

        t_196[k] = f_13 * lk_46[k]
                   + pb_z[k] * mk_154[k];

        t_197[k] = f_15 * lk_84[k]
                   + pa_y[k] * ll_107[k];

        t_198[k] = f_14 * lk_85[k]
                   + pa_y[k] * ll_108[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, pb_z, lk_51, lk_86, \
                         ll_66, ll_110, mk_158, mk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_13 * lk_86[k]
                   + pb_y[k] * mk_158[k];

        t_200[k] = pa_y[k] * ll_110[k];

        t_201[k] = pa_z[k] * ll_66[k];

        t_202[k] = f_13 * lk_51[k]
                   + pb_z[k] * mk_159[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pa_y, pb_y, lk_89, lk_90, lk_91, \
                         lk_92, ll_113, ll_114, ll_115, ll_117, \
                         mk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_16 * lk_89[k]
                   + pa_y[k] * ll_113[k];

        t_204[k] = f_15 * lk_90[k]
                   + pa_y[k] * ll_114[k];

        t_205[k] = f_14 * lk_91[k]
                   + pa_y[k] * ll_115[k];

        t_206[k] = f_13 * lk_92[k]
                   + pb_y[k] * mk_164[k];

        t_207[k] = pa_y[k] * ll_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pa_z, pb_x, lk_173, lk_174, \
                         lk_175, lk_176, ll_73, mk_173, mk_174, mk_175, \
                         mk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_z[k] * ll_73[k];

        t_209[k] = f_22 * lk_173[k]
                   + pb_x[k] * mk_173[k];

        t_210[k] = f_22 * lk_174[k]
                   + pb_x[k] * mk_174[k];

        t_211[k] = f_22 * lk_175[k]
                   + pb_x[k] * mk_175[k];

        t_212[k] = f_22 * lk_176[k]
                   + pb_x[k] * mk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pa_z, pb_x, lk_177, lk_178, ll_81, \
                         ll_125, mk_177, mk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_22 * lk_177[k]
                   + pb_x[k] * mk_177[k];

        t_214[k] = f_22 * lk_178[k]
                   + pb_x[k] * mk_178[k];

        t_215[k] = pa_y[k] * ll_125[k];

        t_216[k] = pa_z[k] * ll_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_y, pb_z, lk_64, lk_102, lk_103, \
                         lk_104, ll_128, ll_129, ll_130, mk_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_13 * lk_64[k]
                   + pb_z[k] * mk_172[k];

        t_218[k] = f_18 * lk_102[k]
                   + pa_y[k] * ll_128[k];

        t_219[k] = f_17 * lk_103[k]
                   + pa_y[k] * ll_129[k];

        t_220[k] = f_16 * lk_104[k]
                   + pa_y[k] * ll_130[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pb_y, lk_105, lk_106, lk_107, \
                         ll_131, ll_132, ll_134, mk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_15 * lk_105[k]
                   + pa_y[k] * ll_131[k];

        t_222[k] = f_14 * lk_106[k]
                   + pa_y[k] * ll_132[k];

        t_223[k] = f_13 * lk_107[k]
                   + pb_y[k] * mk_179[k];

        t_224[k] = pa_y[k] * ll_134[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_z, pb_y, pb_z, kl0_0, kl1_0, lk_72, \
                         ll_90, mi0_140, mi1_140, mk_180, mk_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_20 * kl0_0[k]
                   - f_21 * kl1_0[k]
                   + pa_z[k] * ll_90[k];

        t_226[k] = pb_y[k] * mk_180[k];

        t_227[k] = f_14 * lk_72[k]
                   + pb_z[k] * mk_180[k];

        t_228[k] = f_3 * mi0_140[k]
                   - f_4 * mi1_140[k]
                   + pb_y[k] * mk_181[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, lk_75, lk_185, mi0_141, \
                         mi0_145, mi1_141, mi1_145, mk_182, mk_183, \
                         mk_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_y[k] * mk_182[k];

        t_230[k] = f_22 * lk_185[k]
                   + f_11 * mi0_145[k]
                   - f_12 * mi1_145[k]
                   + pb_x[k] * mk_185[k];

        t_231[k] = f_5 * mi0_141[k]
                   - f_6 * mi1_141[k]
                   + pb_y[k] * mk_183[k];

        t_232[k] = f_14 * lk_75[k]
                   + pb_z[k] * mk_183[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pb_y, pb_z, lk_78, lk_189, mi0_143, \
                         mi0_149, mi1_143, mi1_149, mk_185, mk_186, \
                         mk_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pb_y[k] * mk_185[k];

        t_234[k] = f_22 * lk_189[k]
                   + f_9 * mi0_149[k]
                   - f_10 * mi1_149[k]
                   + pb_x[k] * mk_189[k];

        t_235[k] = f_7 * mi0_143[k]
                   - f_8 * mi1_143[k]
                   + pb_y[k] * mk_186[k];

        t_236[k] = f_14 * lk_78[k]
                   + pb_z[k] * mk_186[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_x, pb_y, lk_194, mi0_145, mi0_154, mi1_145, \
                         mi1_154, mk_188, mk_189, mk_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_3 * mi0_145[k]
                   - f_4 * mi1_145[k]
                   + pb_y[k] * mk_188[k];

        t_238[k] = pb_y[k] * mk_189[k];

        t_239[k] = f_22 * lk_194[k]
                   + f_7 * mi0_154[k]
                   - f_8 * mi1_154[k]
                   + pb_x[k] * mk_194[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_y, pb_z, lk_82, mi0_146, mi0_148, \
                         mi0_149, mi1_146, mi1_148, mi1_149, mk_190, mk_192, \
                         mk_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_9 * mi0_146[k]
                   - f_10 * mi1_146[k]
                   + pb_y[k] * mk_190[k];

        t_241[k] = f_14 * lk_82[k]
                   + pb_z[k] * mk_190[k];

        t_242[k] = f_5 * mi0_148[k]
                   - f_6 * mi1_148[k]
                   + pb_y[k] * mk_192[k];

        t_243[k] = f_3 * mi0_149[k]
                   - f_4 * mi1_149[k]
                   + pb_y[k] * mk_193[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pb_y, pb_z, lk_87, lk_200, mi0_150, \
                         mi0_160, mi1_150, mi1_160, mk_194, mk_195, \
                         mk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pb_y[k] * mk_194[k];

        t_245[k] = f_22 * lk_200[k]
                   + f_5 * mi0_160[k]
                   - f_6 * mi1_160[k]
                   + pb_x[k] * mk_200[k];

        t_246[k] = f_11 * mi0_150[k]
                   - f_12 * mi1_150[k]
                   + pb_y[k] * mk_195[k];

        t_247[k] = f_14 * lk_87[k]
                   + pb_z[k] * mk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_y, mi0_152, mi0_153, mi0_154, mi1_152, \
                         mi1_153, mi1_154, mk_197, mk_198, mk_199, \
                         mk_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * mi0_152[k]
                   - f_8 * mi1_152[k]
                   + pb_y[k] * mk_197[k];

        t_249[k] = f_5 * mi0_153[k]
                   - f_6 * mi1_153[k]
                   + pb_y[k] * mk_198[k];

        t_250[k] = f_3 * mi0_154[k]
                   - f_4 * mi1_154[k]
                   + pb_y[k] * mk_199[k];

        t_251[k] = pb_y[k] * mk_200[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pb_x, lk_207, lk_208, lk_209, lk_210, \
                         mi0_167, mi1_167, mk_207, mk_208, mk_209, \
                         mk_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_22 * lk_207[k]
                   + f_3 * mi0_167[k]
                   - f_4 * mi1_167[k]
                   + pb_x[k] * mk_207[k];

        t_253[k] = f_22 * lk_208[k]
                   + pb_x[k] * mk_208[k];

        t_254[k] = f_22 * lk_209[k]
                   + pb_x[k] * mk_209[k];

        t_255[k] = f_22 * lk_210[k]
                   + pb_x[k] * mk_210[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, lk_211, lk_212, \
                         lk_213, lk_215, mk_207, mk_211, mk_212, mk_213, \
                         mk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_22 * lk_211[k]
                   + pb_x[k] * mk_211[k];

        t_257[k] = f_22 * lk_212[k]
                   + pb_x[k] * mk_212[k];

        t_258[k] = f_22 * lk_213[k]
                   + pb_x[k] * mk_213[k];

        t_259[k] = pb_y[k] * mk_207[k];

        t_260[k] = f_22 * lk_215[k]
                   + pb_x[k] * mk_215[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_y, pb_z, lk_100, mi0_161, mi0_163, \
                         mi0_164, mi1_161, mi1_163, mi1_164, mk_208, mk_210, \
                         mk_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_1 * mi0_161[k]
                   - f_2 * mi1_161[k]
                   + pb_y[k] * mk_208[k];

        t_262[k] = f_14 * lk_100[k]
                   + pb_z[k] * mk_208[k];

        t_263[k] = f_11 * mi0_163[k]
                   - f_12 * mi1_163[k]
                   + pb_y[k] * mk_210[k];

        t_264[k] = f_9 * mi0_164[k]
                   - f_10 * mi1_164[k]
                   + pb_y[k] * mk_211[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_y, mi0_165, mi0_166, mi0_167, mi1_165, \
                         mi1_166, mi1_167, mk_212, mk_213, mk_214, \
                         mk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_7 * mi0_165[k]
                   - f_8 * mi1_165[k]
                   + pb_y[k] * mk_212[k];

        t_266[k] = f_5 * mi0_166[k]
                   - f_6 * mi1_166[k]
                   + pb_y[k] * mk_213[k];

        t_267[k] = f_3 * mi0_167[k]
                   - f_4 * mi1_167[k]
                   + pb_y[k] * mk_214[k];

        t_268[k] = pb_y[k] * mk_215[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pa_y, pb_y, pb_z, kl0_45, kl0_269, \
                         kl1_45, kl1_269, lk_108, ll_135, ll_269, \
                         mk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_23 * kl0_269[k]
                   - f_24 * kl1_269[k]
                   + pa_x[k] * ll_269[k];

        t_270[k] = f_25 * kl0_45[k]
                   - f_26 * kl1_45[k]
                   + pa_y[k] * ll_135[k];

        t_271[k] = f_15 * lk_108[k]
                   + pb_y[k] * mk_216[k];

        t_272[k] = pb_z[k] * mk_216[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pb_x, pb_z, lk_219, mi0_168, mi0_171, mi1_168, \
                         mi1_171, mk_217, mk_218, mk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_18 * lk_219[k]
                   + f_11 * mi0_171[k]
                   - f_12 * mi1_171[k]
                   + pb_x[k] * mk_219[k];

        t_274[k] = pb_z[k] * mk_217[k];

        t_275[k] = f_3 * mi0_168[k]
                   - f_4 * mi1_168[k]
                   + pb_z[k] * mk_218[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_27 = 2.5 / alpha;
    const auto f_28 = 2.5 * beta / (alpha * p);

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

    const auto *kl0_90 = buffer.data(kl0 + 90);
    const auto *kl0_306 = buffer.data(kl0 + 306);

    const auto *kl1_90 = buffer.data(kl1 + 90);
    const auto *kl1_306 = buffer.data(kl1 + 306);

    const auto *lk_108 = buffer.data(lk + 108);
    const auto *lk_110 = buffer.data(lk + 110);
    const auto *lk_111 = buffer.data(lk + 111);
    const auto *lk_113 = buffer.data(lk + 113);
    const auto *lk_114 = buffer.data(lk + 114);
    const auto *lk_115 = buffer.data(lk + 115);
    const auto *lk_117 = buffer.data(lk + 117);
    const auto *lk_118 = buffer.data(lk + 118);
    const auto *lk_119 = buffer.data(lk + 119);
    const auto *lk_120 = buffer.data(lk + 120);
    const auto *lk_122 = buffer.data(lk + 122);
    const auto *lk_123 = buffer.data(lk + 123);
    const auto *lk_124 = buffer.data(lk + 124);
    const auto *lk_125 = buffer.data(lk + 125);
    const auto *lk_126 = buffer.data(lk + 126);
    const auto *lk_128 = buffer.data(lk + 128);
    const auto *lk_136 = buffer.data(lk + 136);
    const auto *lk_137 = buffer.data(lk + 137);
    const auto *lk_138 = buffer.data(lk + 138);
    const auto *lk_139 = buffer.data(lk + 139);
    const auto *lk_140 = buffer.data(lk + 140);
    const auto *lk_141 = buffer.data(lk + 141);
    const auto *lk_143 = buffer.data(lk + 143);
    const auto *lk_146 = buffer.data(lk + 146);
    const auto *lk_147 = buffer.data(lk + 147);
    const auto *lk_149 = buffer.data(lk + 149);
    const auto *lk_150 = buffer.data(lk + 150);
    const auto *lk_153 = buffer.data(lk + 153);
    const auto *lk_154 = buffer.data(lk + 154);
    const auto *lk_158 = buffer.data(lk + 158);
    const auto *lk_159 = buffer.data(lk + 159);
    const auto *lk_164 = buffer.data(lk + 164);
    const auto *lk_172 = buffer.data(lk + 172);
    const auto *lk_179 = buffer.data(lk + 179);
    const auto *lk_180 = buffer.data(lk + 180);
    const auto *lk_181 = buffer.data(lk + 181);
    const auto *lk_182 = buffer.data(lk + 182);
    const auto *lk_183 = buffer.data(lk + 183);
    const auto *lk_185 = buffer.data(lk + 185);
    const auto *lk_186 = buffer.data(lk + 186);
    const auto *lk_188 = buffer.data(lk + 188);
    const auto *lk_189 = buffer.data(lk + 189);
    const auto *lk_190 = buffer.data(lk + 190);
    const auto *lk_192 = buffer.data(lk + 192);
    const auto *lk_193 = buffer.data(lk + 193);
    const auto *lk_194 = buffer.data(lk + 194);
    const auto *lk_195 = buffer.data(lk + 195);
    const auto *lk_197 = buffer.data(lk + 197);
    const auto *lk_198 = buffer.data(lk + 198);
    const auto *lk_199 = buffer.data(lk + 199);
    const auto *lk_200 = buffer.data(lk + 200);
    const auto *lk_208 = buffer.data(lk + 208);
    const auto *lk_210 = buffer.data(lk + 210);
    const auto *lk_211 = buffer.data(lk + 211);
    const auto *lk_212 = buffer.data(lk + 212);
    const auto *lk_213 = buffer.data(lk + 213);
    const auto *lk_214 = buffer.data(lk + 214);
    const auto *lk_215 = buffer.data(lk + 215);
    const auto *lk_222 = buffer.data(lk + 222);
    const auto *lk_226 = buffer.data(lk + 226);
    const auto *lk_231 = buffer.data(lk + 231);
    const auto *lk_237 = buffer.data(lk + 237);
    const auto *lk_244 = buffer.data(lk + 244);
    const auto *lk_246 = buffer.data(lk + 246);
    const auto *lk_247 = buffer.data(lk + 247);
    const auto *lk_248 = buffer.data(lk + 248);
    const auto *lk_249 = buffer.data(lk + 249);
    const auto *lk_250 = buffer.data(lk + 250);
    const auto *lk_251 = buffer.data(lk + 251);
    const auto *lk_281 = buffer.data(lk + 281);
    const auto *lk_282 = buffer.data(lk + 282);
    const auto *lk_283 = buffer.data(lk + 283);
    const auto *lk_284 = buffer.data(lk + 284);
    const auto *lk_285 = buffer.data(lk + 285);
    const auto *lk_286 = buffer.data(lk + 286);
    const auto *lk_287 = buffer.data(lk + 287);
    const auto *lk_316 = buffer.data(lk + 316);
    const auto *lk_317 = buffer.data(lk + 317);
    const auto *lk_318 = buffer.data(lk + 318);
    const auto *lk_319 = buffer.data(lk + 319);
    const auto *lk_320 = buffer.data(lk + 320);
    const auto *lk_321 = buffer.data(lk + 321);
    const auto *lk_322 = buffer.data(lk + 322);
    const auto *lk_329 = buffer.data(lk + 329);

    const auto *ll_135 = buffer.data(ll + 135);
    const auto *ll_136 = buffer.data(ll + 136);
    const auto *ll_138 = buffer.data(ll + 138);
    const auto *ll_140 = buffer.data(ll + 140);
    const auto *ll_141 = buffer.data(ll + 141);
    const auto *ll_144 = buffer.data(ll + 144);
    const auto *ll_145 = buffer.data(ll + 145);
    const auto *ll_147 = buffer.data(ll + 147);
    const auto *ll_149 = buffer.data(ll + 149);
    const auto *ll_150 = buffer.data(ll + 150);
    const auto *ll_152 = buffer.data(ll + 152);
    const auto *ll_153 = buffer.data(ll + 153);
    const auto *ll_155 = buffer.data(ll + 155);
    const auto *ll_156 = buffer.data(ll + 156);
    const auto *ll_158 = buffer.data(ll + 158);
    const auto *ll_159 = buffer.data(ll + 159);
    const auto *ll_160 = buffer.data(ll + 160);
    const auto *ll_162 = buffer.data(ll + 162);
    const auto *ll_163 = buffer.data(ll + 163);
    const auto *ll_171 = buffer.data(ll + 171);
    const auto *ll_173 = buffer.data(ll + 173);
    const auto *ll_174 = buffer.data(ll + 174);
    const auto *ll_175 = buffer.data(ll + 175);
    const auto *ll_176 = buffer.data(ll + 176);
    const auto *ll_177 = buffer.data(ll + 177);
    const auto *ll_179 = buffer.data(ll + 179);
    const auto *ll_225 = buffer.data(ll + 225);
    const auto *ll_227 = buffer.data(ll + 227);
    const auto *ll_228 = buffer.data(ll + 228);
    const auto *ll_230 = buffer.data(ll + 230);
    const auto *ll_231 = buffer.data(ll + 231);
    const auto *ll_234 = buffer.data(ll + 234);
    const auto *ll_235 = buffer.data(ll + 235);
    const auto *ll_237 = buffer.data(ll + 237);
    const auto *ll_239 = buffer.data(ll + 239);
    const auto *ll_240 = buffer.data(ll + 240);
    const auto *ll_242 = buffer.data(ll + 242);
    const auto *ll_243 = buffer.data(ll + 243);
    const auto *ll_245 = buffer.data(ll + 245);
    const auto *ll_246 = buffer.data(ll + 246);
    const auto *ll_248 = buffer.data(ll + 248);
    const auto *ll_249 = buffer.data(ll + 249);
    const auto *ll_250 = buffer.data(ll + 250);
    const auto *ll_252 = buffer.data(ll + 252);
    const auto *ll_260 = buffer.data(ll + 260);
    const auto *ll_261 = buffer.data(ll + 261);
    const auto *ll_263 = buffer.data(ll + 263);
    const auto *ll_264 = buffer.data(ll + 264);
    const auto *ll_265 = buffer.data(ll + 265);
    const auto *ll_266 = buffer.data(ll + 266);
    const auto *ll_267 = buffer.data(ll + 267);
    const auto *ll_269 = buffer.data(ll + 269);
    const auto *ll_306 = buffer.data(ll + 306);

    const auto *mi0_170 = buffer.data(mi0 + 170);
    const auto *mi0_171 = buffer.data(mi0 + 171);
    const auto *mi0_173 = buffer.data(mi0 + 173);
    const auto *mi0_174 = buffer.data(mi0 + 174);
    const auto *mi0_175 = buffer.data(mi0 + 175);
    const auto *mi0_177 = buffer.data(mi0 + 177);
    const auto *mi0_178 = buffer.data(mi0 + 178);
    const auto *mi0_179 = buffer.data(mi0 + 179);
    const auto *mi0_180 = buffer.data(mi0 + 180);
    const auto *mi0_182 = buffer.data(mi0 + 182);
    const auto *mi0_183 = buffer.data(mi0 + 183);
    const auto *mi0_189 = buffer.data(mi0 + 189);
    const auto *mi0_190 = buffer.data(mi0 + 190);
    const auto *mi0_191 = buffer.data(mi0 + 191);
    const auto *mi0_192 = buffer.data(mi0 + 192);
    const auto *mi0_193 = buffer.data(mi0 + 193);
    const auto *mi0_195 = buffer.data(mi0 + 195);
    const auto *mi0_252 = buffer.data(mi0 + 252);
    const auto *mi0_253 = buffer.data(mi0 + 253);
    const auto *mi0_257 = buffer.data(mi0 + 257);

    const auto *mi1_170 = buffer.data(mi1 + 170);
    const auto *mi1_171 = buffer.data(mi1 + 171);
    const auto *mi1_173 = buffer.data(mi1 + 173);
    const auto *mi1_174 = buffer.data(mi1 + 174);
    const auto *mi1_175 = buffer.data(mi1 + 175);
    const auto *mi1_177 = buffer.data(mi1 + 177);
    const auto *mi1_178 = buffer.data(mi1 + 178);
    const auto *mi1_179 = buffer.data(mi1 + 179);
    const auto *mi1_180 = buffer.data(mi1 + 180);
    const auto *mi1_182 = buffer.data(mi1 + 182);
    const auto *mi1_183 = buffer.data(mi1 + 183);
    const auto *mi1_189 = buffer.data(mi1 + 189);
    const auto *mi1_190 = buffer.data(mi1 + 190);
    const auto *mi1_191 = buffer.data(mi1 + 191);
    const auto *mi1_192 = buffer.data(mi1 + 192);
    const auto *mi1_193 = buffer.data(mi1 + 193);
    const auto *mi1_195 = buffer.data(mi1 + 195);
    const auto *mi1_252 = buffer.data(mi1 + 252);
    const auto *mi1_253 = buffer.data(mi1 + 253);
    const auto *mi1_257 = buffer.data(mi1 + 257);

    const auto *mk_219 = buffer.data(mk + 219);
    const auto *mk_221 = buffer.data(mk + 221);
    const auto *mk_222 = buffer.data(mk + 222);
    const auto *mk_223 = buffer.data(mk + 223);
    const auto *mk_225 = buffer.data(mk + 225);
    const auto *mk_226 = buffer.data(mk + 226);
    const auto *mk_227 = buffer.data(mk + 227);
    const auto *mk_228 = buffer.data(mk + 228);
    const auto *mk_230 = buffer.data(mk + 230);
    const auto *mk_231 = buffer.data(mk + 231);
    const auto *mk_232 = buffer.data(mk + 232);
    const auto *mk_233 = buffer.data(mk + 233);
    const auto *mk_234 = buffer.data(mk + 234);
    const auto *mk_236 = buffer.data(mk + 236);
    const auto *mk_237 = buffer.data(mk + 237);
    const auto *mk_244 = buffer.data(mk + 244);
    const auto *mk_245 = buffer.data(mk + 245);
    const auto *mk_246 = buffer.data(mk + 246);
    const auto *mk_247 = buffer.data(mk + 247);
    const auto *mk_248 = buffer.data(mk + 248);
    const auto *mk_249 = buffer.data(mk + 249);
    const auto *mk_250 = buffer.data(mk + 250);
    const auto *mk_251 = buffer.data(mk + 251);
    const auto *mk_252 = buffer.data(mk + 252);
    const auto *mk_254 = buffer.data(mk + 254);
    const auto *mk_255 = buffer.data(mk + 255);
    const auto *mk_257 = buffer.data(mk + 257);
    const auto *mk_258 = buffer.data(mk + 258);
    const auto *mk_261 = buffer.data(mk + 261);
    const auto *mk_262 = buffer.data(mk + 262);
    const auto *mk_266 = buffer.data(mk + 266);
    const auto *mk_267 = buffer.data(mk + 267);
    const auto *mk_272 = buffer.data(mk + 272);
    const auto *mk_280 = buffer.data(mk + 280);
    const auto *mk_281 = buffer.data(mk + 281);
    const auto *mk_282 = buffer.data(mk + 282);
    const auto *mk_283 = buffer.data(mk + 283);
    const auto *mk_284 = buffer.data(mk + 284);
    const auto *mk_285 = buffer.data(mk + 285);
    const auto *mk_286 = buffer.data(mk + 286);
    const auto *mk_287 = buffer.data(mk + 287);
    const auto *mk_288 = buffer.data(mk + 288);
    const auto *mk_290 = buffer.data(mk + 290);
    const auto *mk_291 = buffer.data(mk + 291);
    const auto *mk_293 = buffer.data(mk + 293);
    const auto *mk_294 = buffer.data(mk + 294);
    const auto *mk_297 = buffer.data(mk + 297);
    const auto *mk_298 = buffer.data(mk + 298);
    const auto *mk_302 = buffer.data(mk + 302);
    const auto *mk_303 = buffer.data(mk + 303);
    const auto *mk_308 = buffer.data(mk + 308);
    const auto *mk_316 = buffer.data(mk + 316);
    const auto *mk_317 = buffer.data(mk + 317);
    const auto *mk_318 = buffer.data(mk + 318);
    const auto *mk_319 = buffer.data(mk + 319);
    const auto *mk_320 = buffer.data(mk + 320);
    const auto *mk_321 = buffer.data(mk + 321);
    const auto *mk_322 = buffer.data(mk + 322);
    const auto *mk_323 = buffer.data(mk + 323);
    const auto *mk_324 = buffer.data(mk + 324);
    const auto *mk_325 = buffer.data(mk + 325);
    const auto *mk_326 = buffer.data(mk + 326);
    const auto *mk_327 = buffer.data(mk + 327);
    const auto *mk_329 = buffer.data(mk + 329);

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, pb_y, pb_z, lk_113, lk_222, \
                         mi0_170, mi0_174, mi1_170, mi1_174, mk_219, mk_221, \
                         mk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_18 * lk_222[k]
                   + f_9 * mi0_174[k]
                   - f_10 * mi1_174[k]
                   + pb_x[k] * mk_222[k];

        t_277[k] = pb_z[k] * mk_219[k];

        t_278[k] = f_15 * lk_113[k]
                   + pb_y[k] * mk_221[k];

        t_279[k] = f_5 * mi0_170[k]
                   - f_6 * mi1_170[k]
                   + pb_z[k] * mk_221[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pb_x, pb_z, lk_226, mi0_171, mi0_178, mi1_171, \
                         mi1_178, mk_222, mk_223, mk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_18 * lk_226[k]
                   + f_7 * mi0_178[k]
                   - f_8 * mi1_178[k]
                   + pb_x[k] * mk_226[k];

        t_281[k] = pb_z[k] * mk_222[k];

        t_282[k] = f_3 * mi0_171[k]
                   - f_4 * mi1_171[k]
                   + pb_z[k] * mk_223[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_x, pb_y, pb_z, lk_117, lk_231, \
                         mi0_173, mi0_183, mi1_173, mi1_183, mk_225, mk_226, \
                         mk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_15 * lk_117[k]
                   + pb_y[k] * mk_225[k];

        t_284[k] = f_7 * mi0_173[k]
                   - f_8 * mi1_173[k]
                   + pb_z[k] * mk_225[k];

        t_285[k] = f_18 * lk_231[k]
                   + f_5 * mi0_183[k]
                   - f_6 * mi1_183[k]
                   + pb_x[k] * mk_231[k];

        t_286[k] = pb_z[k] * mk_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_y, pb_z, lk_122, mi0_174, mi0_175, \
                         mi0_177, mi1_174, mi1_175, mi1_177, mk_227, mk_228, \
                         mk_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * mi0_174[k]
                   - f_4 * mi1_174[k]
                   + pb_z[k] * mk_227[k];

        t_288[k] = f_5 * mi0_175[k]
                   - f_6 * mi1_175[k]
                   + pb_z[k] * mk_228[k];

        t_289[k] = f_15 * lk_122[k]
                   + pb_y[k] * mk_230[k];

        t_290[k] = f_9 * mi0_177[k]
                   - f_10 * mi1_177[k]
                   + pb_z[k] * mk_230[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pb_x, pb_z, lk_237, mi0_178, mi0_189, mi1_178, \
                         mi1_189, mk_231, mk_232, mk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_18 * lk_237[k]
                   + f_3 * mi0_189[k]
                   - f_4 * mi1_189[k]
                   + pb_x[k] * mk_237[k];

        t_292[k] = pb_z[k] * mk_231[k];

        t_293[k] = f_3 * mi0_178[k]
                   - f_4 * mi1_178[k]
                   + pb_z[k] * mk_232[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_y, pb_z, lk_128, mi0_179, mi0_180, \
                         mi0_182, mi1_179, mi1_180, mi1_182, mk_233, mk_234, \
                         mk_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_5 * mi0_179[k]
                   - f_6 * mi1_179[k]
                   + pb_z[k] * mk_233[k];

        t_295[k] = f_7 * mi0_180[k]
                   - f_8 * mi1_180[k]
                   + pb_z[k] * mk_234[k];

        t_296[k] = f_15 * lk_128[k]
                   + pb_y[k] * mk_236[k];

        t_297[k] = f_11 * mi0_182[k]
                   - f_12 * mi1_182[k]
                   + pb_z[k] * mk_236[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pb_x, pb_z, lk_244, lk_246, \
                         lk_247, lk_248, mk_237, mk_244, mk_246, mk_247, \
                         mk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_18 * lk_244[k]
                   + pb_x[k] * mk_244[k];

        t_299[k] = pb_z[k] * mk_237[k];

        t_300[k] = f_18 * lk_246[k]
                   + pb_x[k] * mk_246[k];

        t_301[k] = f_18 * lk_247[k]
                   + pb_x[k] * mk_247[k];

        t_302[k] = f_18 * lk_248[k]
                   + pb_x[k] * mk_248[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pb_x, kl0_306, kl1_306, lk_249, \
                         lk_250, lk_251, ll_306, mk_249, mk_250, \
                         mk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_18 * lk_249[k]
                   + pb_x[k] * mk_249[k];

        t_304[k] = f_18 * lk_250[k]
                   + pb_x[k] * mk_250[k];

        t_305[k] = f_18 * lk_251[k]
                   + pb_x[k] * mk_251[k];

        t_306[k] = f_27 * kl0_306[k]
                   - f_28 * kl1_306[k]
                   + pa_x[k] * ll_306[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pb_z, mi0_189, mi0_190, mi0_191, mi1_189, \
                         mi1_190, mi1_191, mk_244, mk_245, mk_246, \
                         mk_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = pb_z[k] * mk_244[k];

        t_308[k] = f_3 * mi0_189[k]
                   - f_4 * mi1_189[k]
                   + pb_z[k] * mk_245[k];

        t_309[k] = f_5 * mi0_190[k]
                   - f_6 * mi1_190[k]
                   + pb_z[k] * mk_246[k];

        t_310[k] = f_7 * mi0_191[k]
                   - f_8 * mi1_191[k]
                   + pb_z[k] * mk_247[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pb_y, pb_z, lk_143, mi0_192, mi0_193, \
                         mi0_195, mi1_192, mi1_193, mi1_195, mk_248, mk_249, \
                         mk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_9 * mi0_192[k]
                   - f_10 * mi1_192[k]
                   + pb_z[k] * mk_248[k];

        t_312[k] = f_11 * mi0_193[k]
                   - f_12 * mi1_193[k]
                   + pb_z[k] * mk_249[k];

        t_313[k] = f_15 * lk_143[k]
                   + pb_y[k] * mk_251[k];

        t_314[k] = f_1 * mi0_195[k]
                   - f_2 * mi1_195[k]
                   + pb_z[k] * mk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, pa_z, pb_y, pb_z, lk_108, lk_146, \
                         ll_135, ll_136, ll_138, mk_252, mk_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_z[k] * ll_135[k];

        t_316[k] = pa_z[k] * ll_136[k];

        t_317[k] = f_13 * lk_108[k]
                   + pb_z[k] * mk_252[k];

        t_318[k] = pa_z[k] * ll_138[k];

        t_319[k] = f_14 * lk_146[k]
                   + pb_y[k] * mk_254[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_z, pb_y, pb_z, lk_110, lk_111, lk_149, \
                         ll_140, ll_141, mk_255, mk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_14 * lk_110[k]
                   + pa_z[k] * ll_140[k];

        t_321[k] = pa_z[k] * ll_141[k];

        t_322[k] = f_13 * lk_111[k]
                   + pb_z[k] * mk_255[k];

        t_323[k] = f_14 * lk_149[k]
                   + pb_y[k] * mk_257[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_z, lk_113, lk_114, lk_115, \
                         ll_144, ll_145, ll_147, mk_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_15 * lk_113[k]
                   + pa_z[k] * ll_144[k];

        t_325[k] = pa_z[k] * ll_145[k];

        t_326[k] = f_13 * lk_114[k]
                   + pb_z[k] * mk_258[k];

        t_327[k] = f_14 * lk_115[k]
                   + pa_z[k] * ll_147[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_z, pb_y, pb_z, lk_117, lk_118, lk_153, \
                         ll_149, ll_150, mk_261, mk_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_14 * lk_153[k]
                   + pb_y[k] * mk_261[k];

        t_329[k] = f_16 * lk_117[k]
                   + pa_z[k] * ll_149[k];

        t_330[k] = pa_z[k] * ll_150[k];

        t_331[k] = f_13 * lk_118[k]
                   + pb_z[k] * mk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pa_z, pb_y, lk_119, lk_120, \
                         lk_122, lk_158, ll_152, ll_153, ll_155, ll_156, \
                         mk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_14 * lk_119[k]
                   + pa_z[k] * ll_152[k];

        t_333[k] = f_15 * lk_120[k]
                   + pa_z[k] * ll_153[k];

        t_334[k] = f_14 * lk_158[k]
                   + pb_y[k] * mk_266[k];

        t_335[k] = f_17 * lk_122[k]
                   + pa_z[k] * ll_155[k];

        t_336[k] = pa_z[k] * ll_156[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pb_z, lk_123, lk_124, lk_125, \
                         lk_126, ll_158, ll_159, ll_160, mk_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_13 * lk_123[k]
                   + pb_z[k] * mk_267[k];

        t_338[k] = f_14 * lk_124[k]
                   + pa_z[k] * ll_158[k];

        t_339[k] = f_15 * lk_125[k]
                   + pa_z[k] * ll_159[k];

        t_340[k] = f_16 * lk_126[k]
                   + pa_z[k] * ll_160[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_z, pb_x, pb_y, lk_128, lk_164, lk_281, \
                         ll_162, ll_163, mk_272, mk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_14 * lk_164[k]
                   + pb_y[k] * mk_272[k];

        t_342[k] = f_18 * lk_128[k]
                   + pa_z[k] * ll_162[k];

        t_343[k] = pa_z[k] * ll_163[k];

        t_344[k] = f_18 * lk_281[k]
                   + pb_x[k] * mk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pb_x, lk_282, lk_283, lk_284, \
                         lk_285, lk_286, mk_282, mk_283, mk_284, mk_285, \
                         mk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * lk_282[k]
                   + pb_x[k] * mk_282[k];

        t_346[k] = f_18 * lk_283[k]
                   + pb_x[k] * mk_283[k];

        t_347[k] = f_18 * lk_284[k]
                   + pb_x[k] * mk_284[k];

        t_348[k] = f_18 * lk_285[k]
                   + pb_x[k] * mk_285[k];

        t_349[k] = f_18 * lk_286[k]
                   + pb_x[k] * mk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_z, pb_x, pb_z, lk_136, lk_137, lk_287, \
                         ll_171, ll_173, mk_280, mk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_18 * lk_287[k]
                   + pb_x[k] * mk_287[k];

        t_351[k] = pa_z[k] * ll_171[k];

        t_352[k] = f_13 * lk_136[k]
                   + pb_z[k] * mk_280[k];

        t_353[k] = f_14 * lk_137[k]
                   + pa_z[k] * ll_173[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_z, lk_138, lk_139, lk_140, lk_141, \
                         ll_174, ll_175, ll_176, ll_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_15 * lk_138[k]
                   + pa_z[k] * ll_174[k];

        t_355[k] = f_16 * lk_139[k]
                   + pa_z[k] * ll_175[k];

        t_356[k] = f_17 * lk_140[k]
                   + pa_z[k] * ll_176[k];

        t_357[k] = f_18 * lk_141[k]
                   + pa_z[k] * ll_177[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, pa_y, pa_z, pb_y, lk_143, lk_179, \
                         lk_180, ll_179, ll_225, ll_227, mk_287, \
                         mk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_14 * lk_179[k]
                   + pb_y[k] * mk_287[k];

        t_359[k] = f_19 * lk_143[k]
                   + pa_z[k] * ll_179[k];

        t_360[k] = pa_y[k] * ll_225[k];

        t_361[k] = f_13 * lk_180[k]
                   + pb_y[k] * mk_288[k];

        t_362[k] = pa_y[k] * ll_227[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pa_y, pb_y, lk_181, lk_182, lk_183, \
                         ll_228, ll_230, ll_231, mk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_14 * lk_181[k]
                   + pa_y[k] * ll_228[k];

        t_364[k] = f_13 * lk_182[k]
                   + pb_y[k] * mk_290[k];

        t_365[k] = pa_y[k] * ll_230[k];

        t_366[k] = f_15 * lk_183[k]
                   + pa_y[k] * ll_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_y, pb_y, pb_z, lk_147, lk_185, lk_186, \
                         ll_234, ll_235, mk_291, mk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_14 * lk_147[k]
                   + pb_z[k] * mk_291[k];

        t_368[k] = f_13 * lk_185[k]
                   + pb_y[k] * mk_293[k];

        t_369[k] = pa_y[k] * ll_234[k];

        t_370[k] = f_16 * lk_186[k]
                   + pa_y[k] * ll_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pa_y, pb_y, pb_z, lk_150, lk_188, lk_189, \
                         ll_237, ll_239, mk_294, mk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_14 * lk_150[k]
                   + pb_z[k] * mk_294[k];

        t_372[k] = f_14 * lk_188[k]
                   + pa_y[k] * ll_237[k];

        t_373[k] = f_13 * lk_189[k]
                   + pb_y[k] * mk_297[k];

        t_374[k] = pa_y[k] * ll_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pa_y, pb_z, lk_154, lk_190, lk_192, \
                         lk_193, ll_240, ll_242, ll_243, mk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_17 * lk_190[k]
                   + pa_y[k] * ll_240[k];

        t_376[k] = f_14 * lk_154[k]
                   + pb_z[k] * mk_298[k];

        t_377[k] = f_15 * lk_192[k]
                   + pa_y[k] * ll_242[k];

        t_378[k] = f_14 * lk_193[k]
                   + pa_y[k] * ll_243[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pa_y, pb_y, pb_z, lk_159, lk_194, lk_195, \
                         ll_245, ll_246, mk_302, mk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_13 * lk_194[k]
                   + pb_y[k] * mk_302[k];

        t_380[k] = pa_y[k] * ll_245[k];

        t_381[k] = f_18 * lk_195[k]
                   + pa_y[k] * ll_246[k];

        t_382[k] = f_14 * lk_159[k]
                   + pb_z[k] * mk_303[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, pa_y, pb_y, lk_197, lk_198, \
                         lk_199, lk_200, ll_248, ll_249, ll_250, ll_252, \
                         mk_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_16 * lk_197[k]
                   + pa_y[k] * ll_248[k];

        t_384[k] = f_15 * lk_198[k]
                   + pa_y[k] * ll_249[k];

        t_385[k] = f_14 * lk_199[k]
                   + pa_y[k] * ll_250[k];

        t_386[k] = f_13 * lk_200[k]
                   + pb_y[k] * mk_308[k];

        t_387[k] = pa_y[k] * ll_252[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, pb_x, lk_316, lk_317, lk_318, \
                         lk_319, lk_320, mk_316, mk_317, mk_318, mk_319, \
                         mk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_18 * lk_316[k]
                   + pb_x[k] * mk_316[k];

        t_389[k] = f_18 * lk_317[k]
                   + pb_x[k] * mk_317[k];

        t_390[k] = f_18 * lk_318[k]
                   + pb_x[k] * mk_318[k];

        t_391[k] = f_18 * lk_319[k]
                   + pb_x[k] * mk_319[k];

        t_392[k] = f_18 * lk_320[k]
                   + pb_x[k] * mk_320[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_y, pb_x, lk_208, lk_321, lk_322, \
                         ll_260, ll_261, mk_321, mk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_18 * lk_321[k]
                   + pb_x[k] * mk_321[k];

        t_394[k] = f_18 * lk_322[k]
                   + pb_x[k] * mk_322[k];

        t_395[k] = pa_y[k] * ll_260[k];

        t_396[k] = f_19 * lk_208[k]
                   + pa_y[k] * ll_261[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pb_z, lk_172, lk_210, lk_211, \
                         lk_212, ll_263, ll_264, ll_265, mk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_14 * lk_172[k]
                   + pb_z[k] * mk_316[k];

        t_398[k] = f_18 * lk_210[k]
                   + pa_y[k] * ll_263[k];

        t_399[k] = f_17 * lk_211[k]
                   + pa_y[k] * ll_264[k];

        t_400[k] = f_16 * lk_212[k]
                   + pa_y[k] * ll_265[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pb_y, lk_213, lk_214, lk_215, \
                         ll_266, ll_267, ll_269, mk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * lk_213[k]
                   + pa_y[k] * ll_266[k];

        t_402[k] = f_14 * lk_214[k]
                   + pa_y[k] * ll_267[k];

        t_403[k] = f_13 * lk_215[k]
                   + pb_y[k] * mk_323[k];

        t_404[k] = pa_y[k] * ll_269[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_z, pb_y, pb_z, kl0_90, kl1_90, lk_180, \
                         ll_225, mi0_252, mi1_252, mk_324, mk_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_25 * kl0_90[k]
                   - f_26 * kl1_90[k]
                   + pa_z[k] * ll_225[k];

        t_406[k] = pb_y[k] * mk_324[k];

        t_407[k] = f_15 * lk_180[k]
                   + pb_z[k] * mk_324[k];

        t_408[k] = f_3 * mi0_252[k]
                   - f_4 * mi1_252[k]
                   + pb_y[k] * mk_325[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, pb_x, pb_y, pb_z, lk_183, lk_329, \
                         mi0_253, mi0_257, mi1_253, mi1_257, mk_326, mk_327, \
                         mk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = pb_y[k] * mk_326[k];

        t_410[k] = f_18 * lk_329[k]
                   + f_11 * mi0_257[k]
                   - f_12 * mi1_257[k]
                   + pb_x[k] * mk_329[k];

        t_411[k] = f_5 * mi0_253[k]
                   - f_6 * mi1_253[k]
                   + pb_y[k] * mk_327[k];

        t_412[k] = f_15 * lk_183[k]
                   + pb_z[k] * mk_327[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_27 = 2.5 / alpha;
    const auto f_28 = 2.5 * beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

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

    const auto *kl0_135 = buffer.data(kl0 + 135);
    const auto *kl0_225 = buffer.data(kl0 + 225);
    const auto *kl0_449 = buffer.data(kl0 + 449);
    const auto *kl0_486 = buffer.data(kl0 + 486);

    const auto *kl1_135 = buffer.data(kl1 + 135);
    const auto *kl1_225 = buffer.data(kl1 + 225);
    const auto *kl1_449 = buffer.data(kl1 + 449);
    const auto *kl1_486 = buffer.data(kl1 + 486);

    const auto *lk_186 = buffer.data(lk + 186);
    const auto *lk_190 = buffer.data(lk + 190);
    const auto *lk_195 = buffer.data(lk + 195);
    const auto *lk_208 = buffer.data(lk + 208);
    const auto *lk_216 = buffer.data(lk + 216);
    const auto *lk_218 = buffer.data(lk + 218);
    const auto *lk_219 = buffer.data(lk + 219);
    const auto *lk_221 = buffer.data(lk + 221);
    const auto *lk_222 = buffer.data(lk + 222);
    const auto *lk_223 = buffer.data(lk + 223);
    const auto *lk_225 = buffer.data(lk + 225);
    const auto *lk_226 = buffer.data(lk + 226);
    const auto *lk_227 = buffer.data(lk + 227);
    const auto *lk_228 = buffer.data(lk + 228);
    const auto *lk_230 = buffer.data(lk + 230);
    const auto *lk_231 = buffer.data(lk + 231);
    const auto *lk_232 = buffer.data(lk + 232);
    const auto *lk_233 = buffer.data(lk + 233);
    const auto *lk_234 = buffer.data(lk + 234);
    const auto *lk_236 = buffer.data(lk + 236);
    const auto *lk_244 = buffer.data(lk + 244);
    const auto *lk_245 = buffer.data(lk + 245);
    const auto *lk_246 = buffer.data(lk + 246);
    const auto *lk_247 = buffer.data(lk + 247);
    const auto *lk_248 = buffer.data(lk + 248);
    const auto *lk_249 = buffer.data(lk + 249);
    const auto *lk_251 = buffer.data(lk + 251);
    const auto *lk_254 = buffer.data(lk + 254);
    const auto *lk_257 = buffer.data(lk + 257);
    const auto *lk_261 = buffer.data(lk + 261);
    const auto *lk_266 = buffer.data(lk + 266);
    const auto *lk_272 = buffer.data(lk + 272);
    const auto *lk_287 = buffer.data(lk + 287);
    const auto *lk_288 = buffer.data(lk + 288);
    const auto *lk_333 = buffer.data(lk + 333);
    const auto *lk_338 = buffer.data(lk + 338);
    const auto *lk_344 = buffer.data(lk + 344);
    const auto *lk_351 = buffer.data(lk + 351);
    const auto *lk_352 = buffer.data(lk + 352);
    const auto *lk_353 = buffer.data(lk + 353);
    const auto *lk_354 = buffer.data(lk + 354);
    const auto *lk_355 = buffer.data(lk + 355);
    const auto *lk_356 = buffer.data(lk + 356);
    const auto *lk_357 = buffer.data(lk + 357);
    const auto *lk_359 = buffer.data(lk + 359);
    const auto *lk_363 = buffer.data(lk + 363);
    const auto *lk_366 = buffer.data(lk + 366);
    const auto *lk_370 = buffer.data(lk + 370);
    const auto *lk_375 = buffer.data(lk + 375);
    const auto *lk_381 = buffer.data(lk + 381);
    const auto *lk_388 = buffer.data(lk + 388);
    const auto *lk_390 = buffer.data(lk + 390);
    const auto *lk_391 = buffer.data(lk + 391);
    const auto *lk_392 = buffer.data(lk + 392);
    const auto *lk_393 = buffer.data(lk + 393);
    const auto *lk_394 = buffer.data(lk + 394);
    const auto *lk_395 = buffer.data(lk + 395);
    const auto *lk_425 = buffer.data(lk + 425);
    const auto *lk_426 = buffer.data(lk + 426);
    const auto *lk_427 = buffer.data(lk + 427);
    const auto *lk_428 = buffer.data(lk + 428);
    const auto *lk_429 = buffer.data(lk + 429);
    const auto *lk_430 = buffer.data(lk + 430);
    const auto *lk_431 = buffer.data(lk + 431);

    const auto *ll_270 = buffer.data(ll + 270);
    const auto *ll_271 = buffer.data(ll + 271);
    const auto *ll_273 = buffer.data(ll + 273);
    const auto *ll_275 = buffer.data(ll + 275);
    const auto *ll_276 = buffer.data(ll + 276);
    const auto *ll_279 = buffer.data(ll + 279);
    const auto *ll_280 = buffer.data(ll + 280);
    const auto *ll_282 = buffer.data(ll + 282);
    const auto *ll_284 = buffer.data(ll + 284);
    const auto *ll_285 = buffer.data(ll + 285);
    const auto *ll_287 = buffer.data(ll + 287);
    const auto *ll_288 = buffer.data(ll + 288);
    const auto *ll_290 = buffer.data(ll + 290);
    const auto *ll_291 = buffer.data(ll + 291);
    const auto *ll_293 = buffer.data(ll + 293);
    const auto *ll_294 = buffer.data(ll + 294);
    const auto *ll_295 = buffer.data(ll + 295);
    const auto *ll_297 = buffer.data(ll + 297);
    const auto *ll_298 = buffer.data(ll + 298);
    const auto *ll_306 = buffer.data(ll + 306);
    const auto *ll_308 = buffer.data(ll + 308);
    const auto *ll_309 = buffer.data(ll + 309);
    const auto *ll_310 = buffer.data(ll + 310);
    const auto *ll_311 = buffer.data(ll + 311);
    const auto *ll_312 = buffer.data(ll + 312);
    const auto *ll_314 = buffer.data(ll + 314);
    const auto *ll_360 = buffer.data(ll + 360);
    const auto *ll_449 = buffer.data(ll + 449);
    const auto *ll_486 = buffer.data(ll + 486);

    const auto *mi0_255 = buffer.data(mi0 + 255);
    const auto *mi0_257 = buffer.data(mi0 + 257);
    const auto *mi0_258 = buffer.data(mi0 + 258);
    const auto *mi0_260 = buffer.data(mi0 + 260);
    const auto *mi0_261 = buffer.data(mi0 + 261);
    const auto *mi0_262 = buffer.data(mi0 + 262);
    const auto *mi0_264 = buffer.data(mi0 + 264);
    const auto *mi0_265 = buffer.data(mi0 + 265);
    const auto *mi0_266 = buffer.data(mi0 + 266);
    const auto *mi0_272 = buffer.data(mi0 + 272);
    const auto *mi0_273 = buffer.data(mi0 + 273);
    const auto *mi0_275 = buffer.data(mi0 + 275);
    const auto *mi0_276 = buffer.data(mi0 + 276);
    const auto *mi0_277 = buffer.data(mi0 + 277);
    const auto *mi0_278 = buffer.data(mi0 + 278);
    const auto *mi0_279 = buffer.data(mi0 + 279);
    const auto *mi0_280 = buffer.data(mi0 + 280);
    const auto *mi0_282 = buffer.data(mi0 + 282);
    const auto *mi0_283 = buffer.data(mi0 + 283);
    const auto *mi0_285 = buffer.data(mi0 + 285);
    const auto *mi0_286 = buffer.data(mi0 + 286);
    const auto *mi0_287 = buffer.data(mi0 + 287);
    const auto *mi0_289 = buffer.data(mi0 + 289);
    const auto *mi0_290 = buffer.data(mi0 + 290);
    const auto *mi0_291 = buffer.data(mi0 + 291);
    const auto *mi0_292 = buffer.data(mi0 + 292);
    const auto *mi0_294 = buffer.data(mi0 + 294);
    const auto *mi0_295 = buffer.data(mi0 + 295);
    const auto *mi0_301 = buffer.data(mi0 + 301);
    const auto *mi0_302 = buffer.data(mi0 + 302);
    const auto *mi0_303 = buffer.data(mi0 + 303);
    const auto *mi0_304 = buffer.data(mi0 + 304);
    const auto *mi0_305 = buffer.data(mi0 + 305);
    const auto *mi0_307 = buffer.data(mi0 + 307);

    const auto *mi1_255 = buffer.data(mi1 + 255);
    const auto *mi1_257 = buffer.data(mi1 + 257);
    const auto *mi1_258 = buffer.data(mi1 + 258);
    const auto *mi1_260 = buffer.data(mi1 + 260);
    const auto *mi1_261 = buffer.data(mi1 + 261);
    const auto *mi1_262 = buffer.data(mi1 + 262);
    const auto *mi1_264 = buffer.data(mi1 + 264);
    const auto *mi1_265 = buffer.data(mi1 + 265);
    const auto *mi1_266 = buffer.data(mi1 + 266);
    const auto *mi1_272 = buffer.data(mi1 + 272);
    const auto *mi1_273 = buffer.data(mi1 + 273);
    const auto *mi1_275 = buffer.data(mi1 + 275);
    const auto *mi1_276 = buffer.data(mi1 + 276);
    const auto *mi1_277 = buffer.data(mi1 + 277);
    const auto *mi1_278 = buffer.data(mi1 + 278);
    const auto *mi1_279 = buffer.data(mi1 + 279);
    const auto *mi1_280 = buffer.data(mi1 + 280);
    const auto *mi1_282 = buffer.data(mi1 + 282);
    const auto *mi1_283 = buffer.data(mi1 + 283);
    const auto *mi1_285 = buffer.data(mi1 + 285);
    const auto *mi1_286 = buffer.data(mi1 + 286);
    const auto *mi1_287 = buffer.data(mi1 + 287);
    const auto *mi1_289 = buffer.data(mi1 + 289);
    const auto *mi1_290 = buffer.data(mi1 + 290);
    const auto *mi1_291 = buffer.data(mi1 + 291);
    const auto *mi1_292 = buffer.data(mi1 + 292);
    const auto *mi1_294 = buffer.data(mi1 + 294);
    const auto *mi1_295 = buffer.data(mi1 + 295);
    const auto *mi1_301 = buffer.data(mi1 + 301);
    const auto *mi1_302 = buffer.data(mi1 + 302);
    const auto *mi1_303 = buffer.data(mi1 + 303);
    const auto *mi1_304 = buffer.data(mi1 + 304);
    const auto *mi1_305 = buffer.data(mi1 + 305);
    const auto *mi1_307 = buffer.data(mi1 + 307);

    const auto *mk_329 = buffer.data(mk + 329);
    const auto *mk_330 = buffer.data(mk + 330);
    const auto *mk_332 = buffer.data(mk + 332);
    const auto *mk_333 = buffer.data(mk + 333);
    const auto *mk_334 = buffer.data(mk + 334);
    const auto *mk_336 = buffer.data(mk + 336);
    const auto *mk_337 = buffer.data(mk + 337);
    const auto *mk_338 = buffer.data(mk + 338);
    const auto *mk_339 = buffer.data(mk + 339);
    const auto *mk_341 = buffer.data(mk + 341);
    const auto *mk_342 = buffer.data(mk + 342);
    const auto *mk_343 = buffer.data(mk + 343);
    const auto *mk_344 = buffer.data(mk + 344);
    const auto *mk_351 = buffer.data(mk + 351);
    const auto *mk_352 = buffer.data(mk + 352);
    const auto *mk_353 = buffer.data(mk + 353);
    const auto *mk_354 = buffer.data(mk + 354);
    const auto *mk_355 = buffer.data(mk + 355);
    const auto *mk_356 = buffer.data(mk + 356);
    const auto *mk_357 = buffer.data(mk + 357);
    const auto *mk_358 = buffer.data(mk + 358);
    const auto *mk_359 = buffer.data(mk + 359);
    const auto *mk_360 = buffer.data(mk + 360);
    const auto *mk_361 = buffer.data(mk + 361);
    const auto *mk_362 = buffer.data(mk + 362);
    const auto *mk_363 = buffer.data(mk + 363);
    const auto *mk_365 = buffer.data(mk + 365);
    const auto *mk_366 = buffer.data(mk + 366);
    const auto *mk_367 = buffer.data(mk + 367);
    const auto *mk_369 = buffer.data(mk + 369);
    const auto *mk_370 = buffer.data(mk + 370);
    const auto *mk_371 = buffer.data(mk + 371);
    const auto *mk_372 = buffer.data(mk + 372);
    const auto *mk_374 = buffer.data(mk + 374);
    const auto *mk_375 = buffer.data(mk + 375);
    const auto *mk_376 = buffer.data(mk + 376);
    const auto *mk_377 = buffer.data(mk + 377);
    const auto *mk_378 = buffer.data(mk + 378);
    const auto *mk_380 = buffer.data(mk + 380);
    const auto *mk_381 = buffer.data(mk + 381);
    const auto *mk_388 = buffer.data(mk + 388);
    const auto *mk_389 = buffer.data(mk + 389);
    const auto *mk_390 = buffer.data(mk + 390);
    const auto *mk_391 = buffer.data(mk + 391);
    const auto *mk_392 = buffer.data(mk + 392);
    const auto *mk_393 = buffer.data(mk + 393);
    const auto *mk_394 = buffer.data(mk + 394);
    const auto *mk_395 = buffer.data(mk + 395);
    const auto *mk_396 = buffer.data(mk + 396);
    const auto *mk_398 = buffer.data(mk + 398);
    const auto *mk_399 = buffer.data(mk + 399);
    const auto *mk_401 = buffer.data(mk + 401);
    const auto *mk_402 = buffer.data(mk + 402);
    const auto *mk_405 = buffer.data(mk + 405);
    const auto *mk_406 = buffer.data(mk + 406);
    const auto *mk_410 = buffer.data(mk + 410);
    const auto *mk_411 = buffer.data(mk + 411);
    const auto *mk_416 = buffer.data(mk + 416);
    const auto *mk_424 = buffer.data(mk + 424);
    const auto *mk_425 = buffer.data(mk + 425);
    const auto *mk_426 = buffer.data(mk + 426);
    const auto *mk_427 = buffer.data(mk + 427);
    const auto *mk_428 = buffer.data(mk + 428);
    const auto *mk_429 = buffer.data(mk + 429);
    const auto *mk_430 = buffer.data(mk + 430);
    const auto *mk_431 = buffer.data(mk + 431);
    const auto *mk_432 = buffer.data(mk + 432);

#pragma omp simd aligned(t_413, t_414, t_415, t_416, pb_x, pb_y, pb_z, lk_186, lk_333, \
                         mi0_255, mi0_261, mi1_255, mi1_261, mk_329, mk_330, \
                         mk_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_y[k] * mk_329[k];

        t_414[k] = f_18 * lk_333[k]
                   + f_9 * mi0_261[k]
                   - f_10 * mi1_261[k]
                   + pb_x[k] * mk_333[k];

        t_415[k] = f_7 * mi0_255[k]
                   - f_8 * mi1_255[k]
                   + pb_y[k] * mk_330[k];

        t_416[k] = f_15 * lk_186[k]
                   + pb_z[k] * mk_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_x, pb_y, lk_338, mi0_257, mi0_266, mi1_257, \
                         mi1_266, mk_332, mk_333, mk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_3 * mi0_257[k]
                   - f_4 * mi1_257[k]
                   + pb_y[k] * mk_332[k];

        t_418[k] = pb_y[k] * mk_333[k];

        t_419[k] = f_18 * lk_338[k]
                   + f_7 * mi0_266[k]
                   - f_8 * mi1_266[k]
                   + pb_x[k] * mk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_y, pb_z, lk_190, mi0_258, mi0_260, \
                         mi0_261, mi1_258, mi1_260, mi1_261, mk_334, mk_336, \
                         mk_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_9 * mi0_258[k]
                   - f_10 * mi1_258[k]
                   + pb_y[k] * mk_334[k];

        t_421[k] = f_15 * lk_190[k]
                   + pb_z[k] * mk_334[k];

        t_422[k] = f_5 * mi0_260[k]
                   - f_6 * mi1_260[k]
                   + pb_y[k] * mk_336[k];

        t_423[k] = f_3 * mi0_261[k]
                   - f_4 * mi1_261[k]
                   + pb_y[k] * mk_337[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pb_x, pb_y, pb_z, lk_195, lk_344, \
                         mi0_262, mi0_272, mi1_262, mi1_272, mk_338, mk_339, \
                         mk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = pb_y[k] * mk_338[k];

        t_425[k] = f_18 * lk_344[k]
                   + f_5 * mi0_272[k]
                   - f_6 * mi1_272[k]
                   + pb_x[k] * mk_344[k];

        t_426[k] = f_11 * mi0_262[k]
                   - f_12 * mi1_262[k]
                   + pb_y[k] * mk_339[k];

        t_427[k] = f_15 * lk_195[k]
                   + pb_z[k] * mk_339[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pb_y, mi0_264, mi0_265, mi0_266, mi1_264, \
                         mi1_265, mi1_266, mk_341, mk_342, mk_343, \
                         mk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_7 * mi0_264[k]
                   - f_8 * mi1_264[k]
                   + pb_y[k] * mk_341[k];

        t_429[k] = f_5 * mi0_265[k]
                   - f_6 * mi1_265[k]
                   + pb_y[k] * mk_342[k];

        t_430[k] = f_3 * mi0_266[k]
                   - f_4 * mi1_266[k]
                   + pb_y[k] * mk_343[k];

        t_431[k] = pb_y[k] * mk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_x, lk_351, lk_352, lk_353, lk_354, \
                         mi0_279, mi1_279, mk_351, mk_352, mk_353, \
                         mk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_18 * lk_351[k]
                   + f_3 * mi0_279[k]
                   - f_4 * mi1_279[k]
                   + pb_x[k] * mk_351[k];

        t_433[k] = f_18 * lk_352[k]
                   + pb_x[k] * mk_352[k];

        t_434[k] = f_18 * lk_353[k]
                   + pb_x[k] * mk_353[k];

        t_435[k] = f_18 * lk_354[k]
                   + pb_x[k] * mk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pb_x, pb_y, lk_355, lk_356, \
                         lk_357, lk_359, mk_351, mk_355, mk_356, mk_357, \
                         mk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_18 * lk_355[k]
                   + pb_x[k] * mk_355[k];

        t_437[k] = f_18 * lk_356[k]
                   + pb_x[k] * mk_356[k];

        t_438[k] = f_18 * lk_357[k]
                   + pb_x[k] * mk_357[k];

        t_439[k] = pb_y[k] * mk_351[k];

        t_440[k] = f_18 * lk_359[k]
                   + pb_x[k] * mk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_y, pb_z, lk_208, mi0_273, mi0_275, \
                         mi0_276, mi1_273, mi1_275, mi1_276, mk_352, mk_354, \
                         mk_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * mi0_273[k]
                   - f_2 * mi1_273[k]
                   + pb_y[k] * mk_352[k];

        t_442[k] = f_15 * lk_208[k]
                   + pb_z[k] * mk_352[k];

        t_443[k] = f_11 * mi0_275[k]
                   - f_12 * mi1_275[k]
                   + pb_y[k] * mk_354[k];

        t_444[k] = f_9 * mi0_276[k]
                   - f_10 * mi1_276[k]
                   + pb_y[k] * mk_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pb_y, mi0_277, mi0_278, mi0_279, mi1_277, \
                         mi1_278, mi1_279, mk_356, mk_357, mk_358, \
                         mk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_7 * mi0_277[k]
                   - f_8 * mi1_277[k]
                   + pb_y[k] * mk_356[k];

        t_446[k] = f_5 * mi0_278[k]
                   - f_6 * mi1_278[k]
                   + pb_y[k] * mk_357[k];

        t_447[k] = f_3 * mi0_279[k]
                   - f_4 * mi1_279[k]
                   + pb_y[k] * mk_358[k];

        t_448[k] = pb_y[k] * mk_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_x, pa_y, pb_y, pb_z, kl0_135, kl0_449, \
                         kl1_135, kl1_449, lk_216, ll_270, ll_449, \
                         mk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_27 * kl0_449[k]
                   - f_28 * kl1_449[k]
                   + pa_x[k] * ll_449[k];

        t_450[k] = f_29 * kl0_135[k]
                   - f_30 * kl1_135[k]
                   + pa_y[k] * ll_270[k];

        t_451[k] = f_16 * lk_216[k]
                   + pb_y[k] * mk_360[k];

        t_452[k] = pb_z[k] * mk_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pb_x, pb_z, lk_363, mi0_280, mi0_283, mi1_280, \
                         mi1_283, mk_361, mk_362, mk_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_17 * lk_363[k]
                   + f_11 * mi0_283[k]
                   - f_12 * mi1_283[k]
                   + pb_x[k] * mk_363[k];

        t_454[k] = pb_z[k] * mk_361[k];

        t_455[k] = f_3 * mi0_280[k]
                   - f_4 * mi1_280[k]
                   + pb_z[k] * mk_362[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_x, pb_y, pb_z, lk_221, lk_366, \
                         mi0_282, mi0_286, mi1_282, mi1_286, mk_363, mk_365, \
                         mk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_17 * lk_366[k]
                   + f_9 * mi0_286[k]
                   - f_10 * mi1_286[k]
                   + pb_x[k] * mk_366[k];

        t_457[k] = pb_z[k] * mk_363[k];

        t_458[k] = f_16 * lk_221[k]
                   + pb_y[k] * mk_365[k];

        t_459[k] = f_5 * mi0_282[k]
                   - f_6 * mi1_282[k]
                   + pb_z[k] * mk_365[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pb_x, pb_z, lk_370, mi0_283, mi0_290, mi1_283, \
                         mi1_290, mk_366, mk_367, mk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_17 * lk_370[k]
                   + f_7 * mi0_290[k]
                   - f_8 * mi1_290[k]
                   + pb_x[k] * mk_370[k];

        t_461[k] = pb_z[k] * mk_366[k];

        t_462[k] = f_3 * mi0_283[k]
                   - f_4 * mi1_283[k]
                   + pb_z[k] * mk_367[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pb_x, pb_y, pb_z, lk_225, lk_375, \
                         mi0_285, mi0_295, mi1_285, mi1_295, mk_369, mk_370, \
                         mk_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_16 * lk_225[k]
                   + pb_y[k] * mk_369[k];

        t_464[k] = f_7 * mi0_285[k]
                   - f_8 * mi1_285[k]
                   + pb_z[k] * mk_369[k];

        t_465[k] = f_17 * lk_375[k]
                   + f_5 * mi0_295[k]
                   - f_6 * mi1_295[k]
                   + pb_x[k] * mk_375[k];

        t_466[k] = pb_z[k] * mk_370[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, pb_y, pb_z, lk_230, mi0_286, mi0_287, \
                         mi0_289, mi1_286, mi1_287, mi1_289, mk_371, mk_372, \
                         mk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_3 * mi0_286[k]
                   - f_4 * mi1_286[k]
                   + pb_z[k] * mk_371[k];

        t_468[k] = f_5 * mi0_287[k]
                   - f_6 * mi1_287[k]
                   + pb_z[k] * mk_372[k];

        t_469[k] = f_16 * lk_230[k]
                   + pb_y[k] * mk_374[k];

        t_470[k] = f_9 * mi0_289[k]
                   - f_10 * mi1_289[k]
                   + pb_z[k] * mk_374[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pb_x, pb_z, lk_381, mi0_290, mi0_301, mi1_290, \
                         mi1_301, mk_375, mk_376, mk_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_17 * lk_381[k]
                   + f_3 * mi0_301[k]
                   - f_4 * mi1_301[k]
                   + pb_x[k] * mk_381[k];

        t_472[k] = pb_z[k] * mk_375[k];

        t_473[k] = f_3 * mi0_290[k]
                   - f_4 * mi1_290[k]
                   + pb_z[k] * mk_376[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pb_z, lk_236, mi0_291, mi0_292, \
                         mi0_294, mi1_291, mi1_292, mi1_294, mk_377, mk_378, \
                         mk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_5 * mi0_291[k]
                   - f_6 * mi1_291[k]
                   + pb_z[k] * mk_377[k];

        t_475[k] = f_7 * mi0_292[k]
                   - f_8 * mi1_292[k]
                   + pb_z[k] * mk_378[k];

        t_476[k] = f_16 * lk_236[k]
                   + pb_y[k] * mk_380[k];

        t_477[k] = f_11 * mi0_294[k]
                   - f_12 * mi1_294[k]
                   + pb_z[k] * mk_380[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pb_x, pb_z, lk_388, lk_390, \
                         lk_391, lk_392, mk_381, mk_388, mk_390, mk_391, \
                         mk_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_17 * lk_388[k]
                   + pb_x[k] * mk_388[k];

        t_479[k] = pb_z[k] * mk_381[k];

        t_480[k] = f_17 * lk_390[k]
                   + pb_x[k] * mk_390[k];

        t_481[k] = f_17 * lk_391[k]
                   + pb_x[k] * mk_391[k];

        t_482[k] = f_17 * lk_392[k]
                   + pb_x[k] * mk_392[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_x, pb_x, kl0_486, kl1_486, lk_393, \
                         lk_394, lk_395, ll_486, mk_393, mk_394, \
                         mk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_17 * lk_393[k]
                   + pb_x[k] * mk_393[k];

        t_484[k] = f_17 * lk_394[k]
                   + pb_x[k] * mk_394[k];

        t_485[k] = f_17 * lk_395[k]
                   + pb_x[k] * mk_395[k];

        t_486[k] = f_31 * kl0_486[k]
                   - f_32 * kl1_486[k]
                   + pa_x[k] * ll_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pb_z, mi0_301, mi0_302, mi0_303, mi1_301, \
                         mi1_302, mi1_303, mk_388, mk_389, mk_390, \
                         mk_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = pb_z[k] * mk_388[k];

        t_488[k] = f_3 * mi0_301[k]
                   - f_4 * mi1_301[k]
                   + pb_z[k] * mk_389[k];

        t_489[k] = f_5 * mi0_302[k]
                   - f_6 * mi1_302[k]
                   + pb_z[k] * mk_390[k];

        t_490[k] = f_7 * mi0_303[k]
                   - f_8 * mi1_303[k]
                   + pb_z[k] * mk_391[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pb_y, pb_z, lk_251, mi0_304, mi0_305, \
                         mi0_307, mi1_304, mi1_305, mi1_307, mk_392, mk_393, \
                         mk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_9 * mi0_304[k]
                   - f_10 * mi1_304[k]
                   + pb_z[k] * mk_392[k];

        t_492[k] = f_11 * mi0_305[k]
                   - f_12 * mi1_305[k]
                   + pb_z[k] * mk_393[k];

        t_493[k] = f_16 * lk_251[k]
                   + pb_y[k] * mk_395[k];

        t_494[k] = f_1 * mi0_307[k]
                   - f_2 * mi1_307[k]
                   + pb_z[k] * mk_395[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, pa_z, pb_y, pb_z, lk_216, lk_254, \
                         ll_270, ll_271, ll_273, mk_396, mk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * ll_270[k];

        t_496[k] = pa_z[k] * ll_271[k];

        t_497[k] = f_13 * lk_216[k]
                   + pb_z[k] * mk_396[k];

        t_498[k] = pa_z[k] * ll_273[k];

        t_499[k] = f_15 * lk_254[k]
                   + pb_y[k] * mk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_z, pb_y, pb_z, lk_218, lk_219, lk_257, \
                         ll_275, ll_276, mk_399, mk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_14 * lk_218[k]
                   + pa_z[k] * ll_275[k];

        t_501[k] = pa_z[k] * ll_276[k];

        t_502[k] = f_13 * lk_219[k]
                   + pb_z[k] * mk_399[k];

        t_503[k] = f_15 * lk_257[k]
                   + pb_y[k] * mk_401[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_z, pb_z, lk_221, lk_222, lk_223, \
                         ll_279, ll_280, ll_282, mk_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_15 * lk_221[k]
                   + pa_z[k] * ll_279[k];

        t_505[k] = pa_z[k] * ll_280[k];

        t_506[k] = f_13 * lk_222[k]
                   + pb_z[k] * mk_402[k];

        t_507[k] = f_14 * lk_223[k]
                   + pa_z[k] * ll_282[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_z, pb_y, pb_z, lk_225, lk_226, lk_261, \
                         ll_284, ll_285, mk_405, mk_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_15 * lk_261[k]
                   + pb_y[k] * mk_405[k];

        t_509[k] = f_16 * lk_225[k]
                   + pa_z[k] * ll_284[k];

        t_510[k] = pa_z[k] * ll_285[k];

        t_511[k] = f_13 * lk_226[k]
                   + pb_z[k] * mk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pa_z, pb_y, lk_227, lk_228, \
                         lk_230, lk_266, ll_287, ll_288, ll_290, ll_291, \
                         mk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_14 * lk_227[k]
                   + pa_z[k] * ll_287[k];

        t_513[k] = f_15 * lk_228[k]
                   + pa_z[k] * ll_288[k];

        t_514[k] = f_15 * lk_266[k]
                   + pb_y[k] * mk_410[k];

        t_515[k] = f_17 * lk_230[k]
                   + pa_z[k] * ll_290[k];

        t_516[k] = pa_z[k] * ll_291[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_z, pb_z, lk_231, lk_232, lk_233, \
                         lk_234, ll_293, ll_294, ll_295, mk_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_13 * lk_231[k]
                   + pb_z[k] * mk_411[k];

        t_518[k] = f_14 * lk_232[k]
                   + pa_z[k] * ll_293[k];

        t_519[k] = f_15 * lk_233[k]
                   + pa_z[k] * ll_294[k];

        t_520[k] = f_16 * lk_234[k]
                   + pa_z[k] * ll_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pa_z, pb_x, pb_y, lk_236, lk_272, lk_425, \
                         ll_297, ll_298, mk_416, mk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_15 * lk_272[k]
                   + pb_y[k] * mk_416[k];

        t_522[k] = f_18 * lk_236[k]
                   + pa_z[k] * ll_297[k];

        t_523[k] = pa_z[k] * ll_298[k];

        t_524[k] = f_17 * lk_425[k]
                   + pb_x[k] * mk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pb_x, lk_426, lk_427, lk_428, \
                         lk_429, lk_430, mk_426, mk_427, mk_428, mk_429, \
                         mk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_17 * lk_426[k]
                   + pb_x[k] * mk_426[k];

        t_526[k] = f_17 * lk_427[k]
                   + pb_x[k] * mk_427[k];

        t_527[k] = f_17 * lk_428[k]
                   + pb_x[k] * mk_428[k];

        t_528[k] = f_17 * lk_429[k]
                   + pb_x[k] * mk_429[k];

        t_529[k] = f_17 * lk_430[k]
                   + pb_x[k] * mk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pa_z, pb_x, pb_z, lk_244, lk_245, lk_431, \
                         ll_306, ll_308, mk_424, mk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_17 * lk_431[k]
                   + pb_x[k] * mk_431[k];

        t_531[k] = pa_z[k] * ll_306[k];

        t_532[k] = f_13 * lk_244[k]
                   + pb_z[k] * mk_424[k];

        t_533[k] = f_14 * lk_245[k]
                   + pa_z[k] * ll_308[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pa_z, lk_246, lk_247, lk_248, lk_249, \
                         ll_309, ll_310, ll_311, ll_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_15 * lk_246[k]
                   + pa_z[k] * ll_309[k];

        t_535[k] = f_16 * lk_247[k]
                   + pa_z[k] * ll_310[k];

        t_536[k] = f_17 * lk_248[k]
                   + pa_z[k] * ll_311[k];

        t_537[k] = f_18 * lk_249[k]
                   + pa_z[k] * ll_312[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_y, pa_z, pb_y, kl0_225, kl1_225, \
                         lk_251, lk_287, lk_288, ll_314, ll_360, mk_431, \
                         mk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_15 * lk_287[k]
                   + pb_y[k] * mk_431[k];

        t_539[k] = f_19 * lk_251[k]
                   + pa_z[k] * ll_314[k];

        t_540[k] = f_20 * kl0_225[k]
                   - f_21 * kl1_225[k]
                   + pa_y[k] * ll_360[k];

        t_541[k] = f_14 * lk_288[k]
                   + pb_y[k] * mk_432[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece4(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

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

    const auto *kl0_138 = buffer.data(kl0 + 138);
    const auto *kl0_141 = buffer.data(kl0 + 141);
    const auto *kl0_145 = buffer.data(kl0 + 145);
    const auto *kl0_150 = buffer.data(kl0 + 150);
    const auto *kl0_156 = buffer.data(kl0 + 156);
    const auto *kl0_225 = buffer.data(kl0 + 225);
    const auto *kl0_230 = buffer.data(kl0 + 230);
    const auto *kl0_234 = buffer.data(kl0 + 234);
    const auto *kl0_239 = buffer.data(kl0 + 239);
    const auto *kl0_245 = buffer.data(kl0 + 245);
    const auto *kl0_252 = buffer.data(kl0 + 252);
    const auto *kl0_576 = buffer.data(kl0 + 576);
    const auto *kl0_578 = buffer.data(kl0 + 578);
    const auto *kl0_579 = buffer.data(kl0 + 579);
    const auto *kl0_580 = buffer.data(kl0 + 580);
    const auto *kl0_581 = buffer.data(kl0 + 581);
    const auto *kl0_582 = buffer.data(kl0 + 582);
    const auto *kl0_584 = buffer.data(kl0 + 584);

    const auto *kl1_138 = buffer.data(kl1 + 138);
    const auto *kl1_141 = buffer.data(kl1 + 141);
    const auto *kl1_145 = buffer.data(kl1 + 145);
    const auto *kl1_150 = buffer.data(kl1 + 150);
    const auto *kl1_156 = buffer.data(kl1 + 156);
    const auto *kl1_225 = buffer.data(kl1 + 225);
    const auto *kl1_230 = buffer.data(kl1 + 230);
    const auto *kl1_234 = buffer.data(kl1 + 234);
    const auto *kl1_239 = buffer.data(kl1 + 239);
    const auto *kl1_245 = buffer.data(kl1 + 245);
    const auto *kl1_252 = buffer.data(kl1 + 252);
    const auto *kl1_576 = buffer.data(kl1 + 576);
    const auto *kl1_578 = buffer.data(kl1 + 578);
    const auto *kl1_579 = buffer.data(kl1 + 579);
    const auto *kl1_580 = buffer.data(kl1 + 580);
    const auto *kl1_581 = buffer.data(kl1 + 581);
    const auto *kl1_582 = buffer.data(kl1 + 582);
    const auto *kl1_584 = buffer.data(kl1 + 584);

    const auto *lk_252 = buffer.data(lk + 252);
    const auto *lk_255 = buffer.data(lk + 255);
    const auto *lk_258 = buffer.data(lk + 258);
    const auto *lk_262 = buffer.data(lk + 262);
    const auto *lk_267 = buffer.data(lk + 267);
    const auto *lk_280 = buffer.data(lk + 280);
    const auto *lk_290 = buffer.data(lk + 290);
    const auto *lk_291 = buffer.data(lk + 291);
    const auto *lk_293 = buffer.data(lk + 293);
    const auto *lk_294 = buffer.data(lk + 294);
    const auto *lk_297 = buffer.data(lk + 297);
    const auto *lk_298 = buffer.data(lk + 298);
    const auto *lk_302 = buffer.data(lk + 302);
    const auto *lk_303 = buffer.data(lk + 303);
    const auto *lk_308 = buffer.data(lk + 308);
    const auto *lk_316 = buffer.data(lk + 316);
    const auto *lk_323 = buffer.data(lk + 323);
    const auto *lk_324 = buffer.data(lk + 324);
    const auto *lk_325 = buffer.data(lk + 325);
    const auto *lk_326 = buffer.data(lk + 326);
    const auto *lk_327 = buffer.data(lk + 327);
    const auto *lk_329 = buffer.data(lk + 329);
    const auto *lk_330 = buffer.data(lk + 330);
    const auto *lk_332 = buffer.data(lk + 332);
    const auto *lk_333 = buffer.data(lk + 333);
    const auto *lk_334 = buffer.data(lk + 334);
    const auto *lk_336 = buffer.data(lk + 336);
    const auto *lk_337 = buffer.data(lk + 337);
    const auto *lk_338 = buffer.data(lk + 338);
    const auto *lk_339 = buffer.data(lk + 339);
    const auto *lk_341 = buffer.data(lk + 341);
    const auto *lk_342 = buffer.data(lk + 342);
    const auto *lk_343 = buffer.data(lk + 343);
    const auto *lk_344 = buffer.data(lk + 344);
    const auto *lk_352 = buffer.data(lk + 352);
    const auto *lk_354 = buffer.data(lk + 354);
    const auto *lk_355 = buffer.data(lk + 355);
    const auto *lk_356 = buffer.data(lk + 356);
    const auto *lk_357 = buffer.data(lk + 357);
    const auto *lk_358 = buffer.data(lk + 358);
    const auto *lk_359 = buffer.data(lk + 359);
    const auto *lk_444 = buffer.data(lk + 444);
    const auto *lk_449 = buffer.data(lk + 449);
    const auto *lk_450 = buffer.data(lk + 450);
    const auto *lk_455 = buffer.data(lk + 455);
    const auto *lk_456 = buffer.data(lk + 456);
    const auto *lk_457 = buffer.data(lk + 457);
    const auto *lk_460 = buffer.data(lk + 460);
    const auto *lk_461 = buffer.data(lk + 461);
    const auto *lk_462 = buffer.data(lk + 462);
    const auto *lk_463 = buffer.data(lk + 463);
    const auto *lk_464 = buffer.data(lk + 464);
    const auto *lk_465 = buffer.data(lk + 465);
    const auto *lk_466 = buffer.data(lk + 466);
    const auto *lk_467 = buffer.data(lk + 467);
    const auto *lk_496 = buffer.data(lk + 496);
    const auto *lk_497 = buffer.data(lk + 497);
    const auto *lk_498 = buffer.data(lk + 498);
    const auto *lk_499 = buffer.data(lk + 499);
    const auto *lk_500 = buffer.data(lk + 500);
    const auto *lk_501 = buffer.data(lk + 501);
    const auto *lk_502 = buffer.data(lk + 502);
    const auto *lk_509 = buffer.data(lk + 509);
    const auto *lk_513 = buffer.data(lk + 513);
    const auto *lk_518 = buffer.data(lk + 518);
    const auto *lk_524 = buffer.data(lk + 524);
    const auto *lk_531 = buffer.data(lk + 531);
    const auto *lk_532 = buffer.data(lk + 532);
    const auto *lk_533 = buffer.data(lk + 533);
    const auto *lk_534 = buffer.data(lk + 534);
    const auto *lk_535 = buffer.data(lk + 535);
    const auto *lk_536 = buffer.data(lk + 536);
    const auto *lk_537 = buffer.data(lk + 537);
    const auto *lk_539 = buffer.data(lk + 539);

    const auto *ll_318 = buffer.data(ll + 318);
    const auto *ll_321 = buffer.data(ll + 321);
    const auto *ll_325 = buffer.data(ll + 325);
    const auto *ll_330 = buffer.data(ll + 330);
    const auto *ll_336 = buffer.data(ll + 336);
    const auto *ll_365 = buffer.data(ll + 365);
    const auto *ll_369 = buffer.data(ll + 369);
    const auto *ll_374 = buffer.data(ll + 374);
    const auto *ll_380 = buffer.data(ll + 380);
    const auto *ll_387 = buffer.data(ll + 387);
    const auto *ll_405 = buffer.data(ll + 405);
    const auto *ll_407 = buffer.data(ll + 407);
    const auto *ll_408 = buffer.data(ll + 408);
    const auto *ll_410 = buffer.data(ll + 410);
    const auto *ll_411 = buffer.data(ll + 411);
    const auto *ll_414 = buffer.data(ll + 414);
    const auto *ll_415 = buffer.data(ll + 415);
    const auto *ll_417 = buffer.data(ll + 417);
    const auto *ll_419 = buffer.data(ll + 419);
    const auto *ll_420 = buffer.data(ll + 420);
    const auto *ll_422 = buffer.data(ll + 422);
    const auto *ll_423 = buffer.data(ll + 423);
    const auto *ll_425 = buffer.data(ll + 425);
    const auto *ll_426 = buffer.data(ll + 426);
    const auto *ll_428 = buffer.data(ll + 428);
    const auto *ll_429 = buffer.data(ll + 429);
    const auto *ll_430 = buffer.data(ll + 430);
    const auto *ll_432 = buffer.data(ll + 432);
    const auto *ll_440 = buffer.data(ll + 440);
    const auto *ll_441 = buffer.data(ll + 441);
    const auto *ll_443 = buffer.data(ll + 443);
    const auto *ll_444 = buffer.data(ll + 444);
    const auto *ll_445 = buffer.data(ll + 445);
    const auto *ll_446 = buffer.data(ll + 446);
    const auto *ll_447 = buffer.data(ll + 447);
    const auto *ll_449 = buffer.data(ll + 449);
    const auto *ll_576 = buffer.data(ll + 576);
    const auto *ll_578 = buffer.data(ll + 578);
    const auto *ll_579 = buffer.data(ll + 579);
    const auto *ll_580 = buffer.data(ll + 580);
    const auto *ll_581 = buffer.data(ll + 581);
    const auto *ll_582 = buffer.data(ll + 582);
    const auto *ll_584 = buffer.data(ll + 584);

    const auto *mi0_348 = buffer.data(mi0 + 348);
    const auto *mi0_353 = buffer.data(mi0 + 353);
    const auto *mi0_354 = buffer.data(mi0 + 354);
    const auto *mi0_359 = buffer.data(mi0 + 359);
    const auto *mi0_360 = buffer.data(mi0 + 360);
    const auto *mi0_361 = buffer.data(mi0 + 361);
    const auto *mi0_392 = buffer.data(mi0 + 392);
    const auto *mi0_393 = buffer.data(mi0 + 393);
    const auto *mi0_395 = buffer.data(mi0 + 395);
    const auto *mi0_397 = buffer.data(mi0 + 397);
    const auto *mi0_398 = buffer.data(mi0 + 398);
    const auto *mi0_400 = buffer.data(mi0 + 400);
    const auto *mi0_401 = buffer.data(mi0 + 401);
    const auto *mi0_402 = buffer.data(mi0 + 402);
    const auto *mi0_404 = buffer.data(mi0 + 404);
    const auto *mi0_405 = buffer.data(mi0 + 405);
    const auto *mi0_406 = buffer.data(mi0 + 406);
    const auto *mi0_412 = buffer.data(mi0 + 412);
    const auto *mi0_419 = buffer.data(mi0 + 419);

    const auto *mi1_348 = buffer.data(mi1 + 348);
    const auto *mi1_353 = buffer.data(mi1 + 353);
    const auto *mi1_354 = buffer.data(mi1 + 354);
    const auto *mi1_359 = buffer.data(mi1 + 359);
    const auto *mi1_360 = buffer.data(mi1 + 360);
    const auto *mi1_361 = buffer.data(mi1 + 361);
    const auto *mi1_392 = buffer.data(mi1 + 392);
    const auto *mi1_393 = buffer.data(mi1 + 393);
    const auto *mi1_395 = buffer.data(mi1 + 395);
    const auto *mi1_397 = buffer.data(mi1 + 397);
    const auto *mi1_398 = buffer.data(mi1 + 398);
    const auto *mi1_400 = buffer.data(mi1 + 400);
    const auto *mi1_401 = buffer.data(mi1 + 401);
    const auto *mi1_402 = buffer.data(mi1 + 402);
    const auto *mi1_404 = buffer.data(mi1 + 404);
    const auto *mi1_405 = buffer.data(mi1 + 405);
    const auto *mi1_406 = buffer.data(mi1 + 406);
    const auto *mi1_412 = buffer.data(mi1 + 412);
    const auto *mi1_419 = buffer.data(mi1 + 419);

    const auto *mk_432 = buffer.data(mk + 432);
    const auto *mk_434 = buffer.data(mk + 434);
    const auto *mk_435 = buffer.data(mk + 435);
    const auto *mk_437 = buffer.data(mk + 437);
    const auto *mk_438 = buffer.data(mk + 438);
    const auto *mk_441 = buffer.data(mk + 441);
    const auto *mk_442 = buffer.data(mk + 442);
    const auto *mk_444 = buffer.data(mk + 444);
    const auto *mk_446 = buffer.data(mk + 446);
    const auto *mk_447 = buffer.data(mk + 447);
    const auto *mk_449 = buffer.data(mk + 449);
    const auto *mk_450 = buffer.data(mk + 450);
    const auto *mk_452 = buffer.data(mk + 452);
    const auto *mk_455 = buffer.data(mk + 455);
    const auto *mk_456 = buffer.data(mk + 456);
    const auto *mk_457 = buffer.data(mk + 457);
    const auto *mk_460 = buffer.data(mk + 460);
    const auto *mk_461 = buffer.data(mk + 461);
    const auto *mk_462 = buffer.data(mk + 462);
    const auto *mk_463 = buffer.data(mk + 463);
    const auto *mk_464 = buffer.data(mk + 464);
    const auto *mk_465 = buffer.data(mk + 465);
    const auto *mk_466 = buffer.data(mk + 466);
    const auto *mk_467 = buffer.data(mk + 467);
    const auto *mk_468 = buffer.data(mk + 468);
    const auto *mk_470 = buffer.data(mk + 470);
    const auto *mk_471 = buffer.data(mk + 471);
    const auto *mk_473 = buffer.data(mk + 473);
    const auto *mk_474 = buffer.data(mk + 474);
    const auto *mk_477 = buffer.data(mk + 477);
    const auto *mk_478 = buffer.data(mk + 478);
    const auto *mk_482 = buffer.data(mk + 482);
    const auto *mk_483 = buffer.data(mk + 483);
    const auto *mk_488 = buffer.data(mk + 488);
    const auto *mk_496 = buffer.data(mk + 496);
    const auto *mk_497 = buffer.data(mk + 497);
    const auto *mk_498 = buffer.data(mk + 498);
    const auto *mk_499 = buffer.data(mk + 499);
    const auto *mk_500 = buffer.data(mk + 500);
    const auto *mk_501 = buffer.data(mk + 501);
    const auto *mk_502 = buffer.data(mk + 502);
    const auto *mk_503 = buffer.data(mk + 503);
    const auto *mk_504 = buffer.data(mk + 504);
    const auto *mk_505 = buffer.data(mk + 505);
    const auto *mk_506 = buffer.data(mk + 506);
    const auto *mk_507 = buffer.data(mk + 507);
    const auto *mk_509 = buffer.data(mk + 509);
    const auto *mk_510 = buffer.data(mk + 510);
    const auto *mk_512 = buffer.data(mk + 512);
    const auto *mk_513 = buffer.data(mk + 513);
    const auto *mk_514 = buffer.data(mk + 514);
    const auto *mk_516 = buffer.data(mk + 516);
    const auto *mk_517 = buffer.data(mk + 517);
    const auto *mk_518 = buffer.data(mk + 518);
    const auto *mk_519 = buffer.data(mk + 519);
    const auto *mk_521 = buffer.data(mk + 521);
    const auto *mk_522 = buffer.data(mk + 522);
    const auto *mk_523 = buffer.data(mk + 523);
    const auto *mk_524 = buffer.data(mk + 524);
    const auto *mk_531 = buffer.data(mk + 531);
    const auto *mk_532 = buffer.data(mk + 532);
    const auto *mk_533 = buffer.data(mk + 533);
    const auto *mk_534 = buffer.data(mk + 534);
    const auto *mk_535 = buffer.data(mk + 535);
    const auto *mk_536 = buffer.data(mk + 536);
    const auto *mk_537 = buffer.data(mk + 537);
    const auto *mk_539 = buffer.data(mk + 539);

#pragma omp simd aligned(t_542, t_543, t_544, pa_z, pb_y, pb_z, kl0_138, kl1_138, lk_252, \
                         lk_290, ll_318, mk_432, mk_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_14 * lk_252[k]
                   + pb_z[k] * mk_432[k];

        t_543[k] = f_20 * kl0_138[k]
                   - f_21 * kl1_138[k]
                   + pa_z[k] * ll_318[k];

        t_544[k] = f_14 * lk_290[k]
                   + pb_y[k] * mk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pa_y, pa_z, pb_z, kl0_141, kl0_230, kl1_141, \
                         kl1_230, lk_255, ll_321, ll_365, mk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_20 * kl0_230[k]
                   - f_21 * kl1_230[k]
                   + pa_y[k] * ll_365[k];

        t_546[k] = f_20 * kl0_141[k]
                   - f_21 * kl1_141[k]
                   + pa_z[k] * ll_321[k];

        t_547[k] = f_14 * lk_255[k]
                   + pb_z[k] * mk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pa_y, pa_z, pb_y, kl0_145, kl0_234, kl1_145, \
                         kl1_234, lk_293, ll_325, ll_369, mk_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_14 * lk_293[k]
                   + pb_y[k] * mk_437[k];

        t_549[k] = f_20 * kl0_234[k]
                   - f_21 * kl1_234[k]
                   + pa_y[k] * ll_369[k];

        t_550[k] = f_20 * kl0_145[k]
                   - f_21 * kl1_145[k]
                   + pa_z[k] * ll_325[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pb_x, pb_y, pb_z, lk_258, lk_297, lk_444, \
                         mi0_348, mi1_348, mk_438, mk_441, mk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_14 * lk_258[k]
                   + pb_z[k] * mk_438[k];

        t_552[k] = f_17 * lk_444[k]
                   + f_7 * mi0_348[k]
                   - f_8 * mi1_348[k]
                   + pb_x[k] * mk_444[k];

        t_553[k] = f_14 * lk_297[k]
                   + pb_y[k] * mk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pa_y, pa_z, pb_z, kl0_150, kl0_239, kl1_150, \
                         kl1_239, lk_262, ll_330, ll_374, mk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_20 * kl0_239[k]
                   - f_21 * kl1_239[k]
                   + pa_y[k] * ll_374[k];

        t_555[k] = f_20 * kl0_150[k]
                   - f_21 * kl1_150[k]
                   + pa_z[k] * ll_330[k];

        t_556[k] = f_14 * lk_262[k]
                   + pb_z[k] * mk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pb_x, pb_y, lk_302, lk_449, lk_450, mi0_353, \
                         mi0_354, mi1_353, mi1_354, mk_446, mk_449, \
                         mk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_17 * lk_449[k]
                   + f_5 * mi0_353[k]
                   - f_6 * mi1_353[k]
                   + pb_x[k] * mk_449[k];

        t_558[k] = f_17 * lk_450[k]
                   + f_5 * mi0_354[k]
                   - f_6 * mi1_354[k]
                   + pb_x[k] * mk_450[k];

        t_559[k] = f_14 * lk_302[k]
                   + pb_y[k] * mk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pa_y, pa_z, pb_z, kl0_156, kl0_245, kl1_156, \
                         kl1_245, lk_267, ll_336, ll_380, mk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_20 * kl0_245[k]
                   - f_21 * kl1_245[k]
                   + pa_y[k] * ll_380[k];

        t_561[k] = f_20 * kl0_156[k]
                   - f_21 * kl1_156[k]
                   + pa_z[k] * ll_336[k];

        t_562[k] = f_14 * lk_267[k]
                   + pb_z[k] * mk_447[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pb_x, lk_455, lk_456, lk_457, mi0_359, mi0_360, \
                         mi0_361, mi1_359, mi1_360, mi1_361, mk_455, mk_456, \
                         mk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_17 * lk_455[k]
                   + f_3 * mi0_359[k]
                   - f_4 * mi1_359[k]
                   + pb_x[k] * mk_455[k];

        t_564[k] = f_17 * lk_456[k]
                   + f_3 * mi0_360[k]
                   - f_4 * mi1_360[k]
                   + pb_x[k] * mk_456[k];

        t_565[k] = f_17 * lk_457[k]
                   + f_3 * mi0_361[k]
                   - f_4 * mi1_361[k]
                   + pb_x[k] * mk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_y, pb_x, pb_y, kl0_252, kl1_252, \
                         lk_308, lk_460, lk_461, ll_387, mk_452, mk_460, \
                         mk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_14 * lk_308[k]
                   + pb_y[k] * mk_452[k];

        t_567[k] = f_20 * kl0_252[k]
                   - f_21 * kl1_252[k]
                   + pa_y[k] * ll_387[k];

        t_568[k] = f_17 * lk_460[k]
                   + pb_x[k] * mk_460[k];

        t_569[k] = f_17 * lk_461[k]
                   + pb_x[k] * mk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pb_x, lk_462, lk_463, lk_464, \
                         lk_465, lk_466, mk_462, mk_463, mk_464, mk_465, \
                         mk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_17 * lk_462[k]
                   + pb_x[k] * mk_462[k];

        t_571[k] = f_17 * lk_463[k]
                   + pb_x[k] * mk_463[k];

        t_572[k] = f_17 * lk_464[k]
                   + pb_x[k] * mk_464[k];

        t_573[k] = f_17 * lk_465[k]
                   + pb_x[k] * mk_465[k];

        t_574[k] = f_17 * lk_466[k]
                   + pb_x[k] * mk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pa_x, pb_x, pb_z, kl0_576, kl1_576, lk_280, \
                         lk_467, ll_576, mk_460, mk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_17 * lk_467[k]
                   + pb_x[k] * mk_467[k];

        t_576[k] = f_31 * kl0_576[k]
                   - f_32 * kl1_576[k]
                   + pa_x[k] * ll_576[k];

        t_577[k] = f_14 * lk_280[k]
                   + pb_z[k] * mk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pa_x, kl0_578, kl0_579, kl0_580, kl1_578, \
                         kl1_579, kl1_580, ll_578, ll_579, ll_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_31 * kl0_578[k]
                   - f_32 * kl1_578[k]
                   + pa_x[k] * ll_578[k];

        t_579[k] = f_31 * kl0_579[k]
                   - f_32 * kl1_579[k]
                   + pa_x[k] * ll_579[k];

        t_580[k] = f_31 * kl0_580[k]
                   - f_32 * kl1_580[k]
                   + pa_x[k] * ll_580[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pa_x, pb_y, kl0_581, kl0_582, kl1_581, kl1_582, \
                         lk_323, ll_581, ll_582, mk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_31 * kl0_581[k]
                   - f_32 * kl1_581[k]
                   + pa_x[k] * ll_581[k];

        t_582[k] = f_31 * kl0_582[k]
                   - f_32 * kl1_582[k]
                   + pa_x[k] * ll_582[k];

        t_583[k] = f_14 * lk_323[k]
                   + pb_y[k] * mk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_x, pa_y, pb_y, kl0_584, kl1_584, \
                         lk_324, ll_405, ll_407, ll_584, mk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_31 * kl0_584[k]
                   - f_32 * kl1_584[k]
                   + pa_x[k] * ll_584[k];

        t_585[k] = pa_y[k] * ll_405[k];

        t_586[k] = f_13 * lk_324[k]
                   + pb_y[k] * mk_468[k];

        t_587[k] = pa_y[k] * ll_407[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_y, pb_y, lk_325, lk_326, lk_327, \
                         ll_408, ll_410, ll_411, mk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_14 * lk_325[k]
                   + pa_y[k] * ll_408[k];

        t_589[k] = f_13 * lk_326[k]
                   + pb_y[k] * mk_470[k];

        t_590[k] = pa_y[k] * ll_410[k];

        t_591[k] = f_15 * lk_327[k]
                   + pa_y[k] * ll_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_y, pb_y, pb_z, lk_291, lk_329, lk_330, \
                         ll_414, ll_415, mk_471, mk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_15 * lk_291[k]
                   + pb_z[k] * mk_471[k];

        t_593[k] = f_13 * lk_329[k]
                   + pb_y[k] * mk_473[k];

        t_594[k] = pa_y[k] * ll_414[k];

        t_595[k] = f_16 * lk_330[k]
                   + pa_y[k] * ll_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_y, pb_y, pb_z, lk_294, lk_332, lk_333, \
                         ll_417, ll_419, mk_474, mk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * lk_294[k]
                   + pb_z[k] * mk_474[k];

        t_597[k] = f_14 * lk_332[k]
                   + pa_y[k] * ll_417[k];

        t_598[k] = f_13 * lk_333[k]
                   + pb_y[k] * mk_477[k];

        t_599[k] = pa_y[k] * ll_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pa_y, pb_z, lk_298, lk_334, lk_336, \
                         lk_337, ll_420, ll_422, ll_423, mk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_17 * lk_334[k]
                   + pa_y[k] * ll_420[k];

        t_601[k] = f_15 * lk_298[k]
                   + pb_z[k] * mk_478[k];

        t_602[k] = f_15 * lk_336[k]
                   + pa_y[k] * ll_422[k];

        t_603[k] = f_14 * lk_337[k]
                   + pa_y[k] * ll_423[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_y, pb_y, pb_z, lk_303, lk_338, lk_339, \
                         ll_425, ll_426, mk_482, mk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_13 * lk_338[k]
                   + pb_y[k] * mk_482[k];

        t_605[k] = pa_y[k] * ll_425[k];

        t_606[k] = f_18 * lk_339[k]
                   + pa_y[k] * ll_426[k];

        t_607[k] = f_15 * lk_303[k]
                   + pb_z[k] * mk_483[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, t_612, pa_y, pb_y, lk_341, lk_342, \
                         lk_343, lk_344, ll_428, ll_429, ll_430, ll_432, \
                         mk_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_16 * lk_341[k]
                   + pa_y[k] * ll_428[k];

        t_609[k] = f_15 * lk_342[k]
                   + pa_y[k] * ll_429[k];

        t_610[k] = f_14 * lk_343[k]
                   + pa_y[k] * ll_430[k];

        t_611[k] = f_13 * lk_344[k]
                   + pb_y[k] * mk_488[k];

        t_612[k] = pa_y[k] * ll_432[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, t_617, pb_x, lk_496, lk_497, lk_498, \
                         lk_499, lk_500, mk_496, mk_497, mk_498, mk_499, \
                         mk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_17 * lk_496[k]
                   + pb_x[k] * mk_496[k];

        t_614[k] = f_17 * lk_497[k]
                   + pb_x[k] * mk_497[k];

        t_615[k] = f_17 * lk_498[k]
                   + pb_x[k] * mk_498[k];

        t_616[k] = f_17 * lk_499[k]
                   + pb_x[k] * mk_499[k];

        t_617[k] = f_17 * lk_500[k]
                   + pb_x[k] * mk_500[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pa_y, pb_x, lk_352, lk_501, lk_502, \
                         ll_440, ll_441, mk_501, mk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_17 * lk_501[k]
                   + pb_x[k] * mk_501[k];

        t_619[k] = f_17 * lk_502[k]
                   + pb_x[k] * mk_502[k];

        t_620[k] = pa_y[k] * ll_440[k];

        t_621[k] = f_19 * lk_352[k]
                   + pa_y[k] * ll_441[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_y, pb_z, lk_316, lk_354, lk_355, \
                         lk_356, ll_443, ll_444, ll_445, mk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_15 * lk_316[k]
                   + pb_z[k] * mk_496[k];

        t_623[k] = f_18 * lk_354[k]
                   + pa_y[k] * ll_443[k];

        t_624[k] = f_17 * lk_355[k]
                   + pa_y[k] * ll_444[k];

        t_625[k] = f_16 * lk_356[k]
                   + pa_y[k] * ll_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pa_y, pb_y, lk_357, lk_358, lk_359, \
                         ll_446, ll_447, ll_449, mk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * lk_357[k]
                   + pa_y[k] * ll_446[k];

        t_627[k] = f_14 * lk_358[k]
                   + pa_y[k] * ll_447[k];

        t_628[k] = f_13 * lk_359[k]
                   + pb_y[k] * mk_503[k];

        t_629[k] = pa_y[k] * ll_449[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pa_z, pb_y, pb_z, kl0_225, kl1_225, \
                         lk_324, ll_405, mi0_392, mi1_392, mk_504, \
                         mk_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_29 * kl0_225[k]
                   - f_30 * kl1_225[k]
                   + pa_z[k] * ll_405[k];

        t_631[k] = pb_y[k] * mk_504[k];

        t_632[k] = f_16 * lk_324[k]
                   + pb_z[k] * mk_504[k];

        t_633[k] = f_3 * mi0_392[k]
                   - f_4 * mi1_392[k]
                   + pb_y[k] * mk_505[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pb_x, pb_y, pb_z, lk_327, lk_509, \
                         mi0_393, mi0_397, mi1_393, mi1_397, mk_506, mk_507, \
                         mk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = pb_y[k] * mk_506[k];

        t_635[k] = f_17 * lk_509[k]
                   + f_11 * mi0_397[k]
                   - f_12 * mi1_397[k]
                   + pb_x[k] * mk_509[k];

        t_636[k] = f_5 * mi0_393[k]
                   - f_6 * mi1_393[k]
                   + pb_y[k] * mk_507[k];

        t_637[k] = f_16 * lk_327[k]
                   + pb_z[k] * mk_507[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pb_x, pb_y, pb_z, lk_330, lk_513, \
                         mi0_395, mi0_401, mi1_395, mi1_401, mk_509, mk_510, \
                         mk_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = pb_y[k] * mk_509[k];

        t_639[k] = f_17 * lk_513[k]
                   + f_9 * mi0_401[k]
                   - f_10 * mi1_401[k]
                   + pb_x[k] * mk_513[k];

        t_640[k] = f_7 * mi0_395[k]
                   - f_8 * mi1_395[k]
                   + pb_y[k] * mk_510[k];

        t_641[k] = f_16 * lk_330[k]
                   + pb_z[k] * mk_510[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pb_x, pb_y, lk_518, mi0_397, mi0_406, mi1_397, \
                         mi1_406, mk_512, mk_513, mk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_3 * mi0_397[k]
                   - f_4 * mi1_397[k]
                   + pb_y[k] * mk_512[k];

        t_643[k] = pb_y[k] * mk_513[k];

        t_644[k] = f_17 * lk_518[k]
                   + f_7 * mi0_406[k]
                   - f_8 * mi1_406[k]
                   + pb_x[k] * mk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, pb_y, pb_z, lk_334, mi0_398, mi0_400, \
                         mi0_401, mi1_398, mi1_400, mi1_401, mk_514, mk_516, \
                         mk_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_9 * mi0_398[k]
                   - f_10 * mi1_398[k]
                   + pb_y[k] * mk_514[k];

        t_646[k] = f_16 * lk_334[k]
                   + pb_z[k] * mk_514[k];

        t_647[k] = f_5 * mi0_400[k]
                   - f_6 * mi1_400[k]
                   + pb_y[k] * mk_516[k];

        t_648[k] = f_3 * mi0_401[k]
                   - f_4 * mi1_401[k]
                   + pb_y[k] * mk_517[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, pb_x, pb_y, pb_z, lk_339, lk_524, \
                         mi0_402, mi0_412, mi1_402, mi1_412, mk_518, mk_519, \
                         mk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pb_y[k] * mk_518[k];

        t_650[k] = f_17 * lk_524[k]
                   + f_5 * mi0_412[k]
                   - f_6 * mi1_412[k]
                   + pb_x[k] * mk_524[k];

        t_651[k] = f_11 * mi0_402[k]
                   - f_12 * mi1_402[k]
                   + pb_y[k] * mk_519[k];

        t_652[k] = f_16 * lk_339[k]
                   + pb_z[k] * mk_519[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, pb_y, mi0_404, mi0_405, mi0_406, mi1_404, \
                         mi1_405, mi1_406, mk_521, mk_522, mk_523, \
                         mk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_7 * mi0_404[k]
                   - f_8 * mi1_404[k]
                   + pb_y[k] * mk_521[k];

        t_654[k] = f_5 * mi0_405[k]
                   - f_6 * mi1_405[k]
                   + pb_y[k] * mk_522[k];

        t_655[k] = f_3 * mi0_406[k]
                   - f_4 * mi1_406[k]
                   + pb_y[k] * mk_523[k];

        t_656[k] = pb_y[k] * mk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pb_x, lk_531, lk_532, lk_533, lk_534, \
                         mi0_419, mi1_419, mk_531, mk_532, mk_533, \
                         mk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_17 * lk_531[k]
                   + f_3 * mi0_419[k]
                   - f_4 * mi1_419[k]
                   + pb_x[k] * mk_531[k];

        t_658[k] = f_17 * lk_532[k]
                   + pb_x[k] * mk_532[k];

        t_659[k] = f_17 * lk_533[k]
                   + pb_x[k] * mk_533[k];

        t_660[k] = f_17 * lk_534[k]
                   + pb_x[k] * mk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pb_x, pb_y, lk_535, lk_536, \
                         lk_537, lk_539, mk_531, mk_535, mk_536, mk_537, \
                         mk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_17 * lk_535[k]
                   + pb_x[k] * mk_535[k];

        t_662[k] = f_17 * lk_536[k]
                   + pb_x[k] * mk_536[k];

        t_663[k] = f_17 * lk_537[k]
                   + pb_x[k] * mk_537[k];

        t_664[k] = pb_y[k] * mk_531[k];

        t_665[k] = f_17 * lk_539[k]
                   + pb_x[k] * mk_539[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece5(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_270 = buffer.data(kl0 + 270);
    const auto *kl0_273 = buffer.data(kl0 + 273);
    const auto *kl0_276 = buffer.data(kl0 + 276);
    const auto *kl0_280 = buffer.data(kl0 + 280);
    const auto *kl0_285 = buffer.data(kl0 + 285);
    const auto *kl0_291 = buffer.data(kl0 + 291);
    const auto *kl0_360 = buffer.data(kl0 + 360);
    const auto *kl0_365 = buffer.data(kl0 + 365);
    const auto *kl0_369 = buffer.data(kl0 + 369);
    const auto *kl0_374 = buffer.data(kl0 + 374);
    const auto *kl0_380 = buffer.data(kl0 + 380);
    const auto *kl0_674 = buffer.data(kl0 + 674);
    const auto *kl0_711 = buffer.data(kl0 + 711);

    const auto *kl1_270 = buffer.data(kl1 + 270);
    const auto *kl1_273 = buffer.data(kl1 + 273);
    const auto *kl1_276 = buffer.data(kl1 + 276);
    const auto *kl1_280 = buffer.data(kl1 + 280);
    const auto *kl1_285 = buffer.data(kl1 + 285);
    const auto *kl1_291 = buffer.data(kl1 + 291);
    const auto *kl1_360 = buffer.data(kl1 + 360);
    const auto *kl1_365 = buffer.data(kl1 + 365);
    const auto *kl1_369 = buffer.data(kl1 + 369);
    const auto *kl1_374 = buffer.data(kl1 + 374);
    const auto *kl1_380 = buffer.data(kl1 + 380);
    const auto *kl1_674 = buffer.data(kl1 + 674);
    const auto *kl1_711 = buffer.data(kl1 + 711);

    const auto *lk_352 = buffer.data(lk + 352);
    const auto *lk_360 = buffer.data(lk + 360);
    const auto *lk_362 = buffer.data(lk + 362);
    const auto *lk_363 = buffer.data(lk + 363);
    const auto *lk_365 = buffer.data(lk + 365);
    const auto *lk_366 = buffer.data(lk + 366);
    const auto *lk_367 = buffer.data(lk + 367);
    const auto *lk_369 = buffer.data(lk + 369);
    const auto *lk_370 = buffer.data(lk + 370);
    const auto *lk_371 = buffer.data(lk + 371);
    const auto *lk_372 = buffer.data(lk + 372);
    const auto *lk_374 = buffer.data(lk + 374);
    const auto *lk_375 = buffer.data(lk + 375);
    const auto *lk_376 = buffer.data(lk + 376);
    const auto *lk_377 = buffer.data(lk + 377);
    const auto *lk_378 = buffer.data(lk + 378);
    const auto *lk_380 = buffer.data(lk + 380);
    const auto *lk_388 = buffer.data(lk + 388);
    const auto *lk_389 = buffer.data(lk + 389);
    const auto *lk_390 = buffer.data(lk + 390);
    const auto *lk_391 = buffer.data(lk + 391);
    const auto *lk_392 = buffer.data(lk + 392);
    const auto *lk_393 = buffer.data(lk + 393);
    const auto *lk_395 = buffer.data(lk + 395);
    const auto *lk_396 = buffer.data(lk + 396);
    const auto *lk_398 = buffer.data(lk + 398);
    const auto *lk_399 = buffer.data(lk + 399);
    const auto *lk_401 = buffer.data(lk + 401);
    const auto *lk_402 = buffer.data(lk + 402);
    const auto *lk_405 = buffer.data(lk + 405);
    const auto *lk_406 = buffer.data(lk + 406);
    const auto *lk_410 = buffer.data(lk + 410);
    const auto *lk_411 = buffer.data(lk + 411);
    const auto *lk_416 = buffer.data(lk + 416);
    const auto *lk_431 = buffer.data(lk + 431);
    const auto *lk_432 = buffer.data(lk + 432);
    const auto *lk_434 = buffer.data(lk + 434);
    const auto *lk_437 = buffer.data(lk + 437);
    const auto *lk_441 = buffer.data(lk + 441);
    const auto *lk_446 = buffer.data(lk + 446);
    const auto *lk_543 = buffer.data(lk + 543);
    const auto *lk_546 = buffer.data(lk + 546);
    const auto *lk_550 = buffer.data(lk + 550);
    const auto *lk_555 = buffer.data(lk + 555);
    const auto *lk_561 = buffer.data(lk + 561);
    const auto *lk_568 = buffer.data(lk + 568);
    const auto *lk_570 = buffer.data(lk + 570);
    const auto *lk_571 = buffer.data(lk + 571);
    const auto *lk_572 = buffer.data(lk + 572);
    const auto *lk_573 = buffer.data(lk + 573);
    const auto *lk_574 = buffer.data(lk + 574);
    const auto *lk_575 = buffer.data(lk + 575);
    const auto *lk_605 = buffer.data(lk + 605);
    const auto *lk_606 = buffer.data(lk + 606);
    const auto *lk_607 = buffer.data(lk + 607);
    const auto *lk_608 = buffer.data(lk + 608);
    const auto *lk_609 = buffer.data(lk + 609);
    const auto *lk_610 = buffer.data(lk + 610);
    const auto *lk_611 = buffer.data(lk + 611);
    const auto *lk_624 = buffer.data(lk + 624);
    const auto *lk_629 = buffer.data(lk + 629);
    const auto *lk_630 = buffer.data(lk + 630);
    const auto *lk_635 = buffer.data(lk + 635);
    const auto *lk_636 = buffer.data(lk + 636);
    const auto *lk_637 = buffer.data(lk + 637);

    const auto *ll_450 = buffer.data(ll + 450);
    const auto *ll_451 = buffer.data(ll + 451);
    const auto *ll_453 = buffer.data(ll + 453);
    const auto *ll_455 = buffer.data(ll + 455);
    const auto *ll_456 = buffer.data(ll + 456);
    const auto *ll_459 = buffer.data(ll + 459);
    const auto *ll_460 = buffer.data(ll + 460);
    const auto *ll_462 = buffer.data(ll + 462);
    const auto *ll_464 = buffer.data(ll + 464);
    const auto *ll_465 = buffer.data(ll + 465);
    const auto *ll_467 = buffer.data(ll + 467);
    const auto *ll_468 = buffer.data(ll + 468);
    const auto *ll_470 = buffer.data(ll + 470);
    const auto *ll_471 = buffer.data(ll + 471);
    const auto *ll_473 = buffer.data(ll + 473);
    const auto *ll_474 = buffer.data(ll + 474);
    const auto *ll_475 = buffer.data(ll + 475);
    const auto *ll_477 = buffer.data(ll + 477);
    const auto *ll_478 = buffer.data(ll + 478);
    const auto *ll_486 = buffer.data(ll + 486);
    const auto *ll_488 = buffer.data(ll + 488);
    const auto *ll_489 = buffer.data(ll + 489);
    const auto *ll_490 = buffer.data(ll + 490);
    const auto *ll_491 = buffer.data(ll + 491);
    const auto *ll_492 = buffer.data(ll + 492);
    const auto *ll_494 = buffer.data(ll + 494);
    const auto *ll_498 = buffer.data(ll + 498);
    const auto *ll_501 = buffer.data(ll + 501);
    const auto *ll_505 = buffer.data(ll + 505);
    const auto *ll_510 = buffer.data(ll + 510);
    const auto *ll_516 = buffer.data(ll + 516);
    const auto *ll_540 = buffer.data(ll + 540);
    const auto *ll_545 = buffer.data(ll + 545);
    const auto *ll_549 = buffer.data(ll + 549);
    const auto *ll_554 = buffer.data(ll + 554);
    const auto *ll_560 = buffer.data(ll + 560);
    const auto *ll_674 = buffer.data(ll + 674);
    const auto *ll_711 = buffer.data(ll + 711);

    const auto *mi0_413 = buffer.data(mi0 + 413);
    const auto *mi0_415 = buffer.data(mi0 + 415);
    const auto *mi0_416 = buffer.data(mi0 + 416);
    const auto *mi0_417 = buffer.data(mi0 + 417);
    const auto *mi0_418 = buffer.data(mi0 + 418);
    const auto *mi0_419 = buffer.data(mi0 + 419);
    const auto *mi0_420 = buffer.data(mi0 + 420);
    const auto *mi0_422 = buffer.data(mi0 + 422);
    const auto *mi0_423 = buffer.data(mi0 + 423);
    const auto *mi0_425 = buffer.data(mi0 + 425);
    const auto *mi0_426 = buffer.data(mi0 + 426);
    const auto *mi0_427 = buffer.data(mi0 + 427);
    const auto *mi0_429 = buffer.data(mi0 + 429);
    const auto *mi0_430 = buffer.data(mi0 + 430);
    const auto *mi0_431 = buffer.data(mi0 + 431);
    const auto *mi0_432 = buffer.data(mi0 + 432);
    const auto *mi0_434 = buffer.data(mi0 + 434);
    const auto *mi0_435 = buffer.data(mi0 + 435);
    const auto *mi0_441 = buffer.data(mi0 + 441);
    const auto *mi0_442 = buffer.data(mi0 + 442);
    const auto *mi0_443 = buffer.data(mi0 + 443);
    const auto *mi0_444 = buffer.data(mi0 + 444);
    const auto *mi0_445 = buffer.data(mi0 + 445);
    const auto *mi0_447 = buffer.data(mi0 + 447);
    const auto *mi0_488 = buffer.data(mi0 + 488);
    const auto *mi0_493 = buffer.data(mi0 + 493);
    const auto *mi0_494 = buffer.data(mi0 + 494);
    const auto *mi0_499 = buffer.data(mi0 + 499);
    const auto *mi0_500 = buffer.data(mi0 + 500);
    const auto *mi0_501 = buffer.data(mi0 + 501);

    const auto *mi1_413 = buffer.data(mi1 + 413);
    const auto *mi1_415 = buffer.data(mi1 + 415);
    const auto *mi1_416 = buffer.data(mi1 + 416);
    const auto *mi1_417 = buffer.data(mi1 + 417);
    const auto *mi1_418 = buffer.data(mi1 + 418);
    const auto *mi1_419 = buffer.data(mi1 + 419);
    const auto *mi1_420 = buffer.data(mi1 + 420);
    const auto *mi1_422 = buffer.data(mi1 + 422);
    const auto *mi1_423 = buffer.data(mi1 + 423);
    const auto *mi1_425 = buffer.data(mi1 + 425);
    const auto *mi1_426 = buffer.data(mi1 + 426);
    const auto *mi1_427 = buffer.data(mi1 + 427);
    const auto *mi1_429 = buffer.data(mi1 + 429);
    const auto *mi1_430 = buffer.data(mi1 + 430);
    const auto *mi1_431 = buffer.data(mi1 + 431);
    const auto *mi1_432 = buffer.data(mi1 + 432);
    const auto *mi1_434 = buffer.data(mi1 + 434);
    const auto *mi1_435 = buffer.data(mi1 + 435);
    const auto *mi1_441 = buffer.data(mi1 + 441);
    const auto *mi1_442 = buffer.data(mi1 + 442);
    const auto *mi1_443 = buffer.data(mi1 + 443);
    const auto *mi1_444 = buffer.data(mi1 + 444);
    const auto *mi1_445 = buffer.data(mi1 + 445);
    const auto *mi1_447 = buffer.data(mi1 + 447);
    const auto *mi1_488 = buffer.data(mi1 + 488);
    const auto *mi1_493 = buffer.data(mi1 + 493);
    const auto *mi1_494 = buffer.data(mi1 + 494);
    const auto *mi1_499 = buffer.data(mi1 + 499);
    const auto *mi1_500 = buffer.data(mi1 + 500);
    const auto *mi1_501 = buffer.data(mi1 + 501);

    const auto *mk_532 = buffer.data(mk + 532);
    const auto *mk_534 = buffer.data(mk + 534);
    const auto *mk_535 = buffer.data(mk + 535);
    const auto *mk_536 = buffer.data(mk + 536);
    const auto *mk_537 = buffer.data(mk + 537);
    const auto *mk_538 = buffer.data(mk + 538);
    const auto *mk_539 = buffer.data(mk + 539);
    const auto *mk_540 = buffer.data(mk + 540);
    const auto *mk_541 = buffer.data(mk + 541);
    const auto *mk_542 = buffer.data(mk + 542);
    const auto *mk_543 = buffer.data(mk + 543);
    const auto *mk_545 = buffer.data(mk + 545);
    const auto *mk_546 = buffer.data(mk + 546);
    const auto *mk_547 = buffer.data(mk + 547);
    const auto *mk_549 = buffer.data(mk + 549);
    const auto *mk_550 = buffer.data(mk + 550);
    const auto *mk_551 = buffer.data(mk + 551);
    const auto *mk_552 = buffer.data(mk + 552);
    const auto *mk_554 = buffer.data(mk + 554);
    const auto *mk_555 = buffer.data(mk + 555);
    const auto *mk_556 = buffer.data(mk + 556);
    const auto *mk_557 = buffer.data(mk + 557);
    const auto *mk_558 = buffer.data(mk + 558);
    const auto *mk_560 = buffer.data(mk + 560);
    const auto *mk_561 = buffer.data(mk + 561);
    const auto *mk_568 = buffer.data(mk + 568);
    const auto *mk_569 = buffer.data(mk + 569);
    const auto *mk_570 = buffer.data(mk + 570);
    const auto *mk_571 = buffer.data(mk + 571);
    const auto *mk_572 = buffer.data(mk + 572);
    const auto *mk_573 = buffer.data(mk + 573);
    const auto *mk_574 = buffer.data(mk + 574);
    const auto *mk_575 = buffer.data(mk + 575);
    const auto *mk_576 = buffer.data(mk + 576);
    const auto *mk_578 = buffer.data(mk + 578);
    const auto *mk_579 = buffer.data(mk + 579);
    const auto *mk_581 = buffer.data(mk + 581);
    const auto *mk_582 = buffer.data(mk + 582);
    const auto *mk_585 = buffer.data(mk + 585);
    const auto *mk_586 = buffer.data(mk + 586);
    const auto *mk_590 = buffer.data(mk + 590);
    const auto *mk_591 = buffer.data(mk + 591);
    const auto *mk_596 = buffer.data(mk + 596);
    const auto *mk_604 = buffer.data(mk + 604);
    const auto *mk_605 = buffer.data(mk + 605);
    const auto *mk_606 = buffer.data(mk + 606);
    const auto *mk_607 = buffer.data(mk + 607);
    const auto *mk_608 = buffer.data(mk + 608);
    const auto *mk_609 = buffer.data(mk + 609);
    const auto *mk_610 = buffer.data(mk + 610);
    const auto *mk_611 = buffer.data(mk + 611);
    const auto *mk_612 = buffer.data(mk + 612);
    const auto *mk_614 = buffer.data(mk + 614);
    const auto *mk_615 = buffer.data(mk + 615);
    const auto *mk_617 = buffer.data(mk + 617);
    const auto *mk_618 = buffer.data(mk + 618);
    const auto *mk_621 = buffer.data(mk + 621);
    const auto *mk_622 = buffer.data(mk + 622);
    const auto *mk_624 = buffer.data(mk + 624);
    const auto *mk_626 = buffer.data(mk + 626);
    const auto *mk_627 = buffer.data(mk + 627);
    const auto *mk_629 = buffer.data(mk + 629);
    const auto *mk_630 = buffer.data(mk + 630);
    const auto *mk_635 = buffer.data(mk + 635);
    const auto *mk_636 = buffer.data(mk + 636);
    const auto *mk_637 = buffer.data(mk + 637);

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_y, pb_z, lk_352, mi0_413, mi0_415, \
                         mi0_416, mi1_413, mi1_415, mi1_416, mk_532, mk_534, \
                         mk_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * mi0_413[k]
                   - f_2 * mi1_413[k]
                   + pb_y[k] * mk_532[k];

        t_667[k] = f_16 * lk_352[k]
                   + pb_z[k] * mk_532[k];

        t_668[k] = f_11 * mi0_415[k]
                   - f_12 * mi1_415[k]
                   + pb_y[k] * mk_534[k];

        t_669[k] = f_9 * mi0_416[k]
                   - f_10 * mi1_416[k]
                   + pb_y[k] * mk_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_y, mi0_417, mi0_418, mi0_419, mi1_417, \
                         mi1_418, mi1_419, mk_536, mk_537, mk_538, \
                         mk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_7 * mi0_417[k]
                   - f_8 * mi1_417[k]
                   + pb_y[k] * mk_536[k];

        t_671[k] = f_5 * mi0_418[k]
                   - f_6 * mi1_418[k]
                   + pb_y[k] * mk_537[k];

        t_672[k] = f_3 * mi0_419[k]
                   - f_4 * mi1_419[k]
                   + pb_y[k] * mk_538[k];

        t_673[k] = pb_y[k] * mk_539[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pa_x, pa_y, pb_y, pb_z, kl0_270, kl0_674, \
                         kl1_270, kl1_674, lk_360, ll_450, ll_674, \
                         mk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_31 * kl0_674[k]
                   - f_32 * kl1_674[k]
                   + pa_x[k] * ll_674[k];

        t_675[k] = f_31 * kl0_270[k]
                   - f_32 * kl1_270[k]
                   + pa_y[k] * ll_450[k];

        t_676[k] = f_17 * lk_360[k]
                   + pb_y[k] * mk_540[k];

        t_677[k] = pb_z[k] * mk_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pb_x, pb_z, lk_543, mi0_420, mi0_423, mi1_420, \
                         mi1_423, mk_541, mk_542, mk_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_16 * lk_543[k]
                   + f_11 * mi0_423[k]
                   - f_12 * mi1_423[k]
                   + pb_x[k] * mk_543[k];

        t_679[k] = pb_z[k] * mk_541[k];

        t_680[k] = f_3 * mi0_420[k]
                   - f_4 * mi1_420[k]
                   + pb_z[k] * mk_542[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, pb_x, pb_y, pb_z, lk_365, lk_546, \
                         mi0_422, mi0_426, mi1_422, mi1_426, mk_543, mk_545, \
                         mk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_16 * lk_546[k]
                   + f_9 * mi0_426[k]
                   - f_10 * mi1_426[k]
                   + pb_x[k] * mk_546[k];

        t_682[k] = pb_z[k] * mk_543[k];

        t_683[k] = f_17 * lk_365[k]
                   + pb_y[k] * mk_545[k];

        t_684[k] = f_5 * mi0_422[k]
                   - f_6 * mi1_422[k]
                   + pb_z[k] * mk_545[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, pb_x, pb_z, lk_550, mi0_423, mi0_430, mi1_423, \
                         mi1_430, mk_546, mk_547, mk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_16 * lk_550[k]
                   + f_7 * mi0_430[k]
                   - f_8 * mi1_430[k]
                   + pb_x[k] * mk_550[k];

        t_686[k] = pb_z[k] * mk_546[k];

        t_687[k] = f_3 * mi0_423[k]
                   - f_4 * mi1_423[k]
                   + pb_z[k] * mk_547[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, t_691, pb_x, pb_y, pb_z, lk_369, lk_555, \
                         mi0_425, mi0_435, mi1_425, mi1_435, mk_549, mk_550, \
                         mk_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_17 * lk_369[k]
                   + pb_y[k] * mk_549[k];

        t_689[k] = f_7 * mi0_425[k]
                   - f_8 * mi1_425[k]
                   + pb_z[k] * mk_549[k];

        t_690[k] = f_16 * lk_555[k]
                   + f_5 * mi0_435[k]
                   - f_6 * mi1_435[k]
                   + pb_x[k] * mk_555[k];

        t_691[k] = pb_z[k] * mk_550[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, pb_y, pb_z, lk_374, mi0_426, mi0_427, \
                         mi0_429, mi1_426, mi1_427, mi1_429, mk_551, mk_552, \
                         mk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_3 * mi0_426[k]
                   - f_4 * mi1_426[k]
                   + pb_z[k] * mk_551[k];

        t_693[k] = f_5 * mi0_427[k]
                   - f_6 * mi1_427[k]
                   + pb_z[k] * mk_552[k];

        t_694[k] = f_17 * lk_374[k]
                   + pb_y[k] * mk_554[k];

        t_695[k] = f_9 * mi0_429[k]
                   - f_10 * mi1_429[k]
                   + pb_z[k] * mk_554[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pb_x, pb_z, lk_561, mi0_430, mi0_441, mi1_430, \
                         mi1_441, mk_555, mk_556, mk_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_16 * lk_561[k]
                   + f_3 * mi0_441[k]
                   - f_4 * mi1_441[k]
                   + pb_x[k] * mk_561[k];

        t_697[k] = pb_z[k] * mk_555[k];

        t_698[k] = f_3 * mi0_430[k]
                   - f_4 * mi1_430[k]
                   + pb_z[k] * mk_556[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pb_y, pb_z, lk_380, mi0_431, mi0_432, \
                         mi0_434, mi1_431, mi1_432, mi1_434, mk_557, mk_558, \
                         mk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_5 * mi0_431[k]
                   - f_6 * mi1_431[k]
                   + pb_z[k] * mk_557[k];

        t_700[k] = f_7 * mi0_432[k]
                   - f_8 * mi1_432[k]
                   + pb_z[k] * mk_558[k];

        t_701[k] = f_17 * lk_380[k]
                   + pb_y[k] * mk_560[k];

        t_702[k] = f_11 * mi0_434[k]
                   - f_12 * mi1_434[k]
                   + pb_z[k] * mk_560[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, pb_x, pb_z, lk_568, lk_570, \
                         lk_571, lk_572, mk_561, mk_568, mk_570, mk_571, \
                         mk_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_16 * lk_568[k]
                   + pb_x[k] * mk_568[k];

        t_704[k] = pb_z[k] * mk_561[k];

        t_705[k] = f_16 * lk_570[k]
                   + pb_x[k] * mk_570[k];

        t_706[k] = f_16 * lk_571[k]
                   + pb_x[k] * mk_571[k];

        t_707[k] = f_16 * lk_572[k]
                   + pb_x[k] * mk_572[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pa_x, pb_x, kl0_711, kl1_711, lk_573, \
                         lk_574, lk_575, ll_711, mk_573, mk_574, \
                         mk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_16 * lk_573[k]
                   + pb_x[k] * mk_573[k];

        t_709[k] = f_16 * lk_574[k]
                   + pb_x[k] * mk_574[k];

        t_710[k] = f_16 * lk_575[k]
                   + pb_x[k] * mk_575[k];

        t_711[k] = f_29 * kl0_711[k]
                   - f_30 * kl1_711[k]
                   + pa_x[k] * ll_711[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pb_z, mi0_441, mi0_442, mi0_443, mi1_441, \
                         mi1_442, mi1_443, mk_568, mk_569, mk_570, \
                         mk_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = pb_z[k] * mk_568[k];

        t_713[k] = f_3 * mi0_441[k]
                   - f_4 * mi1_441[k]
                   + pb_z[k] * mk_569[k];

        t_714[k] = f_5 * mi0_442[k]
                   - f_6 * mi1_442[k]
                   + pb_z[k] * mk_570[k];

        t_715[k] = f_7 * mi0_443[k]
                   - f_8 * mi1_443[k]
                   + pb_z[k] * mk_571[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, pb_y, pb_z, lk_395, mi0_444, mi0_445, \
                         mi0_447, mi1_444, mi1_445, mi1_447, mk_572, mk_573, \
                         mk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_9 * mi0_444[k]
                   - f_10 * mi1_444[k]
                   + pb_z[k] * mk_572[k];

        t_717[k] = f_11 * mi0_445[k]
                   - f_12 * mi1_445[k]
                   + pb_z[k] * mk_573[k];

        t_718[k] = f_17 * lk_395[k]
                   + pb_y[k] * mk_575[k];

        t_719[k] = f_1 * mi0_447[k]
                   - f_2 * mi1_447[k]
                   + pb_z[k] * mk_575[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, t_724, pa_z, pb_y, pb_z, lk_360, lk_398, \
                         ll_450, ll_451, ll_453, mk_576, mk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_z[k] * ll_450[k];

        t_721[k] = pa_z[k] * ll_451[k];

        t_722[k] = f_13 * lk_360[k]
                   + pb_z[k] * mk_576[k];

        t_723[k] = pa_z[k] * ll_453[k];

        t_724[k] = f_16 * lk_398[k]
                   + pb_y[k] * mk_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, pa_z, pb_y, pb_z, lk_362, lk_363, lk_401, \
                         ll_455, ll_456, mk_579, mk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_14 * lk_362[k]
                   + pa_z[k] * ll_455[k];

        t_726[k] = pa_z[k] * ll_456[k];

        t_727[k] = f_13 * lk_363[k]
                   + pb_z[k] * mk_579[k];

        t_728[k] = f_16 * lk_401[k]
                   + pb_y[k] * mk_581[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, pa_z, pb_z, lk_365, lk_366, lk_367, \
                         ll_459, ll_460, ll_462, mk_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_15 * lk_365[k]
                   + pa_z[k] * ll_459[k];

        t_730[k] = pa_z[k] * ll_460[k];

        t_731[k] = f_13 * lk_366[k]
                   + pb_z[k] * mk_582[k];

        t_732[k] = f_14 * lk_367[k]
                   + pa_z[k] * ll_462[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, t_736, pa_z, pb_y, pb_z, lk_369, lk_370, lk_405, \
                         ll_464, ll_465, mk_585, mk_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_16 * lk_405[k]
                   + pb_y[k] * mk_585[k];

        t_734[k] = f_16 * lk_369[k]
                   + pa_z[k] * ll_464[k];

        t_735[k] = pa_z[k] * ll_465[k];

        t_736[k] = f_13 * lk_370[k]
                   + pb_z[k] * mk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, t_740, t_741, pa_z, pb_y, lk_371, lk_372, \
                         lk_374, lk_410, ll_467, ll_468, ll_470, ll_471, \
                         mk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_14 * lk_371[k]
                   + pa_z[k] * ll_467[k];

        t_738[k] = f_15 * lk_372[k]
                   + pa_z[k] * ll_468[k];

        t_739[k] = f_16 * lk_410[k]
                   + pb_y[k] * mk_590[k];

        t_740[k] = f_17 * lk_374[k]
                   + pa_z[k] * ll_470[k];

        t_741[k] = pa_z[k] * ll_471[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pa_z, pb_z, lk_375, lk_376, lk_377, \
                         lk_378, ll_473, ll_474, ll_475, mk_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_13 * lk_375[k]
                   + pb_z[k] * mk_591[k];

        t_743[k] = f_14 * lk_376[k]
                   + pa_z[k] * ll_473[k];

        t_744[k] = f_15 * lk_377[k]
                   + pa_z[k] * ll_474[k];

        t_745[k] = f_16 * lk_378[k]
                   + pa_z[k] * ll_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pa_z, pb_x, pb_y, lk_380, lk_416, lk_605, \
                         ll_477, ll_478, mk_596, mk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_16 * lk_416[k]
                   + pb_y[k] * mk_596[k];

        t_747[k] = f_18 * lk_380[k]
                   + pa_z[k] * ll_477[k];

        t_748[k] = pa_z[k] * ll_478[k];

        t_749[k] = f_16 * lk_605[k]
                   + pb_x[k] * mk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pb_x, lk_606, lk_607, lk_608, \
                         lk_609, lk_610, mk_606, mk_607, mk_608, mk_609, \
                         mk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_16 * lk_606[k]
                   + pb_x[k] * mk_606[k];

        t_751[k] = f_16 * lk_607[k]
                   + pb_x[k] * mk_607[k];

        t_752[k] = f_16 * lk_608[k]
                   + pb_x[k] * mk_608[k];

        t_753[k] = f_16 * lk_609[k]
                   + pb_x[k] * mk_609[k];

        t_754[k] = f_16 * lk_610[k]
                   + pb_x[k] * mk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_z, pb_x, pb_z, lk_388, lk_389, lk_611, \
                         ll_486, ll_488, mk_604, mk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_16 * lk_611[k]
                   + pb_x[k] * mk_611[k];

        t_756[k] = pa_z[k] * ll_486[k];

        t_757[k] = f_13 * lk_388[k]
                   + pb_z[k] * mk_604[k];

        t_758[k] = f_14 * lk_389[k]
                   + pa_z[k] * ll_488[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, pa_z, lk_390, lk_391, lk_392, lk_393, \
                         ll_489, ll_490, ll_491, ll_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_15 * lk_390[k]
                   + pa_z[k] * ll_489[k];

        t_760[k] = f_16 * lk_391[k]
                   + pa_z[k] * ll_490[k];

        t_761[k] = f_17 * lk_392[k]
                   + pa_z[k] * ll_491[k];

        t_762[k] = f_18 * lk_393[k]
                   + pa_z[k] * ll_492[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, t_766, pa_y, pa_z, pb_y, kl0_360, kl1_360, \
                         lk_395, lk_431, lk_432, ll_494, ll_540, mk_611, \
                         mk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_16 * lk_431[k]
                   + pb_y[k] * mk_611[k];

        t_764[k] = f_19 * lk_395[k]
                   + pa_z[k] * ll_494[k];

        t_765[k] = f_25 * kl0_360[k]
                   - f_26 * kl1_360[k]
                   + pa_y[k] * ll_540[k];

        t_766[k] = f_15 * lk_432[k]
                   + pb_y[k] * mk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pa_z, pb_y, pb_z, kl0_273, kl1_273, lk_396, \
                         lk_434, ll_498, mk_612, mk_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_14 * lk_396[k]
                   + pb_z[k] * mk_612[k];

        t_768[k] = f_20 * kl0_273[k]
                   - f_21 * kl1_273[k]
                   + pa_z[k] * ll_498[k];

        t_769[k] = f_15 * lk_434[k]
                   + pb_y[k] * mk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pa_y, pa_z, pb_z, kl0_276, kl0_365, kl1_276, \
                         kl1_365, lk_399, ll_501, ll_545, mk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_25 * kl0_365[k]
                   - f_26 * kl1_365[k]
                   + pa_y[k] * ll_545[k];

        t_771[k] = f_20 * kl0_276[k]
                   - f_21 * kl1_276[k]
                   + pa_z[k] * ll_501[k];

        t_772[k] = f_14 * lk_399[k]
                   + pb_z[k] * mk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pa_y, pa_z, pb_y, kl0_280, kl0_369, kl1_280, \
                         kl1_369, lk_437, ll_505, ll_549, mk_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_15 * lk_437[k]
                   + pb_y[k] * mk_617[k];

        t_774[k] = f_25 * kl0_369[k]
                   - f_26 * kl1_369[k]
                   + pa_y[k] * ll_549[k];

        t_775[k] = f_20 * kl0_280[k]
                   - f_21 * kl1_280[k]
                   + pa_z[k] * ll_505[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pb_x, pb_y, pb_z, lk_402, lk_441, lk_624, \
                         mi0_488, mi1_488, mk_618, mk_621, mk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_14 * lk_402[k]
                   + pb_z[k] * mk_618[k];

        t_777[k] = f_16 * lk_624[k]
                   + f_7 * mi0_488[k]
                   - f_8 * mi1_488[k]
                   + pb_x[k] * mk_624[k];

        t_778[k] = f_15 * lk_441[k]
                   + pb_y[k] * mk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pa_y, pa_z, pb_z, kl0_285, kl0_374, kl1_285, \
                         kl1_374, lk_406, ll_510, ll_554, mk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_25 * kl0_374[k]
                   - f_26 * kl1_374[k]
                   + pa_y[k] * ll_554[k];

        t_780[k] = f_20 * kl0_285[k]
                   - f_21 * kl1_285[k]
                   + pa_z[k] * ll_510[k];

        t_781[k] = f_14 * lk_406[k]
                   + pb_z[k] * mk_622[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pb_x, pb_y, lk_446, lk_629, lk_630, mi0_493, \
                         mi0_494, mi1_493, mi1_494, mk_626, mk_629, \
                         mk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_16 * lk_629[k]
                   + f_5 * mi0_493[k]
                   - f_6 * mi1_493[k]
                   + pb_x[k] * mk_629[k];

        t_783[k] = f_16 * lk_630[k]
                   + f_5 * mi0_494[k]
                   - f_6 * mi1_494[k]
                   + pb_x[k] * mk_630[k];

        t_784[k] = f_15 * lk_446[k]
                   + pb_y[k] * mk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pa_y, pa_z, pb_z, kl0_291, kl0_380, kl1_291, \
                         kl1_380, lk_411, ll_516, ll_560, mk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_25 * kl0_380[k]
                   - f_26 * kl1_380[k]
                   + pa_y[k] * ll_560[k];

        t_786[k] = f_20 * kl0_291[k]
                   - f_21 * kl1_291[k]
                   + pa_z[k] * ll_516[k];

        t_787[k] = f_14 * lk_411[k]
                   + pb_z[k] * mk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pb_x, lk_635, lk_636, lk_637, mi0_499, mi0_500, \
                         mi0_501, mi1_499, mi1_500, mi1_501, mk_635, mk_636, \
                         mk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_16 * lk_635[k]
                   + f_3 * mi0_499[k]
                   - f_4 * mi1_499[k]
                   + pb_x[k] * mk_635[k];

        t_789[k] = f_16 * lk_636[k]
                   + f_3 * mi0_500[k]
                   - f_4 * mi1_500[k]
                   + pb_x[k] * mk_636[k];

        t_790[k] = f_16 * lk_637[k]
                   + f_3 * mi0_501[k]
                   - f_4 * mi1_501[k]
                   + pb_x[k] * mk_637[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece6(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_318 = buffer.data(kl0 + 318);
    const auto *kl0_321 = buffer.data(kl0 + 321);
    const auto *kl0_325 = buffer.data(kl0 + 325);
    const auto *kl0_330 = buffer.data(kl0 + 330);
    const auto *kl0_336 = buffer.data(kl0 + 336);
    const auto *kl0_387 = buffer.data(kl0 + 387);
    const auto *kl0_405 = buffer.data(kl0 + 405);
    const auto *kl0_410 = buffer.data(kl0 + 410);
    const auto *kl0_414 = buffer.data(kl0 + 414);
    const auto *kl0_419 = buffer.data(kl0 + 419);
    const auto *kl0_425 = buffer.data(kl0 + 425);
    const auto *kl0_432 = buffer.data(kl0 + 432);
    const auto *kl0_801 = buffer.data(kl0 + 801);
    const auto *kl0_803 = buffer.data(kl0 + 803);
    const auto *kl0_804 = buffer.data(kl0 + 804);
    const auto *kl0_805 = buffer.data(kl0 + 805);
    const auto *kl0_806 = buffer.data(kl0 + 806);
    const auto *kl0_807 = buffer.data(kl0 + 807);
    const auto *kl0_809 = buffer.data(kl0 + 809);
    const auto *kl0_846 = buffer.data(kl0 + 846);
    const auto *kl0_848 = buffer.data(kl0 + 848);
    const auto *kl0_849 = buffer.data(kl0 + 849);
    const auto *kl0_850 = buffer.data(kl0 + 850);
    const auto *kl0_851 = buffer.data(kl0 + 851);
    const auto *kl0_852 = buffer.data(kl0 + 852);
    const auto *kl0_854 = buffer.data(kl0 + 854);

    const auto *kl1_318 = buffer.data(kl1 + 318);
    const auto *kl1_321 = buffer.data(kl1 + 321);
    const auto *kl1_325 = buffer.data(kl1 + 325);
    const auto *kl1_330 = buffer.data(kl1 + 330);
    const auto *kl1_336 = buffer.data(kl1 + 336);
    const auto *kl1_387 = buffer.data(kl1 + 387);
    const auto *kl1_405 = buffer.data(kl1 + 405);
    const auto *kl1_410 = buffer.data(kl1 + 410);
    const auto *kl1_414 = buffer.data(kl1 + 414);
    const auto *kl1_419 = buffer.data(kl1 + 419);
    const auto *kl1_425 = buffer.data(kl1 + 425);
    const auto *kl1_432 = buffer.data(kl1 + 432);
    const auto *kl1_801 = buffer.data(kl1 + 801);
    const auto *kl1_803 = buffer.data(kl1 + 803);
    const auto *kl1_804 = buffer.data(kl1 + 804);
    const auto *kl1_805 = buffer.data(kl1 + 805);
    const auto *kl1_806 = buffer.data(kl1 + 806);
    const auto *kl1_807 = buffer.data(kl1 + 807);
    const auto *kl1_809 = buffer.data(kl1 + 809);
    const auto *kl1_846 = buffer.data(kl1 + 846);
    const auto *kl1_848 = buffer.data(kl1 + 848);
    const auto *kl1_849 = buffer.data(kl1 + 849);
    const auto *kl1_850 = buffer.data(kl1 + 850);
    const auto *kl1_851 = buffer.data(kl1 + 851);
    const auto *kl1_852 = buffer.data(kl1 + 852);
    const auto *kl1_854 = buffer.data(kl1 + 854);

    const auto *lk_424 = buffer.data(lk + 424);
    const auto *lk_432 = buffer.data(lk + 432);
    const auto *lk_435 = buffer.data(lk + 435);
    const auto *lk_438 = buffer.data(lk + 438);
    const auto *lk_442 = buffer.data(lk + 442);
    const auto *lk_447 = buffer.data(lk + 447);
    const auto *lk_452 = buffer.data(lk + 452);
    const auto *lk_460 = buffer.data(lk + 460);
    const auto *lk_467 = buffer.data(lk + 467);
    const auto *lk_468 = buffer.data(lk + 468);
    const auto *lk_470 = buffer.data(lk + 470);
    const auto *lk_471 = buffer.data(lk + 471);
    const auto *lk_473 = buffer.data(lk + 473);
    const auto *lk_474 = buffer.data(lk + 474);
    const auto *lk_477 = buffer.data(lk + 477);
    const auto *lk_478 = buffer.data(lk + 478);
    const auto *lk_482 = buffer.data(lk + 482);
    const auto *lk_483 = buffer.data(lk + 483);
    const auto *lk_488 = buffer.data(lk + 488);
    const auto *lk_496 = buffer.data(lk + 496);
    const auto *lk_503 = buffer.data(lk + 503);
    const auto *lk_504 = buffer.data(lk + 504);
    const auto *lk_505 = buffer.data(lk + 505);
    const auto *lk_506 = buffer.data(lk + 506);
    const auto *lk_507 = buffer.data(lk + 507);
    const auto *lk_509 = buffer.data(lk + 509);
    const auto *lk_510 = buffer.data(lk + 510);
    const auto *lk_512 = buffer.data(lk + 512);
    const auto *lk_513 = buffer.data(lk + 513);
    const auto *lk_514 = buffer.data(lk + 514);
    const auto *lk_516 = buffer.data(lk + 516);
    const auto *lk_517 = buffer.data(lk + 517);
    const auto *lk_518 = buffer.data(lk + 518);
    const auto *lk_519 = buffer.data(lk + 519);
    const auto *lk_521 = buffer.data(lk + 521);
    const auto *lk_522 = buffer.data(lk + 522);
    const auto *lk_523 = buffer.data(lk + 523);
    const auto *lk_524 = buffer.data(lk + 524);
    const auto *lk_532 = buffer.data(lk + 532);
    const auto *lk_534 = buffer.data(lk + 534);
    const auto *lk_535 = buffer.data(lk + 535);
    const auto *lk_536 = buffer.data(lk + 536);
    const auto *lk_537 = buffer.data(lk + 537);
    const auto *lk_538 = buffer.data(lk + 538);
    const auto *lk_539 = buffer.data(lk + 539);
    const auto *lk_640 = buffer.data(lk + 640);
    const auto *lk_641 = buffer.data(lk + 641);
    const auto *lk_642 = buffer.data(lk + 642);
    const auto *lk_643 = buffer.data(lk + 643);
    const auto *lk_644 = buffer.data(lk + 644);
    const auto *lk_645 = buffer.data(lk + 645);
    const auto *lk_646 = buffer.data(lk + 646);
    const auto *lk_647 = buffer.data(lk + 647);
    const auto *lk_660 = buffer.data(lk + 660);
    const auto *lk_665 = buffer.data(lk + 665);
    const auto *lk_666 = buffer.data(lk + 666);
    const auto *lk_671 = buffer.data(lk + 671);
    const auto *lk_672 = buffer.data(lk + 672);
    const auto *lk_673 = buffer.data(lk + 673);
    const auto *lk_676 = buffer.data(lk + 676);
    const auto *lk_677 = buffer.data(lk + 677);
    const auto *lk_678 = buffer.data(lk + 678);
    const auto *lk_679 = buffer.data(lk + 679);
    const auto *lk_680 = buffer.data(lk + 680);
    const auto *lk_681 = buffer.data(lk + 681);
    const auto *lk_682 = buffer.data(lk + 682);
    const auto *lk_683 = buffer.data(lk + 683);
    const auto *lk_712 = buffer.data(lk + 712);
    const auto *lk_713 = buffer.data(lk + 713);
    const auto *lk_714 = buffer.data(lk + 714);
    const auto *lk_715 = buffer.data(lk + 715);
    const auto *lk_716 = buffer.data(lk + 716);
    const auto *lk_717 = buffer.data(lk + 717);
    const auto *lk_718 = buffer.data(lk + 718);
    const auto *lk_725 = buffer.data(lk + 725);
    const auto *lk_729 = buffer.data(lk + 729);
    const auto *lk_734 = buffer.data(lk + 734);

    const auto *ll_543 = buffer.data(ll + 543);
    const auto *ll_546 = buffer.data(ll + 546);
    const auto *ll_550 = buffer.data(ll + 550);
    const auto *ll_555 = buffer.data(ll + 555);
    const auto *ll_561 = buffer.data(ll + 561);
    const auto *ll_567 = buffer.data(ll + 567);
    const auto *ll_585 = buffer.data(ll + 585);
    const auto *ll_590 = buffer.data(ll + 590);
    const auto *ll_594 = buffer.data(ll + 594);
    const auto *ll_599 = buffer.data(ll + 599);
    const auto *ll_605 = buffer.data(ll + 605);
    const auto *ll_612 = buffer.data(ll + 612);
    const auto *ll_630 = buffer.data(ll + 630);
    const auto *ll_632 = buffer.data(ll + 632);
    const auto *ll_633 = buffer.data(ll + 633);
    const auto *ll_635 = buffer.data(ll + 635);
    const auto *ll_636 = buffer.data(ll + 636);
    const auto *ll_639 = buffer.data(ll + 639);
    const auto *ll_640 = buffer.data(ll + 640);
    const auto *ll_642 = buffer.data(ll + 642);
    const auto *ll_644 = buffer.data(ll + 644);
    const auto *ll_645 = buffer.data(ll + 645);
    const auto *ll_647 = buffer.data(ll + 647);
    const auto *ll_648 = buffer.data(ll + 648);
    const auto *ll_650 = buffer.data(ll + 650);
    const auto *ll_651 = buffer.data(ll + 651);
    const auto *ll_653 = buffer.data(ll + 653);
    const auto *ll_654 = buffer.data(ll + 654);
    const auto *ll_655 = buffer.data(ll + 655);
    const auto *ll_657 = buffer.data(ll + 657);
    const auto *ll_665 = buffer.data(ll + 665);
    const auto *ll_666 = buffer.data(ll + 666);
    const auto *ll_668 = buffer.data(ll + 668);
    const auto *ll_669 = buffer.data(ll + 669);
    const auto *ll_670 = buffer.data(ll + 670);
    const auto *ll_671 = buffer.data(ll + 671);
    const auto *ll_672 = buffer.data(ll + 672);
    const auto *ll_674 = buffer.data(ll + 674);
    const auto *ll_801 = buffer.data(ll + 801);
    const auto *ll_803 = buffer.data(ll + 803);
    const auto *ll_804 = buffer.data(ll + 804);
    const auto *ll_805 = buffer.data(ll + 805);
    const auto *ll_806 = buffer.data(ll + 806);
    const auto *ll_807 = buffer.data(ll + 807);
    const auto *ll_809 = buffer.data(ll + 809);
    const auto *ll_846 = buffer.data(ll + 846);
    const auto *ll_848 = buffer.data(ll + 848);
    const auto *ll_849 = buffer.data(ll + 849);
    const auto *ll_850 = buffer.data(ll + 850);
    const auto *ll_851 = buffer.data(ll + 851);
    const auto *ll_852 = buffer.data(ll + 852);
    const auto *ll_854 = buffer.data(ll + 854);

    const auto *mi0_516 = buffer.data(mi0 + 516);
    const auto *mi0_521 = buffer.data(mi0 + 521);
    const auto *mi0_522 = buffer.data(mi0 + 522);
    const auto *mi0_527 = buffer.data(mi0 + 527);
    const auto *mi0_528 = buffer.data(mi0 + 528);
    const auto *mi0_529 = buffer.data(mi0 + 529);
    const auto *mi0_560 = buffer.data(mi0 + 560);
    const auto *mi0_561 = buffer.data(mi0 + 561);
    const auto *mi0_563 = buffer.data(mi0 + 563);
    const auto *mi0_565 = buffer.data(mi0 + 565);
    const auto *mi0_569 = buffer.data(mi0 + 569);
    const auto *mi0_574 = buffer.data(mi0 + 574);

    const auto *mi1_516 = buffer.data(mi1 + 516);
    const auto *mi1_521 = buffer.data(mi1 + 521);
    const auto *mi1_522 = buffer.data(mi1 + 522);
    const auto *mi1_527 = buffer.data(mi1 + 527);
    const auto *mi1_528 = buffer.data(mi1 + 528);
    const auto *mi1_529 = buffer.data(mi1 + 529);
    const auto *mi1_560 = buffer.data(mi1 + 560);
    const auto *mi1_561 = buffer.data(mi1 + 561);
    const auto *mi1_563 = buffer.data(mi1 + 563);
    const auto *mi1_565 = buffer.data(mi1 + 565);
    const auto *mi1_569 = buffer.data(mi1 + 569);
    const auto *mi1_574 = buffer.data(mi1 + 574);

    const auto *mk_632 = buffer.data(mk + 632);
    const auto *mk_640 = buffer.data(mk + 640);
    const auto *mk_641 = buffer.data(mk + 641);
    const auto *mk_642 = buffer.data(mk + 642);
    const auto *mk_643 = buffer.data(mk + 643);
    const auto *mk_644 = buffer.data(mk + 644);
    const auto *mk_645 = buffer.data(mk + 645);
    const auto *mk_646 = buffer.data(mk + 646);
    const auto *mk_647 = buffer.data(mk + 647);
    const auto *mk_648 = buffer.data(mk + 648);
    const auto *mk_650 = buffer.data(mk + 650);
    const auto *mk_651 = buffer.data(mk + 651);
    const auto *mk_653 = buffer.data(mk + 653);
    const auto *mk_654 = buffer.data(mk + 654);
    const auto *mk_657 = buffer.data(mk + 657);
    const auto *mk_658 = buffer.data(mk + 658);
    const auto *mk_660 = buffer.data(mk + 660);
    const auto *mk_662 = buffer.data(mk + 662);
    const auto *mk_663 = buffer.data(mk + 663);
    const auto *mk_665 = buffer.data(mk + 665);
    const auto *mk_666 = buffer.data(mk + 666);
    const auto *mk_668 = buffer.data(mk + 668);
    const auto *mk_671 = buffer.data(mk + 671);
    const auto *mk_672 = buffer.data(mk + 672);
    const auto *mk_673 = buffer.data(mk + 673);
    const auto *mk_676 = buffer.data(mk + 676);
    const auto *mk_677 = buffer.data(mk + 677);
    const auto *mk_678 = buffer.data(mk + 678);
    const auto *mk_679 = buffer.data(mk + 679);
    const auto *mk_680 = buffer.data(mk + 680);
    const auto *mk_681 = buffer.data(mk + 681);
    const auto *mk_682 = buffer.data(mk + 682);
    const auto *mk_683 = buffer.data(mk + 683);
    const auto *mk_684 = buffer.data(mk + 684);
    const auto *mk_686 = buffer.data(mk + 686);
    const auto *mk_687 = buffer.data(mk + 687);
    const auto *mk_689 = buffer.data(mk + 689);
    const auto *mk_690 = buffer.data(mk + 690);
    const auto *mk_693 = buffer.data(mk + 693);
    const auto *mk_694 = buffer.data(mk + 694);
    const auto *mk_698 = buffer.data(mk + 698);
    const auto *mk_699 = buffer.data(mk + 699);
    const auto *mk_704 = buffer.data(mk + 704);
    const auto *mk_712 = buffer.data(mk + 712);
    const auto *mk_713 = buffer.data(mk + 713);
    const auto *mk_714 = buffer.data(mk + 714);
    const auto *mk_715 = buffer.data(mk + 715);
    const auto *mk_716 = buffer.data(mk + 716);
    const auto *mk_717 = buffer.data(mk + 717);
    const auto *mk_718 = buffer.data(mk + 718);
    const auto *mk_719 = buffer.data(mk + 719);
    const auto *mk_720 = buffer.data(mk + 720);
    const auto *mk_721 = buffer.data(mk + 721);
    const auto *mk_722 = buffer.data(mk + 722);
    const auto *mk_723 = buffer.data(mk + 723);
    const auto *mk_725 = buffer.data(mk + 725);
    const auto *mk_726 = buffer.data(mk + 726);
    const auto *mk_728 = buffer.data(mk + 728);
    const auto *mk_729 = buffer.data(mk + 729);
    const auto *mk_734 = buffer.data(mk + 734);

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pa_y, pb_x, pb_y, kl0_387, kl1_387, \
                         lk_452, lk_640, lk_641, ll_567, mk_632, mk_640, \
                         mk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_15 * lk_452[k]
                   + pb_y[k] * mk_632[k];

        t_792[k] = f_25 * kl0_387[k]
                   - f_26 * kl1_387[k]
                   + pa_y[k] * ll_567[k];

        t_793[k] = f_16 * lk_640[k]
                   + pb_x[k] * mk_640[k];

        t_794[k] = f_16 * lk_641[k]
                   + pb_x[k] * mk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pb_x, lk_642, lk_643, lk_644, \
                         lk_645, lk_646, mk_642, mk_643, mk_644, mk_645, \
                         mk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_16 * lk_642[k]
                   + pb_x[k] * mk_642[k];

        t_796[k] = f_16 * lk_643[k]
                   + pb_x[k] * mk_643[k];

        t_797[k] = f_16 * lk_644[k]
                   + pb_x[k] * mk_644[k];

        t_798[k] = f_16 * lk_645[k]
                   + pb_x[k] * mk_645[k];

        t_799[k] = f_16 * lk_646[k]
                   + pb_x[k] * mk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pa_x, pb_x, pb_z, kl0_801, kl1_801, lk_424, \
                         lk_647, ll_801, mk_640, mk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_16 * lk_647[k]
                   + pb_x[k] * mk_647[k];

        t_801[k] = f_29 * kl0_801[k]
                   - f_30 * kl1_801[k]
                   + pa_x[k] * ll_801[k];

        t_802[k] = f_14 * lk_424[k]
                   + pb_z[k] * mk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pa_x, kl0_803, kl0_804, kl0_805, kl1_803, \
                         kl1_804, kl1_805, ll_803, ll_804, ll_805 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_29 * kl0_803[k]
                   - f_30 * kl1_803[k]
                   + pa_x[k] * ll_803[k];

        t_804[k] = f_29 * kl0_804[k]
                   - f_30 * kl1_804[k]
                   + pa_x[k] * ll_804[k];

        t_805[k] = f_29 * kl0_805[k]
                   - f_30 * kl1_805[k]
                   + pa_x[k] * ll_805[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pa_x, pb_y, kl0_806, kl0_807, kl1_806, kl1_807, \
                         lk_467, ll_806, ll_807, mk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_29 * kl0_806[k]
                   - f_30 * kl1_806[k]
                   + pa_x[k] * ll_806[k];

        t_807[k] = f_29 * kl0_807[k]
                   - f_30 * kl1_807[k]
                   + pa_x[k] * ll_807[k];

        t_808[k] = f_15 * lk_467[k]
                   + pb_y[k] * mk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pa_x, pa_y, pb_y, kl0_405, kl0_809, kl1_405, \
                         kl1_809, lk_468, ll_585, ll_809, mk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_29 * kl0_809[k]
                   - f_30 * kl1_809[k]
                   + pa_x[k] * ll_809[k];

        t_810[k] = f_20 * kl0_405[k]
                   - f_21 * kl1_405[k]
                   + pa_y[k] * ll_585[k];

        t_811[k] = f_14 * lk_468[k]
                   + pb_y[k] * mk_648[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pa_z, pb_y, pb_z, kl0_318, kl1_318, lk_432, \
                         lk_470, ll_543, mk_648, mk_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_15 * lk_432[k]
                   + pb_z[k] * mk_648[k];

        t_813[k] = f_25 * kl0_318[k]
                   - f_26 * kl1_318[k]
                   + pa_z[k] * ll_543[k];

        t_814[k] = f_14 * lk_470[k]
                   + pb_y[k] * mk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pa_y, pa_z, pb_z, kl0_321, kl0_410, kl1_321, \
                         kl1_410, lk_435, ll_546, ll_590, mk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_20 * kl0_410[k]
                   - f_21 * kl1_410[k]
                   + pa_y[k] * ll_590[k];

        t_816[k] = f_25 * kl0_321[k]
                   - f_26 * kl1_321[k]
                   + pa_z[k] * ll_546[k];

        t_817[k] = f_15 * lk_435[k]
                   + pb_z[k] * mk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pa_y, pa_z, pb_y, kl0_325, kl0_414, kl1_325, \
                         kl1_414, lk_473, ll_550, ll_594, mk_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_14 * lk_473[k]
                   + pb_y[k] * mk_653[k];

        t_819[k] = f_20 * kl0_414[k]
                   - f_21 * kl1_414[k]
                   + pa_y[k] * ll_594[k];

        t_820[k] = f_25 * kl0_325[k]
                   - f_26 * kl1_325[k]
                   + pa_z[k] * ll_550[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pb_x, pb_y, pb_z, lk_438, lk_477, lk_660, \
                         mi0_516, mi1_516, mk_654, mk_657, mk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_15 * lk_438[k]
                   + pb_z[k] * mk_654[k];

        t_822[k] = f_16 * lk_660[k]
                   + f_7 * mi0_516[k]
                   - f_8 * mi1_516[k]
                   + pb_x[k] * mk_660[k];

        t_823[k] = f_14 * lk_477[k]
                   + pb_y[k] * mk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pa_y, pa_z, pb_z, kl0_330, kl0_419, kl1_330, \
                         kl1_419, lk_442, ll_555, ll_599, mk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_20 * kl0_419[k]
                   - f_21 * kl1_419[k]
                   + pa_y[k] * ll_599[k];

        t_825[k] = f_25 * kl0_330[k]
                   - f_26 * kl1_330[k]
                   + pa_z[k] * ll_555[k];

        t_826[k] = f_15 * lk_442[k]
                   + pb_z[k] * mk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pb_x, pb_y, lk_482, lk_665, lk_666, mi0_521, \
                         mi0_522, mi1_521, mi1_522, mk_662, mk_665, \
                         mk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_16 * lk_665[k]
                   + f_5 * mi0_521[k]
                   - f_6 * mi1_521[k]
                   + pb_x[k] * mk_665[k];

        t_828[k] = f_16 * lk_666[k]
                   + f_5 * mi0_522[k]
                   - f_6 * mi1_522[k]
                   + pb_x[k] * mk_666[k];

        t_829[k] = f_14 * lk_482[k]
                   + pb_y[k] * mk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pa_y, pa_z, pb_z, kl0_336, kl0_425, kl1_336, \
                         kl1_425, lk_447, ll_561, ll_605, mk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_20 * kl0_425[k]
                   - f_21 * kl1_425[k]
                   + pa_y[k] * ll_605[k];

        t_831[k] = f_25 * kl0_336[k]
                   - f_26 * kl1_336[k]
                   + pa_z[k] * ll_561[k];

        t_832[k] = f_15 * lk_447[k]
                   + pb_z[k] * mk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pb_x, lk_671, lk_672, lk_673, mi0_527, mi0_528, \
                         mi0_529, mi1_527, mi1_528, mi1_529, mk_671, mk_672, \
                         mk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_16 * lk_671[k]
                   + f_3 * mi0_527[k]
                   - f_4 * mi1_527[k]
                   + pb_x[k] * mk_671[k];

        t_834[k] = f_16 * lk_672[k]
                   + f_3 * mi0_528[k]
                   - f_4 * mi1_528[k]
                   + pb_x[k] * mk_672[k];

        t_835[k] = f_16 * lk_673[k]
                   + f_3 * mi0_529[k]
                   - f_4 * mi1_529[k]
                   + pb_x[k] * mk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pa_y, pb_x, pb_y, kl0_432, kl1_432, \
                         lk_488, lk_676, lk_677, ll_612, mk_668, mk_676, \
                         mk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_14 * lk_488[k]
                   + pb_y[k] * mk_668[k];

        t_837[k] = f_20 * kl0_432[k]
                   - f_21 * kl1_432[k]
                   + pa_y[k] * ll_612[k];

        t_838[k] = f_16 * lk_676[k]
                   + pb_x[k] * mk_676[k];

        t_839[k] = f_16 * lk_677[k]
                   + pb_x[k] * mk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pb_x, lk_678, lk_679, lk_680, \
                         lk_681, lk_682, mk_678, mk_679, mk_680, mk_681, \
                         mk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_16 * lk_678[k]
                   + pb_x[k] * mk_678[k];

        t_841[k] = f_16 * lk_679[k]
                   + pb_x[k] * mk_679[k];

        t_842[k] = f_16 * lk_680[k]
                   + pb_x[k] * mk_680[k];

        t_843[k] = f_16 * lk_681[k]
                   + pb_x[k] * mk_681[k];

        t_844[k] = f_16 * lk_682[k]
                   + pb_x[k] * mk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pa_x, pb_x, pb_z, kl0_846, kl1_846, lk_460, \
                         lk_683, ll_846, mk_676, mk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_16 * lk_683[k]
                   + pb_x[k] * mk_683[k];

        t_846[k] = f_29 * kl0_846[k]
                   - f_30 * kl1_846[k]
                   + pa_x[k] * ll_846[k];

        t_847[k] = f_15 * lk_460[k]
                   + pb_z[k] * mk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pa_x, kl0_848, kl0_849, kl0_850, kl1_848, \
                         kl1_849, kl1_850, ll_848, ll_849, ll_850 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_29 * kl0_848[k]
                   - f_30 * kl1_848[k]
                   + pa_x[k] * ll_848[k];

        t_849[k] = f_29 * kl0_849[k]
                   - f_30 * kl1_849[k]
                   + pa_x[k] * ll_849[k];

        t_850[k] = f_29 * kl0_850[k]
                   - f_30 * kl1_850[k]
                   + pa_x[k] * ll_850[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pa_x, pb_y, kl0_851, kl0_852, kl1_851, kl1_852, \
                         lk_503, ll_851, ll_852, mk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_29 * kl0_851[k]
                   - f_30 * kl1_851[k]
                   + pa_x[k] * ll_851[k];

        t_852[k] = f_29 * kl0_852[k]
                   - f_30 * kl1_852[k]
                   + pa_x[k] * ll_852[k];

        t_853[k] = f_14 * lk_503[k]
                   + pb_y[k] * mk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_x, pa_y, pb_y, kl0_854, kl1_854, \
                         lk_504, ll_630, ll_632, ll_854, mk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_29 * kl0_854[k]
                   - f_30 * kl1_854[k]
                   + pa_x[k] * ll_854[k];

        t_855[k] = pa_y[k] * ll_630[k];

        t_856[k] = f_13 * lk_504[k]
                   + pb_y[k] * mk_684[k];

        t_857[k] = pa_y[k] * ll_632[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pa_y, pb_y, lk_505, lk_506, lk_507, \
                         ll_633, ll_635, ll_636, mk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_14 * lk_505[k]
                   + pa_y[k] * ll_633[k];

        t_859[k] = f_13 * lk_506[k]
                   + pb_y[k] * mk_686[k];

        t_860[k] = pa_y[k] * ll_635[k];

        t_861[k] = f_15 * lk_507[k]
                   + pa_y[k] * ll_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pa_y, pb_y, pb_z, lk_471, lk_509, lk_510, \
                         ll_639, ll_640, mk_687, mk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_16 * lk_471[k]
                   + pb_z[k] * mk_687[k];

        t_863[k] = f_13 * lk_509[k]
                   + pb_y[k] * mk_689[k];

        t_864[k] = pa_y[k] * ll_639[k];

        t_865[k] = f_16 * lk_510[k]
                   + pa_y[k] * ll_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pa_y, pb_y, pb_z, lk_474, lk_512, lk_513, \
                         ll_642, ll_644, mk_690, mk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_16 * lk_474[k]
                   + pb_z[k] * mk_690[k];

        t_867[k] = f_14 * lk_512[k]
                   + pa_y[k] * ll_642[k];

        t_868[k] = f_13 * lk_513[k]
                   + pb_y[k] * mk_693[k];

        t_869[k] = pa_y[k] * ll_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, pa_y, pb_z, lk_478, lk_514, lk_516, \
                         lk_517, ll_645, ll_647, ll_648, mk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_17 * lk_514[k]
                   + pa_y[k] * ll_645[k];

        t_871[k] = f_16 * lk_478[k]
                   + pb_z[k] * mk_694[k];

        t_872[k] = f_15 * lk_516[k]
                   + pa_y[k] * ll_647[k];

        t_873[k] = f_14 * lk_517[k]
                   + pa_y[k] * ll_648[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, pa_y, pb_y, pb_z, lk_483, lk_518, lk_519, \
                         ll_650, ll_651, mk_698, mk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_13 * lk_518[k]
                   + pb_y[k] * mk_698[k];

        t_875[k] = pa_y[k] * ll_650[k];

        t_876[k] = f_18 * lk_519[k]
                   + pa_y[k] * ll_651[k];

        t_877[k] = f_16 * lk_483[k]
                   + pb_z[k] * mk_699[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, t_881, t_882, pa_y, pb_y, lk_521, lk_522, \
                         lk_523, lk_524, ll_653, ll_654, ll_655, ll_657, \
                         mk_704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_16 * lk_521[k]
                   + pa_y[k] * ll_653[k];

        t_879[k] = f_15 * lk_522[k]
                   + pa_y[k] * ll_654[k];

        t_880[k] = f_14 * lk_523[k]
                   + pa_y[k] * ll_655[k];

        t_881[k] = f_13 * lk_524[k]
                   + pb_y[k] * mk_704[k];

        t_882[k] = pa_y[k] * ll_657[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, t_887, pb_x, lk_712, lk_713, lk_714, \
                         lk_715, lk_716, mk_712, mk_713, mk_714, mk_715, \
                         mk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_16 * lk_712[k]
                   + pb_x[k] * mk_712[k];

        t_884[k] = f_16 * lk_713[k]
                   + pb_x[k] * mk_713[k];

        t_885[k] = f_16 * lk_714[k]
                   + pb_x[k] * mk_714[k];

        t_886[k] = f_16 * lk_715[k]
                   + pb_x[k] * mk_715[k];

        t_887[k] = f_16 * lk_716[k]
                   + pb_x[k] * mk_716[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, t_891, pa_y, pb_x, lk_532, lk_717, lk_718, \
                         ll_665, ll_666, mk_717, mk_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = f_16 * lk_717[k]
                   + pb_x[k] * mk_717[k];

        t_889[k] = f_16 * lk_718[k]
                   + pb_x[k] * mk_718[k];

        t_890[k] = pa_y[k] * ll_665[k];

        t_891[k] = f_19 * lk_532[k]
                   + pa_y[k] * ll_666[k];
    }

#pragma omp simd aligned(t_892, t_893, t_894, t_895, pa_y, pb_z, lk_496, lk_534, lk_535, \
                         lk_536, ll_668, ll_669, ll_670, mk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_892[k] = f_16 * lk_496[k]
                   + pb_z[k] * mk_712[k];

        t_893[k] = f_18 * lk_534[k]
                   + pa_y[k] * ll_668[k];

        t_894[k] = f_17 * lk_535[k]
                   + pa_y[k] * ll_669[k];

        t_895[k] = f_16 * lk_536[k]
                   + pa_y[k] * ll_670[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, t_899, pa_y, pb_y, lk_537, lk_538, lk_539, \
                         ll_671, ll_672, ll_674, mk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * lk_537[k]
                   + pa_y[k] * ll_671[k];

        t_897[k] = f_14 * lk_538[k]
                   + pa_y[k] * ll_672[k];

        t_898[k] = f_13 * lk_539[k]
                   + pb_y[k] * mk_719[k];

        t_899[k] = pa_y[k] * ll_674[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, pa_z, pb_y, pb_z, kl0_405, kl1_405, \
                         lk_504, ll_630, mi0_560, mi1_560, mk_720, \
                         mk_721 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_31 * kl0_405[k]
                   - f_32 * kl1_405[k]
                   + pa_z[k] * ll_630[k];

        t_901[k] = pb_y[k] * mk_720[k];

        t_902[k] = f_17 * lk_504[k]
                   + pb_z[k] * mk_720[k];

        t_903[k] = f_3 * mi0_560[k]
                   - f_4 * mi1_560[k]
                   + pb_y[k] * mk_721[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, t_907, pb_x, pb_y, pb_z, lk_507, lk_725, \
                         mi0_561, mi0_565, mi1_561, mi1_565, mk_722, mk_723, \
                         mk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = pb_y[k] * mk_722[k];

        t_905[k] = f_16 * lk_725[k]
                   + f_11 * mi0_565[k]
                   - f_12 * mi1_565[k]
                   + pb_x[k] * mk_725[k];

        t_906[k] = f_5 * mi0_561[k]
                   - f_6 * mi1_561[k]
                   + pb_y[k] * mk_723[k];

        t_907[k] = f_17 * lk_507[k]
                   + pb_z[k] * mk_723[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, t_911, pb_x, pb_y, pb_z, lk_510, lk_729, \
                         mi0_563, mi0_569, mi1_563, mi1_569, mk_725, mk_726, \
                         mk_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = pb_y[k] * mk_725[k];

        t_909[k] = f_16 * lk_729[k]
                   + f_9 * mi0_569[k]
                   - f_10 * mi1_569[k]
                   + pb_x[k] * mk_729[k];

        t_910[k] = f_7 * mi0_563[k]
                   - f_8 * mi1_563[k]
                   + pb_y[k] * mk_726[k];

        t_911[k] = f_17 * lk_510[k]
                   + pb_z[k] * mk_726[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pb_x, pb_y, lk_734, mi0_565, mi0_574, mi1_565, \
                         mi1_574, mk_728, mk_729, mk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_3 * mi0_565[k]
                   - f_4 * mi1_565[k]
                   + pb_y[k] * mk_728[k];

        t_913[k] = pb_y[k] * mk_729[k];

        t_914[k] = f_16 * lk_734[k]
                   + f_7 * mi0_574[k]
                   - f_8 * mi1_574[k]
                   + pb_x[k] * mk_734[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece7(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_27 = 2.5 / alpha;
    const auto f_28 = 2.5 * beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_450 = buffer.data(kl0 + 450);
    const auto *kl0_453 = buffer.data(kl0 + 453);
    const auto *kl0_456 = buffer.data(kl0 + 456);
    const auto *kl0_540 = buffer.data(kl0 + 540);
    const auto *kl0_545 = buffer.data(kl0 + 545);
    const auto *kl0_944 = buffer.data(kl0 + 944);
    const auto *kl0_981 = buffer.data(kl0 + 981);

    const auto *kl1_450 = buffer.data(kl1 + 450);
    const auto *kl1_453 = buffer.data(kl1 + 453);
    const auto *kl1_456 = buffer.data(kl1 + 456);
    const auto *kl1_540 = buffer.data(kl1 + 540);
    const auto *kl1_545 = buffer.data(kl1 + 545);
    const auto *kl1_944 = buffer.data(kl1 + 944);
    const auto *kl1_981 = buffer.data(kl1 + 981);

    const auto *lk_514 = buffer.data(lk + 514);
    const auto *lk_519 = buffer.data(lk + 519);
    const auto *lk_532 = buffer.data(lk + 532);
    const auto *lk_540 = buffer.data(lk + 540);
    const auto *lk_542 = buffer.data(lk + 542);
    const auto *lk_543 = buffer.data(lk + 543);
    const auto *lk_545 = buffer.data(lk + 545);
    const auto *lk_546 = buffer.data(lk + 546);
    const auto *lk_547 = buffer.data(lk + 547);
    const auto *lk_549 = buffer.data(lk + 549);
    const auto *lk_550 = buffer.data(lk + 550);
    const auto *lk_551 = buffer.data(lk + 551);
    const auto *lk_552 = buffer.data(lk + 552);
    const auto *lk_554 = buffer.data(lk + 554);
    const auto *lk_555 = buffer.data(lk + 555);
    const auto *lk_556 = buffer.data(lk + 556);
    const auto *lk_557 = buffer.data(lk + 557);
    const auto *lk_558 = buffer.data(lk + 558);
    const auto *lk_560 = buffer.data(lk + 560);
    const auto *lk_568 = buffer.data(lk + 568);
    const auto *lk_569 = buffer.data(lk + 569);
    const auto *lk_570 = buffer.data(lk + 570);
    const auto *lk_571 = buffer.data(lk + 571);
    const auto *lk_572 = buffer.data(lk + 572);
    const auto *lk_573 = buffer.data(lk + 573);
    const auto *lk_575 = buffer.data(lk + 575);
    const auto *lk_576 = buffer.data(lk + 576);
    const auto *lk_578 = buffer.data(lk + 578);
    const auto *lk_579 = buffer.data(lk + 579);
    const auto *lk_581 = buffer.data(lk + 581);
    const auto *lk_585 = buffer.data(lk + 585);
    const auto *lk_590 = buffer.data(lk + 590);
    const auto *lk_596 = buffer.data(lk + 596);
    const auto *lk_611 = buffer.data(lk + 611);
    const auto *lk_612 = buffer.data(lk + 612);
    const auto *lk_614 = buffer.data(lk + 614);
    const auto *lk_740 = buffer.data(lk + 740);
    const auto *lk_747 = buffer.data(lk + 747);
    const auto *lk_748 = buffer.data(lk + 748);
    const auto *lk_749 = buffer.data(lk + 749);
    const auto *lk_750 = buffer.data(lk + 750);
    const auto *lk_751 = buffer.data(lk + 751);
    const auto *lk_752 = buffer.data(lk + 752);
    const auto *lk_753 = buffer.data(lk + 753);
    const auto *lk_755 = buffer.data(lk + 755);
    const auto *lk_759 = buffer.data(lk + 759);
    const auto *lk_762 = buffer.data(lk + 762);
    const auto *lk_766 = buffer.data(lk + 766);
    const auto *lk_771 = buffer.data(lk + 771);
    const auto *lk_777 = buffer.data(lk + 777);
    const auto *lk_784 = buffer.data(lk + 784);
    const auto *lk_786 = buffer.data(lk + 786);
    const auto *lk_787 = buffer.data(lk + 787);
    const auto *lk_788 = buffer.data(lk + 788);
    const auto *lk_789 = buffer.data(lk + 789);
    const auto *lk_790 = buffer.data(lk + 790);
    const auto *lk_791 = buffer.data(lk + 791);
    const auto *lk_821 = buffer.data(lk + 821);
    const auto *lk_822 = buffer.data(lk + 822);
    const auto *lk_823 = buffer.data(lk + 823);
    const auto *lk_824 = buffer.data(lk + 824);
    const auto *lk_825 = buffer.data(lk + 825);
    const auto *lk_826 = buffer.data(lk + 826);
    const auto *lk_827 = buffer.data(lk + 827);

    const auto *ll_675 = buffer.data(ll + 675);
    const auto *ll_676 = buffer.data(ll + 676);
    const auto *ll_678 = buffer.data(ll + 678);
    const auto *ll_680 = buffer.data(ll + 680);
    const auto *ll_681 = buffer.data(ll + 681);
    const auto *ll_684 = buffer.data(ll + 684);
    const auto *ll_685 = buffer.data(ll + 685);
    const auto *ll_687 = buffer.data(ll + 687);
    const auto *ll_689 = buffer.data(ll + 689);
    const auto *ll_690 = buffer.data(ll + 690);
    const auto *ll_692 = buffer.data(ll + 692);
    const auto *ll_693 = buffer.data(ll + 693);
    const auto *ll_695 = buffer.data(ll + 695);
    const auto *ll_696 = buffer.data(ll + 696);
    const auto *ll_698 = buffer.data(ll + 698);
    const auto *ll_699 = buffer.data(ll + 699);
    const auto *ll_700 = buffer.data(ll + 700);
    const auto *ll_702 = buffer.data(ll + 702);
    const auto *ll_703 = buffer.data(ll + 703);
    const auto *ll_711 = buffer.data(ll + 711);
    const auto *ll_713 = buffer.data(ll + 713);
    const auto *ll_714 = buffer.data(ll + 714);
    const auto *ll_715 = buffer.data(ll + 715);
    const auto *ll_716 = buffer.data(ll + 716);
    const auto *ll_717 = buffer.data(ll + 717);
    const auto *ll_719 = buffer.data(ll + 719);
    const auto *ll_723 = buffer.data(ll + 723);
    const auto *ll_726 = buffer.data(ll + 726);
    const auto *ll_765 = buffer.data(ll + 765);
    const auto *ll_770 = buffer.data(ll + 770);
    const auto *ll_944 = buffer.data(ll + 944);
    const auto *ll_981 = buffer.data(ll + 981);

    const auto *mi0_566 = buffer.data(mi0 + 566);
    const auto *mi0_568 = buffer.data(mi0 + 568);
    const auto *mi0_569 = buffer.data(mi0 + 569);
    const auto *mi0_570 = buffer.data(mi0 + 570);
    const auto *mi0_572 = buffer.data(mi0 + 572);
    const auto *mi0_573 = buffer.data(mi0 + 573);
    const auto *mi0_574 = buffer.data(mi0 + 574);
    const auto *mi0_580 = buffer.data(mi0 + 580);
    const auto *mi0_581 = buffer.data(mi0 + 581);
    const auto *mi0_583 = buffer.data(mi0 + 583);
    const auto *mi0_584 = buffer.data(mi0 + 584);
    const auto *mi0_585 = buffer.data(mi0 + 585);
    const auto *mi0_586 = buffer.data(mi0 + 586);
    const auto *mi0_587 = buffer.data(mi0 + 587);
    const auto *mi0_588 = buffer.data(mi0 + 588);
    const auto *mi0_590 = buffer.data(mi0 + 590);
    const auto *mi0_591 = buffer.data(mi0 + 591);
    const auto *mi0_593 = buffer.data(mi0 + 593);
    const auto *mi0_594 = buffer.data(mi0 + 594);
    const auto *mi0_595 = buffer.data(mi0 + 595);
    const auto *mi0_597 = buffer.data(mi0 + 597);
    const auto *mi0_598 = buffer.data(mi0 + 598);
    const auto *mi0_599 = buffer.data(mi0 + 599);
    const auto *mi0_600 = buffer.data(mi0 + 600);
    const auto *mi0_602 = buffer.data(mi0 + 602);
    const auto *mi0_603 = buffer.data(mi0 + 603);
    const auto *mi0_609 = buffer.data(mi0 + 609);
    const auto *mi0_610 = buffer.data(mi0 + 610);
    const auto *mi0_611 = buffer.data(mi0 + 611);
    const auto *mi0_612 = buffer.data(mi0 + 612);
    const auto *mi0_613 = buffer.data(mi0 + 613);
    const auto *mi0_615 = buffer.data(mi0 + 615);

    const auto *mi1_566 = buffer.data(mi1 + 566);
    const auto *mi1_568 = buffer.data(mi1 + 568);
    const auto *mi1_569 = buffer.data(mi1 + 569);
    const auto *mi1_570 = buffer.data(mi1 + 570);
    const auto *mi1_572 = buffer.data(mi1 + 572);
    const auto *mi1_573 = buffer.data(mi1 + 573);
    const auto *mi1_574 = buffer.data(mi1 + 574);
    const auto *mi1_580 = buffer.data(mi1 + 580);
    const auto *mi1_581 = buffer.data(mi1 + 581);
    const auto *mi1_583 = buffer.data(mi1 + 583);
    const auto *mi1_584 = buffer.data(mi1 + 584);
    const auto *mi1_585 = buffer.data(mi1 + 585);
    const auto *mi1_586 = buffer.data(mi1 + 586);
    const auto *mi1_587 = buffer.data(mi1 + 587);
    const auto *mi1_588 = buffer.data(mi1 + 588);
    const auto *mi1_590 = buffer.data(mi1 + 590);
    const auto *mi1_591 = buffer.data(mi1 + 591);
    const auto *mi1_593 = buffer.data(mi1 + 593);
    const auto *mi1_594 = buffer.data(mi1 + 594);
    const auto *mi1_595 = buffer.data(mi1 + 595);
    const auto *mi1_597 = buffer.data(mi1 + 597);
    const auto *mi1_598 = buffer.data(mi1 + 598);
    const auto *mi1_599 = buffer.data(mi1 + 599);
    const auto *mi1_600 = buffer.data(mi1 + 600);
    const auto *mi1_602 = buffer.data(mi1 + 602);
    const auto *mi1_603 = buffer.data(mi1 + 603);
    const auto *mi1_609 = buffer.data(mi1 + 609);
    const auto *mi1_610 = buffer.data(mi1 + 610);
    const auto *mi1_611 = buffer.data(mi1 + 611);
    const auto *mi1_612 = buffer.data(mi1 + 612);
    const auto *mi1_613 = buffer.data(mi1 + 613);
    const auto *mi1_615 = buffer.data(mi1 + 615);

    const auto *mk_730 = buffer.data(mk + 730);
    const auto *mk_732 = buffer.data(mk + 732);
    const auto *mk_733 = buffer.data(mk + 733);
    const auto *mk_734 = buffer.data(mk + 734);
    const auto *mk_735 = buffer.data(mk + 735);
    const auto *mk_737 = buffer.data(mk + 737);
    const auto *mk_738 = buffer.data(mk + 738);
    const auto *mk_739 = buffer.data(mk + 739);
    const auto *mk_740 = buffer.data(mk + 740);
    const auto *mk_747 = buffer.data(mk + 747);
    const auto *mk_748 = buffer.data(mk + 748);
    const auto *mk_749 = buffer.data(mk + 749);
    const auto *mk_750 = buffer.data(mk + 750);
    const auto *mk_751 = buffer.data(mk + 751);
    const auto *mk_752 = buffer.data(mk + 752);
    const auto *mk_753 = buffer.data(mk + 753);
    const auto *mk_754 = buffer.data(mk + 754);
    const auto *mk_755 = buffer.data(mk + 755);
    const auto *mk_756 = buffer.data(mk + 756);
    const auto *mk_757 = buffer.data(mk + 757);
    const auto *mk_758 = buffer.data(mk + 758);
    const auto *mk_759 = buffer.data(mk + 759);
    const auto *mk_761 = buffer.data(mk + 761);
    const auto *mk_762 = buffer.data(mk + 762);
    const auto *mk_763 = buffer.data(mk + 763);
    const auto *mk_765 = buffer.data(mk + 765);
    const auto *mk_766 = buffer.data(mk + 766);
    const auto *mk_767 = buffer.data(mk + 767);
    const auto *mk_768 = buffer.data(mk + 768);
    const auto *mk_770 = buffer.data(mk + 770);
    const auto *mk_771 = buffer.data(mk + 771);
    const auto *mk_772 = buffer.data(mk + 772);
    const auto *mk_773 = buffer.data(mk + 773);
    const auto *mk_774 = buffer.data(mk + 774);
    const auto *mk_776 = buffer.data(mk + 776);
    const auto *mk_777 = buffer.data(mk + 777);
    const auto *mk_784 = buffer.data(mk + 784);
    const auto *mk_785 = buffer.data(mk + 785);
    const auto *mk_786 = buffer.data(mk + 786);
    const auto *mk_787 = buffer.data(mk + 787);
    const auto *mk_788 = buffer.data(mk + 788);
    const auto *mk_789 = buffer.data(mk + 789);
    const auto *mk_790 = buffer.data(mk + 790);
    const auto *mk_791 = buffer.data(mk + 791);
    const auto *mk_792 = buffer.data(mk + 792);
    const auto *mk_794 = buffer.data(mk + 794);
    const auto *mk_795 = buffer.data(mk + 795);
    const auto *mk_797 = buffer.data(mk + 797);
    const auto *mk_798 = buffer.data(mk + 798);
    const auto *mk_801 = buffer.data(mk + 801);
    const auto *mk_802 = buffer.data(mk + 802);
    const auto *mk_806 = buffer.data(mk + 806);
    const auto *mk_807 = buffer.data(mk + 807);
    const auto *mk_812 = buffer.data(mk + 812);
    const auto *mk_820 = buffer.data(mk + 820);
    const auto *mk_821 = buffer.data(mk + 821);
    const auto *mk_822 = buffer.data(mk + 822);
    const auto *mk_823 = buffer.data(mk + 823);
    const auto *mk_824 = buffer.data(mk + 824);
    const auto *mk_825 = buffer.data(mk + 825);
    const auto *mk_826 = buffer.data(mk + 826);
    const auto *mk_827 = buffer.data(mk + 827);
    const auto *mk_828 = buffer.data(mk + 828);
    const auto *mk_830 = buffer.data(mk + 830);
    const auto *mk_831 = buffer.data(mk + 831);

#pragma omp simd aligned(t_915, t_916, t_917, t_918, pb_y, pb_z, lk_514, mi0_566, mi0_568, \
                         mi0_569, mi1_566, mi1_568, mi1_569, mk_730, mk_732, \
                         mk_733 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_9 * mi0_566[k]
                   - f_10 * mi1_566[k]
                   + pb_y[k] * mk_730[k];

        t_916[k] = f_17 * lk_514[k]
                   + pb_z[k] * mk_730[k];

        t_917[k] = f_5 * mi0_568[k]
                   - f_6 * mi1_568[k]
                   + pb_y[k] * mk_732[k];

        t_918[k] = f_3 * mi0_569[k]
                   - f_4 * mi1_569[k]
                   + pb_y[k] * mk_733[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, pb_x, pb_y, pb_z, lk_519, lk_740, \
                         mi0_570, mi0_580, mi1_570, mi1_580, mk_734, mk_735, \
                         mk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = pb_y[k] * mk_734[k];

        t_920[k] = f_16 * lk_740[k]
                   + f_5 * mi0_580[k]
                   - f_6 * mi1_580[k]
                   + pb_x[k] * mk_740[k];

        t_921[k] = f_11 * mi0_570[k]
                   - f_12 * mi1_570[k]
                   + pb_y[k] * mk_735[k];

        t_922[k] = f_17 * lk_519[k]
                   + pb_z[k] * mk_735[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, pb_y, mi0_572, mi0_573, mi0_574, mi1_572, \
                         mi1_573, mi1_574, mk_737, mk_738, mk_739, \
                         mk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = f_7 * mi0_572[k]
                   - f_8 * mi1_572[k]
                   + pb_y[k] * mk_737[k];

        t_924[k] = f_5 * mi0_573[k]
                   - f_6 * mi1_573[k]
                   + pb_y[k] * mk_738[k];

        t_925[k] = f_3 * mi0_574[k]
                   - f_4 * mi1_574[k]
                   + pb_y[k] * mk_739[k];

        t_926[k] = pb_y[k] * mk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pb_x, lk_747, lk_748, lk_749, lk_750, \
                         mi0_587, mi1_587, mk_747, mk_748, mk_749, \
                         mk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_16 * lk_747[k]
                   + f_3 * mi0_587[k]
                   - f_4 * mi1_587[k]
                   + pb_x[k] * mk_747[k];

        t_928[k] = f_16 * lk_748[k]
                   + pb_x[k] * mk_748[k];

        t_929[k] = f_16 * lk_749[k]
                   + pb_x[k] * mk_749[k];

        t_930[k] = f_16 * lk_750[k]
                   + pb_x[k] * mk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pb_x, pb_y, lk_751, lk_752, \
                         lk_753, lk_755, mk_747, mk_751, mk_752, mk_753, \
                         mk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_16 * lk_751[k]
                   + pb_x[k] * mk_751[k];

        t_932[k] = f_16 * lk_752[k]
                   + pb_x[k] * mk_752[k];

        t_933[k] = f_16 * lk_753[k]
                   + pb_x[k] * mk_753[k];

        t_934[k] = pb_y[k] * mk_747[k];

        t_935[k] = f_16 * lk_755[k]
                   + pb_x[k] * mk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pb_y, pb_z, lk_532, mi0_581, mi0_583, \
                         mi0_584, mi1_581, mi1_583, mi1_584, mk_748, mk_750, \
                         mk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * mi0_581[k]
                   - f_2 * mi1_581[k]
                   + pb_y[k] * mk_748[k];

        t_937[k] = f_17 * lk_532[k]
                   + pb_z[k] * mk_748[k];

        t_938[k] = f_11 * mi0_583[k]
                   - f_12 * mi1_583[k]
                   + pb_y[k] * mk_750[k];

        t_939[k] = f_9 * mi0_584[k]
                   - f_10 * mi1_584[k]
                   + pb_y[k] * mk_751[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pb_y, mi0_585, mi0_586, mi0_587, mi1_585, \
                         mi1_586, mi1_587, mk_752, mk_753, mk_754, \
                         mk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_7 * mi0_585[k]
                   - f_8 * mi1_585[k]
                   + pb_y[k] * mk_752[k];

        t_941[k] = f_5 * mi0_586[k]
                   - f_6 * mi1_586[k]
                   + pb_y[k] * mk_753[k];

        t_942[k] = f_3 * mi0_587[k]
                   - f_4 * mi1_587[k]
                   + pb_y[k] * mk_754[k];

        t_943[k] = pb_y[k] * mk_755[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, t_947, pa_x, pa_y, pb_y, pb_z, kl0_450, kl0_944, \
                         kl1_450, kl1_944, lk_540, ll_675, ll_944, \
                         mk_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_29 * kl0_944[k]
                   - f_30 * kl1_944[k]
                   + pa_x[k] * ll_944[k];

        t_945[k] = f_27 * kl0_450[k]
                   - f_28 * kl1_450[k]
                   + pa_y[k] * ll_675[k];

        t_946[k] = f_18 * lk_540[k]
                   + pb_y[k] * mk_756[k];

        t_947[k] = pb_z[k] * mk_756[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, pb_x, pb_z, lk_759, mi0_588, mi0_591, mi1_588, \
                         mi1_591, mk_757, mk_758, mk_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_15 * lk_759[k]
                   + f_11 * mi0_591[k]
                   - f_12 * mi1_591[k]
                   + pb_x[k] * mk_759[k];

        t_949[k] = pb_z[k] * mk_757[k];

        t_950[k] = f_3 * mi0_588[k]
                   - f_4 * mi1_588[k]
                   + pb_z[k] * mk_758[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pb_x, pb_y, pb_z, lk_545, lk_762, \
                         mi0_590, mi0_594, mi1_590, mi1_594, mk_759, mk_761, \
                         mk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_15 * lk_762[k]
                   + f_9 * mi0_594[k]
                   - f_10 * mi1_594[k]
                   + pb_x[k] * mk_762[k];

        t_952[k] = pb_z[k] * mk_759[k];

        t_953[k] = f_18 * lk_545[k]
                   + pb_y[k] * mk_761[k];

        t_954[k] = f_5 * mi0_590[k]
                   - f_6 * mi1_590[k]
                   + pb_z[k] * mk_761[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, pb_x, pb_z, lk_766, mi0_591, mi0_598, mi1_591, \
                         mi1_598, mk_762, mk_763, mk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_15 * lk_766[k]
                   + f_7 * mi0_598[k]
                   - f_8 * mi1_598[k]
                   + pb_x[k] * mk_766[k];

        t_956[k] = pb_z[k] * mk_762[k];

        t_957[k] = f_3 * mi0_591[k]
                   - f_4 * mi1_591[k]
                   + pb_z[k] * mk_763[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, pb_x, pb_y, pb_z, lk_549, lk_771, \
                         mi0_593, mi0_603, mi1_593, mi1_603, mk_765, mk_766, \
                         mk_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_18 * lk_549[k]
                   + pb_y[k] * mk_765[k];

        t_959[k] = f_7 * mi0_593[k]
                   - f_8 * mi1_593[k]
                   + pb_z[k] * mk_765[k];

        t_960[k] = f_15 * lk_771[k]
                   + f_5 * mi0_603[k]
                   - f_6 * mi1_603[k]
                   + pb_x[k] * mk_771[k];

        t_961[k] = pb_z[k] * mk_766[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, t_965, pb_y, pb_z, lk_554, mi0_594, mi0_595, \
                         mi0_597, mi1_594, mi1_595, mi1_597, mk_767, mk_768, \
                         mk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_3 * mi0_594[k]
                   - f_4 * mi1_594[k]
                   + pb_z[k] * mk_767[k];

        t_963[k] = f_5 * mi0_595[k]
                   - f_6 * mi1_595[k]
                   + pb_z[k] * mk_768[k];

        t_964[k] = f_18 * lk_554[k]
                   + pb_y[k] * mk_770[k];

        t_965[k] = f_9 * mi0_597[k]
                   - f_10 * mi1_597[k]
                   + pb_z[k] * mk_770[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pb_x, pb_z, lk_777, mi0_598, mi0_609, mi1_598, \
                         mi1_609, mk_771, mk_772, mk_777 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_15 * lk_777[k]
                   + f_3 * mi0_609[k]
                   - f_4 * mi1_609[k]
                   + pb_x[k] * mk_777[k];

        t_967[k] = pb_z[k] * mk_771[k];

        t_968[k] = f_3 * mi0_598[k]
                   - f_4 * mi1_598[k]
                   + pb_z[k] * mk_772[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, pb_y, pb_z, lk_560, mi0_599, mi0_600, \
                         mi0_602, mi1_599, mi1_600, mi1_602, mk_773, mk_774, \
                         mk_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_5 * mi0_599[k]
                   - f_6 * mi1_599[k]
                   + pb_z[k] * mk_773[k];

        t_970[k] = f_7 * mi0_600[k]
                   - f_8 * mi1_600[k]
                   + pb_z[k] * mk_774[k];

        t_971[k] = f_18 * lk_560[k]
                   + pb_y[k] * mk_776[k];

        t_972[k] = f_11 * mi0_602[k]
                   - f_12 * mi1_602[k]
                   + pb_z[k] * mk_776[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, t_977, pb_x, pb_z, lk_784, lk_786, \
                         lk_787, lk_788, mk_777, mk_784, mk_786, mk_787, \
                         mk_788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_15 * lk_784[k]
                   + pb_x[k] * mk_784[k];

        t_974[k] = pb_z[k] * mk_777[k];

        t_975[k] = f_15 * lk_786[k]
                   + pb_x[k] * mk_786[k];

        t_976[k] = f_15 * lk_787[k]
                   + pb_x[k] * mk_787[k];

        t_977[k] = f_15 * lk_788[k]
                   + pb_x[k] * mk_788[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pa_x, pb_x, kl0_981, kl1_981, lk_789, \
                         lk_790, lk_791, ll_981, mk_789, mk_790, \
                         mk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_15 * lk_789[k]
                   + pb_x[k] * mk_789[k];

        t_979[k] = f_15 * lk_790[k]
                   + pb_x[k] * mk_790[k];

        t_980[k] = f_15 * lk_791[k]
                   + pb_x[k] * mk_791[k];

        t_981[k] = f_25 * kl0_981[k]
                   - f_26 * kl1_981[k]
                   + pa_x[k] * ll_981[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, t_985, pb_z, mi0_609, mi0_610, mi0_611, mi1_609, \
                         mi1_610, mi1_611, mk_784, mk_785, mk_786, \
                         mk_787 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = pb_z[k] * mk_784[k];

        t_983[k] = f_3 * mi0_609[k]
                   - f_4 * mi1_609[k]
                   + pb_z[k] * mk_785[k];

        t_984[k] = f_5 * mi0_610[k]
                   - f_6 * mi1_610[k]
                   + pb_z[k] * mk_786[k];

        t_985[k] = f_7 * mi0_611[k]
                   - f_8 * mi1_611[k]
                   + pb_z[k] * mk_787[k];
    }

#pragma omp simd aligned(t_986, t_987, t_988, t_989, pb_y, pb_z, lk_575, mi0_612, mi0_613, \
                         mi0_615, mi1_612, mi1_613, mi1_615, mk_788, mk_789, \
                         mk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_986[k] = f_9 * mi0_612[k]
                   - f_10 * mi1_612[k]
                   + pb_z[k] * mk_788[k];

        t_987[k] = f_11 * mi0_613[k]
                   - f_12 * mi1_613[k]
                   + pb_z[k] * mk_789[k];

        t_988[k] = f_18 * lk_575[k]
                   + pb_y[k] * mk_791[k];

        t_989[k] = f_1 * mi0_615[k]
                   - f_2 * mi1_615[k]
                   + pb_z[k] * mk_791[k];
    }

#pragma omp simd aligned(t_990, t_991, t_992, t_993, t_994, pa_z, pb_y, pb_z, lk_540, lk_578, \
                         ll_675, ll_676, ll_678, mk_792, mk_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_990[k] = pa_z[k] * ll_675[k];

        t_991[k] = pa_z[k] * ll_676[k];

        t_992[k] = f_13 * lk_540[k]
                   + pb_z[k] * mk_792[k];

        t_993[k] = pa_z[k] * ll_678[k];

        t_994[k] = f_17 * lk_578[k]
                   + pb_y[k] * mk_794[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pa_z, pb_y, pb_z, lk_542, lk_543, lk_581, \
                         ll_680, ll_681, mk_795, mk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_14 * lk_542[k]
                   + pa_z[k] * ll_680[k];

        t_996[k] = pa_z[k] * ll_681[k];

        t_997[k] = f_13 * lk_543[k]
                   + pb_z[k] * mk_795[k];

        t_998[k] = f_17 * lk_581[k]
                   + pb_y[k] * mk_797[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, t_1002, pa_z, pb_z, lk_545, lk_546, lk_547, \
                         ll_684, ll_685, ll_687, mk_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_15 * lk_545[k]
                   + pa_z[k] * ll_684[k];

        t_1000[k] = pa_z[k] * ll_685[k];

        t_1001[k] = f_13 * lk_546[k]
                    + pb_z[k] * mk_798[k];

        t_1002[k] = f_14 * lk_547[k]
                    + pa_z[k] * ll_687[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, t_1006, pa_z, pb_y, pb_z, lk_549, lk_550, \
                         lk_585, ll_689, ll_690, mk_801, mk_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_17 * lk_585[k]
                    + pb_y[k] * mk_801[k];

        t_1004[k] = f_16 * lk_549[k]
                    + pa_z[k] * ll_689[k];

        t_1005[k] = pa_z[k] * ll_690[k];

        t_1006[k] = f_13 * lk_550[k]
                    + pb_z[k] * mk_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, t_1010, t_1011, pa_z, pb_y, lk_551, lk_552, \
                         lk_554, lk_590, ll_692, ll_693, ll_695, ll_696, \
                         mk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = f_14 * lk_551[k]
                    + pa_z[k] * ll_692[k];

        t_1008[k] = f_15 * lk_552[k]
                    + pa_z[k] * ll_693[k];

        t_1009[k] = f_17 * lk_590[k]
                    + pb_y[k] * mk_806[k];

        t_1010[k] = f_17 * lk_554[k]
                    + pa_z[k] * ll_695[k];

        t_1011[k] = pa_z[k] * ll_696[k];
    }

#pragma omp simd aligned(t_1012, t_1013, t_1014, t_1015, pa_z, pb_z, lk_555, lk_556, lk_557, \
                         lk_558, ll_698, ll_699, ll_700, mk_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1012[k] = f_13 * lk_555[k]
                    + pb_z[k] * mk_807[k];

        t_1013[k] = f_14 * lk_556[k]
                    + pa_z[k] * ll_698[k];

        t_1014[k] = f_15 * lk_557[k]
                    + pa_z[k] * ll_699[k];

        t_1015[k] = f_16 * lk_558[k]
                    + pa_z[k] * ll_700[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pa_z, pb_x, pb_y, lk_560, lk_596, \
                         lk_821, ll_702, ll_703, mk_812, mk_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_17 * lk_596[k]
                    + pb_y[k] * mk_812[k];

        t_1017[k] = f_18 * lk_560[k]
                    + pa_z[k] * ll_702[k];

        t_1018[k] = pa_z[k] * ll_703[k];

        t_1019[k] = f_15 * lk_821[k]
                    + pb_x[k] * mk_821[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, pb_x, lk_822, lk_823, lk_824, \
                         lk_825, lk_826, mk_822, mk_823, mk_824, mk_825, \
                         mk_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_15 * lk_822[k]
                    + pb_x[k] * mk_822[k];

        t_1021[k] = f_15 * lk_823[k]
                    + pb_x[k] * mk_823[k];

        t_1022[k] = f_15 * lk_824[k]
                    + pb_x[k] * mk_824[k];

        t_1023[k] = f_15 * lk_825[k]
                    + pb_x[k] * mk_825[k];

        t_1024[k] = f_15 * lk_826[k]
                    + pb_x[k] * mk_826[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pa_z, pb_x, pb_z, lk_568, lk_569, \
                         lk_827, ll_711, ll_713, mk_820, mk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_15 * lk_827[k]
                    + pb_x[k] * mk_827[k];

        t_1026[k] = pa_z[k] * ll_711[k];

        t_1027[k] = f_13 * lk_568[k]
                    + pb_z[k] * mk_820[k];

        t_1028[k] = f_14 * lk_569[k]
                    + pa_z[k] * ll_713[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pa_z, lk_570, lk_571, lk_572, lk_573, \
                         ll_714, ll_715, ll_716, ll_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = f_15 * lk_570[k]
                    + pa_z[k] * ll_714[k];

        t_1030[k] = f_16 * lk_571[k]
                    + pa_z[k] * ll_715[k];

        t_1031[k] = f_17 * lk_572[k]
                    + pa_z[k] * ll_716[k];

        t_1032[k] = f_18 * lk_573[k]
                    + pa_z[k] * ll_717[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, pa_y, pa_z, pb_y, kl0_540, kl1_540, \
                         lk_575, lk_611, lk_612, ll_719, ll_765, mk_827, \
                         mk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_17 * lk_611[k]
                    + pb_y[k] * mk_827[k];

        t_1034[k] = f_19 * lk_575[k]
                    + pa_z[k] * ll_719[k];

        t_1035[k] = f_29 * kl0_540[k]
                    - f_30 * kl1_540[k]
                    + pa_y[k] * ll_765[k];

        t_1036[k] = f_16 * lk_612[k]
                    + pb_y[k] * mk_828[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, pa_z, pb_y, pb_z, kl0_453, kl1_453, lk_576, \
                         lk_614, ll_723, mk_828, mk_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_14 * lk_576[k]
                    + pb_z[k] * mk_828[k];

        t_1038[k] = f_20 * kl0_453[k]
                    - f_21 * kl1_453[k]
                    + pa_z[k] * ll_723[k];

        t_1039[k] = f_16 * lk_614[k]
                    + pb_y[k] * mk_830[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pa_y, pa_z, pb_z, kl0_456, kl0_545, kl1_456, \
                         kl1_545, lk_579, ll_726, ll_770, mk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_29 * kl0_545[k]
                    - f_30 * kl1_545[k]
                    + pa_y[k] * ll_770[k];

        t_1041[k] = f_20 * kl0_456[k]
                    - f_21 * kl1_456[k]
                    + pa_z[k] * ll_726[k];

        t_1042[k] = f_14 * lk_579[k]
                    + pb_z[k] * mk_831[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece8(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_460 = buffer.data(kl0 + 460);
    const auto *kl0_465 = buffer.data(kl0 + 465);
    const auto *kl0_471 = buffer.data(kl0 + 471);
    const auto *kl0_498 = buffer.data(kl0 + 498);
    const auto *kl0_501 = buffer.data(kl0 + 501);
    const auto *kl0_505 = buffer.data(kl0 + 505);
    const auto *kl0_510 = buffer.data(kl0 + 510);
    const auto *kl0_516 = buffer.data(kl0 + 516);
    const auto *kl0_543 = buffer.data(kl0 + 543);
    const auto *kl0_546 = buffer.data(kl0 + 546);
    const auto *kl0_549 = buffer.data(kl0 + 549);
    const auto *kl0_550 = buffer.data(kl0 + 550);
    const auto *kl0_554 = buffer.data(kl0 + 554);
    const auto *kl0_555 = buffer.data(kl0 + 555);
    const auto *kl0_560 = buffer.data(kl0 + 560);
    const auto *kl0_561 = buffer.data(kl0 + 561);
    const auto *kl0_567 = buffer.data(kl0 + 567);
    const auto *kl0_585 = buffer.data(kl0 + 585);
    const auto *kl0_590 = buffer.data(kl0 + 590);
    const auto *kl0_594 = buffer.data(kl0 + 594);
    const auto *kl0_599 = buffer.data(kl0 + 599);
    const auto *kl0_605 = buffer.data(kl0 + 605);
    const auto *kl0_612 = buffer.data(kl0 + 612);
    const auto *kl0_630 = buffer.data(kl0 + 630);
    const auto *kl0_635 = buffer.data(kl0 + 635);
    const auto *kl0_639 = buffer.data(kl0 + 639);
    const auto *kl0_644 = buffer.data(kl0 + 644);
    const auto *kl0_650 = buffer.data(kl0 + 650);
    const auto *kl0_1071 = buffer.data(kl0 + 1071);
    const auto *kl0_1073 = buffer.data(kl0 + 1073);
    const auto *kl0_1074 = buffer.data(kl0 + 1074);
    const auto *kl0_1075 = buffer.data(kl0 + 1075);
    const auto *kl0_1076 = buffer.data(kl0 + 1076);
    const auto *kl0_1077 = buffer.data(kl0 + 1077);
    const auto *kl0_1079 = buffer.data(kl0 + 1079);
    const auto *kl0_1116 = buffer.data(kl0 + 1116);
    const auto *kl0_1118 = buffer.data(kl0 + 1118);
    const auto *kl0_1119 = buffer.data(kl0 + 1119);
    const auto *kl0_1120 = buffer.data(kl0 + 1120);
    const auto *kl0_1121 = buffer.data(kl0 + 1121);
    const auto *kl0_1122 = buffer.data(kl0 + 1122);
    const auto *kl0_1124 = buffer.data(kl0 + 1124);

    const auto *kl1_460 = buffer.data(kl1 + 460);
    const auto *kl1_465 = buffer.data(kl1 + 465);
    const auto *kl1_471 = buffer.data(kl1 + 471);
    const auto *kl1_498 = buffer.data(kl1 + 498);
    const auto *kl1_501 = buffer.data(kl1 + 501);
    const auto *kl1_505 = buffer.data(kl1 + 505);
    const auto *kl1_510 = buffer.data(kl1 + 510);
    const auto *kl1_516 = buffer.data(kl1 + 516);
    const auto *kl1_543 = buffer.data(kl1 + 543);
    const auto *kl1_546 = buffer.data(kl1 + 546);
    const auto *kl1_549 = buffer.data(kl1 + 549);
    const auto *kl1_550 = buffer.data(kl1 + 550);
    const auto *kl1_554 = buffer.data(kl1 + 554);
    const auto *kl1_555 = buffer.data(kl1 + 555);
    const auto *kl1_560 = buffer.data(kl1 + 560);
    const auto *kl1_561 = buffer.data(kl1 + 561);
    const auto *kl1_567 = buffer.data(kl1 + 567);
    const auto *kl1_585 = buffer.data(kl1 + 585);
    const auto *kl1_590 = buffer.data(kl1 + 590);
    const auto *kl1_594 = buffer.data(kl1 + 594);
    const auto *kl1_599 = buffer.data(kl1 + 599);
    const auto *kl1_605 = buffer.data(kl1 + 605);
    const auto *kl1_612 = buffer.data(kl1 + 612);
    const auto *kl1_630 = buffer.data(kl1 + 630);
    const auto *kl1_635 = buffer.data(kl1 + 635);
    const auto *kl1_639 = buffer.data(kl1 + 639);
    const auto *kl1_644 = buffer.data(kl1 + 644);
    const auto *kl1_650 = buffer.data(kl1 + 650);
    const auto *kl1_1071 = buffer.data(kl1 + 1071);
    const auto *kl1_1073 = buffer.data(kl1 + 1073);
    const auto *kl1_1074 = buffer.data(kl1 + 1074);
    const auto *kl1_1075 = buffer.data(kl1 + 1075);
    const auto *kl1_1076 = buffer.data(kl1 + 1076);
    const auto *kl1_1077 = buffer.data(kl1 + 1077);
    const auto *kl1_1079 = buffer.data(kl1 + 1079);
    const auto *kl1_1116 = buffer.data(kl1 + 1116);
    const auto *kl1_1118 = buffer.data(kl1 + 1118);
    const auto *kl1_1119 = buffer.data(kl1 + 1119);
    const auto *kl1_1120 = buffer.data(kl1 + 1120);
    const auto *kl1_1121 = buffer.data(kl1 + 1121);
    const auto *kl1_1122 = buffer.data(kl1 + 1122);
    const auto *kl1_1124 = buffer.data(kl1 + 1124);

    const auto *lk_582 = buffer.data(lk + 582);
    const auto *lk_586 = buffer.data(lk + 586);
    const auto *lk_591 = buffer.data(lk + 591);
    const auto *lk_604 = buffer.data(lk + 604);
    const auto *lk_612 = buffer.data(lk + 612);
    const auto *lk_615 = buffer.data(lk + 615);
    const auto *lk_617 = buffer.data(lk + 617);
    const auto *lk_618 = buffer.data(lk + 618);
    const auto *lk_621 = buffer.data(lk + 621);
    const auto *lk_622 = buffer.data(lk + 622);
    const auto *lk_626 = buffer.data(lk + 626);
    const auto *lk_627 = buffer.data(lk + 627);
    const auto *lk_632 = buffer.data(lk + 632);
    const auto *lk_640 = buffer.data(lk + 640);
    const auto *lk_647 = buffer.data(lk + 647);
    const auto *lk_648 = buffer.data(lk + 648);
    const auto *lk_650 = buffer.data(lk + 650);
    const auto *lk_651 = buffer.data(lk + 651);
    const auto *lk_653 = buffer.data(lk + 653);
    const auto *lk_654 = buffer.data(lk + 654);
    const auto *lk_657 = buffer.data(lk + 657);
    const auto *lk_658 = buffer.data(lk + 658);
    const auto *lk_662 = buffer.data(lk + 662);
    const auto *lk_663 = buffer.data(lk + 663);
    const auto *lk_668 = buffer.data(lk + 668);
    const auto *lk_683 = buffer.data(lk + 683);
    const auto *lk_684 = buffer.data(lk + 684);
    const auto *lk_686 = buffer.data(lk + 686);
    const auto *lk_689 = buffer.data(lk + 689);
    const auto *lk_693 = buffer.data(lk + 693);
    const auto *lk_698 = buffer.data(lk + 698);
    const auto *lk_840 = buffer.data(lk + 840);
    const auto *lk_845 = buffer.data(lk + 845);
    const auto *lk_846 = buffer.data(lk + 846);
    const auto *lk_851 = buffer.data(lk + 851);
    const auto *lk_852 = buffer.data(lk + 852);
    const auto *lk_853 = buffer.data(lk + 853);
    const auto *lk_856 = buffer.data(lk + 856);
    const auto *lk_857 = buffer.data(lk + 857);
    const auto *lk_858 = buffer.data(lk + 858);
    const auto *lk_859 = buffer.data(lk + 859);
    const auto *lk_860 = buffer.data(lk + 860);
    const auto *lk_861 = buffer.data(lk + 861);
    const auto *lk_862 = buffer.data(lk + 862);
    const auto *lk_863 = buffer.data(lk + 863);
    const auto *lk_876 = buffer.data(lk + 876);
    const auto *lk_881 = buffer.data(lk + 881);
    const auto *lk_882 = buffer.data(lk + 882);
    const auto *lk_887 = buffer.data(lk + 887);
    const auto *lk_888 = buffer.data(lk + 888);
    const auto *lk_889 = buffer.data(lk + 889);
    const auto *lk_892 = buffer.data(lk + 892);
    const auto *lk_893 = buffer.data(lk + 893);
    const auto *lk_894 = buffer.data(lk + 894);
    const auto *lk_895 = buffer.data(lk + 895);
    const auto *lk_896 = buffer.data(lk + 896);
    const auto *lk_897 = buffer.data(lk + 897);
    const auto *lk_898 = buffer.data(lk + 898);
    const auto *lk_899 = buffer.data(lk + 899);
    const auto *lk_912 = buffer.data(lk + 912);
    const auto *lk_917 = buffer.data(lk + 917);
    const auto *lk_918 = buffer.data(lk + 918);
    const auto *lk_923 = buffer.data(lk + 923);
    const auto *lk_924 = buffer.data(lk + 924);
    const auto *lk_925 = buffer.data(lk + 925);

    const auto *ll_730 = buffer.data(ll + 730);
    const auto *ll_735 = buffer.data(ll + 735);
    const auto *ll_741 = buffer.data(ll + 741);
    const auto *ll_768 = buffer.data(ll + 768);
    const auto *ll_771 = buffer.data(ll + 771);
    const auto *ll_774 = buffer.data(ll + 774);
    const auto *ll_775 = buffer.data(ll + 775);
    const auto *ll_779 = buffer.data(ll + 779);
    const auto *ll_780 = buffer.data(ll + 780);
    const auto *ll_785 = buffer.data(ll + 785);
    const auto *ll_786 = buffer.data(ll + 786);
    const auto *ll_792 = buffer.data(ll + 792);
    const auto *ll_810 = buffer.data(ll + 810);
    const auto *ll_813 = buffer.data(ll + 813);
    const auto *ll_815 = buffer.data(ll + 815);
    const auto *ll_816 = buffer.data(ll + 816);
    const auto *ll_819 = buffer.data(ll + 819);
    const auto *ll_820 = buffer.data(ll + 820);
    const auto *ll_824 = buffer.data(ll + 824);
    const auto *ll_825 = buffer.data(ll + 825);
    const auto *ll_830 = buffer.data(ll + 830);
    const auto *ll_831 = buffer.data(ll + 831);
    const auto *ll_837 = buffer.data(ll + 837);
    const auto *ll_855 = buffer.data(ll + 855);
    const auto *ll_860 = buffer.data(ll + 860);
    const auto *ll_864 = buffer.data(ll + 864);
    const auto *ll_869 = buffer.data(ll + 869);
    const auto *ll_875 = buffer.data(ll + 875);
    const auto *ll_1071 = buffer.data(ll + 1071);
    const auto *ll_1073 = buffer.data(ll + 1073);
    const auto *ll_1074 = buffer.data(ll + 1074);
    const auto *ll_1075 = buffer.data(ll + 1075);
    const auto *ll_1076 = buffer.data(ll + 1076);
    const auto *ll_1077 = buffer.data(ll + 1077);
    const auto *ll_1079 = buffer.data(ll + 1079);
    const auto *ll_1116 = buffer.data(ll + 1116);
    const auto *ll_1118 = buffer.data(ll + 1118);
    const auto *ll_1119 = buffer.data(ll + 1119);
    const auto *ll_1120 = buffer.data(ll + 1120);
    const auto *ll_1121 = buffer.data(ll + 1121);
    const auto *ll_1122 = buffer.data(ll + 1122);
    const auto *ll_1124 = buffer.data(ll + 1124);

    const auto *mi0_656 = buffer.data(mi0 + 656);
    const auto *mi0_661 = buffer.data(mi0 + 661);
    const auto *mi0_662 = buffer.data(mi0 + 662);
    const auto *mi0_667 = buffer.data(mi0 + 667);
    const auto *mi0_668 = buffer.data(mi0 + 668);
    const auto *mi0_669 = buffer.data(mi0 + 669);
    const auto *mi0_684 = buffer.data(mi0 + 684);
    const auto *mi0_689 = buffer.data(mi0 + 689);
    const auto *mi0_690 = buffer.data(mi0 + 690);
    const auto *mi0_695 = buffer.data(mi0 + 695);
    const auto *mi0_696 = buffer.data(mi0 + 696);
    const auto *mi0_697 = buffer.data(mi0 + 697);
    const auto *mi0_712 = buffer.data(mi0 + 712);
    const auto *mi0_717 = buffer.data(mi0 + 717);
    const auto *mi0_718 = buffer.data(mi0 + 718);
    const auto *mi0_723 = buffer.data(mi0 + 723);
    const auto *mi0_724 = buffer.data(mi0 + 724);
    const auto *mi0_725 = buffer.data(mi0 + 725);

    const auto *mi1_656 = buffer.data(mi1 + 656);
    const auto *mi1_661 = buffer.data(mi1 + 661);
    const auto *mi1_662 = buffer.data(mi1 + 662);
    const auto *mi1_667 = buffer.data(mi1 + 667);
    const auto *mi1_668 = buffer.data(mi1 + 668);
    const auto *mi1_669 = buffer.data(mi1 + 669);
    const auto *mi1_684 = buffer.data(mi1 + 684);
    const auto *mi1_689 = buffer.data(mi1 + 689);
    const auto *mi1_690 = buffer.data(mi1 + 690);
    const auto *mi1_695 = buffer.data(mi1 + 695);
    const auto *mi1_696 = buffer.data(mi1 + 696);
    const auto *mi1_697 = buffer.data(mi1 + 697);
    const auto *mi1_712 = buffer.data(mi1 + 712);
    const auto *mi1_717 = buffer.data(mi1 + 717);
    const auto *mi1_718 = buffer.data(mi1 + 718);
    const auto *mi1_723 = buffer.data(mi1 + 723);
    const auto *mi1_724 = buffer.data(mi1 + 724);
    const auto *mi1_725 = buffer.data(mi1 + 725);

    const auto *mk_833 = buffer.data(mk + 833);
    const auto *mk_834 = buffer.data(mk + 834);
    const auto *mk_837 = buffer.data(mk + 837);
    const auto *mk_838 = buffer.data(mk + 838);
    const auto *mk_840 = buffer.data(mk + 840);
    const auto *mk_842 = buffer.data(mk + 842);
    const auto *mk_843 = buffer.data(mk + 843);
    const auto *mk_845 = buffer.data(mk + 845);
    const auto *mk_846 = buffer.data(mk + 846);
    const auto *mk_848 = buffer.data(mk + 848);
    const auto *mk_851 = buffer.data(mk + 851);
    const auto *mk_852 = buffer.data(mk + 852);
    const auto *mk_853 = buffer.data(mk + 853);
    const auto *mk_856 = buffer.data(mk + 856);
    const auto *mk_857 = buffer.data(mk + 857);
    const auto *mk_858 = buffer.data(mk + 858);
    const auto *mk_859 = buffer.data(mk + 859);
    const auto *mk_860 = buffer.data(mk + 860);
    const auto *mk_861 = buffer.data(mk + 861);
    const auto *mk_862 = buffer.data(mk + 862);
    const auto *mk_863 = buffer.data(mk + 863);
    const auto *mk_864 = buffer.data(mk + 864);
    const auto *mk_866 = buffer.data(mk + 866);
    const auto *mk_867 = buffer.data(mk + 867);
    const auto *mk_869 = buffer.data(mk + 869);
    const auto *mk_870 = buffer.data(mk + 870);
    const auto *mk_873 = buffer.data(mk + 873);
    const auto *mk_874 = buffer.data(mk + 874);
    const auto *mk_876 = buffer.data(mk + 876);
    const auto *mk_878 = buffer.data(mk + 878);
    const auto *mk_879 = buffer.data(mk + 879);
    const auto *mk_881 = buffer.data(mk + 881);
    const auto *mk_882 = buffer.data(mk + 882);
    const auto *mk_884 = buffer.data(mk + 884);
    const auto *mk_887 = buffer.data(mk + 887);
    const auto *mk_888 = buffer.data(mk + 888);
    const auto *mk_889 = buffer.data(mk + 889);
    const auto *mk_892 = buffer.data(mk + 892);
    const auto *mk_893 = buffer.data(mk + 893);
    const auto *mk_894 = buffer.data(mk + 894);
    const auto *mk_895 = buffer.data(mk + 895);
    const auto *mk_896 = buffer.data(mk + 896);
    const auto *mk_897 = buffer.data(mk + 897);
    const auto *mk_898 = buffer.data(mk + 898);
    const auto *mk_899 = buffer.data(mk + 899);
    const auto *mk_900 = buffer.data(mk + 900);
    const auto *mk_902 = buffer.data(mk + 902);
    const auto *mk_903 = buffer.data(mk + 903);
    const auto *mk_905 = buffer.data(mk + 905);
    const auto *mk_906 = buffer.data(mk + 906);
    const auto *mk_909 = buffer.data(mk + 909);
    const auto *mk_910 = buffer.data(mk + 910);
    const auto *mk_912 = buffer.data(mk + 912);
    const auto *mk_914 = buffer.data(mk + 914);
    const auto *mk_915 = buffer.data(mk + 915);
    const auto *mk_917 = buffer.data(mk + 917);
    const auto *mk_918 = buffer.data(mk + 918);
    const auto *mk_923 = buffer.data(mk + 923);
    const auto *mk_924 = buffer.data(mk + 924);
    const auto *mk_925 = buffer.data(mk + 925);

#pragma omp simd aligned(t_1043, t_1044, t_1045, pa_y, pa_z, pb_y, kl0_460, kl0_549, kl1_460, \
                         kl1_549, lk_617, ll_730, ll_774, mk_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_16 * lk_617[k]
                    + pb_y[k] * mk_833[k];

        t_1044[k] = f_29 * kl0_549[k]
                    - f_30 * kl1_549[k]
                    + pa_y[k] * ll_774[k];

        t_1045[k] = f_20 * kl0_460[k]
                    - f_21 * kl1_460[k]
                    + pa_z[k] * ll_730[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pb_x, pb_y, pb_z, lk_582, lk_621, lk_840, \
                         mi0_656, mi1_656, mk_834, mk_837, mk_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_14 * lk_582[k]
                    + pb_z[k] * mk_834[k];

        t_1047[k] = f_15 * lk_840[k]
                    + f_7 * mi0_656[k]
                    - f_8 * mi1_656[k]
                    + pb_x[k] * mk_840[k];

        t_1048[k] = f_16 * lk_621[k]
                    + pb_y[k] * mk_837[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pa_y, pa_z, pb_z, kl0_465, kl0_554, kl1_465, \
                         kl1_554, lk_586, ll_735, ll_779, mk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_29 * kl0_554[k]
                    - f_30 * kl1_554[k]
                    + pa_y[k] * ll_779[k];

        t_1050[k] = f_20 * kl0_465[k]
                    - f_21 * kl1_465[k]
                    + pa_z[k] * ll_735[k];

        t_1051[k] = f_14 * lk_586[k]
                    + pb_z[k] * mk_838[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pb_x, pb_y, lk_626, lk_845, lk_846, mi0_661, \
                         mi0_662, mi1_661, mi1_662, mk_842, mk_845, \
                         mk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_15 * lk_845[k]
                    + f_5 * mi0_661[k]
                    - f_6 * mi1_661[k]
                    + pb_x[k] * mk_845[k];

        t_1053[k] = f_15 * lk_846[k]
                    + f_5 * mi0_662[k]
                    - f_6 * mi1_662[k]
                    + pb_x[k] * mk_846[k];

        t_1054[k] = f_16 * lk_626[k]
                    + pb_y[k] * mk_842[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pa_y, pa_z, pb_z, kl0_471, kl0_560, kl1_471, \
                         kl1_560, lk_591, ll_741, ll_785, mk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_29 * kl0_560[k]
                    - f_30 * kl1_560[k]
                    + pa_y[k] * ll_785[k];

        t_1056[k] = f_20 * kl0_471[k]
                    - f_21 * kl1_471[k]
                    + pa_z[k] * ll_741[k];

        t_1057[k] = f_14 * lk_591[k]
                    + pb_z[k] * mk_843[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pb_x, lk_851, lk_852, lk_853, mi0_667, \
                         mi0_668, mi0_669, mi1_667, mi1_668, mi1_669, mk_851, mk_852, \
                         mk_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_15 * lk_851[k]
                    + f_3 * mi0_667[k]
                    - f_4 * mi1_667[k]
                    + pb_x[k] * mk_851[k];

        t_1059[k] = f_15 * lk_852[k]
                    + f_3 * mi0_668[k]
                    - f_4 * mi1_668[k]
                    + pb_x[k] * mk_852[k];

        t_1060[k] = f_15 * lk_853[k]
                    + f_3 * mi0_669[k]
                    - f_4 * mi1_669[k]
                    + pb_x[k] * mk_853[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pa_y, pb_x, pb_y, kl0_567, kl1_567, \
                         lk_632, lk_856, lk_857, ll_792, mk_848, mk_856, \
                         mk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_16 * lk_632[k]
                    + pb_y[k] * mk_848[k];

        t_1062[k] = f_29 * kl0_567[k]
                    - f_30 * kl1_567[k]
                    + pa_y[k] * ll_792[k];

        t_1063[k] = f_15 * lk_856[k]
                    + pb_x[k] * mk_856[k];

        t_1064[k] = f_15 * lk_857[k]
                    + pb_x[k] * mk_857[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pb_x, lk_858, lk_859, lk_860, \
                         lk_861, lk_862, mk_858, mk_859, mk_860, mk_861, \
                         mk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_15 * lk_858[k]
                    + pb_x[k] * mk_858[k];

        t_1066[k] = f_15 * lk_859[k]
                    + pb_x[k] * mk_859[k];

        t_1067[k] = f_15 * lk_860[k]
                    + pb_x[k] * mk_860[k];

        t_1068[k] = f_15 * lk_861[k]
                    + pb_x[k] * mk_861[k];

        t_1069[k] = f_15 * lk_862[k]
                    + pb_x[k] * mk_862[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pa_x, pb_x, pb_z, kl0_1071, kl1_1071, lk_604, \
                         lk_863, ll_1071, mk_856, mk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_15 * lk_863[k]
                    + pb_x[k] * mk_863[k];

        t_1071[k] = f_25 * kl0_1071[k]
                    - f_26 * kl1_1071[k]
                    + pa_x[k] * ll_1071[k];

        t_1072[k] = f_14 * lk_604[k]
                    + pb_z[k] * mk_856[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pa_x, kl0_1073, kl0_1074, kl0_1075, kl1_1073, \
                         kl1_1074, kl1_1075, ll_1073, ll_1074, \
                         ll_1075 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_25 * kl0_1073[k]
                    - f_26 * kl1_1073[k]
                    + pa_x[k] * ll_1073[k];

        t_1074[k] = f_25 * kl0_1074[k]
                    - f_26 * kl1_1074[k]
                    + pa_x[k] * ll_1074[k];

        t_1075[k] = f_25 * kl0_1075[k]
                    - f_26 * kl1_1075[k]
                    + pa_x[k] * ll_1075[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pa_x, pb_y, kl0_1076, kl0_1077, kl1_1076, \
                         kl1_1077, lk_647, ll_1076, ll_1077, mk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_25 * kl0_1076[k]
                    - f_26 * kl1_1076[k]
                    + pa_x[k] * ll_1076[k];

        t_1077[k] = f_25 * kl0_1077[k]
                    - f_26 * kl1_1077[k]
                    + pa_x[k] * ll_1077[k];

        t_1078[k] = f_16 * lk_647[k]
                    + pb_y[k] * mk_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, pa_x, pa_y, pb_y, kl0_585, kl0_1079, kl1_585, \
                         kl1_1079, lk_648, ll_810, ll_1079, mk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_25 * kl0_1079[k]
                    - f_26 * kl1_1079[k]
                    + pa_x[k] * ll_1079[k];

        t_1080[k] = f_25 * kl0_585[k]
                    - f_26 * kl1_585[k]
                    + pa_y[k] * ll_810[k];

        t_1081[k] = f_15 * lk_648[k]
                    + pb_y[k] * mk_864[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pa_z, pb_y, pb_z, kl0_498, kl1_498, lk_612, \
                         lk_650, ll_768, mk_864, mk_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_15 * lk_612[k]
                    + pb_z[k] * mk_864[k];

        t_1083[k] = f_25 * kl0_498[k]
                    - f_26 * kl1_498[k]
                    + pa_z[k] * ll_768[k];

        t_1084[k] = f_15 * lk_650[k]
                    + pb_y[k] * mk_866[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pa_y, pa_z, pb_z, kl0_501, kl0_590, kl1_501, \
                         kl1_590, lk_615, ll_771, ll_815, mk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_25 * kl0_590[k]
                    - f_26 * kl1_590[k]
                    + pa_y[k] * ll_815[k];

        t_1086[k] = f_25 * kl0_501[k]
                    - f_26 * kl1_501[k]
                    + pa_z[k] * ll_771[k];

        t_1087[k] = f_15 * lk_615[k]
                    + pb_z[k] * mk_867[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pa_y, pa_z, pb_y, kl0_505, kl0_594, kl1_505, \
                         kl1_594, lk_653, ll_775, ll_819, mk_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_15 * lk_653[k]
                    + pb_y[k] * mk_869[k];

        t_1089[k] = f_25 * kl0_594[k]
                    - f_26 * kl1_594[k]
                    + pa_y[k] * ll_819[k];

        t_1090[k] = f_25 * kl0_505[k]
                    - f_26 * kl1_505[k]
                    + pa_z[k] * ll_775[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pb_x, pb_y, pb_z, lk_618, lk_657, lk_876, \
                         mi0_684, mi1_684, mk_870, mk_873, mk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_15 * lk_618[k]
                    + pb_z[k] * mk_870[k];

        t_1092[k] = f_15 * lk_876[k]
                    + f_7 * mi0_684[k]
                    - f_8 * mi1_684[k]
                    + pb_x[k] * mk_876[k];

        t_1093[k] = f_15 * lk_657[k]
                    + pb_y[k] * mk_873[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pa_y, pa_z, pb_z, kl0_510, kl0_599, kl1_510, \
                         kl1_599, lk_622, ll_780, ll_824, mk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_25 * kl0_599[k]
                    - f_26 * kl1_599[k]
                    + pa_y[k] * ll_824[k];

        t_1095[k] = f_25 * kl0_510[k]
                    - f_26 * kl1_510[k]
                    + pa_z[k] * ll_780[k];

        t_1096[k] = f_15 * lk_622[k]
                    + pb_z[k] * mk_874[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pb_x, pb_y, lk_662, lk_881, lk_882, mi0_689, \
                         mi0_690, mi1_689, mi1_690, mk_878, mk_881, \
                         mk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_15 * lk_881[k]
                    + f_5 * mi0_689[k]
                    - f_6 * mi1_689[k]
                    + pb_x[k] * mk_881[k];

        t_1098[k] = f_15 * lk_882[k]
                    + f_5 * mi0_690[k]
                    - f_6 * mi1_690[k]
                    + pb_x[k] * mk_882[k];

        t_1099[k] = f_15 * lk_662[k]
                    + pb_y[k] * mk_878[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pa_y, pa_z, pb_z, kl0_516, kl0_605, kl1_516, \
                         kl1_605, lk_627, ll_786, ll_830, mk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_25 * kl0_605[k]
                    - f_26 * kl1_605[k]
                    + pa_y[k] * ll_830[k];

        t_1101[k] = f_25 * kl0_516[k]
                    - f_26 * kl1_516[k]
                    + pa_z[k] * ll_786[k];

        t_1102[k] = f_15 * lk_627[k]
                    + pb_z[k] * mk_879[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pb_x, lk_887, lk_888, lk_889, mi0_695, \
                         mi0_696, mi0_697, mi1_695, mi1_696, mi1_697, mk_887, mk_888, \
                         mk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_15 * lk_887[k]
                    + f_3 * mi0_695[k]
                    - f_4 * mi1_695[k]
                    + pb_x[k] * mk_887[k];

        t_1104[k] = f_15 * lk_888[k]
                    + f_3 * mi0_696[k]
                    - f_4 * mi1_696[k]
                    + pb_x[k] * mk_888[k];

        t_1105[k] = f_15 * lk_889[k]
                    + f_3 * mi0_697[k]
                    - f_4 * mi1_697[k]
                    + pb_x[k] * mk_889[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pa_y, pb_x, pb_y, kl0_612, kl1_612, \
                         lk_668, lk_892, lk_893, ll_837, mk_884, mk_892, \
                         mk_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_15 * lk_668[k]
                    + pb_y[k] * mk_884[k];

        t_1107[k] = f_25 * kl0_612[k]
                    - f_26 * kl1_612[k]
                    + pa_y[k] * ll_837[k];

        t_1108[k] = f_15 * lk_892[k]
                    + pb_x[k] * mk_892[k];

        t_1109[k] = f_15 * lk_893[k]
                    + pb_x[k] * mk_893[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pb_x, lk_894, lk_895, lk_896, \
                         lk_897, lk_898, mk_894, mk_895, mk_896, mk_897, \
                         mk_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_15 * lk_894[k]
                    + pb_x[k] * mk_894[k];

        t_1111[k] = f_15 * lk_895[k]
                    + pb_x[k] * mk_895[k];

        t_1112[k] = f_15 * lk_896[k]
                    + pb_x[k] * mk_896[k];

        t_1113[k] = f_15 * lk_897[k]
                    + pb_x[k] * mk_897[k];

        t_1114[k] = f_15 * lk_898[k]
                    + pb_x[k] * mk_898[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pa_x, pb_x, pb_z, kl0_1116, kl1_1116, lk_640, \
                         lk_899, ll_1116, mk_892, mk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_15 * lk_899[k]
                    + pb_x[k] * mk_899[k];

        t_1116[k] = f_25 * kl0_1116[k]
                    - f_26 * kl1_1116[k]
                    + pa_x[k] * ll_1116[k];

        t_1117[k] = f_15 * lk_640[k]
                    + pb_z[k] * mk_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pa_x, kl0_1118, kl0_1119, kl0_1120, kl1_1118, \
                         kl1_1119, kl1_1120, ll_1118, ll_1119, \
                         ll_1120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_25 * kl0_1118[k]
                    - f_26 * kl1_1118[k]
                    + pa_x[k] * ll_1118[k];

        t_1119[k] = f_25 * kl0_1119[k]
                    - f_26 * kl1_1119[k]
                    + pa_x[k] * ll_1119[k];

        t_1120[k] = f_25 * kl0_1120[k]
                    - f_26 * kl1_1120[k]
                    + pa_x[k] * ll_1120[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pa_x, pb_y, kl0_1121, kl0_1122, kl1_1121, \
                         kl1_1122, lk_683, ll_1121, ll_1122, mk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_25 * kl0_1121[k]
                    - f_26 * kl1_1121[k]
                    + pa_x[k] * ll_1121[k];

        t_1122[k] = f_25 * kl0_1122[k]
                    - f_26 * kl1_1122[k]
                    + pa_x[k] * ll_1122[k];

        t_1123[k] = f_15 * lk_683[k]
                    + pb_y[k] * mk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pa_x, pa_y, pb_y, kl0_630, kl0_1124, kl1_630, \
                         kl1_1124, lk_684, ll_855, ll_1124, mk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_25 * kl0_1124[k]
                    - f_26 * kl1_1124[k]
                    + pa_x[k] * ll_1124[k];

        t_1125[k] = f_20 * kl0_630[k]
                    - f_21 * kl1_630[k]
                    + pa_y[k] * ll_855[k];

        t_1126[k] = f_14 * lk_684[k]
                    + pb_y[k] * mk_900[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pa_z, pb_y, pb_z, kl0_543, kl1_543, lk_648, \
                         lk_686, ll_813, mk_900, mk_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_16 * lk_648[k]
                    + pb_z[k] * mk_900[k];

        t_1128[k] = f_29 * kl0_543[k]
                    - f_30 * kl1_543[k]
                    + pa_z[k] * ll_813[k];

        t_1129[k] = f_14 * lk_686[k]
                    + pb_y[k] * mk_902[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pa_y, pa_z, pb_z, kl0_546, kl0_635, kl1_546, \
                         kl1_635, lk_651, ll_816, ll_860, mk_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_20 * kl0_635[k]
                    - f_21 * kl1_635[k]
                    + pa_y[k] * ll_860[k];

        t_1131[k] = f_29 * kl0_546[k]
                    - f_30 * kl1_546[k]
                    + pa_z[k] * ll_816[k];

        t_1132[k] = f_16 * lk_651[k]
                    + pb_z[k] * mk_903[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pa_y, pa_z, pb_y, kl0_550, kl0_639, kl1_550, \
                         kl1_639, lk_689, ll_820, ll_864, mk_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_14 * lk_689[k]
                    + pb_y[k] * mk_905[k];

        t_1134[k] = f_20 * kl0_639[k]
                    - f_21 * kl1_639[k]
                    + pa_y[k] * ll_864[k];

        t_1135[k] = f_29 * kl0_550[k]
                    - f_30 * kl1_550[k]
                    + pa_z[k] * ll_820[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pb_x, pb_y, pb_z, lk_654, lk_693, lk_912, \
                         mi0_712, mi1_712, mk_906, mk_909, mk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_16 * lk_654[k]
                    + pb_z[k] * mk_906[k];

        t_1137[k] = f_15 * lk_912[k]
                    + f_7 * mi0_712[k]
                    - f_8 * mi1_712[k]
                    + pb_x[k] * mk_912[k];

        t_1138[k] = f_14 * lk_693[k]
                    + pb_y[k] * mk_909[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, pa_y, pa_z, pb_z, kl0_555, kl0_644, kl1_555, \
                         kl1_644, lk_658, ll_825, ll_869, mk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_20 * kl0_644[k]
                    - f_21 * kl1_644[k]
                    + pa_y[k] * ll_869[k];

        t_1140[k] = f_29 * kl0_555[k]
                    - f_30 * kl1_555[k]
                    + pa_z[k] * ll_825[k];

        t_1141[k] = f_16 * lk_658[k]
                    + pb_z[k] * mk_910[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pb_x, pb_y, lk_698, lk_917, lk_918, mi0_717, \
                         mi0_718, mi1_717, mi1_718, mk_914, mk_917, \
                         mk_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_15 * lk_917[k]
                    + f_5 * mi0_717[k]
                    - f_6 * mi1_717[k]
                    + pb_x[k] * mk_917[k];

        t_1143[k] = f_15 * lk_918[k]
                    + f_5 * mi0_718[k]
                    - f_6 * mi1_718[k]
                    + pb_x[k] * mk_918[k];

        t_1144[k] = f_14 * lk_698[k]
                    + pb_y[k] * mk_914[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pa_y, pa_z, pb_z, kl0_561, kl0_650, kl1_561, \
                         kl1_650, lk_663, ll_831, ll_875, mk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_20 * kl0_650[k]
                    - f_21 * kl1_650[k]
                    + pa_y[k] * ll_875[k];

        t_1146[k] = f_29 * kl0_561[k]
                    - f_30 * kl1_561[k]
                    + pa_z[k] * ll_831[k];

        t_1147[k] = f_16 * lk_663[k]
                    + pb_z[k] * mk_915[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pb_x, lk_923, lk_924, lk_925, mi0_723, \
                         mi0_724, mi0_725, mi1_723, mi1_724, mi1_725, mk_923, mk_924, \
                         mk_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_15 * lk_923[k]
                    + f_3 * mi0_723[k]
                    - f_4 * mi1_723[k]
                    + pb_x[k] * mk_923[k];

        t_1149[k] = f_15 * lk_924[k]
                    + f_3 * mi0_724[k]
                    - f_4 * mi1_724[k]
                    + pb_x[k] * mk_924[k];

        t_1150[k] = f_15 * lk_925[k]
                    + f_3 * mi0_725[k]
                    - f_4 * mi1_725[k]
                    + pb_x[k] * mk_925[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece9(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kl0,
                                            const size_t kl1, const size_t lk, const size_t ll,
                                            const size_t mi0, const size_t mi1, const size_t mk,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;
    const auto f_23 = 3.0 / alpha;
    const auto f_24 = 3.0 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_27 = 2.5 / alpha;
    const auto f_28 = 2.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_630 = buffer.data(kl0 + 630);
    const auto *kl0_657 = buffer.data(kl0 + 657);
    const auto *kl0_675 = buffer.data(kl0 + 675);
    const auto *kl0_1161 = buffer.data(kl0 + 1161);
    const auto *kl0_1163 = buffer.data(kl0 + 1163);
    const auto *kl0_1164 = buffer.data(kl0 + 1164);
    const auto *kl0_1165 = buffer.data(kl0 + 1165);
    const auto *kl0_1166 = buffer.data(kl0 + 1166);
    const auto *kl0_1167 = buffer.data(kl0 + 1167);
    const auto *kl0_1169 = buffer.data(kl0 + 1169);
    const auto *kl0_1259 = buffer.data(kl0 + 1259);

    const auto *kl1_630 = buffer.data(kl1 + 630);
    const auto *kl1_657 = buffer.data(kl1 + 657);
    const auto *kl1_675 = buffer.data(kl1 + 675);
    const auto *kl1_1161 = buffer.data(kl1 + 1161);
    const auto *kl1_1163 = buffer.data(kl1 + 1163);
    const auto *kl1_1164 = buffer.data(kl1 + 1164);
    const auto *kl1_1165 = buffer.data(kl1 + 1165);
    const auto *kl1_1166 = buffer.data(kl1 + 1166);
    const auto *kl1_1167 = buffer.data(kl1 + 1167);
    const auto *kl1_1169 = buffer.data(kl1 + 1169);
    const auto *kl1_1259 = buffer.data(kl1 + 1259);

    const auto *lk_676 = buffer.data(lk + 676);
    const auto *lk_687 = buffer.data(lk + 687);
    const auto *lk_690 = buffer.data(lk + 690);
    const auto *lk_694 = buffer.data(lk + 694);
    const auto *lk_699 = buffer.data(lk + 699);
    const auto *lk_704 = buffer.data(lk + 704);
    const auto *lk_712 = buffer.data(lk + 712);
    const auto *lk_719 = buffer.data(lk + 719);
    const auto *lk_720 = buffer.data(lk + 720);
    const auto *lk_721 = buffer.data(lk + 721);
    const auto *lk_722 = buffer.data(lk + 722);
    const auto *lk_723 = buffer.data(lk + 723);
    const auto *lk_725 = buffer.data(lk + 725);
    const auto *lk_726 = buffer.data(lk + 726);
    const auto *lk_728 = buffer.data(lk + 728);
    const auto *lk_729 = buffer.data(lk + 729);
    const auto *lk_730 = buffer.data(lk + 730);
    const auto *lk_732 = buffer.data(lk + 732);
    const auto *lk_733 = buffer.data(lk + 733);
    const auto *lk_734 = buffer.data(lk + 734);
    const auto *lk_735 = buffer.data(lk + 735);
    const auto *lk_737 = buffer.data(lk + 737);
    const auto *lk_738 = buffer.data(lk + 738);
    const auto *lk_739 = buffer.data(lk + 739);
    const auto *lk_740 = buffer.data(lk + 740);
    const auto *lk_748 = buffer.data(lk + 748);
    const auto *lk_750 = buffer.data(lk + 750);
    const auto *lk_751 = buffer.data(lk + 751);
    const auto *lk_752 = buffer.data(lk + 752);
    const auto *lk_753 = buffer.data(lk + 753);
    const auto *lk_754 = buffer.data(lk + 754);
    const auto *lk_755 = buffer.data(lk + 755);
    const auto *lk_756 = buffer.data(lk + 756);
    const auto *lk_761 = buffer.data(lk + 761);
    const auto *lk_765 = buffer.data(lk + 765);
    const auto *lk_928 = buffer.data(lk + 928);
    const auto *lk_929 = buffer.data(lk + 929);
    const auto *lk_930 = buffer.data(lk + 930);
    const auto *lk_931 = buffer.data(lk + 931);
    const auto *lk_932 = buffer.data(lk + 932);
    const auto *lk_933 = buffer.data(lk + 933);
    const auto *lk_934 = buffer.data(lk + 934);
    const auto *lk_935 = buffer.data(lk + 935);
    const auto *lk_964 = buffer.data(lk + 964);
    const auto *lk_965 = buffer.data(lk + 965);
    const auto *lk_966 = buffer.data(lk + 966);
    const auto *lk_967 = buffer.data(lk + 967);
    const auto *lk_968 = buffer.data(lk + 968);
    const auto *lk_969 = buffer.data(lk + 969);
    const auto *lk_970 = buffer.data(lk + 970);
    const auto *lk_977 = buffer.data(lk + 977);
    const auto *lk_981 = buffer.data(lk + 981);
    const auto *lk_986 = buffer.data(lk + 986);
    const auto *lk_992 = buffer.data(lk + 992);
    const auto *lk_999 = buffer.data(lk + 999);
    const auto *lk_1000 = buffer.data(lk + 1000);
    const auto *lk_1001 = buffer.data(lk + 1001);
    const auto *lk_1002 = buffer.data(lk + 1002);
    const auto *lk_1003 = buffer.data(lk + 1003);
    const auto *lk_1004 = buffer.data(lk + 1004);
    const auto *lk_1005 = buffer.data(lk + 1005);
    const auto *lk_1007 = buffer.data(lk + 1007);
    const auto *lk_1011 = buffer.data(lk + 1011);
    const auto *lk_1014 = buffer.data(lk + 1014);
    const auto *lk_1018 = buffer.data(lk + 1018);
    const auto *lk_1023 = buffer.data(lk + 1023);

    const auto *ll_882 = buffer.data(ll + 882);
    const auto *ll_900 = buffer.data(ll + 900);
    const auto *ll_902 = buffer.data(ll + 902);
    const auto *ll_903 = buffer.data(ll + 903);
    const auto *ll_905 = buffer.data(ll + 905);
    const auto *ll_906 = buffer.data(ll + 906);
    const auto *ll_909 = buffer.data(ll + 909);
    const auto *ll_910 = buffer.data(ll + 910);
    const auto *ll_912 = buffer.data(ll + 912);
    const auto *ll_914 = buffer.data(ll + 914);
    const auto *ll_915 = buffer.data(ll + 915);
    const auto *ll_917 = buffer.data(ll + 917);
    const auto *ll_918 = buffer.data(ll + 918);
    const auto *ll_920 = buffer.data(ll + 920);
    const auto *ll_921 = buffer.data(ll + 921);
    const auto *ll_923 = buffer.data(ll + 923);
    const auto *ll_924 = buffer.data(ll + 924);
    const auto *ll_925 = buffer.data(ll + 925);
    const auto *ll_927 = buffer.data(ll + 927);
    const auto *ll_935 = buffer.data(ll + 935);
    const auto *ll_936 = buffer.data(ll + 936);
    const auto *ll_938 = buffer.data(ll + 938);
    const auto *ll_939 = buffer.data(ll + 939);
    const auto *ll_940 = buffer.data(ll + 940);
    const auto *ll_941 = buffer.data(ll + 941);
    const auto *ll_942 = buffer.data(ll + 942);
    const auto *ll_944 = buffer.data(ll + 944);
    const auto *ll_945 = buffer.data(ll + 945);
    const auto *ll_1161 = buffer.data(ll + 1161);
    const auto *ll_1163 = buffer.data(ll + 1163);
    const auto *ll_1164 = buffer.data(ll + 1164);
    const auto *ll_1165 = buffer.data(ll + 1165);
    const auto *ll_1166 = buffer.data(ll + 1166);
    const auto *ll_1167 = buffer.data(ll + 1167);
    const auto *ll_1169 = buffer.data(ll + 1169);
    const auto *ll_1259 = buffer.data(ll + 1259);

    const auto *mi0_756 = buffer.data(mi0 + 756);
    const auto *mi0_757 = buffer.data(mi0 + 757);
    const auto *mi0_759 = buffer.data(mi0 + 759);
    const auto *mi0_761 = buffer.data(mi0 + 761);
    const auto *mi0_762 = buffer.data(mi0 + 762);
    const auto *mi0_764 = buffer.data(mi0 + 764);
    const auto *mi0_765 = buffer.data(mi0 + 765);
    const auto *mi0_766 = buffer.data(mi0 + 766);
    const auto *mi0_768 = buffer.data(mi0 + 768);
    const auto *mi0_769 = buffer.data(mi0 + 769);
    const auto *mi0_770 = buffer.data(mi0 + 770);
    const auto *mi0_776 = buffer.data(mi0 + 776);
    const auto *mi0_777 = buffer.data(mi0 + 777);
    const auto *mi0_779 = buffer.data(mi0 + 779);
    const auto *mi0_780 = buffer.data(mi0 + 780);
    const auto *mi0_781 = buffer.data(mi0 + 781);
    const auto *mi0_782 = buffer.data(mi0 + 782);
    const auto *mi0_783 = buffer.data(mi0 + 783);
    const auto *mi0_784 = buffer.data(mi0 + 784);
    const auto *mi0_786 = buffer.data(mi0 + 786);
    const auto *mi0_787 = buffer.data(mi0 + 787);
    const auto *mi0_789 = buffer.data(mi0 + 789);
    const auto *mi0_790 = buffer.data(mi0 + 790);
    const auto *mi0_794 = buffer.data(mi0 + 794);
    const auto *mi0_799 = buffer.data(mi0 + 799);

    const auto *mi1_756 = buffer.data(mi1 + 756);
    const auto *mi1_757 = buffer.data(mi1 + 757);
    const auto *mi1_759 = buffer.data(mi1 + 759);
    const auto *mi1_761 = buffer.data(mi1 + 761);
    const auto *mi1_762 = buffer.data(mi1 + 762);
    const auto *mi1_764 = buffer.data(mi1 + 764);
    const auto *mi1_765 = buffer.data(mi1 + 765);
    const auto *mi1_766 = buffer.data(mi1 + 766);
    const auto *mi1_768 = buffer.data(mi1 + 768);
    const auto *mi1_769 = buffer.data(mi1 + 769);
    const auto *mi1_770 = buffer.data(mi1 + 770);
    const auto *mi1_776 = buffer.data(mi1 + 776);
    const auto *mi1_777 = buffer.data(mi1 + 777);
    const auto *mi1_779 = buffer.data(mi1 + 779);
    const auto *mi1_780 = buffer.data(mi1 + 780);
    const auto *mi1_781 = buffer.data(mi1 + 781);
    const auto *mi1_782 = buffer.data(mi1 + 782);
    const auto *mi1_783 = buffer.data(mi1 + 783);
    const auto *mi1_784 = buffer.data(mi1 + 784);
    const auto *mi1_786 = buffer.data(mi1 + 786);
    const auto *mi1_787 = buffer.data(mi1 + 787);
    const auto *mi1_789 = buffer.data(mi1 + 789);
    const auto *mi1_790 = buffer.data(mi1 + 790);
    const auto *mi1_794 = buffer.data(mi1 + 794);
    const auto *mi1_799 = buffer.data(mi1 + 799);

    const auto *mk_920 = buffer.data(mk + 920);
    const auto *mk_928 = buffer.data(mk + 928);
    const auto *mk_929 = buffer.data(mk + 929);
    const auto *mk_930 = buffer.data(mk + 930);
    const auto *mk_931 = buffer.data(mk + 931);
    const auto *mk_932 = buffer.data(mk + 932);
    const auto *mk_933 = buffer.data(mk + 933);
    const auto *mk_934 = buffer.data(mk + 934);
    const auto *mk_935 = buffer.data(mk + 935);
    const auto *mk_936 = buffer.data(mk + 936);
    const auto *mk_938 = buffer.data(mk + 938);
    const auto *mk_939 = buffer.data(mk + 939);
    const auto *mk_941 = buffer.data(mk + 941);
    const auto *mk_942 = buffer.data(mk + 942);
    const auto *mk_945 = buffer.data(mk + 945);
    const auto *mk_946 = buffer.data(mk + 946);
    const auto *mk_950 = buffer.data(mk + 950);
    const auto *mk_951 = buffer.data(mk + 951);
    const auto *mk_956 = buffer.data(mk + 956);
    const auto *mk_964 = buffer.data(mk + 964);
    const auto *mk_965 = buffer.data(mk + 965);
    const auto *mk_966 = buffer.data(mk + 966);
    const auto *mk_967 = buffer.data(mk + 967);
    const auto *mk_968 = buffer.data(mk + 968);
    const auto *mk_969 = buffer.data(mk + 969);
    const auto *mk_970 = buffer.data(mk + 970);
    const auto *mk_971 = buffer.data(mk + 971);
    const auto *mk_972 = buffer.data(mk + 972);
    const auto *mk_973 = buffer.data(mk + 973);
    const auto *mk_974 = buffer.data(mk + 974);
    const auto *mk_975 = buffer.data(mk + 975);
    const auto *mk_977 = buffer.data(mk + 977);
    const auto *mk_978 = buffer.data(mk + 978);
    const auto *mk_980 = buffer.data(mk + 980);
    const auto *mk_981 = buffer.data(mk + 981);
    const auto *mk_982 = buffer.data(mk + 982);
    const auto *mk_984 = buffer.data(mk + 984);
    const auto *mk_985 = buffer.data(mk + 985);
    const auto *mk_986 = buffer.data(mk + 986);
    const auto *mk_987 = buffer.data(mk + 987);
    const auto *mk_989 = buffer.data(mk + 989);
    const auto *mk_990 = buffer.data(mk + 990);
    const auto *mk_991 = buffer.data(mk + 991);
    const auto *mk_992 = buffer.data(mk + 992);
    const auto *mk_999 = buffer.data(mk + 999);
    const auto *mk_1000 = buffer.data(mk + 1000);
    const auto *mk_1001 = buffer.data(mk + 1001);
    const auto *mk_1002 = buffer.data(mk + 1002);
    const auto *mk_1003 = buffer.data(mk + 1003);
    const auto *mk_1004 = buffer.data(mk + 1004);
    const auto *mk_1005 = buffer.data(mk + 1005);
    const auto *mk_1006 = buffer.data(mk + 1006);
    const auto *mk_1007 = buffer.data(mk + 1007);
    const auto *mk_1008 = buffer.data(mk + 1008);
    const auto *mk_1009 = buffer.data(mk + 1009);
    const auto *mk_1010 = buffer.data(mk + 1010);
    const auto *mk_1011 = buffer.data(mk + 1011);
    const auto *mk_1013 = buffer.data(mk + 1013);
    const auto *mk_1014 = buffer.data(mk + 1014);
    const auto *mk_1015 = buffer.data(mk + 1015);
    const auto *mk_1017 = buffer.data(mk + 1017);
    const auto *mk_1018 = buffer.data(mk + 1018);
    const auto *mk_1023 = buffer.data(mk + 1023);

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pa_y, pb_x, pb_y, kl0_657, kl1_657, \
                         lk_704, lk_928, lk_929, ll_882, mk_920, mk_928, \
                         mk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_14 * lk_704[k]
                    + pb_y[k] * mk_920[k];

        t_1152[k] = f_20 * kl0_657[k]
                    - f_21 * kl1_657[k]
                    + pa_y[k] * ll_882[k];

        t_1153[k] = f_15 * lk_928[k]
                    + pb_x[k] * mk_928[k];

        t_1154[k] = f_15 * lk_929[k]
                    + pb_x[k] * mk_929[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pb_x, lk_930, lk_931, lk_932, \
                         lk_933, lk_934, mk_930, mk_931, mk_932, mk_933, \
                         mk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_15 * lk_930[k]
                    + pb_x[k] * mk_930[k];

        t_1156[k] = f_15 * lk_931[k]
                    + pb_x[k] * mk_931[k];

        t_1157[k] = f_15 * lk_932[k]
                    + pb_x[k] * mk_932[k];

        t_1158[k] = f_15 * lk_933[k]
                    + pb_x[k] * mk_933[k];

        t_1159[k] = f_15 * lk_934[k]
                    + pb_x[k] * mk_934[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pa_x, pb_x, pb_z, kl0_1161, kl1_1161, lk_676, \
                         lk_935, ll_1161, mk_928, mk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * lk_935[k]
                    + pb_x[k] * mk_935[k];

        t_1161[k] = f_25 * kl0_1161[k]
                    - f_26 * kl1_1161[k]
                    + pa_x[k] * ll_1161[k];

        t_1162[k] = f_16 * lk_676[k]
                    + pb_z[k] * mk_928[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pa_x, kl0_1163, kl0_1164, kl0_1165, kl1_1163, \
                         kl1_1164, kl1_1165, ll_1163, ll_1164, \
                         ll_1165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_25 * kl0_1163[k]
                    - f_26 * kl1_1163[k]
                    + pa_x[k] * ll_1163[k];

        t_1164[k] = f_25 * kl0_1164[k]
                    - f_26 * kl1_1164[k]
                    + pa_x[k] * ll_1164[k];

        t_1165[k] = f_25 * kl0_1165[k]
                    - f_26 * kl1_1165[k]
                    + pa_x[k] * ll_1165[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pa_x, pb_y, kl0_1166, kl0_1167, kl1_1166, \
                         kl1_1167, lk_719, ll_1166, ll_1167, mk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_25 * kl0_1166[k]
                    - f_26 * kl1_1166[k]
                    + pa_x[k] * ll_1166[k];

        t_1167[k] = f_25 * kl0_1167[k]
                    - f_26 * kl1_1167[k]
                    + pa_x[k] * ll_1167[k];

        t_1168[k] = f_14 * lk_719[k]
                    + pb_y[k] * mk_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pa_x, pa_y, pb_y, kl0_1169, kl1_1169, \
                         lk_720, ll_900, ll_902, ll_1169, mk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_25 * kl0_1169[k]
                    - f_26 * kl1_1169[k]
                    + pa_x[k] * ll_1169[k];

        t_1170[k] = pa_y[k] * ll_900[k];

        t_1171[k] = f_13 * lk_720[k]
                    + pb_y[k] * mk_936[k];

        t_1172[k] = pa_y[k] * ll_902[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, pa_y, pb_y, lk_721, lk_722, lk_723, \
                         ll_903, ll_905, ll_906, mk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_14 * lk_721[k]
                    + pa_y[k] * ll_903[k];

        t_1174[k] = f_13 * lk_722[k]
                    + pb_y[k] * mk_938[k];

        t_1175[k] = pa_y[k] * ll_905[k];

        t_1176[k] = f_15 * lk_723[k]
                    + pa_y[k] * ll_906[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, pa_y, pb_y, pb_z, lk_687, lk_725, \
                         lk_726, ll_909, ll_910, mk_939, mk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_17 * lk_687[k]
                    + pb_z[k] * mk_939[k];

        t_1178[k] = f_13 * lk_725[k]
                    + pb_y[k] * mk_941[k];

        t_1179[k] = pa_y[k] * ll_909[k];

        t_1180[k] = f_16 * lk_726[k]
                    + pa_y[k] * ll_910[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pa_y, pb_y, pb_z, lk_690, lk_728, \
                         lk_729, ll_912, ll_914, mk_942, mk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_17 * lk_690[k]
                    + pb_z[k] * mk_942[k];

        t_1182[k] = f_14 * lk_728[k]
                    + pa_y[k] * ll_912[k];

        t_1183[k] = f_13 * lk_729[k]
                    + pb_y[k] * mk_945[k];

        t_1184[k] = pa_y[k] * ll_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, t_1188, pa_y, pb_z, lk_694, lk_730, lk_732, \
                         lk_733, ll_915, ll_917, ll_918, mk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_17 * lk_730[k]
                    + pa_y[k] * ll_915[k];

        t_1186[k] = f_17 * lk_694[k]
                    + pb_z[k] * mk_946[k];

        t_1187[k] = f_15 * lk_732[k]
                    + pa_y[k] * ll_917[k];

        t_1188[k] = f_14 * lk_733[k]
                    + pa_y[k] * ll_918[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, pa_y, pb_y, pb_z, lk_699, lk_734, \
                         lk_735, ll_920, ll_921, mk_950, mk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_13 * lk_734[k]
                    + pb_y[k] * mk_950[k];

        t_1190[k] = pa_y[k] * ll_920[k];

        t_1191[k] = f_18 * lk_735[k]
                    + pa_y[k] * ll_921[k];

        t_1192[k] = f_17 * lk_699[k]
                    + pb_z[k] * mk_951[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, t_1196, t_1197, pa_y, pb_y, lk_737, lk_738, \
                         lk_739, lk_740, ll_923, ll_924, ll_925, ll_927, \
                         mk_956 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = f_16 * lk_737[k]
                    + pa_y[k] * ll_923[k];

        t_1194[k] = f_15 * lk_738[k]
                    + pa_y[k] * ll_924[k];

        t_1195[k] = f_14 * lk_739[k]
                    + pa_y[k] * ll_925[k];

        t_1196[k] = f_13 * lk_740[k]
                    + pb_y[k] * mk_956[k];

        t_1197[k] = pa_y[k] * ll_927[k];
    }

#pragma omp simd aligned(t_1198, t_1199, t_1200, t_1201, t_1202, pb_x, lk_964, lk_965, lk_966, \
                         lk_967, lk_968, mk_964, mk_965, mk_966, mk_967, \
                         mk_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1198[k] = f_15 * lk_964[k]
                    + pb_x[k] * mk_964[k];

        t_1199[k] = f_15 * lk_965[k]
                    + pb_x[k] * mk_965[k];

        t_1200[k] = f_15 * lk_966[k]
                    + pb_x[k] * mk_966[k];

        t_1201[k] = f_15 * lk_967[k]
                    + pb_x[k] * mk_967[k];

        t_1202[k] = f_15 * lk_968[k]
                    + pb_x[k] * mk_968[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, pa_y, pb_x, lk_748, lk_969, lk_970, \
                         ll_935, ll_936, mk_969, mk_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_15 * lk_969[k]
                    + pb_x[k] * mk_969[k];

        t_1204[k] = f_15 * lk_970[k]
                    + pb_x[k] * mk_970[k];

        t_1205[k] = pa_y[k] * ll_935[k];

        t_1206[k] = f_19 * lk_748[k]
                    + pa_y[k] * ll_936[k];
    }

#pragma omp simd aligned(t_1207, t_1208, t_1209, t_1210, pa_y, pb_z, lk_712, lk_750, lk_751, \
                         lk_752, ll_938, ll_939, ll_940, mk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1207[k] = f_17 * lk_712[k]
                    + pb_z[k] * mk_964[k];

        t_1208[k] = f_18 * lk_750[k]
                    + pa_y[k] * ll_938[k];

        t_1209[k] = f_17 * lk_751[k]
                    + pa_y[k] * ll_939[k];

        t_1210[k] = f_16 * lk_752[k]
                    + pa_y[k] * ll_940[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, t_1214, pa_y, pb_y, lk_753, lk_754, lk_755, \
                         ll_941, ll_942, ll_944, mk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_15 * lk_753[k]
                    + pa_y[k] * ll_941[k];

        t_1212[k] = f_14 * lk_754[k]
                    + pa_y[k] * ll_942[k];

        t_1213[k] = f_13 * lk_755[k]
                    + pb_y[k] * mk_971[k];

        t_1214[k] = pa_y[k] * ll_944[k];
    }

#pragma omp simd aligned(t_1215, t_1216, t_1217, t_1218, pa_z, pb_y, pb_z, kl0_630, kl1_630, \
                         lk_720, ll_900, mi0_756, mi1_756, mk_972, \
                         mk_973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1215[k] = f_27 * kl0_630[k]
                    - f_28 * kl1_630[k]
                    + pa_z[k] * ll_900[k];

        t_1216[k] = pb_y[k] * mk_972[k];

        t_1217[k] = f_18 * lk_720[k]
                    + pb_z[k] * mk_972[k];

        t_1218[k] = f_3 * mi0_756[k]
                    - f_4 * mi1_756[k]
                    + pb_y[k] * mk_973[k];
    }

#pragma omp simd aligned(t_1219, t_1220, t_1221, t_1222, pb_x, pb_y, pb_z, lk_723, lk_977, \
                         mi0_757, mi0_761, mi1_757, mi1_761, mk_974, mk_975, \
                         mk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1219[k] = pb_y[k] * mk_974[k];

        t_1220[k] = f_15 * lk_977[k]
                    + f_11 * mi0_761[k]
                    - f_12 * mi1_761[k]
                    + pb_x[k] * mk_977[k];

        t_1221[k] = f_5 * mi0_757[k]
                    - f_6 * mi1_757[k]
                    + pb_y[k] * mk_975[k];

        t_1222[k] = f_18 * lk_723[k]
                    + pb_z[k] * mk_975[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, t_1226, pb_x, pb_y, pb_z, lk_726, lk_981, \
                         mi0_759, mi0_765, mi1_759, mi1_765, mk_977, mk_978, \
                         mk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = pb_y[k] * mk_977[k];

        t_1224[k] = f_15 * lk_981[k]
                    + f_9 * mi0_765[k]
                    - f_10 * mi1_765[k]
                    + pb_x[k] * mk_981[k];

        t_1225[k] = f_7 * mi0_759[k]
                    - f_8 * mi1_759[k]
                    + pb_y[k] * mk_978[k];

        t_1226[k] = f_18 * lk_726[k]
                    + pb_z[k] * mk_978[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pb_x, pb_y, lk_986, mi0_761, mi0_770, \
                         mi1_761, mi1_770, mk_980, mk_981, mk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_3 * mi0_761[k]
                    - f_4 * mi1_761[k]
                    + pb_y[k] * mk_980[k];

        t_1228[k] = pb_y[k] * mk_981[k];

        t_1229[k] = f_15 * lk_986[k]
                    + f_7 * mi0_770[k]
                    - f_8 * mi1_770[k]
                    + pb_x[k] * mk_986[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, t_1233, pb_y, pb_z, lk_730, mi0_762, mi0_764, \
                         mi0_765, mi1_762, mi1_764, mi1_765, mk_982, mk_984, \
                         mk_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_9 * mi0_762[k]
                    - f_10 * mi1_762[k]
                    + pb_y[k] * mk_982[k];

        t_1231[k] = f_18 * lk_730[k]
                    + pb_z[k] * mk_982[k];

        t_1232[k] = f_5 * mi0_764[k]
                    - f_6 * mi1_764[k]
                    + pb_y[k] * mk_984[k];

        t_1233[k] = f_3 * mi0_765[k]
                    - f_4 * mi1_765[k]
                    + pb_y[k] * mk_985[k];
    }

#pragma omp simd aligned(t_1234, t_1235, t_1236, t_1237, pb_x, pb_y, pb_z, lk_735, lk_992, \
                         mi0_766, mi0_776, mi1_766, mi1_776, mk_986, mk_987, \
                         mk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1234[k] = pb_y[k] * mk_986[k];

        t_1235[k] = f_15 * lk_992[k]
                    + f_5 * mi0_776[k]
                    - f_6 * mi1_776[k]
                    + pb_x[k] * mk_992[k];

        t_1236[k] = f_11 * mi0_766[k]
                    - f_12 * mi1_766[k]
                    + pb_y[k] * mk_987[k];

        t_1237[k] = f_18 * lk_735[k]
                    + pb_z[k] * mk_987[k];
    }

#pragma omp simd aligned(t_1238, t_1239, t_1240, t_1241, pb_y, mi0_768, mi0_769, mi0_770, \
                         mi1_768, mi1_769, mi1_770, mk_989, mk_990, mk_991, \
                         mk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = f_7 * mi0_768[k]
                    - f_8 * mi1_768[k]
                    + pb_y[k] * mk_989[k];

        t_1239[k] = f_5 * mi0_769[k]
                    - f_6 * mi1_769[k]
                    + pb_y[k] * mk_990[k];

        t_1240[k] = f_3 * mi0_770[k]
                    - f_4 * mi1_770[k]
                    + pb_y[k] * mk_991[k];

        t_1241[k] = pb_y[k] * mk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pb_x, lk_999, lk_1000, lk_1001, \
                         lk_1002, mi0_783, mi1_783, mk_999, mk_1000, mk_1001, \
                         mk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_15 * lk_999[k]
                    + f_3 * mi0_783[k]
                    - f_4 * mi1_783[k]
                    + pb_x[k] * mk_999[k];

        t_1243[k] = f_15 * lk_1000[k]
                    + pb_x[k] * mk_1000[k];

        t_1244[k] = f_15 * lk_1001[k]
                    + pb_x[k] * mk_1001[k];

        t_1245[k] = f_15 * lk_1002[k]
                    + pb_x[k] * mk_1002[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, pb_x, pb_y, lk_1003, lk_1004, \
                         lk_1005, lk_1007, mk_999, mk_1003, mk_1004, mk_1005, \
                         mk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_15 * lk_1003[k]
                    + pb_x[k] * mk_1003[k];

        t_1247[k] = f_15 * lk_1004[k]
                    + pb_x[k] * mk_1004[k];

        t_1248[k] = f_15 * lk_1005[k]
                    + pb_x[k] * mk_1005[k];

        t_1249[k] = pb_y[k] * mk_999[k];

        t_1250[k] = f_15 * lk_1007[k]
                    + pb_x[k] * mk_1007[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, pb_y, pb_z, lk_748, mi0_777, mi0_779, \
                         mi0_780, mi1_777, mi1_779, mi1_780, mk_1000, mk_1002, \
                         mk_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_1 * mi0_777[k]
                    - f_2 * mi1_777[k]
                    + pb_y[k] * mk_1000[k];

        t_1252[k] = f_18 * lk_748[k]
                    + pb_z[k] * mk_1000[k];

        t_1253[k] = f_11 * mi0_779[k]
                    - f_12 * mi1_779[k]
                    + pb_y[k] * mk_1002[k];

        t_1254[k] = f_9 * mi0_780[k]
                    - f_10 * mi1_780[k]
                    + pb_y[k] * mk_1003[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, t_1258, pb_y, mi0_781, mi0_782, mi0_783, \
                         mi1_781, mi1_782, mi1_783, mk_1004, mk_1005, mk_1006, \
                         mk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = f_7 * mi0_781[k]
                    - f_8 * mi1_781[k]
                    + pb_y[k] * mk_1004[k];

        t_1256[k] = f_5 * mi0_782[k]
                    - f_6 * mi1_782[k]
                    + pb_y[k] * mk_1005[k];

        t_1257[k] = f_3 * mi0_783[k]
                    - f_4 * mi1_783[k]
                    + pb_y[k] * mk_1006[k];

        t_1258[k] = pb_y[k] * mk_1007[k];
    }

#pragma omp simd aligned(t_1259, t_1260, t_1261, t_1262, pa_x, pa_y, pb_y, pb_z, kl0_675, \
                         kl0_1259, kl1_675, kl1_1259, lk_756, ll_945, ll_1259, \
                         mk_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_25 * kl0_1259[k]
                    - f_26 * kl1_1259[k]
                    + pa_x[k] * ll_1259[k];

        t_1260[k] = f_23 * kl0_675[k]
                    - f_24 * kl1_675[k]
                    + pa_y[k] * ll_945[k];

        t_1261[k] = f_22 * lk_756[k]
                    + pb_y[k] * mk_1008[k];

        t_1262[k] = pb_z[k] * mk_1008[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, pb_x, pb_z, lk_1011, mi0_784, mi0_787, \
                         mi1_784, mi1_787, mk_1009, mk_1010, mk_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_14 * lk_1011[k]
                    + f_11 * mi0_787[k]
                    - f_12 * mi1_787[k]
                    + pb_x[k] * mk_1011[k];

        t_1264[k] = pb_z[k] * mk_1009[k];

        t_1265[k] = f_3 * mi0_784[k]
                    - f_4 * mi1_784[k]
                    + pb_z[k] * mk_1010[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, t_1269, pb_x, pb_y, pb_z, lk_761, lk_1014, \
                         mi0_786, mi0_790, mi1_786, mi1_790, mk_1011, mk_1013, \
                         mk_1014 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = f_14 * lk_1014[k]
                    + f_9 * mi0_790[k]
                    - f_10 * mi1_790[k]
                    + pb_x[k] * mk_1014[k];

        t_1267[k] = pb_z[k] * mk_1011[k];

        t_1268[k] = f_22 * lk_761[k]
                    + pb_y[k] * mk_1013[k];

        t_1269[k] = f_5 * mi0_786[k]
                    - f_6 * mi1_786[k]
                    + pb_z[k] * mk_1013[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, pb_x, pb_z, lk_1018, mi0_787, mi0_794, \
                         mi1_787, mi1_794, mk_1014, mk_1015, mk_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_14 * lk_1018[k]
                    + f_7 * mi0_794[k]
                    - f_8 * mi1_794[k]
                    + pb_x[k] * mk_1018[k];

        t_1271[k] = pb_z[k] * mk_1014[k];

        t_1272[k] = f_3 * mi0_787[k]
                    - f_4 * mi1_787[k]
                    + pb_z[k] * mk_1015[k];
    }

#pragma omp simd aligned(t_1273, t_1274, t_1275, t_1276, pb_x, pb_y, pb_z, lk_765, lk_1023, \
                         mi0_789, mi0_799, mi1_789, mi1_799, mk_1017, mk_1018, \
                         mk_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1273[k] = f_22 * lk_765[k]
                    + pb_y[k] * mk_1017[k];

        t_1274[k] = f_7 * mi0_789[k]
                    - f_8 * mi1_789[k]
                    + pb_z[k] * mk_1017[k];

        t_1275[k] = f_14 * lk_1023[k]
                    + f_5 * mi0_799[k]
                    - f_6 * mi1_799[k]
                    + pb_x[k] * mk_1023[k];

        t_1276[k] = pb_z[k] * mk_1018[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece10(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kl0,
                                             const size_t kl1, const size_t lk, const size_t ll,
                                             const size_t mi0, const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_678 = buffer.data(kl0 + 678);
    const auto *kl0_681 = buffer.data(kl0 + 681);
    const auto *kl0_685 = buffer.data(kl0 + 685);
    const auto *kl0_690 = buffer.data(kl0 + 690);
    const auto *kl0_696 = buffer.data(kl0 + 696);
    const auto *kl0_723 = buffer.data(kl0 + 723);
    const auto *kl0_765 = buffer.data(kl0 + 765);
    const auto *kl0_770 = buffer.data(kl0 + 770);
    const auto *kl0_774 = buffer.data(kl0 + 774);
    const auto *kl0_779 = buffer.data(kl0 + 779);
    const auto *kl0_785 = buffer.data(kl0 + 785);
    const auto *kl0_792 = buffer.data(kl0 + 792);
    const auto *kl0_810 = buffer.data(kl0 + 810);
    const auto *kl0_1296 = buffer.data(kl0 + 1296);
    const auto *kl0_1386 = buffer.data(kl0 + 1386);
    const auto *kl0_1388 = buffer.data(kl0 + 1388);
    const auto *kl0_1389 = buffer.data(kl0 + 1389);
    const auto *kl0_1390 = buffer.data(kl0 + 1390);
    const auto *kl0_1391 = buffer.data(kl0 + 1391);
    const auto *kl0_1392 = buffer.data(kl0 + 1392);
    const auto *kl0_1394 = buffer.data(kl0 + 1394);

    const auto *kl1_678 = buffer.data(kl1 + 678);
    const auto *kl1_681 = buffer.data(kl1 + 681);
    const auto *kl1_685 = buffer.data(kl1 + 685);
    const auto *kl1_690 = buffer.data(kl1 + 690);
    const auto *kl1_696 = buffer.data(kl1 + 696);
    const auto *kl1_723 = buffer.data(kl1 + 723);
    const auto *kl1_765 = buffer.data(kl1 + 765);
    const auto *kl1_770 = buffer.data(kl1 + 770);
    const auto *kl1_774 = buffer.data(kl1 + 774);
    const auto *kl1_779 = buffer.data(kl1 + 779);
    const auto *kl1_785 = buffer.data(kl1 + 785);
    const auto *kl1_792 = buffer.data(kl1 + 792);
    const auto *kl1_810 = buffer.data(kl1 + 810);
    const auto *kl1_1296 = buffer.data(kl1 + 1296);
    const auto *kl1_1386 = buffer.data(kl1 + 1386);
    const auto *kl1_1388 = buffer.data(kl1 + 1388);
    const auto *kl1_1389 = buffer.data(kl1 + 1389);
    const auto *kl1_1390 = buffer.data(kl1 + 1390);
    const auto *kl1_1391 = buffer.data(kl1 + 1391);
    const auto *kl1_1392 = buffer.data(kl1 + 1392);
    const auto *kl1_1394 = buffer.data(kl1 + 1394);

    const auto *lk_756 = buffer.data(lk + 756);
    const auto *lk_758 = buffer.data(lk + 758);
    const auto *lk_759 = buffer.data(lk + 759);
    const auto *lk_761 = buffer.data(lk + 761);
    const auto *lk_762 = buffer.data(lk + 762);
    const auto *lk_763 = buffer.data(lk + 763);
    const auto *lk_765 = buffer.data(lk + 765);
    const auto *lk_766 = buffer.data(lk + 766);
    const auto *lk_767 = buffer.data(lk + 767);
    const auto *lk_768 = buffer.data(lk + 768);
    const auto *lk_770 = buffer.data(lk + 770);
    const auto *lk_771 = buffer.data(lk + 771);
    const auto *lk_772 = buffer.data(lk + 772);
    const auto *lk_773 = buffer.data(lk + 773);
    const auto *lk_774 = buffer.data(lk + 774);
    const auto *lk_776 = buffer.data(lk + 776);
    const auto *lk_784 = buffer.data(lk + 784);
    const auto *lk_785 = buffer.data(lk + 785);
    const auto *lk_786 = buffer.data(lk + 786);
    const auto *lk_787 = buffer.data(lk + 787);
    const auto *lk_788 = buffer.data(lk + 788);
    const auto *lk_789 = buffer.data(lk + 789);
    const auto *lk_791 = buffer.data(lk + 791);
    const auto *lk_792 = buffer.data(lk + 792);
    const auto *lk_794 = buffer.data(lk + 794);
    const auto *lk_795 = buffer.data(lk + 795);
    const auto *lk_797 = buffer.data(lk + 797);
    const auto *lk_798 = buffer.data(lk + 798);
    const auto *lk_801 = buffer.data(lk + 801);
    const auto *lk_802 = buffer.data(lk + 802);
    const auto *lk_806 = buffer.data(lk + 806);
    const auto *lk_807 = buffer.data(lk + 807);
    const auto *lk_812 = buffer.data(lk + 812);
    const auto *lk_820 = buffer.data(lk + 820);
    const auto *lk_827 = buffer.data(lk + 827);
    const auto *lk_828 = buffer.data(lk + 828);
    const auto *lk_830 = buffer.data(lk + 830);
    const auto *lk_833 = buffer.data(lk + 833);
    const auto *lk_837 = buffer.data(lk + 837);
    const auto *lk_842 = buffer.data(lk + 842);
    const auto *lk_848 = buffer.data(lk + 848);
    const auto *lk_863 = buffer.data(lk + 863);
    const auto *lk_864 = buffer.data(lk + 864);
    const auto *lk_866 = buffer.data(lk + 866);
    const auto *lk_1029 = buffer.data(lk + 1029);
    const auto *lk_1036 = buffer.data(lk + 1036);
    const auto *lk_1038 = buffer.data(lk + 1038);
    const auto *lk_1039 = buffer.data(lk + 1039);
    const auto *lk_1040 = buffer.data(lk + 1040);
    const auto *lk_1041 = buffer.data(lk + 1041);
    const auto *lk_1042 = buffer.data(lk + 1042);
    const auto *lk_1043 = buffer.data(lk + 1043);
    const auto *lk_1073 = buffer.data(lk + 1073);
    const auto *lk_1074 = buffer.data(lk + 1074);
    const auto *lk_1075 = buffer.data(lk + 1075);
    const auto *lk_1076 = buffer.data(lk + 1076);
    const auto *lk_1077 = buffer.data(lk + 1077);
    const auto *lk_1078 = buffer.data(lk + 1078);
    const auto *lk_1079 = buffer.data(lk + 1079);
    const auto *lk_1092 = buffer.data(lk + 1092);
    const auto *lk_1097 = buffer.data(lk + 1097);
    const auto *lk_1098 = buffer.data(lk + 1098);
    const auto *lk_1103 = buffer.data(lk + 1103);
    const auto *lk_1104 = buffer.data(lk + 1104);
    const auto *lk_1105 = buffer.data(lk + 1105);
    const auto *lk_1108 = buffer.data(lk + 1108);
    const auto *lk_1109 = buffer.data(lk + 1109);
    const auto *lk_1110 = buffer.data(lk + 1110);
    const auto *lk_1111 = buffer.data(lk + 1111);
    const auto *lk_1112 = buffer.data(lk + 1112);
    const auto *lk_1113 = buffer.data(lk + 1113);
    const auto *lk_1114 = buffer.data(lk + 1114);
    const auto *lk_1115 = buffer.data(lk + 1115);

    const auto *ll_945 = buffer.data(ll + 945);
    const auto *ll_946 = buffer.data(ll + 946);
    const auto *ll_948 = buffer.data(ll + 948);
    const auto *ll_950 = buffer.data(ll + 950);
    const auto *ll_951 = buffer.data(ll + 951);
    const auto *ll_954 = buffer.data(ll + 954);
    const auto *ll_955 = buffer.data(ll + 955);
    const auto *ll_957 = buffer.data(ll + 957);
    const auto *ll_959 = buffer.data(ll + 959);
    const auto *ll_960 = buffer.data(ll + 960);
    const auto *ll_962 = buffer.data(ll + 962);
    const auto *ll_963 = buffer.data(ll + 963);
    const auto *ll_965 = buffer.data(ll + 965);
    const auto *ll_966 = buffer.data(ll + 966);
    const auto *ll_968 = buffer.data(ll + 968);
    const auto *ll_969 = buffer.data(ll + 969);
    const auto *ll_970 = buffer.data(ll + 970);
    const auto *ll_972 = buffer.data(ll + 972);
    const auto *ll_973 = buffer.data(ll + 973);
    const auto *ll_981 = buffer.data(ll + 981);
    const auto *ll_983 = buffer.data(ll + 983);
    const auto *ll_984 = buffer.data(ll + 984);
    const auto *ll_985 = buffer.data(ll + 985);
    const auto *ll_986 = buffer.data(ll + 986);
    const auto *ll_987 = buffer.data(ll + 987);
    const auto *ll_989 = buffer.data(ll + 989);
    const auto *ll_993 = buffer.data(ll + 993);
    const auto *ll_996 = buffer.data(ll + 996);
    const auto *ll_1000 = buffer.data(ll + 1000);
    const auto *ll_1005 = buffer.data(ll + 1005);
    const auto *ll_1011 = buffer.data(ll + 1011);
    const auto *ll_1035 = buffer.data(ll + 1035);
    const auto *ll_1038 = buffer.data(ll + 1038);
    const auto *ll_1040 = buffer.data(ll + 1040);
    const auto *ll_1044 = buffer.data(ll + 1044);
    const auto *ll_1049 = buffer.data(ll + 1049);
    const auto *ll_1055 = buffer.data(ll + 1055);
    const auto *ll_1062 = buffer.data(ll + 1062);
    const auto *ll_1080 = buffer.data(ll + 1080);
    const auto *ll_1296 = buffer.data(ll + 1296);
    const auto *ll_1386 = buffer.data(ll + 1386);
    const auto *ll_1388 = buffer.data(ll + 1388);
    const auto *ll_1389 = buffer.data(ll + 1389);
    const auto *ll_1390 = buffer.data(ll + 1390);
    const auto *ll_1391 = buffer.data(ll + 1391);
    const auto *ll_1392 = buffer.data(ll + 1392);
    const auto *ll_1394 = buffer.data(ll + 1394);

    const auto *mi0_790 = buffer.data(mi0 + 790);
    const auto *mi0_791 = buffer.data(mi0 + 791);
    const auto *mi0_793 = buffer.data(mi0 + 793);
    const auto *mi0_794 = buffer.data(mi0 + 794);
    const auto *mi0_795 = buffer.data(mi0 + 795);
    const auto *mi0_796 = buffer.data(mi0 + 796);
    const auto *mi0_798 = buffer.data(mi0 + 798);
    const auto *mi0_805 = buffer.data(mi0 + 805);
    const auto *mi0_806 = buffer.data(mi0 + 806);
    const auto *mi0_807 = buffer.data(mi0 + 807);
    const auto *mi0_808 = buffer.data(mi0 + 808);
    const auto *mi0_809 = buffer.data(mi0 + 809);
    const auto *mi0_811 = buffer.data(mi0 + 811);
    const auto *mi0_852 = buffer.data(mi0 + 852);
    const auto *mi0_857 = buffer.data(mi0 + 857);
    const auto *mi0_858 = buffer.data(mi0 + 858);
    const auto *mi0_863 = buffer.data(mi0 + 863);
    const auto *mi0_864 = buffer.data(mi0 + 864);
    const auto *mi0_865 = buffer.data(mi0 + 865);

    const auto *mi1_790 = buffer.data(mi1 + 790);
    const auto *mi1_791 = buffer.data(mi1 + 791);
    const auto *mi1_793 = buffer.data(mi1 + 793);
    const auto *mi1_794 = buffer.data(mi1 + 794);
    const auto *mi1_795 = buffer.data(mi1 + 795);
    const auto *mi1_796 = buffer.data(mi1 + 796);
    const auto *mi1_798 = buffer.data(mi1 + 798);
    const auto *mi1_805 = buffer.data(mi1 + 805);
    const auto *mi1_806 = buffer.data(mi1 + 806);
    const auto *mi1_807 = buffer.data(mi1 + 807);
    const auto *mi1_808 = buffer.data(mi1 + 808);
    const auto *mi1_809 = buffer.data(mi1 + 809);
    const auto *mi1_811 = buffer.data(mi1 + 811);
    const auto *mi1_852 = buffer.data(mi1 + 852);
    const auto *mi1_857 = buffer.data(mi1 + 857);
    const auto *mi1_858 = buffer.data(mi1 + 858);
    const auto *mi1_863 = buffer.data(mi1 + 863);
    const auto *mi1_864 = buffer.data(mi1 + 864);
    const auto *mi1_865 = buffer.data(mi1 + 865);

    const auto *mk_1019 = buffer.data(mk + 1019);
    const auto *mk_1020 = buffer.data(mk + 1020);
    const auto *mk_1022 = buffer.data(mk + 1022);
    const auto *mk_1023 = buffer.data(mk + 1023);
    const auto *mk_1024 = buffer.data(mk + 1024);
    const auto *mk_1025 = buffer.data(mk + 1025);
    const auto *mk_1026 = buffer.data(mk + 1026);
    const auto *mk_1028 = buffer.data(mk + 1028);
    const auto *mk_1029 = buffer.data(mk + 1029);
    const auto *mk_1036 = buffer.data(mk + 1036);
    const auto *mk_1037 = buffer.data(mk + 1037);
    const auto *mk_1038 = buffer.data(mk + 1038);
    const auto *mk_1039 = buffer.data(mk + 1039);
    const auto *mk_1040 = buffer.data(mk + 1040);
    const auto *mk_1041 = buffer.data(mk + 1041);
    const auto *mk_1042 = buffer.data(mk + 1042);
    const auto *mk_1043 = buffer.data(mk + 1043);
    const auto *mk_1044 = buffer.data(mk + 1044);
    const auto *mk_1046 = buffer.data(mk + 1046);
    const auto *mk_1047 = buffer.data(mk + 1047);
    const auto *mk_1049 = buffer.data(mk + 1049);
    const auto *mk_1050 = buffer.data(mk + 1050);
    const auto *mk_1053 = buffer.data(mk + 1053);
    const auto *mk_1054 = buffer.data(mk + 1054);
    const auto *mk_1058 = buffer.data(mk + 1058);
    const auto *mk_1059 = buffer.data(mk + 1059);
    const auto *mk_1064 = buffer.data(mk + 1064);
    const auto *mk_1072 = buffer.data(mk + 1072);
    const auto *mk_1073 = buffer.data(mk + 1073);
    const auto *mk_1074 = buffer.data(mk + 1074);
    const auto *mk_1075 = buffer.data(mk + 1075);
    const auto *mk_1076 = buffer.data(mk + 1076);
    const auto *mk_1077 = buffer.data(mk + 1077);
    const auto *mk_1078 = buffer.data(mk + 1078);
    const auto *mk_1079 = buffer.data(mk + 1079);
    const auto *mk_1080 = buffer.data(mk + 1080);
    const auto *mk_1082 = buffer.data(mk + 1082);
    const auto *mk_1083 = buffer.data(mk + 1083);
    const auto *mk_1085 = buffer.data(mk + 1085);
    const auto *mk_1086 = buffer.data(mk + 1086);
    const auto *mk_1089 = buffer.data(mk + 1089);
    const auto *mk_1090 = buffer.data(mk + 1090);
    const auto *mk_1092 = buffer.data(mk + 1092);
    const auto *mk_1094 = buffer.data(mk + 1094);
    const auto *mk_1095 = buffer.data(mk + 1095);
    const auto *mk_1097 = buffer.data(mk + 1097);
    const auto *mk_1098 = buffer.data(mk + 1098);
    const auto *mk_1100 = buffer.data(mk + 1100);
    const auto *mk_1103 = buffer.data(mk + 1103);
    const auto *mk_1104 = buffer.data(mk + 1104);
    const auto *mk_1105 = buffer.data(mk + 1105);
    const auto *mk_1108 = buffer.data(mk + 1108);
    const auto *mk_1109 = buffer.data(mk + 1109);
    const auto *mk_1110 = buffer.data(mk + 1110);
    const auto *mk_1111 = buffer.data(mk + 1111);
    const auto *mk_1112 = buffer.data(mk + 1112);
    const auto *mk_1113 = buffer.data(mk + 1113);
    const auto *mk_1114 = buffer.data(mk + 1114);
    const auto *mk_1115 = buffer.data(mk + 1115);
    const auto *mk_1116 = buffer.data(mk + 1116);
    const auto *mk_1118 = buffer.data(mk + 1118);

#pragma omp simd aligned(t_1277, t_1278, t_1279, t_1280, pb_y, pb_z, lk_770, mi0_790, mi0_791, \
                         mi0_793, mi1_790, mi1_791, mi1_793, mk_1019, mk_1020, \
                         mk_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = f_3 * mi0_790[k]
                    - f_4 * mi1_790[k]
                    + pb_z[k] * mk_1019[k];

        t_1278[k] = f_5 * mi0_791[k]
                    - f_6 * mi1_791[k]
                    + pb_z[k] * mk_1020[k];

        t_1279[k] = f_22 * lk_770[k]
                    + pb_y[k] * mk_1022[k];

        t_1280[k] = f_9 * mi0_793[k]
                    - f_10 * mi1_793[k]
                    + pb_z[k] * mk_1022[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, pb_x, pb_z, lk_1029, mi0_794, mi0_805, \
                         mi1_794, mi1_805, mk_1023, mk_1024, mk_1029 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_14 * lk_1029[k]
                    + f_3 * mi0_805[k]
                    - f_4 * mi1_805[k]
                    + pb_x[k] * mk_1029[k];

        t_1282[k] = pb_z[k] * mk_1023[k];

        t_1283[k] = f_3 * mi0_794[k]
                    - f_4 * mi1_794[k]
                    + pb_z[k] * mk_1024[k];
    }

#pragma omp simd aligned(t_1284, t_1285, t_1286, t_1287, pb_y, pb_z, lk_776, mi0_795, mi0_796, \
                         mi0_798, mi1_795, mi1_796, mi1_798, mk_1025, mk_1026, \
                         mk_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1284[k] = f_5 * mi0_795[k]
                    - f_6 * mi1_795[k]
                    + pb_z[k] * mk_1025[k];

        t_1285[k] = f_7 * mi0_796[k]
                    - f_8 * mi1_796[k]
                    + pb_z[k] * mk_1026[k];

        t_1286[k] = f_22 * lk_776[k]
                    + pb_y[k] * mk_1028[k];

        t_1287[k] = f_11 * mi0_798[k]
                    - f_12 * mi1_798[k]
                    + pb_z[k] * mk_1028[k];
    }

#pragma omp simd aligned(t_1288, t_1289, t_1290, t_1291, t_1292, pb_x, pb_z, lk_1036, lk_1038, \
                         lk_1039, lk_1040, mk_1029, mk_1036, mk_1038, mk_1039, \
                         mk_1040 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1288[k] = f_14 * lk_1036[k]
                    + pb_x[k] * mk_1036[k];

        t_1289[k] = pb_z[k] * mk_1029[k];

        t_1290[k] = f_14 * lk_1038[k]
                    + pb_x[k] * mk_1038[k];

        t_1291[k] = f_14 * lk_1039[k]
                    + pb_x[k] * mk_1039[k];

        t_1292[k] = f_14 * lk_1040[k]
                    + pb_x[k] * mk_1040[k];
    }

#pragma omp simd aligned(t_1293, t_1294, t_1295, t_1296, pa_x, pb_x, kl0_1296, kl1_1296, \
                         lk_1041, lk_1042, lk_1043, ll_1296, mk_1041, mk_1042, \
                         mk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1293[k] = f_14 * lk_1041[k]
                    + pb_x[k] * mk_1041[k];

        t_1294[k] = f_14 * lk_1042[k]
                    + pb_x[k] * mk_1042[k];

        t_1295[k] = f_14 * lk_1043[k]
                    + pb_x[k] * mk_1043[k];

        t_1296[k] = f_20 * kl0_1296[k]
                    - f_21 * kl1_1296[k]
                    + pa_x[k] * ll_1296[k];
    }

#pragma omp simd aligned(t_1297, t_1298, t_1299, t_1300, pb_z, mi0_805, mi0_806, mi0_807, \
                         mi1_805, mi1_806, mi1_807, mk_1036, mk_1037, mk_1038, \
                         mk_1039 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1297[k] = pb_z[k] * mk_1036[k];

        t_1298[k] = f_3 * mi0_805[k]
                    - f_4 * mi1_805[k]
                    + pb_z[k] * mk_1037[k];

        t_1299[k] = f_5 * mi0_806[k]
                    - f_6 * mi1_806[k]
                    + pb_z[k] * mk_1038[k];

        t_1300[k] = f_7 * mi0_807[k]
                    - f_8 * mi1_807[k]
                    + pb_z[k] * mk_1039[k];
    }

#pragma omp simd aligned(t_1301, t_1302, t_1303, t_1304, pb_y, pb_z, lk_791, mi0_808, mi0_809, \
                         mi0_811, mi1_808, mi1_809, mi1_811, mk_1040, mk_1041, \
                         mk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1301[k] = f_9 * mi0_808[k]
                    - f_10 * mi1_808[k]
                    + pb_z[k] * mk_1040[k];

        t_1302[k] = f_11 * mi0_809[k]
                    - f_12 * mi1_809[k]
                    + pb_z[k] * mk_1041[k];

        t_1303[k] = f_22 * lk_791[k]
                    + pb_y[k] * mk_1043[k];

        t_1304[k] = f_1 * mi0_811[k]
                    - f_2 * mi1_811[k]
                    + pb_z[k] * mk_1043[k];
    }

#pragma omp simd aligned(t_1305, t_1306, t_1307, t_1308, t_1309, pa_z, pb_y, pb_z, lk_756, \
                         lk_794, ll_945, ll_946, ll_948, mk_1044, \
                         mk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1305[k] = pa_z[k] * ll_945[k];

        t_1306[k] = pa_z[k] * ll_946[k];

        t_1307[k] = f_13 * lk_756[k]
                    + pb_z[k] * mk_1044[k];

        t_1308[k] = pa_z[k] * ll_948[k];

        t_1309[k] = f_18 * lk_794[k]
                    + pb_y[k] * mk_1046[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, t_1313, pa_z, pb_y, pb_z, lk_758, lk_759, \
                         lk_797, ll_950, ll_951, mk_1047, mk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_14 * lk_758[k]
                    + pa_z[k] * ll_950[k];

        t_1311[k] = pa_z[k] * ll_951[k];

        t_1312[k] = f_13 * lk_759[k]
                    + pb_z[k] * mk_1047[k];

        t_1313[k] = f_18 * lk_797[k]
                    + pb_y[k] * mk_1049[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, t_1317, pa_z, pb_z, lk_761, lk_762, lk_763, \
                         ll_954, ll_955, ll_957, mk_1050 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_15 * lk_761[k]
                    + pa_z[k] * ll_954[k];

        t_1315[k] = pa_z[k] * ll_955[k];

        t_1316[k] = f_13 * lk_762[k]
                    + pb_z[k] * mk_1050[k];

        t_1317[k] = f_14 * lk_763[k]
                    + pa_z[k] * ll_957[k];
    }

#pragma omp simd aligned(t_1318, t_1319, t_1320, t_1321, pa_z, pb_y, pb_z, lk_765, lk_766, \
                         lk_801, ll_959, ll_960, mk_1053, mk_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1318[k] = f_18 * lk_801[k]
                    + pb_y[k] * mk_1053[k];

        t_1319[k] = f_16 * lk_765[k]
                    + pa_z[k] * ll_959[k];

        t_1320[k] = pa_z[k] * ll_960[k];

        t_1321[k] = f_13 * lk_766[k]
                    + pb_z[k] * mk_1054[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, t_1325, t_1326, pa_z, pb_y, lk_767, lk_768, \
                         lk_770, lk_806, ll_962, ll_963, ll_965, ll_966, \
                         mk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = f_14 * lk_767[k]
                    + pa_z[k] * ll_962[k];

        t_1323[k] = f_15 * lk_768[k]
                    + pa_z[k] * ll_963[k];

        t_1324[k] = f_18 * lk_806[k]
                    + pb_y[k] * mk_1058[k];

        t_1325[k] = f_17 * lk_770[k]
                    + pa_z[k] * ll_965[k];

        t_1326[k] = pa_z[k] * ll_966[k];
    }

#pragma omp simd aligned(t_1327, t_1328, t_1329, t_1330, pa_z, pb_z, lk_771, lk_772, lk_773, \
                         lk_774, ll_968, ll_969, ll_970, mk_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1327[k] = f_13 * lk_771[k]
                    + pb_z[k] * mk_1059[k];

        t_1328[k] = f_14 * lk_772[k]
                    + pa_z[k] * ll_968[k];

        t_1329[k] = f_15 * lk_773[k]
                    + pa_z[k] * ll_969[k];

        t_1330[k] = f_16 * lk_774[k]
                    + pa_z[k] * ll_970[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, pa_z, pb_x, pb_y, lk_776, lk_812, \
                         lk_1073, ll_972, ll_973, mk_1064, mk_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_18 * lk_812[k]
                    + pb_y[k] * mk_1064[k];

        t_1332[k] = f_18 * lk_776[k]
                    + pa_z[k] * ll_972[k];

        t_1333[k] = pa_z[k] * ll_973[k];

        t_1334[k] = f_14 * lk_1073[k]
                    + pb_x[k] * mk_1073[k];
    }

#pragma omp simd aligned(t_1335, t_1336, t_1337, t_1338, t_1339, pb_x, lk_1074, lk_1075, \
                         lk_1076, lk_1077, lk_1078, mk_1074, mk_1075, mk_1076, mk_1077, \
                         mk_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1335[k] = f_14 * lk_1074[k]
                    + pb_x[k] * mk_1074[k];

        t_1336[k] = f_14 * lk_1075[k]
                    + pb_x[k] * mk_1075[k];

        t_1337[k] = f_14 * lk_1076[k]
                    + pb_x[k] * mk_1076[k];

        t_1338[k] = f_14 * lk_1077[k]
                    + pb_x[k] * mk_1077[k];

        t_1339[k] = f_14 * lk_1078[k]
                    + pb_x[k] * mk_1078[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, t_1343, pa_z, pb_x, pb_z, lk_784, lk_785, \
                         lk_1079, ll_981, ll_983, mk_1072, mk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_14 * lk_1079[k]
                    + pb_x[k] * mk_1079[k];

        t_1341[k] = pa_z[k] * ll_981[k];

        t_1342[k] = f_13 * lk_784[k]
                    + pb_z[k] * mk_1072[k];

        t_1343[k] = f_14 * lk_785[k]
                    + pa_z[k] * ll_983[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, t_1347, pa_z, lk_786, lk_787, lk_788, lk_789, \
                         ll_984, ll_985, ll_986, ll_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = f_15 * lk_786[k]
                    + pa_z[k] * ll_984[k];

        t_1345[k] = f_16 * lk_787[k]
                    + pa_z[k] * ll_985[k];

        t_1346[k] = f_17 * lk_788[k]
                    + pa_z[k] * ll_986[k];

        t_1347[k] = f_18 * lk_789[k]
                    + pa_z[k] * ll_987[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, t_1351, pa_y, pa_z, pb_y, kl0_765, kl1_765, \
                         lk_791, lk_827, lk_828, ll_989, ll_1035, mk_1079, \
                         mk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_18 * lk_827[k]
                    + pb_y[k] * mk_1079[k];

        t_1349[k] = f_19 * lk_791[k]
                    + pa_z[k] * ll_989[k];

        t_1350[k] = f_31 * kl0_765[k]
                    - f_32 * kl1_765[k]
                    + pa_y[k] * ll_1035[k];

        t_1351[k] = f_17 * lk_828[k]
                    + pb_y[k] * mk_1080[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pa_z, pb_y, pb_z, kl0_678, kl1_678, lk_792, \
                         lk_830, ll_993, mk_1080, mk_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_14 * lk_792[k]
                    + pb_z[k] * mk_1080[k];

        t_1353[k] = f_20 * kl0_678[k]
                    - f_21 * kl1_678[k]
                    + pa_z[k] * ll_993[k];

        t_1354[k] = f_17 * lk_830[k]
                    + pb_y[k] * mk_1082[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pa_y, pa_z, pb_z, kl0_681, kl0_770, kl1_681, \
                         kl1_770, lk_795, ll_996, ll_1040, mk_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = f_31 * kl0_770[k]
                    - f_32 * kl1_770[k]
                    + pa_y[k] * ll_1040[k];

        t_1356[k] = f_20 * kl0_681[k]
                    - f_21 * kl1_681[k]
                    + pa_z[k] * ll_996[k];

        t_1357[k] = f_14 * lk_795[k]
                    + pb_z[k] * mk_1083[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, pa_y, pa_z, pb_y, kl0_685, kl0_774, kl1_685, \
                         kl1_774, lk_833, ll_1000, ll_1044, mk_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_17 * lk_833[k]
                    + pb_y[k] * mk_1085[k];

        t_1359[k] = f_31 * kl0_774[k]
                    - f_32 * kl1_774[k]
                    + pa_y[k] * ll_1044[k];

        t_1360[k] = f_20 * kl0_685[k]
                    - f_21 * kl1_685[k]
                    + pa_z[k] * ll_1000[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, pb_x, pb_y, pb_z, lk_798, lk_837, lk_1092, \
                         mi0_852, mi1_852, mk_1086, mk_1089, mk_1092 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = f_14 * lk_798[k]
                    + pb_z[k] * mk_1086[k];

        t_1362[k] = f_14 * lk_1092[k]
                    + f_7 * mi0_852[k]
                    - f_8 * mi1_852[k]
                    + pb_x[k] * mk_1092[k];

        t_1363[k] = f_17 * lk_837[k]
                    + pb_y[k] * mk_1089[k];
    }

#pragma omp simd aligned(t_1364, t_1365, t_1366, pa_y, pa_z, pb_z, kl0_690, kl0_779, kl1_690, \
                         kl1_779, lk_802, ll_1005, ll_1049, mk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = f_31 * kl0_779[k]
                    - f_32 * kl1_779[k]
                    + pa_y[k] * ll_1049[k];

        t_1365[k] = f_20 * kl0_690[k]
                    - f_21 * kl1_690[k]
                    + pa_z[k] * ll_1005[k];

        t_1366[k] = f_14 * lk_802[k]
                    + pb_z[k] * mk_1090[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, pb_x, pb_y, lk_842, lk_1097, lk_1098, \
                         mi0_857, mi0_858, mi1_857, mi1_858, mk_1094, mk_1097, \
                         mk_1098 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = f_14 * lk_1097[k]
                    + f_5 * mi0_857[k]
                    - f_6 * mi1_857[k]
                    + pb_x[k] * mk_1097[k];

        t_1368[k] = f_14 * lk_1098[k]
                    + f_5 * mi0_858[k]
                    - f_6 * mi1_858[k]
                    + pb_x[k] * mk_1098[k];

        t_1369[k] = f_17 * lk_842[k]
                    + pb_y[k] * mk_1094[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, pa_y, pa_z, pb_z, kl0_696, kl0_785, kl1_696, \
                         kl1_785, lk_807, ll_1011, ll_1055, mk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = f_31 * kl0_785[k]
                    - f_32 * kl1_785[k]
                    + pa_y[k] * ll_1055[k];

        t_1371[k] = f_20 * kl0_696[k]
                    - f_21 * kl1_696[k]
                    + pa_z[k] * ll_1011[k];

        t_1372[k] = f_14 * lk_807[k]
                    + pb_z[k] * mk_1095[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pb_x, lk_1103, lk_1104, lk_1105, mi0_863, \
                         mi0_864, mi0_865, mi1_863, mi1_864, mi1_865, mk_1103, mk_1104, \
                         mk_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = f_14 * lk_1103[k]
                    + f_3 * mi0_863[k]
                    - f_4 * mi1_863[k]
                    + pb_x[k] * mk_1103[k];

        t_1374[k] = f_14 * lk_1104[k]
                    + f_3 * mi0_864[k]
                    - f_4 * mi1_864[k]
                    + pb_x[k] * mk_1104[k];

        t_1375[k] = f_14 * lk_1105[k]
                    + f_3 * mi0_865[k]
                    - f_4 * mi1_865[k]
                    + pb_x[k] * mk_1105[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, t_1379, pa_y, pb_x, pb_y, kl0_792, kl1_792, \
                         lk_848, lk_1108, lk_1109, ll_1062, mk_1100, mk_1108, \
                         mk_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_17 * lk_848[k]
                    + pb_y[k] * mk_1100[k];

        t_1377[k] = f_31 * kl0_792[k]
                    - f_32 * kl1_792[k]
                    + pa_y[k] * ll_1062[k];

        t_1378[k] = f_14 * lk_1108[k]
                    + pb_x[k] * mk_1108[k];

        t_1379[k] = f_14 * lk_1109[k]
                    + pb_x[k] * mk_1109[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, t_1384, pb_x, lk_1110, lk_1111, \
                         lk_1112, lk_1113, lk_1114, mk_1110, mk_1111, mk_1112, mk_1113, \
                         mk_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_14 * lk_1110[k]
                    + pb_x[k] * mk_1110[k];

        t_1381[k] = f_14 * lk_1111[k]
                    + pb_x[k] * mk_1111[k];

        t_1382[k] = f_14 * lk_1112[k]
                    + pb_x[k] * mk_1112[k];

        t_1383[k] = f_14 * lk_1113[k]
                    + pb_x[k] * mk_1113[k];

        t_1384[k] = f_14 * lk_1114[k]
                    + pb_x[k] * mk_1114[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, pa_x, pb_x, pb_z, kl0_1386, kl1_1386, lk_820, \
                         lk_1115, ll_1386, mk_1108, mk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_14 * lk_1115[k]
                    + pb_x[k] * mk_1115[k];

        t_1386[k] = f_20 * kl0_1386[k]
                    - f_21 * kl1_1386[k]
                    + pa_x[k] * ll_1386[k];

        t_1387[k] = f_14 * lk_820[k]
                    + pb_z[k] * mk_1108[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, pa_x, kl0_1388, kl0_1389, kl0_1390, kl1_1388, \
                         kl1_1389, kl1_1390, ll_1388, ll_1389, \
                         ll_1390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_20 * kl0_1388[k]
                    - f_21 * kl1_1388[k]
                    + pa_x[k] * ll_1388[k];

        t_1389[k] = f_20 * kl0_1389[k]
                    - f_21 * kl1_1389[k]
                    + pa_x[k] * ll_1389[k];

        t_1390[k] = f_20 * kl0_1390[k]
                    - f_21 * kl1_1390[k]
                    + pa_x[k] * ll_1390[k];
    }

#pragma omp simd aligned(t_1391, t_1392, t_1393, pa_x, pb_y, kl0_1391, kl0_1392, kl1_1391, \
                         kl1_1392, lk_863, ll_1391, ll_1392, mk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1391[k] = f_20 * kl0_1391[k]
                    - f_21 * kl1_1391[k]
                    + pa_x[k] * ll_1391[k];

        t_1392[k] = f_20 * kl0_1392[k]
                    - f_21 * kl1_1392[k]
                    + pa_x[k] * ll_1392[k];

        t_1393[k] = f_17 * lk_863[k]
                    + pb_y[k] * mk_1115[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, pa_x, pa_y, pb_y, kl0_810, kl0_1394, kl1_810, \
                         kl1_1394, lk_864, ll_1080, ll_1394, mk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_20 * kl0_1394[k]
                    - f_21 * kl1_1394[k]
                    + pa_x[k] * ll_1394[k];

        t_1395[k] = f_29 * kl0_810[k]
                    - f_30 * kl1_810[k]
                    + pa_y[k] * ll_1080[k];

        t_1396[k] = f_16 * lk_864[k]
                    + pb_y[k] * mk_1116[k];
    }

#pragma omp simd aligned(t_1397, t_1398, t_1399, pa_z, pb_y, pb_z, kl0_723, kl1_723, lk_828, \
                         lk_866, ll_1038, mk_1116, mk_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1397[k] = f_15 * lk_828[k]
                    + pb_z[k] * mk_1116[k];

        t_1398[k] = f_25 * kl0_723[k]
                    - f_26 * kl1_723[k]
                    + pa_z[k] * ll_1038[k];

        t_1399[k] = f_16 * lk_866[k]
                    + pb_y[k] * mk_1118[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece11(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kl0,
                                             const size_t kl1, const size_t lk, const size_t ll,
                                             const size_t mi0, const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_726 = buffer.data(kl0 + 726);
    const auto *kl0_730 = buffer.data(kl0 + 730);
    const auto *kl0_735 = buffer.data(kl0 + 735);
    const auto *kl0_741 = buffer.data(kl0 + 741);
    const auto *kl0_768 = buffer.data(kl0 + 768);
    const auto *kl0_771 = buffer.data(kl0 + 771);
    const auto *kl0_775 = buffer.data(kl0 + 775);
    const auto *kl0_780 = buffer.data(kl0 + 780);
    const auto *kl0_786 = buffer.data(kl0 + 786);
    const auto *kl0_813 = buffer.data(kl0 + 813);
    const auto *kl0_815 = buffer.data(kl0 + 815);
    const auto *kl0_816 = buffer.data(kl0 + 816);
    const auto *kl0_819 = buffer.data(kl0 + 819);
    const auto *kl0_820 = buffer.data(kl0 + 820);
    const auto *kl0_824 = buffer.data(kl0 + 824);
    const auto *kl0_825 = buffer.data(kl0 + 825);
    const auto *kl0_830 = buffer.data(kl0 + 830);
    const auto *kl0_831 = buffer.data(kl0 + 831);
    const auto *kl0_837 = buffer.data(kl0 + 837);
    const auto *kl0_855 = buffer.data(kl0 + 855);
    const auto *kl0_860 = buffer.data(kl0 + 860);
    const auto *kl0_864 = buffer.data(kl0 + 864);
    const auto *kl0_869 = buffer.data(kl0 + 869);
    const auto *kl0_875 = buffer.data(kl0 + 875);
    const auto *kl0_882 = buffer.data(kl0 + 882);
    const auto *kl0_900 = buffer.data(kl0 + 900);
    const auto *kl0_905 = buffer.data(kl0 + 905);
    const auto *kl0_909 = buffer.data(kl0 + 909);
    const auto *kl0_914 = buffer.data(kl0 + 914);
    const auto *kl0_920 = buffer.data(kl0 + 920);
    const auto *kl0_1431 = buffer.data(kl0 + 1431);
    const auto *kl0_1433 = buffer.data(kl0 + 1433);
    const auto *kl0_1434 = buffer.data(kl0 + 1434);
    const auto *kl0_1435 = buffer.data(kl0 + 1435);
    const auto *kl0_1436 = buffer.data(kl0 + 1436);
    const auto *kl0_1437 = buffer.data(kl0 + 1437);
    const auto *kl0_1439 = buffer.data(kl0 + 1439);
    const auto *kl0_1476 = buffer.data(kl0 + 1476);
    const auto *kl0_1478 = buffer.data(kl0 + 1478);
    const auto *kl0_1479 = buffer.data(kl0 + 1479);
    const auto *kl0_1480 = buffer.data(kl0 + 1480);
    const auto *kl0_1481 = buffer.data(kl0 + 1481);
    const auto *kl0_1482 = buffer.data(kl0 + 1482);
    const auto *kl0_1484 = buffer.data(kl0 + 1484);

    const auto *kl1_726 = buffer.data(kl1 + 726);
    const auto *kl1_730 = buffer.data(kl1 + 730);
    const auto *kl1_735 = buffer.data(kl1 + 735);
    const auto *kl1_741 = buffer.data(kl1 + 741);
    const auto *kl1_768 = buffer.data(kl1 + 768);
    const auto *kl1_771 = buffer.data(kl1 + 771);
    const auto *kl1_775 = buffer.data(kl1 + 775);
    const auto *kl1_780 = buffer.data(kl1 + 780);
    const auto *kl1_786 = buffer.data(kl1 + 786);
    const auto *kl1_813 = buffer.data(kl1 + 813);
    const auto *kl1_815 = buffer.data(kl1 + 815);
    const auto *kl1_816 = buffer.data(kl1 + 816);
    const auto *kl1_819 = buffer.data(kl1 + 819);
    const auto *kl1_820 = buffer.data(kl1 + 820);
    const auto *kl1_824 = buffer.data(kl1 + 824);
    const auto *kl1_825 = buffer.data(kl1 + 825);
    const auto *kl1_830 = buffer.data(kl1 + 830);
    const auto *kl1_831 = buffer.data(kl1 + 831);
    const auto *kl1_837 = buffer.data(kl1 + 837);
    const auto *kl1_855 = buffer.data(kl1 + 855);
    const auto *kl1_860 = buffer.data(kl1 + 860);
    const auto *kl1_864 = buffer.data(kl1 + 864);
    const auto *kl1_869 = buffer.data(kl1 + 869);
    const auto *kl1_875 = buffer.data(kl1 + 875);
    const auto *kl1_882 = buffer.data(kl1 + 882);
    const auto *kl1_900 = buffer.data(kl1 + 900);
    const auto *kl1_905 = buffer.data(kl1 + 905);
    const auto *kl1_909 = buffer.data(kl1 + 909);
    const auto *kl1_914 = buffer.data(kl1 + 914);
    const auto *kl1_920 = buffer.data(kl1 + 920);
    const auto *kl1_1431 = buffer.data(kl1 + 1431);
    const auto *kl1_1433 = buffer.data(kl1 + 1433);
    const auto *kl1_1434 = buffer.data(kl1 + 1434);
    const auto *kl1_1435 = buffer.data(kl1 + 1435);
    const auto *kl1_1436 = buffer.data(kl1 + 1436);
    const auto *kl1_1437 = buffer.data(kl1 + 1437);
    const auto *kl1_1439 = buffer.data(kl1 + 1439);
    const auto *kl1_1476 = buffer.data(kl1 + 1476);
    const auto *kl1_1478 = buffer.data(kl1 + 1478);
    const auto *kl1_1479 = buffer.data(kl1 + 1479);
    const auto *kl1_1480 = buffer.data(kl1 + 1480);
    const auto *kl1_1481 = buffer.data(kl1 + 1481);
    const auto *kl1_1482 = buffer.data(kl1 + 1482);
    const auto *kl1_1484 = buffer.data(kl1 + 1484);

    const auto *lk_831 = buffer.data(lk + 831);
    const auto *lk_834 = buffer.data(lk + 834);
    const auto *lk_838 = buffer.data(lk + 838);
    const auto *lk_843 = buffer.data(lk + 843);
    const auto *lk_856 = buffer.data(lk + 856);
    const auto *lk_864 = buffer.data(lk + 864);
    const auto *lk_867 = buffer.data(lk + 867);
    const auto *lk_869 = buffer.data(lk + 869);
    const auto *lk_870 = buffer.data(lk + 870);
    const auto *lk_873 = buffer.data(lk + 873);
    const auto *lk_874 = buffer.data(lk + 874);
    const auto *lk_878 = buffer.data(lk + 878);
    const auto *lk_879 = buffer.data(lk + 879);
    const auto *lk_884 = buffer.data(lk + 884);
    const auto *lk_892 = buffer.data(lk + 892);
    const auto *lk_899 = buffer.data(lk + 899);
    const auto *lk_900 = buffer.data(lk + 900);
    const auto *lk_902 = buffer.data(lk + 902);
    const auto *lk_903 = buffer.data(lk + 903);
    const auto *lk_905 = buffer.data(lk + 905);
    const auto *lk_906 = buffer.data(lk + 906);
    const auto *lk_909 = buffer.data(lk + 909);
    const auto *lk_910 = buffer.data(lk + 910);
    const auto *lk_914 = buffer.data(lk + 914);
    const auto *lk_915 = buffer.data(lk + 915);
    const auto *lk_920 = buffer.data(lk + 920);
    const auto *lk_935 = buffer.data(lk + 935);
    const auto *lk_936 = buffer.data(lk + 936);
    const auto *lk_938 = buffer.data(lk + 938);
    const auto *lk_941 = buffer.data(lk + 941);
    const auto *lk_945 = buffer.data(lk + 945);
    const auto *lk_950 = buffer.data(lk + 950);
    const auto *lk_1128 = buffer.data(lk + 1128);
    const auto *lk_1133 = buffer.data(lk + 1133);
    const auto *lk_1134 = buffer.data(lk + 1134);
    const auto *lk_1139 = buffer.data(lk + 1139);
    const auto *lk_1140 = buffer.data(lk + 1140);
    const auto *lk_1141 = buffer.data(lk + 1141);
    const auto *lk_1144 = buffer.data(lk + 1144);
    const auto *lk_1145 = buffer.data(lk + 1145);
    const auto *lk_1146 = buffer.data(lk + 1146);
    const auto *lk_1147 = buffer.data(lk + 1147);
    const auto *lk_1148 = buffer.data(lk + 1148);
    const auto *lk_1149 = buffer.data(lk + 1149);
    const auto *lk_1150 = buffer.data(lk + 1150);
    const auto *lk_1151 = buffer.data(lk + 1151);
    const auto *lk_1164 = buffer.data(lk + 1164);
    const auto *lk_1169 = buffer.data(lk + 1169);
    const auto *lk_1170 = buffer.data(lk + 1170);
    const auto *lk_1175 = buffer.data(lk + 1175);
    const auto *lk_1176 = buffer.data(lk + 1176);
    const auto *lk_1177 = buffer.data(lk + 1177);
    const auto *lk_1180 = buffer.data(lk + 1180);
    const auto *lk_1181 = buffer.data(lk + 1181);
    const auto *lk_1182 = buffer.data(lk + 1182);
    const auto *lk_1183 = buffer.data(lk + 1183);
    const auto *lk_1184 = buffer.data(lk + 1184);
    const auto *lk_1185 = buffer.data(lk + 1185);
    const auto *lk_1186 = buffer.data(lk + 1186);
    const auto *lk_1187 = buffer.data(lk + 1187);
    const auto *lk_1200 = buffer.data(lk + 1200);
    const auto *lk_1205 = buffer.data(lk + 1205);
    const auto *lk_1206 = buffer.data(lk + 1206);

    const auto *ll_1041 = buffer.data(ll + 1041);
    const auto *ll_1045 = buffer.data(ll + 1045);
    const auto *ll_1050 = buffer.data(ll + 1050);
    const auto *ll_1056 = buffer.data(ll + 1056);
    const auto *ll_1083 = buffer.data(ll + 1083);
    const auto *ll_1085 = buffer.data(ll + 1085);
    const auto *ll_1086 = buffer.data(ll + 1086);
    const auto *ll_1089 = buffer.data(ll + 1089);
    const auto *ll_1090 = buffer.data(ll + 1090);
    const auto *ll_1094 = buffer.data(ll + 1094);
    const auto *ll_1095 = buffer.data(ll + 1095);
    const auto *ll_1100 = buffer.data(ll + 1100);
    const auto *ll_1101 = buffer.data(ll + 1101);
    const auto *ll_1107 = buffer.data(ll + 1107);
    const auto *ll_1125 = buffer.data(ll + 1125);
    const auto *ll_1128 = buffer.data(ll + 1128);
    const auto *ll_1130 = buffer.data(ll + 1130);
    const auto *ll_1131 = buffer.data(ll + 1131);
    const auto *ll_1134 = buffer.data(ll + 1134);
    const auto *ll_1135 = buffer.data(ll + 1135);
    const auto *ll_1139 = buffer.data(ll + 1139);
    const auto *ll_1140 = buffer.data(ll + 1140);
    const auto *ll_1145 = buffer.data(ll + 1145);
    const auto *ll_1146 = buffer.data(ll + 1146);
    const auto *ll_1152 = buffer.data(ll + 1152);
    const auto *ll_1170 = buffer.data(ll + 1170);
    const auto *ll_1175 = buffer.data(ll + 1175);
    const auto *ll_1179 = buffer.data(ll + 1179);
    const auto *ll_1184 = buffer.data(ll + 1184);
    const auto *ll_1190 = buffer.data(ll + 1190);
    const auto *ll_1431 = buffer.data(ll + 1431);
    const auto *ll_1433 = buffer.data(ll + 1433);
    const auto *ll_1434 = buffer.data(ll + 1434);
    const auto *ll_1435 = buffer.data(ll + 1435);
    const auto *ll_1436 = buffer.data(ll + 1436);
    const auto *ll_1437 = buffer.data(ll + 1437);
    const auto *ll_1439 = buffer.data(ll + 1439);
    const auto *ll_1476 = buffer.data(ll + 1476);
    const auto *ll_1478 = buffer.data(ll + 1478);
    const auto *ll_1479 = buffer.data(ll + 1479);
    const auto *ll_1480 = buffer.data(ll + 1480);
    const auto *ll_1481 = buffer.data(ll + 1481);
    const auto *ll_1482 = buffer.data(ll + 1482);
    const auto *ll_1484 = buffer.data(ll + 1484);

    const auto *mi0_880 = buffer.data(mi0 + 880);
    const auto *mi0_885 = buffer.data(mi0 + 885);
    const auto *mi0_886 = buffer.data(mi0 + 886);
    const auto *mi0_891 = buffer.data(mi0 + 891);
    const auto *mi0_892 = buffer.data(mi0 + 892);
    const auto *mi0_893 = buffer.data(mi0 + 893);
    const auto *mi0_908 = buffer.data(mi0 + 908);
    const auto *mi0_913 = buffer.data(mi0 + 913);
    const auto *mi0_914 = buffer.data(mi0 + 914);
    const auto *mi0_919 = buffer.data(mi0 + 919);
    const auto *mi0_920 = buffer.data(mi0 + 920);
    const auto *mi0_921 = buffer.data(mi0 + 921);
    const auto *mi0_936 = buffer.data(mi0 + 936);
    const auto *mi0_941 = buffer.data(mi0 + 941);
    const auto *mi0_942 = buffer.data(mi0 + 942);

    const auto *mi1_880 = buffer.data(mi1 + 880);
    const auto *mi1_885 = buffer.data(mi1 + 885);
    const auto *mi1_886 = buffer.data(mi1 + 886);
    const auto *mi1_891 = buffer.data(mi1 + 891);
    const auto *mi1_892 = buffer.data(mi1 + 892);
    const auto *mi1_893 = buffer.data(mi1 + 893);
    const auto *mi1_908 = buffer.data(mi1 + 908);
    const auto *mi1_913 = buffer.data(mi1 + 913);
    const auto *mi1_914 = buffer.data(mi1 + 914);
    const auto *mi1_919 = buffer.data(mi1 + 919);
    const auto *mi1_920 = buffer.data(mi1 + 920);
    const auto *mi1_921 = buffer.data(mi1 + 921);
    const auto *mi1_936 = buffer.data(mi1 + 936);
    const auto *mi1_941 = buffer.data(mi1 + 941);
    const auto *mi1_942 = buffer.data(mi1 + 942);

    const auto *mk_1119 = buffer.data(mk + 1119);
    const auto *mk_1121 = buffer.data(mk + 1121);
    const auto *mk_1122 = buffer.data(mk + 1122);
    const auto *mk_1125 = buffer.data(mk + 1125);
    const auto *mk_1126 = buffer.data(mk + 1126);
    const auto *mk_1128 = buffer.data(mk + 1128);
    const auto *mk_1130 = buffer.data(mk + 1130);
    const auto *mk_1131 = buffer.data(mk + 1131);
    const auto *mk_1133 = buffer.data(mk + 1133);
    const auto *mk_1134 = buffer.data(mk + 1134);
    const auto *mk_1136 = buffer.data(mk + 1136);
    const auto *mk_1139 = buffer.data(mk + 1139);
    const auto *mk_1140 = buffer.data(mk + 1140);
    const auto *mk_1141 = buffer.data(mk + 1141);
    const auto *mk_1144 = buffer.data(mk + 1144);
    const auto *mk_1145 = buffer.data(mk + 1145);
    const auto *mk_1146 = buffer.data(mk + 1146);
    const auto *mk_1147 = buffer.data(mk + 1147);
    const auto *mk_1148 = buffer.data(mk + 1148);
    const auto *mk_1149 = buffer.data(mk + 1149);
    const auto *mk_1150 = buffer.data(mk + 1150);
    const auto *mk_1151 = buffer.data(mk + 1151);
    const auto *mk_1152 = buffer.data(mk + 1152);
    const auto *mk_1154 = buffer.data(mk + 1154);
    const auto *mk_1155 = buffer.data(mk + 1155);
    const auto *mk_1157 = buffer.data(mk + 1157);
    const auto *mk_1158 = buffer.data(mk + 1158);
    const auto *mk_1161 = buffer.data(mk + 1161);
    const auto *mk_1162 = buffer.data(mk + 1162);
    const auto *mk_1164 = buffer.data(mk + 1164);
    const auto *mk_1166 = buffer.data(mk + 1166);
    const auto *mk_1167 = buffer.data(mk + 1167);
    const auto *mk_1169 = buffer.data(mk + 1169);
    const auto *mk_1170 = buffer.data(mk + 1170);
    const auto *mk_1172 = buffer.data(mk + 1172);
    const auto *mk_1175 = buffer.data(mk + 1175);
    const auto *mk_1176 = buffer.data(mk + 1176);
    const auto *mk_1177 = buffer.data(mk + 1177);
    const auto *mk_1180 = buffer.data(mk + 1180);
    const auto *mk_1181 = buffer.data(mk + 1181);
    const auto *mk_1182 = buffer.data(mk + 1182);
    const auto *mk_1183 = buffer.data(mk + 1183);
    const auto *mk_1184 = buffer.data(mk + 1184);
    const auto *mk_1185 = buffer.data(mk + 1185);
    const auto *mk_1186 = buffer.data(mk + 1186);
    const auto *mk_1187 = buffer.data(mk + 1187);
    const auto *mk_1188 = buffer.data(mk + 1188);
    const auto *mk_1190 = buffer.data(mk + 1190);
    const auto *mk_1191 = buffer.data(mk + 1191);
    const auto *mk_1193 = buffer.data(mk + 1193);
    const auto *mk_1194 = buffer.data(mk + 1194);
    const auto *mk_1197 = buffer.data(mk + 1197);
    const auto *mk_1198 = buffer.data(mk + 1198);
    const auto *mk_1200 = buffer.data(mk + 1200);
    const auto *mk_1202 = buffer.data(mk + 1202);
    const auto *mk_1203 = buffer.data(mk + 1203);
    const auto *mk_1205 = buffer.data(mk + 1205);
    const auto *mk_1206 = buffer.data(mk + 1206);

#pragma omp simd aligned(t_1400, t_1401, t_1402, pa_y, pa_z, pb_z, kl0_726, kl0_815, kl1_726, \
                         kl1_815, lk_831, ll_1041, ll_1085, mk_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1400[k] = f_29 * kl0_815[k]
                    - f_30 * kl1_815[k]
                    + pa_y[k] * ll_1085[k];

        t_1401[k] = f_25 * kl0_726[k]
                    - f_26 * kl1_726[k]
                    + pa_z[k] * ll_1041[k];

        t_1402[k] = f_15 * lk_831[k]
                    + pb_z[k] * mk_1119[k];
    }

#pragma omp simd aligned(t_1403, t_1404, t_1405, pa_y, pa_z, pb_y, kl0_730, kl0_819, kl1_730, \
                         kl1_819, lk_869, ll_1045, ll_1089, mk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1403[k] = f_16 * lk_869[k]
                    + pb_y[k] * mk_1121[k];

        t_1404[k] = f_29 * kl0_819[k]
                    - f_30 * kl1_819[k]
                    + pa_y[k] * ll_1089[k];

        t_1405[k] = f_25 * kl0_730[k]
                    - f_26 * kl1_730[k]
                    + pa_z[k] * ll_1045[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pb_x, pb_y, pb_z, lk_834, lk_873, lk_1128, \
                         mi0_880, mi1_880, mk_1122, mk_1125, mk_1128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_15 * lk_834[k]
                    + pb_z[k] * mk_1122[k];

        t_1407[k] = f_14 * lk_1128[k]
                    + f_7 * mi0_880[k]
                    - f_8 * mi1_880[k]
                    + pb_x[k] * mk_1128[k];

        t_1408[k] = f_16 * lk_873[k]
                    + pb_y[k] * mk_1125[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pa_y, pa_z, pb_z, kl0_735, kl0_824, kl1_735, \
                         kl1_824, lk_838, ll_1050, ll_1094, mk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = f_29 * kl0_824[k]
                    - f_30 * kl1_824[k]
                    + pa_y[k] * ll_1094[k];

        t_1410[k] = f_25 * kl0_735[k]
                    - f_26 * kl1_735[k]
                    + pa_z[k] * ll_1050[k];

        t_1411[k] = f_15 * lk_838[k]
                    + pb_z[k] * mk_1126[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pb_x, pb_y, lk_878, lk_1133, lk_1134, \
                         mi0_885, mi0_886, mi1_885, mi1_886, mk_1130, mk_1133, \
                         mk_1134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = f_14 * lk_1133[k]
                    + f_5 * mi0_885[k]
                    - f_6 * mi1_885[k]
                    + pb_x[k] * mk_1133[k];

        t_1413[k] = f_14 * lk_1134[k]
                    + f_5 * mi0_886[k]
                    - f_6 * mi1_886[k]
                    + pb_x[k] * mk_1134[k];

        t_1414[k] = f_16 * lk_878[k]
                    + pb_y[k] * mk_1130[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pa_y, pa_z, pb_z, kl0_741, kl0_830, kl1_741, \
                         kl1_830, lk_843, ll_1056, ll_1100, mk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = f_29 * kl0_830[k]
                    - f_30 * kl1_830[k]
                    + pa_y[k] * ll_1100[k];

        t_1416[k] = f_25 * kl0_741[k]
                    - f_26 * kl1_741[k]
                    + pa_z[k] * ll_1056[k];

        t_1417[k] = f_15 * lk_843[k]
                    + pb_z[k] * mk_1131[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pb_x, lk_1139, lk_1140, lk_1141, mi0_891, \
                         mi0_892, mi0_893, mi1_891, mi1_892, mi1_893, mk_1139, mk_1140, \
                         mk_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = f_14 * lk_1139[k]
                    + f_3 * mi0_891[k]
                    - f_4 * mi1_891[k]
                    + pb_x[k] * mk_1139[k];

        t_1419[k] = f_14 * lk_1140[k]
                    + f_3 * mi0_892[k]
                    - f_4 * mi1_892[k]
                    + pb_x[k] * mk_1140[k];

        t_1420[k] = f_14 * lk_1141[k]
                    + f_3 * mi0_893[k]
                    - f_4 * mi1_893[k]
                    + pb_x[k] * mk_1141[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, t_1424, pa_y, pb_x, pb_y, kl0_837, kl1_837, \
                         lk_884, lk_1144, lk_1145, ll_1107, mk_1136, mk_1144, \
                         mk_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_16 * lk_884[k]
                    + pb_y[k] * mk_1136[k];

        t_1422[k] = f_29 * kl0_837[k]
                    - f_30 * kl1_837[k]
                    + pa_y[k] * ll_1107[k];

        t_1423[k] = f_14 * lk_1144[k]
                    + pb_x[k] * mk_1144[k];

        t_1424[k] = f_14 * lk_1145[k]
                    + pb_x[k] * mk_1145[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, t_1428, t_1429, pb_x, lk_1146, lk_1147, \
                         lk_1148, lk_1149, lk_1150, mk_1146, mk_1147, mk_1148, mk_1149, \
                         mk_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_14 * lk_1146[k]
                    + pb_x[k] * mk_1146[k];

        t_1426[k] = f_14 * lk_1147[k]
                    + pb_x[k] * mk_1147[k];

        t_1427[k] = f_14 * lk_1148[k]
                    + pb_x[k] * mk_1148[k];

        t_1428[k] = f_14 * lk_1149[k]
                    + pb_x[k] * mk_1149[k];

        t_1429[k] = f_14 * lk_1150[k]
                    + pb_x[k] * mk_1150[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, pa_x, pb_x, pb_z, kl0_1431, kl1_1431, lk_856, \
                         lk_1151, ll_1431, mk_1144, mk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_14 * lk_1151[k]
                    + pb_x[k] * mk_1151[k];

        t_1431[k] = f_20 * kl0_1431[k]
                    - f_21 * kl1_1431[k]
                    + pa_x[k] * ll_1431[k];

        t_1432[k] = f_15 * lk_856[k]
                    + pb_z[k] * mk_1144[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pa_x, kl0_1433, kl0_1434, kl0_1435, kl1_1433, \
                         kl1_1434, kl1_1435, ll_1433, ll_1434, \
                         ll_1435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_20 * kl0_1433[k]
                    - f_21 * kl1_1433[k]
                    + pa_x[k] * ll_1433[k];

        t_1434[k] = f_20 * kl0_1434[k]
                    - f_21 * kl1_1434[k]
                    + pa_x[k] * ll_1434[k];

        t_1435[k] = f_20 * kl0_1435[k]
                    - f_21 * kl1_1435[k]
                    + pa_x[k] * ll_1435[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pa_x, pb_y, kl0_1436, kl0_1437, kl1_1436, \
                         kl1_1437, lk_899, ll_1436, ll_1437, mk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_20 * kl0_1436[k]
                    - f_21 * kl1_1436[k]
                    + pa_x[k] * ll_1436[k];

        t_1437[k] = f_20 * kl0_1437[k]
                    - f_21 * kl1_1437[k]
                    + pa_x[k] * ll_1437[k];

        t_1438[k] = f_16 * lk_899[k]
                    + pb_y[k] * mk_1151[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, pa_x, pa_y, pb_y, kl0_855, kl0_1439, kl1_855, \
                         kl1_1439, lk_900, ll_1125, ll_1439, mk_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_20 * kl0_1439[k]
                    - f_21 * kl1_1439[k]
                    + pa_x[k] * ll_1439[k];

        t_1440[k] = f_25 * kl0_855[k]
                    - f_26 * kl1_855[k]
                    + pa_y[k] * ll_1125[k];

        t_1441[k] = f_15 * lk_900[k]
                    + pb_y[k] * mk_1152[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pa_z, pb_y, pb_z, kl0_768, kl1_768, lk_864, \
                         lk_902, ll_1083, mk_1152, mk_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_16 * lk_864[k]
                    + pb_z[k] * mk_1152[k];

        t_1443[k] = f_29 * kl0_768[k]
                    - f_30 * kl1_768[k]
                    + pa_z[k] * ll_1083[k];

        t_1444[k] = f_15 * lk_902[k]
                    + pb_y[k] * mk_1154[k];
    }

#pragma omp simd aligned(t_1445, t_1446, t_1447, pa_y, pa_z, pb_z, kl0_771, kl0_860, kl1_771, \
                         kl1_860, lk_867, ll_1086, ll_1130, mk_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = f_25 * kl0_860[k]
                    - f_26 * kl1_860[k]
                    + pa_y[k] * ll_1130[k];

        t_1446[k] = f_29 * kl0_771[k]
                    - f_30 * kl1_771[k]
                    + pa_z[k] * ll_1086[k];

        t_1447[k] = f_16 * lk_867[k]
                    + pb_z[k] * mk_1155[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pa_y, pa_z, pb_y, kl0_775, kl0_864, kl1_775, \
                         kl1_864, lk_905, ll_1090, ll_1134, mk_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_15 * lk_905[k]
                    + pb_y[k] * mk_1157[k];

        t_1449[k] = f_25 * kl0_864[k]
                    - f_26 * kl1_864[k]
                    + pa_y[k] * ll_1134[k];

        t_1450[k] = f_29 * kl0_775[k]
                    - f_30 * kl1_775[k]
                    + pa_z[k] * ll_1090[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pb_x, pb_y, pb_z, lk_870, lk_909, lk_1164, \
                         mi0_908, mi1_908, mk_1158, mk_1161, mk_1164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_16 * lk_870[k]
                    + pb_z[k] * mk_1158[k];

        t_1452[k] = f_14 * lk_1164[k]
                    + f_7 * mi0_908[k]
                    - f_8 * mi1_908[k]
                    + pb_x[k] * mk_1164[k];

        t_1453[k] = f_15 * lk_909[k]
                    + pb_y[k] * mk_1161[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pa_y, pa_z, pb_z, kl0_780, kl0_869, kl1_780, \
                         kl1_869, lk_874, ll_1095, ll_1139, mk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = f_25 * kl0_869[k]
                    - f_26 * kl1_869[k]
                    + pa_y[k] * ll_1139[k];

        t_1455[k] = f_29 * kl0_780[k]
                    - f_30 * kl1_780[k]
                    + pa_z[k] * ll_1095[k];

        t_1456[k] = f_16 * lk_874[k]
                    + pb_z[k] * mk_1162[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pb_x, pb_y, lk_914, lk_1169, lk_1170, \
                         mi0_913, mi0_914, mi1_913, mi1_914, mk_1166, mk_1169, \
                         mk_1170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = f_14 * lk_1169[k]
                    + f_5 * mi0_913[k]
                    - f_6 * mi1_913[k]
                    + pb_x[k] * mk_1169[k];

        t_1458[k] = f_14 * lk_1170[k]
                    + f_5 * mi0_914[k]
                    - f_6 * mi1_914[k]
                    + pb_x[k] * mk_1170[k];

        t_1459[k] = f_15 * lk_914[k]
                    + pb_y[k] * mk_1166[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pa_y, pa_z, pb_z, kl0_786, kl0_875, kl1_786, \
                         kl1_875, lk_879, ll_1101, ll_1145, mk_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = f_25 * kl0_875[k]
                    - f_26 * kl1_875[k]
                    + pa_y[k] * ll_1145[k];

        t_1461[k] = f_29 * kl0_786[k]
                    - f_30 * kl1_786[k]
                    + pa_z[k] * ll_1101[k];

        t_1462[k] = f_16 * lk_879[k]
                    + pb_z[k] * mk_1167[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pb_x, lk_1175, lk_1176, lk_1177, mi0_919, \
                         mi0_920, mi0_921, mi1_919, mi1_920, mi1_921, mk_1175, mk_1176, \
                         mk_1177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = f_14 * lk_1175[k]
                    + f_3 * mi0_919[k]
                    - f_4 * mi1_919[k]
                    + pb_x[k] * mk_1175[k];

        t_1464[k] = f_14 * lk_1176[k]
                    + f_3 * mi0_920[k]
                    - f_4 * mi1_920[k]
                    + pb_x[k] * mk_1176[k];

        t_1465[k] = f_14 * lk_1177[k]
                    + f_3 * mi0_921[k]
                    - f_4 * mi1_921[k]
                    + pb_x[k] * mk_1177[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pa_y, pb_x, pb_y, kl0_882, kl1_882, \
                         lk_920, lk_1180, lk_1181, ll_1152, mk_1172, mk_1180, \
                         mk_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_15 * lk_920[k]
                    + pb_y[k] * mk_1172[k];

        t_1467[k] = f_25 * kl0_882[k]
                    - f_26 * kl1_882[k]
                    + pa_y[k] * ll_1152[k];

        t_1468[k] = f_14 * lk_1180[k]
                    + pb_x[k] * mk_1180[k];

        t_1469[k] = f_14 * lk_1181[k]
                    + pb_x[k] * mk_1181[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, t_1474, pb_x, lk_1182, lk_1183, \
                         lk_1184, lk_1185, lk_1186, mk_1182, mk_1183, mk_1184, mk_1185, \
                         mk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_14 * lk_1182[k]
                    + pb_x[k] * mk_1182[k];

        t_1471[k] = f_14 * lk_1183[k]
                    + pb_x[k] * mk_1183[k];

        t_1472[k] = f_14 * lk_1184[k]
                    + pb_x[k] * mk_1184[k];

        t_1473[k] = f_14 * lk_1185[k]
                    + pb_x[k] * mk_1185[k];

        t_1474[k] = f_14 * lk_1186[k]
                    + pb_x[k] * mk_1186[k];
    }

#pragma omp simd aligned(t_1475, t_1476, t_1477, pa_x, pb_x, pb_z, kl0_1476, kl1_1476, lk_892, \
                         lk_1187, ll_1476, mk_1180, mk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1475[k] = f_14 * lk_1187[k]
                    + pb_x[k] * mk_1187[k];

        t_1476[k] = f_20 * kl0_1476[k]
                    - f_21 * kl1_1476[k]
                    + pa_x[k] * ll_1476[k];

        t_1477[k] = f_16 * lk_892[k]
                    + pb_z[k] * mk_1180[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pa_x, kl0_1478, kl0_1479, kl0_1480, kl1_1478, \
                         kl1_1479, kl1_1480, ll_1478, ll_1479, \
                         ll_1480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_20 * kl0_1478[k]
                    - f_21 * kl1_1478[k]
                    + pa_x[k] * ll_1478[k];

        t_1479[k] = f_20 * kl0_1479[k]
                    - f_21 * kl1_1479[k]
                    + pa_x[k] * ll_1479[k];

        t_1480[k] = f_20 * kl0_1480[k]
                    - f_21 * kl1_1480[k]
                    + pa_x[k] * ll_1480[k];
    }

#pragma omp simd aligned(t_1481, t_1482, t_1483, pa_x, pb_y, kl0_1481, kl0_1482, kl1_1481, \
                         kl1_1482, lk_935, ll_1481, ll_1482, mk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_20 * kl0_1481[k]
                    - f_21 * kl1_1481[k]
                    + pa_x[k] * ll_1481[k];

        t_1482[k] = f_20 * kl0_1482[k]
                    - f_21 * kl1_1482[k]
                    + pa_x[k] * ll_1482[k];

        t_1483[k] = f_15 * lk_935[k]
                    + pb_y[k] * mk_1187[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, pa_x, pa_y, pb_y, kl0_900, kl0_1484, kl1_900, \
                         kl1_1484, lk_936, ll_1170, ll_1484, mk_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_20 * kl0_1484[k]
                    - f_21 * kl1_1484[k]
                    + pa_x[k] * ll_1484[k];

        t_1485[k] = f_20 * kl0_900[k]
                    - f_21 * kl1_900[k]
                    + pa_y[k] * ll_1170[k];

        t_1486[k] = f_14 * lk_936[k]
                    + pb_y[k] * mk_1188[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pa_z, pb_y, pb_z, kl0_813, kl1_813, lk_900, \
                         lk_938, ll_1128, mk_1188, mk_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_17 * lk_900[k]
                    + pb_z[k] * mk_1188[k];

        t_1488[k] = f_31 * kl0_813[k]
                    - f_32 * kl1_813[k]
                    + pa_z[k] * ll_1128[k];

        t_1489[k] = f_14 * lk_938[k]
                    + pb_y[k] * mk_1190[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pa_y, pa_z, pb_z, kl0_816, kl0_905, kl1_816, \
                         kl1_905, lk_903, ll_1131, ll_1175, mk_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = f_20 * kl0_905[k]
                    - f_21 * kl1_905[k]
                    + pa_y[k] * ll_1175[k];

        t_1491[k] = f_31 * kl0_816[k]
                    - f_32 * kl1_816[k]
                    + pa_z[k] * ll_1131[k];

        t_1492[k] = f_17 * lk_903[k]
                    + pb_z[k] * mk_1191[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pa_y, pa_z, pb_y, kl0_820, kl0_909, kl1_820, \
                         kl1_909, lk_941, ll_1135, ll_1179, mk_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_14 * lk_941[k]
                    + pb_y[k] * mk_1193[k];

        t_1494[k] = f_20 * kl0_909[k]
                    - f_21 * kl1_909[k]
                    + pa_y[k] * ll_1179[k];

        t_1495[k] = f_31 * kl0_820[k]
                    - f_32 * kl1_820[k]
                    + pa_z[k] * ll_1135[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, pb_x, pb_y, pb_z, lk_906, lk_945, lk_1200, \
                         mi0_936, mi1_936, mk_1194, mk_1197, mk_1200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_17 * lk_906[k]
                    + pb_z[k] * mk_1194[k];

        t_1497[k] = f_14 * lk_1200[k]
                    + f_7 * mi0_936[k]
                    - f_8 * mi1_936[k]
                    + pb_x[k] * mk_1200[k];

        t_1498[k] = f_14 * lk_945[k]
                    + pb_y[k] * mk_1197[k];
    }

#pragma omp simd aligned(t_1499, t_1500, t_1501, pa_y, pa_z, pb_z, kl0_825, kl0_914, kl1_825, \
                         kl1_914, lk_910, ll_1140, ll_1184, mk_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1499[k] = f_20 * kl0_914[k]
                    - f_21 * kl1_914[k]
                    + pa_y[k] * ll_1184[k];

        t_1500[k] = f_31 * kl0_825[k]
                    - f_32 * kl1_825[k]
                    + pa_z[k] * ll_1140[k];

        t_1501[k] = f_17 * lk_910[k]
                    + pb_z[k] * mk_1198[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, pb_x, pb_y, lk_950, lk_1205, lk_1206, \
                         mi0_941, mi0_942, mi1_941, mi1_942, mk_1202, mk_1205, \
                         mk_1206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = f_14 * lk_1205[k]
                    + f_5 * mi0_941[k]
                    - f_6 * mi1_941[k]
                    + pb_x[k] * mk_1205[k];

        t_1503[k] = f_14 * lk_1206[k]
                    + f_5 * mi0_942[k]
                    - f_6 * mi1_942[k]
                    + pb_x[k] * mk_1206[k];

        t_1504[k] = f_14 * lk_950[k]
                    + pb_y[k] * mk_1202[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pa_y, pa_z, pb_z, kl0_831, kl0_920, kl1_831, \
                         kl1_920, lk_915, ll_1146, ll_1190, mk_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = f_20 * kl0_920[k]
                    - f_21 * kl1_920[k]
                    + pa_y[k] * ll_1190[k];

        t_1506[k] = f_31 * kl0_831[k]
                    - f_32 * kl1_831[k]
                    + pa_z[k] * ll_1146[k];

        t_1507[k] = f_17 * lk_915[k]
                    + pb_z[k] * mk_1203[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece12(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kl0,
                                             const size_t kl1, const size_t lk, const size_t ll,
                                             const size_t mi0, const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;
    const auto f_23 = 3.0 / alpha;
    const auto f_24 = 3.0 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_900 = buffer.data(kl0 + 900);
    const auto *kl0_927 = buffer.data(kl0 + 927);
    const auto *kl0_1521 = buffer.data(kl0 + 1521);
    const auto *kl0_1523 = buffer.data(kl0 + 1523);
    const auto *kl0_1524 = buffer.data(kl0 + 1524);
    const auto *kl0_1525 = buffer.data(kl0 + 1525);
    const auto *kl0_1526 = buffer.data(kl0 + 1526);
    const auto *kl0_1527 = buffer.data(kl0 + 1527);
    const auto *kl0_1529 = buffer.data(kl0 + 1529);
    const auto *kl0_1619 = buffer.data(kl0 + 1619);

    const auto *kl1_900 = buffer.data(kl1 + 900);
    const auto *kl1_927 = buffer.data(kl1 + 927);
    const auto *kl1_1521 = buffer.data(kl1 + 1521);
    const auto *kl1_1523 = buffer.data(kl1 + 1523);
    const auto *kl1_1524 = buffer.data(kl1 + 1524);
    const auto *kl1_1525 = buffer.data(kl1 + 1525);
    const auto *kl1_1526 = buffer.data(kl1 + 1526);
    const auto *kl1_1527 = buffer.data(kl1 + 1527);
    const auto *kl1_1529 = buffer.data(kl1 + 1529);
    const auto *kl1_1619 = buffer.data(kl1 + 1619);

    const auto *lk_928 = buffer.data(lk + 928);
    const auto *lk_939 = buffer.data(lk + 939);
    const auto *lk_942 = buffer.data(lk + 942);
    const auto *lk_946 = buffer.data(lk + 946);
    const auto *lk_951 = buffer.data(lk + 951);
    const auto *lk_956 = buffer.data(lk + 956);
    const auto *lk_964 = buffer.data(lk + 964);
    const auto *lk_971 = buffer.data(lk + 971);
    const auto *lk_972 = buffer.data(lk + 972);
    const auto *lk_973 = buffer.data(lk + 973);
    const auto *lk_974 = buffer.data(lk + 974);
    const auto *lk_975 = buffer.data(lk + 975);
    const auto *lk_977 = buffer.data(lk + 977);
    const auto *lk_978 = buffer.data(lk + 978);
    const auto *lk_980 = buffer.data(lk + 980);
    const auto *lk_981 = buffer.data(lk + 981);
    const auto *lk_982 = buffer.data(lk + 982);
    const auto *lk_984 = buffer.data(lk + 984);
    const auto *lk_985 = buffer.data(lk + 985);
    const auto *lk_986 = buffer.data(lk + 986);
    const auto *lk_987 = buffer.data(lk + 987);
    const auto *lk_989 = buffer.data(lk + 989);
    const auto *lk_990 = buffer.data(lk + 990);
    const auto *lk_991 = buffer.data(lk + 991);
    const auto *lk_992 = buffer.data(lk + 992);
    const auto *lk_1000 = buffer.data(lk + 1000);
    const auto *lk_1002 = buffer.data(lk + 1002);
    const auto *lk_1003 = buffer.data(lk + 1003);
    const auto *lk_1004 = buffer.data(lk + 1004);
    const auto *lk_1005 = buffer.data(lk + 1005);
    const auto *lk_1006 = buffer.data(lk + 1006);
    const auto *lk_1007 = buffer.data(lk + 1007);
    const auto *lk_1008 = buffer.data(lk + 1008);
    const auto *lk_1013 = buffer.data(lk + 1013);
    const auto *lk_1017 = buffer.data(lk + 1017);
    const auto *lk_1022 = buffer.data(lk + 1022);
    const auto *lk_1211 = buffer.data(lk + 1211);
    const auto *lk_1212 = buffer.data(lk + 1212);
    const auto *lk_1213 = buffer.data(lk + 1213);
    const auto *lk_1216 = buffer.data(lk + 1216);
    const auto *lk_1217 = buffer.data(lk + 1217);
    const auto *lk_1218 = buffer.data(lk + 1218);
    const auto *lk_1219 = buffer.data(lk + 1219);
    const auto *lk_1220 = buffer.data(lk + 1220);
    const auto *lk_1221 = buffer.data(lk + 1221);
    const auto *lk_1222 = buffer.data(lk + 1222);
    const auto *lk_1223 = buffer.data(lk + 1223);
    const auto *lk_1252 = buffer.data(lk + 1252);
    const auto *lk_1253 = buffer.data(lk + 1253);
    const auto *lk_1254 = buffer.data(lk + 1254);
    const auto *lk_1255 = buffer.data(lk + 1255);
    const auto *lk_1256 = buffer.data(lk + 1256);
    const auto *lk_1257 = buffer.data(lk + 1257);
    const auto *lk_1258 = buffer.data(lk + 1258);
    const auto *lk_1265 = buffer.data(lk + 1265);
    const auto *lk_1269 = buffer.data(lk + 1269);
    const auto *lk_1274 = buffer.data(lk + 1274);
    const auto *lk_1280 = buffer.data(lk + 1280);
    const auto *lk_1287 = buffer.data(lk + 1287);
    const auto *lk_1288 = buffer.data(lk + 1288);
    const auto *lk_1289 = buffer.data(lk + 1289);
    const auto *lk_1290 = buffer.data(lk + 1290);
    const auto *lk_1291 = buffer.data(lk + 1291);
    const auto *lk_1292 = buffer.data(lk + 1292);
    const auto *lk_1293 = buffer.data(lk + 1293);
    const auto *lk_1295 = buffer.data(lk + 1295);
    const auto *lk_1296 = buffer.data(lk + 1296);
    const auto *lk_1299 = buffer.data(lk + 1299);
    const auto *lk_1301 = buffer.data(lk + 1301);
    const auto *lk_1302 = buffer.data(lk + 1302);
    const auto *lk_1305 = buffer.data(lk + 1305);
    const auto *lk_1306 = buffer.data(lk + 1306);
    const auto *lk_1308 = buffer.data(lk + 1308);
    const auto *lk_1310 = buffer.data(lk + 1310);
    const auto *lk_1311 = buffer.data(lk + 1311);
    const auto *lk_1313 = buffer.data(lk + 1313);
    const auto *lk_1314 = buffer.data(lk + 1314);

    const auto *ll_1197 = buffer.data(ll + 1197);
    const auto *ll_1215 = buffer.data(ll + 1215);
    const auto *ll_1217 = buffer.data(ll + 1217);
    const auto *ll_1218 = buffer.data(ll + 1218);
    const auto *ll_1220 = buffer.data(ll + 1220);
    const auto *ll_1221 = buffer.data(ll + 1221);
    const auto *ll_1224 = buffer.data(ll + 1224);
    const auto *ll_1225 = buffer.data(ll + 1225);
    const auto *ll_1227 = buffer.data(ll + 1227);
    const auto *ll_1229 = buffer.data(ll + 1229);
    const auto *ll_1230 = buffer.data(ll + 1230);
    const auto *ll_1232 = buffer.data(ll + 1232);
    const auto *ll_1233 = buffer.data(ll + 1233);
    const auto *ll_1235 = buffer.data(ll + 1235);
    const auto *ll_1236 = buffer.data(ll + 1236);
    const auto *ll_1238 = buffer.data(ll + 1238);
    const auto *ll_1239 = buffer.data(ll + 1239);
    const auto *ll_1240 = buffer.data(ll + 1240);
    const auto *ll_1242 = buffer.data(ll + 1242);
    const auto *ll_1250 = buffer.data(ll + 1250);
    const auto *ll_1251 = buffer.data(ll + 1251);
    const auto *ll_1253 = buffer.data(ll + 1253);
    const auto *ll_1254 = buffer.data(ll + 1254);
    const auto *ll_1255 = buffer.data(ll + 1255);
    const auto *ll_1256 = buffer.data(ll + 1256);
    const auto *ll_1257 = buffer.data(ll + 1257);
    const auto *ll_1259 = buffer.data(ll + 1259);
    const auto *ll_1521 = buffer.data(ll + 1521);
    const auto *ll_1523 = buffer.data(ll + 1523);
    const auto *ll_1524 = buffer.data(ll + 1524);
    const auto *ll_1525 = buffer.data(ll + 1525);
    const auto *ll_1526 = buffer.data(ll + 1526);
    const auto *ll_1527 = buffer.data(ll + 1527);
    const auto *ll_1529 = buffer.data(ll + 1529);
    const auto *ll_1619 = buffer.data(ll + 1619);
    const auto *ll_1620 = buffer.data(ll + 1620);
    const auto *ll_1623 = buffer.data(ll + 1623);
    const auto *ll_1625 = buffer.data(ll + 1625);
    const auto *ll_1626 = buffer.data(ll + 1626);
    const auto *ll_1629 = buffer.data(ll + 1629);
    const auto *ll_1630 = buffer.data(ll + 1630);
    const auto *ll_1632 = buffer.data(ll + 1632);
    const auto *ll_1634 = buffer.data(ll + 1634);
    const auto *ll_1635 = buffer.data(ll + 1635);
    const auto *ll_1637 = buffer.data(ll + 1637);
    const auto *ll_1638 = buffer.data(ll + 1638);

    const auto *mi0_947 = buffer.data(mi0 + 947);
    const auto *mi0_948 = buffer.data(mi0 + 948);
    const auto *mi0_949 = buffer.data(mi0 + 949);
    const auto *mi0_980 = buffer.data(mi0 + 980);
    const auto *mi0_981 = buffer.data(mi0 + 981);
    const auto *mi0_983 = buffer.data(mi0 + 983);
    const auto *mi0_985 = buffer.data(mi0 + 985);
    const auto *mi0_986 = buffer.data(mi0 + 986);
    const auto *mi0_988 = buffer.data(mi0 + 988);
    const auto *mi0_989 = buffer.data(mi0 + 989);
    const auto *mi0_990 = buffer.data(mi0 + 990);
    const auto *mi0_992 = buffer.data(mi0 + 992);
    const auto *mi0_993 = buffer.data(mi0 + 993);
    const auto *mi0_994 = buffer.data(mi0 + 994);
    const auto *mi0_1000 = buffer.data(mi0 + 1000);
    const auto *mi0_1001 = buffer.data(mi0 + 1001);
    const auto *mi0_1003 = buffer.data(mi0 + 1003);
    const auto *mi0_1004 = buffer.data(mi0 + 1004);
    const auto *mi0_1005 = buffer.data(mi0 + 1005);
    const auto *mi0_1006 = buffer.data(mi0 + 1006);
    const auto *mi0_1007 = buffer.data(mi0 + 1007);

    const auto *mi1_947 = buffer.data(mi1 + 947);
    const auto *mi1_948 = buffer.data(mi1 + 948);
    const auto *mi1_949 = buffer.data(mi1 + 949);
    const auto *mi1_980 = buffer.data(mi1 + 980);
    const auto *mi1_981 = buffer.data(mi1 + 981);
    const auto *mi1_983 = buffer.data(mi1 + 983);
    const auto *mi1_985 = buffer.data(mi1 + 985);
    const auto *mi1_986 = buffer.data(mi1 + 986);
    const auto *mi1_988 = buffer.data(mi1 + 988);
    const auto *mi1_989 = buffer.data(mi1 + 989);
    const auto *mi1_990 = buffer.data(mi1 + 990);
    const auto *mi1_992 = buffer.data(mi1 + 992);
    const auto *mi1_993 = buffer.data(mi1 + 993);
    const auto *mi1_994 = buffer.data(mi1 + 994);
    const auto *mi1_1000 = buffer.data(mi1 + 1000);
    const auto *mi1_1001 = buffer.data(mi1 + 1001);
    const auto *mi1_1003 = buffer.data(mi1 + 1003);
    const auto *mi1_1004 = buffer.data(mi1 + 1004);
    const auto *mi1_1005 = buffer.data(mi1 + 1005);
    const auto *mi1_1006 = buffer.data(mi1 + 1006);
    const auto *mi1_1007 = buffer.data(mi1 + 1007);

    const auto *mk_1208 = buffer.data(mk + 1208);
    const auto *mk_1211 = buffer.data(mk + 1211);
    const auto *mk_1212 = buffer.data(mk + 1212);
    const auto *mk_1213 = buffer.data(mk + 1213);
    const auto *mk_1216 = buffer.data(mk + 1216);
    const auto *mk_1217 = buffer.data(mk + 1217);
    const auto *mk_1218 = buffer.data(mk + 1218);
    const auto *mk_1219 = buffer.data(mk + 1219);
    const auto *mk_1220 = buffer.data(mk + 1220);
    const auto *mk_1221 = buffer.data(mk + 1221);
    const auto *mk_1222 = buffer.data(mk + 1222);
    const auto *mk_1223 = buffer.data(mk + 1223);
    const auto *mk_1224 = buffer.data(mk + 1224);
    const auto *mk_1226 = buffer.data(mk + 1226);
    const auto *mk_1227 = buffer.data(mk + 1227);
    const auto *mk_1229 = buffer.data(mk + 1229);
    const auto *mk_1230 = buffer.data(mk + 1230);
    const auto *mk_1233 = buffer.data(mk + 1233);
    const auto *mk_1234 = buffer.data(mk + 1234);
    const auto *mk_1238 = buffer.data(mk + 1238);
    const auto *mk_1239 = buffer.data(mk + 1239);
    const auto *mk_1244 = buffer.data(mk + 1244);
    const auto *mk_1252 = buffer.data(mk + 1252);
    const auto *mk_1253 = buffer.data(mk + 1253);
    const auto *mk_1254 = buffer.data(mk + 1254);
    const auto *mk_1255 = buffer.data(mk + 1255);
    const auto *mk_1256 = buffer.data(mk + 1256);
    const auto *mk_1257 = buffer.data(mk + 1257);
    const auto *mk_1258 = buffer.data(mk + 1258);
    const auto *mk_1259 = buffer.data(mk + 1259);
    const auto *mk_1260 = buffer.data(mk + 1260);
    const auto *mk_1261 = buffer.data(mk + 1261);
    const auto *mk_1262 = buffer.data(mk + 1262);
    const auto *mk_1263 = buffer.data(mk + 1263);
    const auto *mk_1265 = buffer.data(mk + 1265);
    const auto *mk_1266 = buffer.data(mk + 1266);
    const auto *mk_1268 = buffer.data(mk + 1268);
    const auto *mk_1269 = buffer.data(mk + 1269);
    const auto *mk_1270 = buffer.data(mk + 1270);
    const auto *mk_1272 = buffer.data(mk + 1272);
    const auto *mk_1273 = buffer.data(mk + 1273);
    const auto *mk_1274 = buffer.data(mk + 1274);
    const auto *mk_1275 = buffer.data(mk + 1275);
    const auto *mk_1277 = buffer.data(mk + 1277);
    const auto *mk_1278 = buffer.data(mk + 1278);
    const auto *mk_1279 = buffer.data(mk + 1279);
    const auto *mk_1280 = buffer.data(mk + 1280);
    const auto *mk_1287 = buffer.data(mk + 1287);
    const auto *mk_1288 = buffer.data(mk + 1288);
    const auto *mk_1289 = buffer.data(mk + 1289);
    const auto *mk_1290 = buffer.data(mk + 1290);
    const auto *mk_1291 = buffer.data(mk + 1291);
    const auto *mk_1292 = buffer.data(mk + 1292);
    const auto *mk_1293 = buffer.data(mk + 1293);
    const auto *mk_1294 = buffer.data(mk + 1294);
    const auto *mk_1295 = buffer.data(mk + 1295);
    const auto *mk_1296 = buffer.data(mk + 1296);
    const auto *mk_1297 = buffer.data(mk + 1297);
    const auto *mk_1299 = buffer.data(mk + 1299);
    const auto *mk_1301 = buffer.data(mk + 1301);
    const auto *mk_1302 = buffer.data(mk + 1302);
    const auto *mk_1305 = buffer.data(mk + 1305);
    const auto *mk_1306 = buffer.data(mk + 1306);
    const auto *mk_1310 = buffer.data(mk + 1310);

#pragma omp simd aligned(t_1508, t_1509, t_1510, pb_x, lk_1211, lk_1212, lk_1213, mi0_947, \
                         mi0_948, mi0_949, mi1_947, mi1_948, mi1_949, mk_1211, mk_1212, \
                         mk_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = f_14 * lk_1211[k]
                    + f_3 * mi0_947[k]
                    - f_4 * mi1_947[k]
                    + pb_x[k] * mk_1211[k];

        t_1509[k] = f_14 * lk_1212[k]
                    + f_3 * mi0_948[k]
                    - f_4 * mi1_948[k]
                    + pb_x[k] * mk_1212[k];

        t_1510[k] = f_14 * lk_1213[k]
                    + f_3 * mi0_949[k]
                    - f_4 * mi1_949[k]
                    + pb_x[k] * mk_1213[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pa_y, pb_x, pb_y, kl0_927, kl1_927, \
                         lk_956, lk_1216, lk_1217, ll_1197, mk_1208, mk_1216, \
                         mk_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = f_14 * lk_956[k]
                    + pb_y[k] * mk_1208[k];

        t_1512[k] = f_20 * kl0_927[k]
                    - f_21 * kl1_927[k]
                    + pa_y[k] * ll_1197[k];

        t_1513[k] = f_14 * lk_1216[k]
                    + pb_x[k] * mk_1216[k];

        t_1514[k] = f_14 * lk_1217[k]
                    + pb_x[k] * mk_1217[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, t_1518, t_1519, pb_x, lk_1218, lk_1219, \
                         lk_1220, lk_1221, lk_1222, mk_1218, mk_1219, mk_1220, mk_1221, \
                         mk_1222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_14 * lk_1218[k]
                    + pb_x[k] * mk_1218[k];

        t_1516[k] = f_14 * lk_1219[k]
                    + pb_x[k] * mk_1219[k];

        t_1517[k] = f_14 * lk_1220[k]
                    + pb_x[k] * mk_1220[k];

        t_1518[k] = f_14 * lk_1221[k]
                    + pb_x[k] * mk_1221[k];

        t_1519[k] = f_14 * lk_1222[k]
                    + pb_x[k] * mk_1222[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, pa_x, pb_x, pb_z, kl0_1521, kl1_1521, lk_928, \
                         lk_1223, ll_1521, mk_1216, mk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_14 * lk_1223[k]
                    + pb_x[k] * mk_1223[k];

        t_1521[k] = f_20 * kl0_1521[k]
                    - f_21 * kl1_1521[k]
                    + pa_x[k] * ll_1521[k];

        t_1522[k] = f_17 * lk_928[k]
                    + pb_z[k] * mk_1216[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pa_x, kl0_1523, kl0_1524, kl0_1525, kl1_1523, \
                         kl1_1524, kl1_1525, ll_1523, ll_1524, \
                         ll_1525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_20 * kl0_1523[k]
                    - f_21 * kl1_1523[k]
                    + pa_x[k] * ll_1523[k];

        t_1524[k] = f_20 * kl0_1524[k]
                    - f_21 * kl1_1524[k]
                    + pa_x[k] * ll_1524[k];

        t_1525[k] = f_20 * kl0_1525[k]
                    - f_21 * kl1_1525[k]
                    + pa_x[k] * ll_1525[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pa_x, pb_y, kl0_1526, kl0_1527, kl1_1526, \
                         kl1_1527, lk_971, ll_1526, ll_1527, mk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_20 * kl0_1526[k]
                    - f_21 * kl1_1526[k]
                    + pa_x[k] * ll_1526[k];

        t_1527[k] = f_20 * kl0_1527[k]
                    - f_21 * kl1_1527[k]
                    + pa_x[k] * ll_1527[k];

        t_1528[k] = f_14 * lk_971[k]
                    + pb_y[k] * mk_1223[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, pa_x, pa_y, pb_y, kl0_1529, kl1_1529, \
                         lk_972, ll_1215, ll_1217, ll_1529, mk_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_20 * kl0_1529[k]
                    - f_21 * kl1_1529[k]
                    + pa_x[k] * ll_1529[k];

        t_1530[k] = pa_y[k] * ll_1215[k];

        t_1531[k] = f_13 * lk_972[k]
                    + pb_y[k] * mk_1224[k];

        t_1532[k] = pa_y[k] * ll_1217[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, t_1536, pa_y, pb_y, lk_973, lk_974, lk_975, \
                         ll_1218, ll_1220, ll_1221, mk_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = f_14 * lk_973[k]
                    + pa_y[k] * ll_1218[k];

        t_1534[k] = f_13 * lk_974[k]
                    + pb_y[k] * mk_1226[k];

        t_1535[k] = pa_y[k] * ll_1220[k];

        t_1536[k] = f_15 * lk_975[k]
                    + pa_y[k] * ll_1221[k];
    }

#pragma omp simd aligned(t_1537, t_1538, t_1539, t_1540, pa_y, pb_y, pb_z, lk_939, lk_977, \
                         lk_978, ll_1224, ll_1225, mk_1227, mk_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1537[k] = f_18 * lk_939[k]
                    + pb_z[k] * mk_1227[k];

        t_1538[k] = f_13 * lk_977[k]
                    + pb_y[k] * mk_1229[k];

        t_1539[k] = pa_y[k] * ll_1224[k];

        t_1540[k] = f_16 * lk_978[k]
                    + pa_y[k] * ll_1225[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, t_1544, pa_y, pb_y, pb_z, lk_942, lk_980, \
                         lk_981, ll_1227, ll_1229, mk_1230, mk_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_18 * lk_942[k]
                    + pb_z[k] * mk_1230[k];

        t_1542[k] = f_14 * lk_980[k]
                    + pa_y[k] * ll_1227[k];

        t_1543[k] = f_13 * lk_981[k]
                    + pb_y[k] * mk_1233[k];

        t_1544[k] = pa_y[k] * ll_1229[k];
    }

#pragma omp simd aligned(t_1545, t_1546, t_1547, t_1548, pa_y, pb_z, lk_946, lk_982, lk_984, \
                         lk_985, ll_1230, ll_1232, ll_1233, mk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1545[k] = f_17 * lk_982[k]
                    + pa_y[k] * ll_1230[k];

        t_1546[k] = f_18 * lk_946[k]
                    + pb_z[k] * mk_1234[k];

        t_1547[k] = f_15 * lk_984[k]
                    + pa_y[k] * ll_1232[k];

        t_1548[k] = f_14 * lk_985[k]
                    + pa_y[k] * ll_1233[k];
    }

#pragma omp simd aligned(t_1549, t_1550, t_1551, t_1552, pa_y, pb_y, pb_z, lk_951, lk_986, \
                         lk_987, ll_1235, ll_1236, mk_1238, mk_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1549[k] = f_13 * lk_986[k]
                    + pb_y[k] * mk_1238[k];

        t_1550[k] = pa_y[k] * ll_1235[k];

        t_1551[k] = f_18 * lk_987[k]
                    + pa_y[k] * ll_1236[k];

        t_1552[k] = f_18 * lk_951[k]
                    + pb_z[k] * mk_1239[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, t_1556, t_1557, pa_y, pb_y, lk_989, lk_990, \
                         lk_991, lk_992, ll_1238, ll_1239, ll_1240, ll_1242, \
                         mk_1244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = f_16 * lk_989[k]
                    + pa_y[k] * ll_1238[k];

        t_1554[k] = f_15 * lk_990[k]
                    + pa_y[k] * ll_1239[k];

        t_1555[k] = f_14 * lk_991[k]
                    + pa_y[k] * ll_1240[k];

        t_1556[k] = f_13 * lk_992[k]
                    + pb_y[k] * mk_1244[k];

        t_1557[k] = pa_y[k] * ll_1242[k];
    }

#pragma omp simd aligned(t_1558, t_1559, t_1560, t_1561, t_1562, pb_x, lk_1252, lk_1253, \
                         lk_1254, lk_1255, lk_1256, mk_1252, mk_1253, mk_1254, mk_1255, \
                         mk_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1558[k] = f_14 * lk_1252[k]
                    + pb_x[k] * mk_1252[k];

        t_1559[k] = f_14 * lk_1253[k]
                    + pb_x[k] * mk_1253[k];

        t_1560[k] = f_14 * lk_1254[k]
                    + pb_x[k] * mk_1254[k];

        t_1561[k] = f_14 * lk_1255[k]
                    + pb_x[k] * mk_1255[k];

        t_1562[k] = f_14 * lk_1256[k]
                    + pb_x[k] * mk_1256[k];
    }

#pragma omp simd aligned(t_1563, t_1564, t_1565, t_1566, pa_y, pb_x, lk_1000, lk_1257, \
                         lk_1258, ll_1250, ll_1251, mk_1257, mk_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1563[k] = f_14 * lk_1257[k]
                    + pb_x[k] * mk_1257[k];

        t_1564[k] = f_14 * lk_1258[k]
                    + pb_x[k] * mk_1258[k];

        t_1565[k] = pa_y[k] * ll_1250[k];

        t_1566[k] = f_19 * lk_1000[k]
                    + pa_y[k] * ll_1251[k];
    }

#pragma omp simd aligned(t_1567, t_1568, t_1569, t_1570, pa_y, pb_z, lk_964, lk_1002, lk_1003, \
                         lk_1004, ll_1253, ll_1254, ll_1255, mk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1567[k] = f_18 * lk_964[k]
                    + pb_z[k] * mk_1252[k];

        t_1568[k] = f_18 * lk_1002[k]
                    + pa_y[k] * ll_1253[k];

        t_1569[k] = f_17 * lk_1003[k]
                    + pa_y[k] * ll_1254[k];

        t_1570[k] = f_16 * lk_1004[k]
                    + pa_y[k] * ll_1255[k];
    }

#pragma omp simd aligned(t_1571, t_1572, t_1573, t_1574, pa_y, pb_y, lk_1005, lk_1006, \
                         lk_1007, ll_1256, ll_1257, ll_1259, mk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1571[k] = f_15 * lk_1005[k]
                    + pa_y[k] * ll_1256[k];

        t_1572[k] = f_14 * lk_1006[k]
                    + pa_y[k] * ll_1257[k];

        t_1573[k] = f_13 * lk_1007[k]
                    + pb_y[k] * mk_1259[k];

        t_1574[k] = pa_y[k] * ll_1259[k];
    }

#pragma omp simd aligned(t_1575, t_1576, t_1577, t_1578, pa_z, pb_y, pb_z, kl0_900, kl1_900, \
                         lk_972, ll_1215, mi0_980, mi1_980, mk_1260, \
                         mk_1261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1575[k] = f_23 * kl0_900[k]
                    - f_24 * kl1_900[k]
                    + pa_z[k] * ll_1215[k];

        t_1576[k] = pb_y[k] * mk_1260[k];

        t_1577[k] = f_22 * lk_972[k]
                    + pb_z[k] * mk_1260[k];

        t_1578[k] = f_3 * mi0_980[k]
                    - f_4 * mi1_980[k]
                    + pb_y[k] * mk_1261[k];
    }

#pragma omp simd aligned(t_1579, t_1580, t_1581, t_1582, pb_x, pb_y, pb_z, lk_975, lk_1265, \
                         mi0_981, mi0_985, mi1_981, mi1_985, mk_1262, mk_1263, \
                         mk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1579[k] = pb_y[k] * mk_1262[k];

        t_1580[k] = f_14 * lk_1265[k]
                    + f_11 * mi0_985[k]
                    - f_12 * mi1_985[k]
                    + pb_x[k] * mk_1265[k];

        t_1581[k] = f_5 * mi0_981[k]
                    - f_6 * mi1_981[k]
                    + pb_y[k] * mk_1263[k];

        t_1582[k] = f_22 * lk_975[k]
                    + pb_z[k] * mk_1263[k];
    }

#pragma omp simd aligned(t_1583, t_1584, t_1585, t_1586, pb_x, pb_y, pb_z, lk_978, lk_1269, \
                         mi0_983, mi0_989, mi1_983, mi1_989, mk_1265, mk_1266, \
                         mk_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1583[k] = pb_y[k] * mk_1265[k];

        t_1584[k] = f_14 * lk_1269[k]
                    + f_9 * mi0_989[k]
                    - f_10 * mi1_989[k]
                    + pb_x[k] * mk_1269[k];

        t_1585[k] = f_7 * mi0_983[k]
                    - f_8 * mi1_983[k]
                    + pb_y[k] * mk_1266[k];

        t_1586[k] = f_22 * lk_978[k]
                    + pb_z[k] * mk_1266[k];
    }

#pragma omp simd aligned(t_1587, t_1588, t_1589, pb_x, pb_y, lk_1274, mi0_985, mi0_994, \
                         mi1_985, mi1_994, mk_1268, mk_1269, mk_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1587[k] = f_3 * mi0_985[k]
                    - f_4 * mi1_985[k]
                    + pb_y[k] * mk_1268[k];

        t_1588[k] = pb_y[k] * mk_1269[k];

        t_1589[k] = f_14 * lk_1274[k]
                    + f_7 * mi0_994[k]
                    - f_8 * mi1_994[k]
                    + pb_x[k] * mk_1274[k];
    }

#pragma omp simd aligned(t_1590, t_1591, t_1592, t_1593, pb_y, pb_z, lk_982, mi0_986, mi0_988, \
                         mi0_989, mi1_986, mi1_988, mi1_989, mk_1270, mk_1272, \
                         mk_1273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1590[k] = f_9 * mi0_986[k]
                    - f_10 * mi1_986[k]
                    + pb_y[k] * mk_1270[k];

        t_1591[k] = f_22 * lk_982[k]
                    + pb_z[k] * mk_1270[k];

        t_1592[k] = f_5 * mi0_988[k]
                    - f_6 * mi1_988[k]
                    + pb_y[k] * mk_1272[k];

        t_1593[k] = f_3 * mi0_989[k]
                    - f_4 * mi1_989[k]
                    + pb_y[k] * mk_1273[k];
    }

#pragma omp simd aligned(t_1594, t_1595, t_1596, t_1597, pb_x, pb_y, pb_z, lk_987, lk_1280, \
                         mi0_990, mi0_1000, mi1_990, mi1_1000, mk_1274, mk_1275, \
                         mk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1594[k] = pb_y[k] * mk_1274[k];

        t_1595[k] = f_14 * lk_1280[k]
                    + f_5 * mi0_1000[k]
                    - f_6 * mi1_1000[k]
                    + pb_x[k] * mk_1280[k];

        t_1596[k] = f_11 * mi0_990[k]
                    - f_12 * mi1_990[k]
                    + pb_y[k] * mk_1275[k];

        t_1597[k] = f_22 * lk_987[k]
                    + pb_z[k] * mk_1275[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, t_1601, pb_y, mi0_992, mi0_993, mi0_994, \
                         mi1_992, mi1_993, mi1_994, mk_1277, mk_1278, mk_1279, \
                         mk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = f_7 * mi0_992[k]
                    - f_8 * mi1_992[k]
                    + pb_y[k] * mk_1277[k];

        t_1599[k] = f_5 * mi0_993[k]
                    - f_6 * mi1_993[k]
                    + pb_y[k] * mk_1278[k];

        t_1600[k] = f_3 * mi0_994[k]
                    - f_4 * mi1_994[k]
                    + pb_y[k] * mk_1279[k];

        t_1601[k] = pb_y[k] * mk_1280[k];
    }

#pragma omp simd aligned(t_1602, t_1603, t_1604, t_1605, pb_x, lk_1287, lk_1288, lk_1289, \
                         lk_1290, mi0_1007, mi1_1007, mk_1287, mk_1288, mk_1289, \
                         mk_1290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1602[k] = f_14 * lk_1287[k]
                    + f_3 * mi0_1007[k]
                    - f_4 * mi1_1007[k]
                    + pb_x[k] * mk_1287[k];

        t_1603[k] = f_14 * lk_1288[k]
                    + pb_x[k] * mk_1288[k];

        t_1604[k] = f_14 * lk_1289[k]
                    + pb_x[k] * mk_1289[k];

        t_1605[k] = f_14 * lk_1290[k]
                    + pb_x[k] * mk_1290[k];
    }

#pragma omp simd aligned(t_1606, t_1607, t_1608, t_1609, t_1610, pb_x, pb_y, lk_1291, lk_1292, \
                         lk_1293, lk_1295, mk_1287, mk_1291, mk_1292, mk_1293, \
                         mk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1606[k] = f_14 * lk_1291[k]
                    + pb_x[k] * mk_1291[k];

        t_1607[k] = f_14 * lk_1292[k]
                    + pb_x[k] * mk_1292[k];

        t_1608[k] = f_14 * lk_1293[k]
                    + pb_x[k] * mk_1293[k];

        t_1609[k] = pb_y[k] * mk_1287[k];

        t_1610[k] = f_14 * lk_1295[k]
                    + pb_x[k] * mk_1295[k];
    }

#pragma omp simd aligned(t_1611, t_1612, t_1613, t_1614, pb_y, pb_z, lk_1000, mi0_1001, \
                         mi0_1003, mi0_1004, mi1_1001, mi1_1003, mi1_1004, mk_1288, mk_1290, \
                         mk_1291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1611[k] = f_1 * mi0_1001[k]
                    - f_2 * mi1_1001[k]
                    + pb_y[k] * mk_1288[k];

        t_1612[k] = f_22 * lk_1000[k]
                    + pb_z[k] * mk_1288[k];

        t_1613[k] = f_11 * mi0_1003[k]
                    - f_12 * mi1_1003[k]
                    + pb_y[k] * mk_1290[k];

        t_1614[k] = f_9 * mi0_1004[k]
                    - f_10 * mi1_1004[k]
                    + pb_y[k] * mk_1291[k];
    }

#pragma omp simd aligned(t_1615, t_1616, t_1617, t_1618, pb_y, mi0_1005, mi0_1006, mi0_1007, \
                         mi1_1005, mi1_1006, mi1_1007, mk_1292, mk_1293, mk_1294, \
                         mk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1615[k] = f_7 * mi0_1005[k]
                    - f_8 * mi1_1005[k]
                    + pb_y[k] * mk_1292[k];

        t_1616[k] = f_5 * mi0_1006[k]
                    - f_6 * mi1_1006[k]
                    + pb_y[k] * mk_1293[k];

        t_1617[k] = f_3 * mi0_1007[k]
                    - f_4 * mi1_1007[k]
                    + pb_y[k] * mk_1294[k];

        t_1618[k] = pb_y[k] * mk_1295[k];
    }

#pragma omp simd aligned(t_1619, t_1620, t_1621, t_1622, pa_x, pb_y, pb_z, kl0_1619, kl1_1619, \
                         lk_1008, lk_1296, ll_1619, ll_1620, mk_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1619[k] = f_20 * kl0_1619[k]
                    - f_21 * kl1_1619[k]
                    + pa_x[k] * ll_1619[k];

        t_1620[k] = f_19 * lk_1296[k]
                    + pa_x[k] * ll_1620[k];

        t_1621[k] = f_19 * lk_1008[k]
                    + pb_y[k] * mk_1296[k];

        t_1622[k] = pb_z[k] * mk_1296[k];
    }

#pragma omp simd aligned(t_1623, t_1624, t_1625, t_1626, t_1627, pa_x, pb_z, lk_1299, lk_1301, \
                         lk_1302, ll_1623, ll_1625, ll_1626, mk_1297, \
                         mk_1299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1623[k] = f_18 * lk_1299[k]
                    + pa_x[k] * ll_1623[k];

        t_1624[k] = pb_z[k] * mk_1297[k];

        t_1625[k] = f_18 * lk_1301[k]
                    + pa_x[k] * ll_1625[k];

        t_1626[k] = f_17 * lk_1302[k]
                    + pa_x[k] * ll_1626[k];

        t_1627[k] = pb_z[k] * mk_1299[k];
    }

#pragma omp simd aligned(t_1628, t_1629, t_1630, t_1631, pa_x, pb_y, pb_z, lk_1013, lk_1305, \
                         lk_1306, ll_1629, ll_1630, mk_1301, mk_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1628[k] = f_19 * lk_1013[k]
                    + pb_y[k] * mk_1301[k];

        t_1629[k] = f_17 * lk_1305[k]
                    + pa_x[k] * ll_1629[k];

        t_1630[k] = f_16 * lk_1306[k]
                    + pa_x[k] * ll_1630[k];

        t_1631[k] = pb_z[k] * mk_1302[k];
    }

#pragma omp simd aligned(t_1632, t_1633, t_1634, t_1635, pa_x, pb_y, lk_1017, lk_1308, \
                         lk_1310, lk_1311, ll_1632, ll_1634, ll_1635, \
                         mk_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1632[k] = f_16 * lk_1308[k]
                    + pa_x[k] * ll_1632[k];

        t_1633[k] = f_19 * lk_1017[k]
                    + pb_y[k] * mk_1305[k];

        t_1634[k] = f_16 * lk_1310[k]
                    + pa_x[k] * ll_1634[k];

        t_1635[k] = f_15 * lk_1311[k]
                    + pa_x[k] * ll_1635[k];
    }

#pragma omp simd aligned(t_1636, t_1637, t_1638, t_1639, pa_x, pb_y, pb_z, lk_1022, lk_1313, \
                         lk_1314, ll_1637, ll_1638, mk_1306, mk_1310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1636[k] = pb_z[k] * mk_1306[k];

        t_1637[k] = f_15 * lk_1313[k]
                    + pa_x[k] * ll_1637[k];

        t_1638[k] = f_15 * lk_1314[k]
                    + pa_x[k] * ll_1638[k];

        t_1639[k] = f_19 * lk_1022[k]
                    + pb_y[k] * mk_1310[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece13(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t lk,
                                             const size_t ll, const size_t mk,
                                             const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_22 = 3.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lk_1008 = buffer.data(lk + 1008);
    const auto *lk_1011 = buffer.data(lk + 1011);
    const auto *lk_1014 = buffer.data(lk + 1014);
    const auto *lk_1018 = buffer.data(lk + 1018);
    const auto *lk_1023 = buffer.data(lk + 1023);
    const auto *lk_1028 = buffer.data(lk + 1028);
    const auto *lk_1044 = buffer.data(lk + 1044);
    const auto *lk_1046 = buffer.data(lk + 1046);
    const auto *lk_1047 = buffer.data(lk + 1047);
    const auto *lk_1049 = buffer.data(lk + 1049);
    const auto *lk_1050 = buffer.data(lk + 1050);
    const auto *lk_1053 = buffer.data(lk + 1053);
    const auto *lk_1054 = buffer.data(lk + 1054);
    const auto *lk_1058 = buffer.data(lk + 1058);
    const auto *lk_1059 = buffer.data(lk + 1059);
    const auto *lk_1064 = buffer.data(lk + 1064);
    const auto *lk_1080 = buffer.data(lk + 1080);
    const auto *lk_1082 = buffer.data(lk + 1082);
    const auto *lk_1083 = buffer.data(lk + 1083);
    const auto *lk_1085 = buffer.data(lk + 1085);
    const auto *lk_1086 = buffer.data(lk + 1086);
    const auto *lk_1089 = buffer.data(lk + 1089);
    const auto *lk_1090 = buffer.data(lk + 1090);
    const auto *lk_1094 = buffer.data(lk + 1094);
    const auto *lk_1095 = buffer.data(lk + 1095);
    const auto *lk_1100 = buffer.data(lk + 1100);
    const auto *lk_1116 = buffer.data(lk + 1116);
    const auto *lk_1118 = buffer.data(lk + 1118);
    const auto *lk_1121 = buffer.data(lk + 1121);
    const auto *lk_1125 = buffer.data(lk + 1125);
    const auto *lk_1130 = buffer.data(lk + 1130);
    const auto *lk_1136 = buffer.data(lk + 1136);
    const auto *lk_1316 = buffer.data(lk + 1316);
    const auto *lk_1317 = buffer.data(lk + 1317);
    const auto *lk_1319 = buffer.data(lk + 1319);
    const auto *lk_1320 = buffer.data(lk + 1320);
    const auto *lk_1321 = buffer.data(lk + 1321);
    const auto *lk_1323 = buffer.data(lk + 1323);
    const auto *lk_1324 = buffer.data(lk + 1324);
    const auto *lk_1326 = buffer.data(lk + 1326);
    const auto *lk_1327 = buffer.data(lk + 1327);
    const auto *lk_1328 = buffer.data(lk + 1328);
    const auto *lk_1329 = buffer.data(lk + 1329);
    const auto *lk_1330 = buffer.data(lk + 1330);
    const auto *lk_1331 = buffer.data(lk + 1331);
    const auto *lk_1337 = buffer.data(lk + 1337);
    const auto *lk_1341 = buffer.data(lk + 1341);
    const auto *lk_1344 = buffer.data(lk + 1344);
    const auto *lk_1346 = buffer.data(lk + 1346);
    const auto *lk_1349 = buffer.data(lk + 1349);
    const auto *lk_1350 = buffer.data(lk + 1350);
    const auto *lk_1352 = buffer.data(lk + 1352);
    const auto *lk_1355 = buffer.data(lk + 1355);
    const auto *lk_1356 = buffer.data(lk + 1356);
    const auto *lk_1357 = buffer.data(lk + 1357);
    const auto *lk_1359 = buffer.data(lk + 1359);
    const auto *lk_1361 = buffer.data(lk + 1361);
    const auto *lk_1362 = buffer.data(lk + 1362);
    const auto *lk_1363 = buffer.data(lk + 1363);
    const auto *lk_1364 = buffer.data(lk + 1364);
    const auto *lk_1365 = buffer.data(lk + 1365);
    const auto *lk_1366 = buffer.data(lk + 1366);
    const auto *lk_1367 = buffer.data(lk + 1367);
    const auto *lk_1368 = buffer.data(lk + 1368);
    const auto *lk_1371 = buffer.data(lk + 1371);
    const auto *lk_1373 = buffer.data(lk + 1373);
    const auto *lk_1374 = buffer.data(lk + 1374);
    const auto *lk_1377 = buffer.data(lk + 1377);
    const auto *lk_1378 = buffer.data(lk + 1378);
    const auto *lk_1380 = buffer.data(lk + 1380);
    const auto *lk_1382 = buffer.data(lk + 1382);
    const auto *lk_1383 = buffer.data(lk + 1383);
    const auto *lk_1385 = buffer.data(lk + 1385);
    const auto *lk_1386 = buffer.data(lk + 1386);
    const auto *lk_1388 = buffer.data(lk + 1388);
    const auto *lk_1389 = buffer.data(lk + 1389);
    const auto *lk_1391 = buffer.data(lk + 1391);
    const auto *lk_1392 = buffer.data(lk + 1392);
    const auto *lk_1393 = buffer.data(lk + 1393);
    const auto *lk_1395 = buffer.data(lk + 1395);
    const auto *lk_1396 = buffer.data(lk + 1396);
    const auto *lk_1397 = buffer.data(lk + 1397);
    const auto *lk_1398 = buffer.data(lk + 1398);
    const auto *lk_1399 = buffer.data(lk + 1399);
    const auto *lk_1400 = buffer.data(lk + 1400);
    const auto *lk_1401 = buffer.data(lk + 1401);
    const auto *lk_1402 = buffer.data(lk + 1402);
    const auto *lk_1403 = buffer.data(lk + 1403);
    const auto *lk_1404 = buffer.data(lk + 1404);
    const auto *lk_1407 = buffer.data(lk + 1407);
    const auto *lk_1409 = buffer.data(lk + 1409);
    const auto *lk_1410 = buffer.data(lk + 1410);
    const auto *lk_1413 = buffer.data(lk + 1413);
    const auto *lk_1414 = buffer.data(lk + 1414);
    const auto *lk_1416 = buffer.data(lk + 1416);
    const auto *lk_1418 = buffer.data(lk + 1418);
    const auto *lk_1419 = buffer.data(lk + 1419);
    const auto *lk_1421 = buffer.data(lk + 1421);
    const auto *lk_1422 = buffer.data(lk + 1422);
    const auto *lk_1424 = buffer.data(lk + 1424);
    const auto *lk_1425 = buffer.data(lk + 1425);
    const auto *lk_1427 = buffer.data(lk + 1427);
    const auto *lk_1428 = buffer.data(lk + 1428);
    const auto *lk_1429 = buffer.data(lk + 1429);
    const auto *lk_1431 = buffer.data(lk + 1431);
    const auto *lk_1432 = buffer.data(lk + 1432);
    const auto *lk_1433 = buffer.data(lk + 1433);
    const auto *lk_1434 = buffer.data(lk + 1434);
    const auto *lk_1435 = buffer.data(lk + 1435);
    const auto *lk_1436 = buffer.data(lk + 1436);
    const auto *lk_1437 = buffer.data(lk + 1437);

    const auto *ll_1260 = buffer.data(ll + 1260);
    const auto *ll_1261 = buffer.data(ll + 1261);
    const auto *ll_1263 = buffer.data(ll + 1263);
    const auto *ll_1266 = buffer.data(ll + 1266);
    const auto *ll_1270 = buffer.data(ll + 1270);
    const auto *ll_1275 = buffer.data(ll + 1275);
    const auto *ll_1281 = buffer.data(ll + 1281);
    const auto *ll_1288 = buffer.data(ll + 1288);
    const auto *ll_1640 = buffer.data(ll + 1640);
    const auto *ll_1641 = buffer.data(ll + 1641);
    const auto *ll_1643 = buffer.data(ll + 1643);
    const auto *ll_1644 = buffer.data(ll + 1644);
    const auto *ll_1645 = buffer.data(ll + 1645);
    const auto *ll_1647 = buffer.data(ll + 1647);
    const auto *ll_1656 = buffer.data(ll + 1656);
    const auto *ll_1658 = buffer.data(ll + 1658);
    const auto *ll_1659 = buffer.data(ll + 1659);
    const auto *ll_1660 = buffer.data(ll + 1660);
    const auto *ll_1661 = buffer.data(ll + 1661);
    const auto *ll_1662 = buffer.data(ll + 1662);
    const auto *ll_1663 = buffer.data(ll + 1663);
    const auto *ll_1664 = buffer.data(ll + 1664);
    const auto *ll_1670 = buffer.data(ll + 1670);
    const auto *ll_1674 = buffer.data(ll + 1674);
    const auto *ll_1677 = buffer.data(ll + 1677);
    const auto *ll_1679 = buffer.data(ll + 1679);
    const auto *ll_1682 = buffer.data(ll + 1682);
    const auto *ll_1683 = buffer.data(ll + 1683);
    const auto *ll_1685 = buffer.data(ll + 1685);
    const auto *ll_1688 = buffer.data(ll + 1688);
    const auto *ll_1689 = buffer.data(ll + 1689);
    const auto *ll_1690 = buffer.data(ll + 1690);
    const auto *ll_1692 = buffer.data(ll + 1692);
    const auto *ll_1701 = buffer.data(ll + 1701);
    const auto *ll_1702 = buffer.data(ll + 1702);
    const auto *ll_1703 = buffer.data(ll + 1703);
    const auto *ll_1704 = buffer.data(ll + 1704);
    const auto *ll_1705 = buffer.data(ll + 1705);
    const auto *ll_1706 = buffer.data(ll + 1706);
    const auto *ll_1707 = buffer.data(ll + 1707);
    const auto *ll_1708 = buffer.data(ll + 1708);
    const auto *ll_1709 = buffer.data(ll + 1709);
    const auto *ll_1710 = buffer.data(ll + 1710);
    const auto *ll_1713 = buffer.data(ll + 1713);
    const auto *ll_1715 = buffer.data(ll + 1715);
    const auto *ll_1716 = buffer.data(ll + 1716);
    const auto *ll_1719 = buffer.data(ll + 1719);
    const auto *ll_1720 = buffer.data(ll + 1720);
    const auto *ll_1722 = buffer.data(ll + 1722);
    const auto *ll_1724 = buffer.data(ll + 1724);
    const auto *ll_1725 = buffer.data(ll + 1725);
    const auto *ll_1727 = buffer.data(ll + 1727);
    const auto *ll_1728 = buffer.data(ll + 1728);
    const auto *ll_1730 = buffer.data(ll + 1730);
    const auto *ll_1731 = buffer.data(ll + 1731);
    const auto *ll_1733 = buffer.data(ll + 1733);
    const auto *ll_1734 = buffer.data(ll + 1734);
    const auto *ll_1735 = buffer.data(ll + 1735);
    const auto *ll_1737 = buffer.data(ll + 1737);
    const auto *ll_1746 = buffer.data(ll + 1746);
    const auto *ll_1747 = buffer.data(ll + 1747);
    const auto *ll_1748 = buffer.data(ll + 1748);
    const auto *ll_1749 = buffer.data(ll + 1749);
    const auto *ll_1750 = buffer.data(ll + 1750);
    const auto *ll_1751 = buffer.data(ll + 1751);
    const auto *ll_1752 = buffer.data(ll + 1752);
    const auto *ll_1753 = buffer.data(ll + 1753);
    const auto *ll_1754 = buffer.data(ll + 1754);
    const auto *ll_1755 = buffer.data(ll + 1755);
    const auto *ll_1758 = buffer.data(ll + 1758);
    const auto *ll_1760 = buffer.data(ll + 1760);
    const auto *ll_1761 = buffer.data(ll + 1761);
    const auto *ll_1764 = buffer.data(ll + 1764);
    const auto *ll_1765 = buffer.data(ll + 1765);
    const auto *ll_1767 = buffer.data(ll + 1767);
    const auto *ll_1769 = buffer.data(ll + 1769);
    const auto *ll_1770 = buffer.data(ll + 1770);
    const auto *ll_1772 = buffer.data(ll + 1772);
    const auto *ll_1773 = buffer.data(ll + 1773);
    const auto *ll_1775 = buffer.data(ll + 1775);
    const auto *ll_1776 = buffer.data(ll + 1776);
    const auto *ll_1778 = buffer.data(ll + 1778);
    const auto *ll_1779 = buffer.data(ll + 1779);
    const auto *ll_1780 = buffer.data(ll + 1780);
    const auto *ll_1782 = buffer.data(ll + 1782);

    const auto *mk_1311 = buffer.data(mk + 1311);
    const auto *mk_1316 = buffer.data(mk + 1316);
    const auto *mk_1317 = buffer.data(mk + 1317);
    const auto *mk_1324 = buffer.data(mk + 1324);
    const auto *mk_1326 = buffer.data(mk + 1326);
    const auto *mk_1327 = buffer.data(mk + 1327);
    const auto *mk_1328 = buffer.data(mk + 1328);
    const auto *mk_1329 = buffer.data(mk + 1329);
    const auto *mk_1330 = buffer.data(mk + 1330);
    const auto *mk_1331 = buffer.data(mk + 1331);
    const auto *mk_1332 = buffer.data(mk + 1332);
    const auto *mk_1334 = buffer.data(mk + 1334);
    const auto *mk_1335 = buffer.data(mk + 1335);
    const auto *mk_1337 = buffer.data(mk + 1337);
    const auto *mk_1338 = buffer.data(mk + 1338);
    const auto *mk_1341 = buffer.data(mk + 1341);
    const auto *mk_1342 = buffer.data(mk + 1342);
    const auto *mk_1346 = buffer.data(mk + 1346);
    const auto *mk_1347 = buffer.data(mk + 1347);
    const auto *mk_1352 = buffer.data(mk + 1352);
    const auto *mk_1361 = buffer.data(mk + 1361);
    const auto *mk_1362 = buffer.data(mk + 1362);
    const auto *mk_1363 = buffer.data(mk + 1363);
    const auto *mk_1364 = buffer.data(mk + 1364);
    const auto *mk_1365 = buffer.data(mk + 1365);
    const auto *mk_1366 = buffer.data(mk + 1366);
    const auto *mk_1367 = buffer.data(mk + 1367);
    const auto *mk_1368 = buffer.data(mk + 1368);
    const auto *mk_1370 = buffer.data(mk + 1370);
    const auto *mk_1371 = buffer.data(mk + 1371);
    const auto *mk_1373 = buffer.data(mk + 1373);
    const auto *mk_1374 = buffer.data(mk + 1374);
    const auto *mk_1377 = buffer.data(mk + 1377);
    const auto *mk_1378 = buffer.data(mk + 1378);
    const auto *mk_1382 = buffer.data(mk + 1382);
    const auto *mk_1383 = buffer.data(mk + 1383);
    const auto *mk_1388 = buffer.data(mk + 1388);
    const auto *mk_1396 = buffer.data(mk + 1396);
    const auto *mk_1397 = buffer.data(mk + 1397);
    const auto *mk_1398 = buffer.data(mk + 1398);
    const auto *mk_1399 = buffer.data(mk + 1399);
    const auto *mk_1400 = buffer.data(mk + 1400);
    const auto *mk_1401 = buffer.data(mk + 1401);
    const auto *mk_1402 = buffer.data(mk + 1402);
    const auto *mk_1403 = buffer.data(mk + 1403);
    const auto *mk_1404 = buffer.data(mk + 1404);
    const auto *mk_1406 = buffer.data(mk + 1406);
    const auto *mk_1407 = buffer.data(mk + 1407);
    const auto *mk_1409 = buffer.data(mk + 1409);
    const auto *mk_1410 = buffer.data(mk + 1410);
    const auto *mk_1413 = buffer.data(mk + 1413);
    const auto *mk_1414 = buffer.data(mk + 1414);
    const auto *mk_1418 = buffer.data(mk + 1418);
    const auto *mk_1419 = buffer.data(mk + 1419);
    const auto *mk_1424 = buffer.data(mk + 1424);
    const auto *mk_1432 = buffer.data(mk + 1432);
    const auto *mk_1433 = buffer.data(mk + 1433);
    const auto *mk_1434 = buffer.data(mk + 1434);
    const auto *mk_1435 = buffer.data(mk + 1435);
    const auto *mk_1436 = buffer.data(mk + 1436);
    const auto *mk_1437 = buffer.data(mk + 1437);

#pragma omp simd aligned(t_1640, t_1641, t_1642, t_1643, t_1644, pa_x, pb_z, lk_1316, lk_1317, \
                         lk_1319, lk_1320, ll_1640, ll_1641, ll_1643, ll_1644, \
                         mk_1311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1640[k] = f_15 * lk_1316[k]
                    + pa_x[k] * ll_1640[k];

        t_1641[k] = f_14 * lk_1317[k]
                    + pa_x[k] * ll_1641[k];

        t_1642[k] = pb_z[k] * mk_1311[k];

        t_1643[k] = f_14 * lk_1319[k]
                    + pa_x[k] * ll_1643[k];

        t_1644[k] = f_14 * lk_1320[k]
                    + pa_x[k] * ll_1644[k];
    }

#pragma omp simd aligned(t_1645, t_1646, t_1647, t_1648, pa_x, pb_x, pb_y, lk_1028, lk_1321, \
                         lk_1323, lk_1324, ll_1645, ll_1647, mk_1316, \
                         mk_1324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1645[k] = f_14 * lk_1321[k]
                    + pa_x[k] * ll_1645[k];

        t_1646[k] = f_19 * lk_1028[k]
                    + pb_y[k] * mk_1316[k];

        t_1647[k] = f_14 * lk_1323[k]
                    + pa_x[k] * ll_1647[k];

        t_1648[k] = f_13 * lk_1324[k]
                    + pb_x[k] * mk_1324[k];
    }

#pragma omp simd aligned(t_1649, t_1650, t_1651, t_1652, t_1653, pb_x, pb_z, lk_1326, lk_1327, \
                         lk_1328, lk_1329, mk_1317, mk_1326, mk_1327, mk_1328, \
                         mk_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1649[k] = pb_z[k] * mk_1317[k];

        t_1650[k] = f_13 * lk_1326[k]
                    + pb_x[k] * mk_1326[k];

        t_1651[k] = f_13 * lk_1327[k]
                    + pb_x[k] * mk_1327[k];

        t_1652[k] = f_13 * lk_1328[k]
                    + pb_x[k] * mk_1328[k];

        t_1653[k] = f_13 * lk_1329[k]
                    + pb_x[k] * mk_1329[k];
    }

#pragma omp simd aligned(t_1654, t_1655, t_1656, t_1657, t_1658, pa_x, pb_x, pb_z, lk_1330, \
                         lk_1331, ll_1656, ll_1658, mk_1324, mk_1330, \
                         mk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1654[k] = f_13 * lk_1330[k]
                    + pb_x[k] * mk_1330[k];

        t_1655[k] = f_13 * lk_1331[k]
                    + pb_x[k] * mk_1331[k];

        t_1656[k] = pa_x[k] * ll_1656[k];

        t_1657[k] = pb_z[k] * mk_1324[k];

        t_1658[k] = pa_x[k] * ll_1658[k];
    }

#pragma omp simd aligned(t_1659, t_1660, t_1661, t_1662, t_1663, t_1664, t_1665, pa_x, pa_z, \
                         ll_1260, ll_1659, ll_1660, ll_1661, ll_1662, ll_1663, \
                         ll_1664 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1659[k] = pa_x[k] * ll_1659[k];

        t_1660[k] = pa_x[k] * ll_1660[k];

        t_1661[k] = pa_x[k] * ll_1661[k];

        t_1662[k] = pa_x[k] * ll_1662[k];

        t_1663[k] = pa_x[k] * ll_1663[k];

        t_1664[k] = pa_x[k] * ll_1664[k];

        t_1665[k] = pa_z[k] * ll_1260[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, t_1669, pa_z, pb_y, pb_z, lk_1008, lk_1046, \
                         ll_1261, ll_1263, mk_1332, mk_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = pa_z[k] * ll_1261[k];

        t_1667[k] = f_13 * lk_1008[k]
                    + pb_z[k] * mk_1332[k];

        t_1668[k] = pa_z[k] * ll_1263[k];

        t_1669[k] = f_22 * lk_1046[k]
                    + pb_y[k] * mk_1334[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, t_1673, pa_x, pa_z, pb_y, pb_z, lk_1011, \
                         lk_1049, lk_1337, ll_1266, ll_1670, mk_1335, \
                         mk_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_18 * lk_1337[k]
                    + pa_x[k] * ll_1670[k];

        t_1671[k] = pa_z[k] * ll_1266[k];

        t_1672[k] = f_13 * lk_1011[k]
                    + pb_z[k] * mk_1335[k];

        t_1673[k] = f_22 * lk_1049[k]
                    + pb_y[k] * mk_1337[k];
    }

#pragma omp simd aligned(t_1674, t_1675, t_1676, t_1677, pa_x, pa_z, pb_z, lk_1014, lk_1341, \
                         lk_1344, ll_1270, ll_1674, ll_1677, mk_1338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1674[k] = f_17 * lk_1341[k]
                    + pa_x[k] * ll_1674[k];

        t_1675[k] = pa_z[k] * ll_1270[k];

        t_1676[k] = f_13 * lk_1014[k]
                    + pb_z[k] * mk_1338[k];

        t_1677[k] = f_16 * lk_1344[k]
                    + pa_x[k] * ll_1677[k];
    }

#pragma omp simd aligned(t_1678, t_1679, t_1680, t_1681, pa_x, pa_z, pb_y, pb_z, lk_1018, \
                         lk_1053, lk_1346, ll_1275, ll_1679, mk_1341, \
                         mk_1342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1678[k] = f_22 * lk_1053[k]
                    + pb_y[k] * mk_1341[k];

        t_1679[k] = f_16 * lk_1346[k]
                    + pa_x[k] * ll_1679[k];

        t_1680[k] = pa_z[k] * ll_1275[k];

        t_1681[k] = f_13 * lk_1018[k]
                    + pb_z[k] * mk_1342[k];
    }

#pragma omp simd aligned(t_1682, t_1683, t_1684, t_1685, pa_x, pb_y, lk_1058, lk_1349, \
                         lk_1350, lk_1352, ll_1682, ll_1683, ll_1685, \
                         mk_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1682[k] = f_15 * lk_1349[k]
                    + pa_x[k] * ll_1682[k];

        t_1683[k] = f_15 * lk_1350[k]
                    + pa_x[k] * ll_1683[k];

        t_1684[k] = f_22 * lk_1058[k]
                    + pb_y[k] * mk_1346[k];

        t_1685[k] = f_15 * lk_1352[k]
                    + pa_x[k] * ll_1685[k];
    }

#pragma omp simd aligned(t_1686, t_1687, t_1688, t_1689, pa_x, pa_z, pb_z, lk_1023, lk_1355, \
                         lk_1356, ll_1281, ll_1688, ll_1689, mk_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1686[k] = pa_z[k] * ll_1281[k];

        t_1687[k] = f_13 * lk_1023[k]
                    + pb_z[k] * mk_1347[k];

        t_1688[k] = f_14 * lk_1355[k]
                    + pa_x[k] * ll_1688[k];

        t_1689[k] = f_14 * lk_1356[k]
                    + pa_x[k] * ll_1689[k];
    }

#pragma omp simd aligned(t_1690, t_1691, t_1692, t_1693, pa_x, pa_z, pb_y, lk_1064, lk_1357, \
                         lk_1359, ll_1288, ll_1690, ll_1692, mk_1352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1690[k] = f_14 * lk_1357[k]
                    + pa_x[k] * ll_1690[k];

        t_1691[k] = f_22 * lk_1064[k]
                    + pb_y[k] * mk_1352[k];

        t_1692[k] = f_14 * lk_1359[k]
                    + pa_x[k] * ll_1692[k];

        t_1693[k] = pa_z[k] * ll_1288[k];
    }

#pragma omp simd aligned(t_1694, t_1695, t_1696, t_1697, t_1698, pb_x, lk_1361, lk_1362, \
                         lk_1363, lk_1364, lk_1365, mk_1361, mk_1362, mk_1363, mk_1364, \
                         mk_1365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1694[k] = f_13 * lk_1361[k]
                    + pb_x[k] * mk_1361[k];

        t_1695[k] = f_13 * lk_1362[k]
                    + pb_x[k] * mk_1362[k];

        t_1696[k] = f_13 * lk_1363[k]
                    + pb_x[k] * mk_1363[k];

        t_1697[k] = f_13 * lk_1364[k]
                    + pb_x[k] * mk_1364[k];

        t_1698[k] = f_13 * lk_1365[k]
                    + pb_x[k] * mk_1365[k];
    }

#pragma omp simd aligned(t_1699, t_1700, t_1701, t_1702, t_1703, t_1704, pa_x, pb_x, lk_1366, \
                         lk_1367, ll_1701, ll_1702, ll_1703, ll_1704, mk_1366, \
                         mk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1699[k] = f_13 * lk_1366[k]
                    + pb_x[k] * mk_1366[k];

        t_1700[k] = f_13 * lk_1367[k]
                    + pb_x[k] * mk_1367[k];

        t_1701[k] = pa_x[k] * ll_1701[k];

        t_1702[k] = pa_x[k] * ll_1702[k];

        t_1703[k] = pa_x[k] * ll_1703[k];

        t_1704[k] = pa_x[k] * ll_1704[k];
    }

#pragma omp simd aligned(t_1705, t_1706, t_1707, t_1708, t_1709, t_1710, pa_x, lk_1368, \
                         ll_1705, ll_1706, ll_1707, ll_1708, ll_1709, \
                         ll_1710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1705[k] = pa_x[k] * ll_1705[k];

        t_1706[k] = pa_x[k] * ll_1706[k];

        t_1707[k] = pa_x[k] * ll_1707[k];

        t_1708[k] = pa_x[k] * ll_1708[k];

        t_1709[k] = pa_x[k] * ll_1709[k];

        t_1710[k] = f_19 * lk_1368[k]
                    + pa_x[k] * ll_1710[k];
    }

#pragma omp simd aligned(t_1711, t_1712, t_1713, t_1714, pa_x, pb_y, pb_z, lk_1044, lk_1080, \
                         lk_1082, lk_1371, ll_1713, mk_1368, mk_1370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1711[k] = f_18 * lk_1080[k]
                    + pb_y[k] * mk_1368[k];

        t_1712[k] = f_14 * lk_1044[k]
                    + pb_z[k] * mk_1368[k];

        t_1713[k] = f_18 * lk_1371[k]
                    + pa_x[k] * ll_1713[k];

        t_1714[k] = f_18 * lk_1082[k]
                    + pb_y[k] * mk_1370[k];
    }

#pragma omp simd aligned(t_1715, t_1716, t_1717, t_1718, pa_x, pb_y, pb_z, lk_1047, lk_1085, \
                         lk_1373, lk_1374, ll_1715, ll_1716, mk_1371, \
                         mk_1373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1715[k] = f_18 * lk_1373[k]
                    + pa_x[k] * ll_1715[k];

        t_1716[k] = f_17 * lk_1374[k]
                    + pa_x[k] * ll_1716[k];

        t_1717[k] = f_14 * lk_1047[k]
                    + pb_z[k] * mk_1371[k];

        t_1718[k] = f_18 * lk_1085[k]
                    + pb_y[k] * mk_1373[k];
    }

#pragma omp simd aligned(t_1719, t_1720, t_1721, t_1722, pa_x, pb_z, lk_1050, lk_1377, \
                         lk_1378, lk_1380, ll_1719, ll_1720, ll_1722, \
                         mk_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1719[k] = f_17 * lk_1377[k]
                    + pa_x[k] * ll_1719[k];

        t_1720[k] = f_16 * lk_1378[k]
                    + pa_x[k] * ll_1720[k];

        t_1721[k] = f_14 * lk_1050[k]
                    + pb_z[k] * mk_1374[k];

        t_1722[k] = f_16 * lk_1380[k]
                    + pa_x[k] * ll_1722[k];
    }

#pragma omp simd aligned(t_1723, t_1724, t_1725, t_1726, pa_x, pb_y, pb_z, lk_1054, lk_1089, \
                         lk_1382, lk_1383, ll_1724, ll_1725, mk_1377, \
                         mk_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1723[k] = f_18 * lk_1089[k]
                    + pb_y[k] * mk_1377[k];

        t_1724[k] = f_16 * lk_1382[k]
                    + pa_x[k] * ll_1724[k];

        t_1725[k] = f_15 * lk_1383[k]
                    + pa_x[k] * ll_1725[k];

        t_1726[k] = f_14 * lk_1054[k]
                    + pb_z[k] * mk_1378[k];
    }

#pragma omp simd aligned(t_1727, t_1728, t_1729, t_1730, pa_x, pb_y, lk_1094, lk_1385, \
                         lk_1386, lk_1388, ll_1727, ll_1728, ll_1730, \
                         mk_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1727[k] = f_15 * lk_1385[k]
                    + pa_x[k] * ll_1727[k];

        t_1728[k] = f_15 * lk_1386[k]
                    + pa_x[k] * ll_1728[k];

        t_1729[k] = f_18 * lk_1094[k]
                    + pb_y[k] * mk_1382[k];

        t_1730[k] = f_15 * lk_1388[k]
                    + pa_x[k] * ll_1730[k];
    }

#pragma omp simd aligned(t_1731, t_1732, t_1733, t_1734, pa_x, pb_z, lk_1059, lk_1389, \
                         lk_1391, lk_1392, ll_1731, ll_1733, ll_1734, \
                         mk_1383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1731[k] = f_14 * lk_1389[k]
                    + pa_x[k] * ll_1731[k];

        t_1732[k] = f_14 * lk_1059[k]
                    + pb_z[k] * mk_1383[k];

        t_1733[k] = f_14 * lk_1391[k]
                    + pa_x[k] * ll_1733[k];

        t_1734[k] = f_14 * lk_1392[k]
                    + pa_x[k] * ll_1734[k];
    }

#pragma omp simd aligned(t_1735, t_1736, t_1737, t_1738, pa_x, pb_x, pb_y, lk_1100, lk_1393, \
                         lk_1395, lk_1396, ll_1735, ll_1737, mk_1388, \
                         mk_1396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1735[k] = f_14 * lk_1393[k]
                    + pa_x[k] * ll_1735[k];

        t_1736[k] = f_18 * lk_1100[k]
                    + pb_y[k] * mk_1388[k];

        t_1737[k] = f_14 * lk_1395[k]
                    + pa_x[k] * ll_1737[k];

        t_1738[k] = f_13 * lk_1396[k]
                    + pb_x[k] * mk_1396[k];
    }

#pragma omp simd aligned(t_1739, t_1740, t_1741, t_1742, t_1743, pb_x, lk_1397, lk_1398, \
                         lk_1399, lk_1400, lk_1401, mk_1397, mk_1398, mk_1399, mk_1400, \
                         mk_1401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1739[k] = f_13 * lk_1397[k]
                    + pb_x[k] * mk_1397[k];

        t_1740[k] = f_13 * lk_1398[k]
                    + pb_x[k] * mk_1398[k];

        t_1741[k] = f_13 * lk_1399[k]
                    + pb_x[k] * mk_1399[k];

        t_1742[k] = f_13 * lk_1400[k]
                    + pb_x[k] * mk_1400[k];

        t_1743[k] = f_13 * lk_1401[k]
                    + pb_x[k] * mk_1401[k];
    }

#pragma omp simd aligned(t_1744, t_1745, t_1746, t_1747, t_1748, t_1749, pa_x, pb_x, lk_1402, \
                         lk_1403, ll_1746, ll_1747, ll_1748, ll_1749, mk_1402, \
                         mk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1744[k] = f_13 * lk_1402[k]
                    + pb_x[k] * mk_1402[k];

        t_1745[k] = f_13 * lk_1403[k]
                    + pb_x[k] * mk_1403[k];

        t_1746[k] = pa_x[k] * ll_1746[k];

        t_1747[k] = pa_x[k] * ll_1747[k];

        t_1748[k] = pa_x[k] * ll_1748[k];

        t_1749[k] = pa_x[k] * ll_1749[k];
    }

#pragma omp simd aligned(t_1750, t_1751, t_1752, t_1753, t_1754, t_1755, pa_x, lk_1404, \
                         ll_1750, ll_1751, ll_1752, ll_1753, ll_1754, \
                         ll_1755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1750[k] = pa_x[k] * ll_1750[k];

        t_1751[k] = pa_x[k] * ll_1751[k];

        t_1752[k] = pa_x[k] * ll_1752[k];

        t_1753[k] = pa_x[k] * ll_1753[k];

        t_1754[k] = pa_x[k] * ll_1754[k];

        t_1755[k] = f_19 * lk_1404[k]
                    + pa_x[k] * ll_1755[k];
    }

#pragma omp simd aligned(t_1756, t_1757, t_1758, t_1759, pa_x, pb_y, pb_z, lk_1080, lk_1116, \
                         lk_1118, lk_1407, ll_1758, mk_1404, mk_1406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1756[k] = f_17 * lk_1116[k]
                    + pb_y[k] * mk_1404[k];

        t_1757[k] = f_15 * lk_1080[k]
                    + pb_z[k] * mk_1404[k];

        t_1758[k] = f_18 * lk_1407[k]
                    + pa_x[k] * ll_1758[k];

        t_1759[k] = f_17 * lk_1118[k]
                    + pb_y[k] * mk_1406[k];
    }

#pragma omp simd aligned(t_1760, t_1761, t_1762, t_1763, pa_x, pb_y, pb_z, lk_1083, lk_1121, \
                         lk_1409, lk_1410, ll_1760, ll_1761, mk_1407, \
                         mk_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1760[k] = f_18 * lk_1409[k]
                    + pa_x[k] * ll_1760[k];

        t_1761[k] = f_17 * lk_1410[k]
                    + pa_x[k] * ll_1761[k];

        t_1762[k] = f_15 * lk_1083[k]
                    + pb_z[k] * mk_1407[k];

        t_1763[k] = f_17 * lk_1121[k]
                    + pb_y[k] * mk_1409[k];
    }

#pragma omp simd aligned(t_1764, t_1765, t_1766, t_1767, pa_x, pb_z, lk_1086, lk_1413, \
                         lk_1414, lk_1416, ll_1764, ll_1765, ll_1767, \
                         mk_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1764[k] = f_17 * lk_1413[k]
                    + pa_x[k] * ll_1764[k];

        t_1765[k] = f_16 * lk_1414[k]
                    + pa_x[k] * ll_1765[k];

        t_1766[k] = f_15 * lk_1086[k]
                    + pb_z[k] * mk_1410[k];

        t_1767[k] = f_16 * lk_1416[k]
                    + pa_x[k] * ll_1767[k];
    }

#pragma omp simd aligned(t_1768, t_1769, t_1770, t_1771, pa_x, pb_y, pb_z, lk_1090, lk_1125, \
                         lk_1418, lk_1419, ll_1769, ll_1770, mk_1413, \
                         mk_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1768[k] = f_17 * lk_1125[k]
                    + pb_y[k] * mk_1413[k];

        t_1769[k] = f_16 * lk_1418[k]
                    + pa_x[k] * ll_1769[k];

        t_1770[k] = f_15 * lk_1419[k]
                    + pa_x[k] * ll_1770[k];

        t_1771[k] = f_15 * lk_1090[k]
                    + pb_z[k] * mk_1414[k];
    }

#pragma omp simd aligned(t_1772, t_1773, t_1774, t_1775, pa_x, pb_y, lk_1130, lk_1421, \
                         lk_1422, lk_1424, ll_1772, ll_1773, ll_1775, \
                         mk_1418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1772[k] = f_15 * lk_1421[k]
                    + pa_x[k] * ll_1772[k];

        t_1773[k] = f_15 * lk_1422[k]
                    + pa_x[k] * ll_1773[k];

        t_1774[k] = f_17 * lk_1130[k]
                    + pb_y[k] * mk_1418[k];

        t_1775[k] = f_15 * lk_1424[k]
                    + pa_x[k] * ll_1775[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, t_1779, pa_x, pb_z, lk_1095, lk_1425, \
                         lk_1427, lk_1428, ll_1776, ll_1778, ll_1779, \
                         mk_1419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_14 * lk_1425[k]
                    + pa_x[k] * ll_1776[k];

        t_1777[k] = f_15 * lk_1095[k]
                    + pb_z[k] * mk_1419[k];

        t_1778[k] = f_14 * lk_1427[k]
                    + pa_x[k] * ll_1778[k];

        t_1779[k] = f_14 * lk_1428[k]
                    + pa_x[k] * ll_1779[k];
    }

#pragma omp simd aligned(t_1780, t_1781, t_1782, t_1783, pa_x, pb_x, pb_y, lk_1136, lk_1429, \
                         lk_1431, lk_1432, ll_1780, ll_1782, mk_1424, \
                         mk_1432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1780[k] = f_14 * lk_1429[k]
                    + pa_x[k] * ll_1780[k];

        t_1781[k] = f_17 * lk_1136[k]
                    + pb_y[k] * mk_1424[k];

        t_1782[k] = f_14 * lk_1431[k]
                    + pa_x[k] * ll_1782[k];

        t_1783[k] = f_13 * lk_1432[k]
                    + pb_x[k] * mk_1432[k];
    }

#pragma omp simd aligned(t_1784, t_1785, t_1786, t_1787, t_1788, pb_x, lk_1433, lk_1434, \
                         lk_1435, lk_1436, lk_1437, mk_1433, mk_1434, mk_1435, mk_1436, \
                         mk_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1784[k] = f_13 * lk_1433[k]
                    + pb_x[k] * mk_1433[k];

        t_1785[k] = f_13 * lk_1434[k]
                    + pb_x[k] * mk_1434[k];

        t_1786[k] = f_13 * lk_1435[k]
                    + pb_x[k] * mk_1435[k];

        t_1787[k] = f_13 * lk_1436[k]
                    + pb_x[k] * mk_1436[k];

        t_1788[k] = f_13 * lk_1437[k]
                    + pb_x[k] * mk_1437[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece14(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t lk,
                                             const size_t ll, const size_t mk,
                                             const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_13 = 0.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lk_1116 = buffer.data(lk + 1116);
    const auto *lk_1119 = buffer.data(lk + 1119);
    const auto *lk_1122 = buffer.data(lk + 1122);
    const auto *lk_1126 = buffer.data(lk + 1126);
    const auto *lk_1131 = buffer.data(lk + 1131);
    const auto *lk_1152 = buffer.data(lk + 1152);
    const auto *lk_1154 = buffer.data(lk + 1154);
    const auto *lk_1155 = buffer.data(lk + 1155);
    const auto *lk_1157 = buffer.data(lk + 1157);
    const auto *lk_1158 = buffer.data(lk + 1158);
    const auto *lk_1161 = buffer.data(lk + 1161);
    const auto *lk_1162 = buffer.data(lk + 1162);
    const auto *lk_1166 = buffer.data(lk + 1166);
    const auto *lk_1167 = buffer.data(lk + 1167);
    const auto *lk_1172 = buffer.data(lk + 1172);
    const auto *lk_1188 = buffer.data(lk + 1188);
    const auto *lk_1190 = buffer.data(lk + 1190);
    const auto *lk_1191 = buffer.data(lk + 1191);
    const auto *lk_1193 = buffer.data(lk + 1193);
    const auto *lk_1194 = buffer.data(lk + 1194);
    const auto *lk_1197 = buffer.data(lk + 1197);
    const auto *lk_1198 = buffer.data(lk + 1198);
    const auto *lk_1202 = buffer.data(lk + 1202);
    const auto *lk_1203 = buffer.data(lk + 1203);
    const auto *lk_1208 = buffer.data(lk + 1208);
    const auto *lk_1224 = buffer.data(lk + 1224);
    const auto *lk_1226 = buffer.data(lk + 1226);
    const auto *lk_1229 = buffer.data(lk + 1229);
    const auto *lk_1233 = buffer.data(lk + 1233);
    const auto *lk_1238 = buffer.data(lk + 1238);
    const auto *lk_1244 = buffer.data(lk + 1244);
    const auto *lk_1438 = buffer.data(lk + 1438);
    const auto *lk_1439 = buffer.data(lk + 1439);
    const auto *lk_1440 = buffer.data(lk + 1440);
    const auto *lk_1443 = buffer.data(lk + 1443);
    const auto *lk_1445 = buffer.data(lk + 1445);
    const auto *lk_1446 = buffer.data(lk + 1446);
    const auto *lk_1449 = buffer.data(lk + 1449);
    const auto *lk_1450 = buffer.data(lk + 1450);
    const auto *lk_1452 = buffer.data(lk + 1452);
    const auto *lk_1454 = buffer.data(lk + 1454);
    const auto *lk_1455 = buffer.data(lk + 1455);
    const auto *lk_1457 = buffer.data(lk + 1457);
    const auto *lk_1458 = buffer.data(lk + 1458);
    const auto *lk_1460 = buffer.data(lk + 1460);
    const auto *lk_1461 = buffer.data(lk + 1461);
    const auto *lk_1463 = buffer.data(lk + 1463);
    const auto *lk_1464 = buffer.data(lk + 1464);
    const auto *lk_1465 = buffer.data(lk + 1465);
    const auto *lk_1467 = buffer.data(lk + 1467);
    const auto *lk_1468 = buffer.data(lk + 1468);
    const auto *lk_1469 = buffer.data(lk + 1469);
    const auto *lk_1470 = buffer.data(lk + 1470);
    const auto *lk_1471 = buffer.data(lk + 1471);
    const auto *lk_1472 = buffer.data(lk + 1472);
    const auto *lk_1473 = buffer.data(lk + 1473);
    const auto *lk_1474 = buffer.data(lk + 1474);
    const auto *lk_1475 = buffer.data(lk + 1475);
    const auto *lk_1476 = buffer.data(lk + 1476);
    const auto *lk_1479 = buffer.data(lk + 1479);
    const auto *lk_1481 = buffer.data(lk + 1481);
    const auto *lk_1482 = buffer.data(lk + 1482);
    const auto *lk_1485 = buffer.data(lk + 1485);
    const auto *lk_1486 = buffer.data(lk + 1486);
    const auto *lk_1488 = buffer.data(lk + 1488);
    const auto *lk_1490 = buffer.data(lk + 1490);
    const auto *lk_1491 = buffer.data(lk + 1491);
    const auto *lk_1493 = buffer.data(lk + 1493);
    const auto *lk_1494 = buffer.data(lk + 1494);
    const auto *lk_1496 = buffer.data(lk + 1496);
    const auto *lk_1497 = buffer.data(lk + 1497);
    const auto *lk_1499 = buffer.data(lk + 1499);
    const auto *lk_1500 = buffer.data(lk + 1500);
    const auto *lk_1501 = buffer.data(lk + 1501);
    const auto *lk_1503 = buffer.data(lk + 1503);
    const auto *lk_1504 = buffer.data(lk + 1504);
    const auto *lk_1505 = buffer.data(lk + 1505);
    const auto *lk_1506 = buffer.data(lk + 1506);
    const auto *lk_1507 = buffer.data(lk + 1507);
    const auto *lk_1508 = buffer.data(lk + 1508);
    const auto *lk_1509 = buffer.data(lk + 1509);
    const auto *lk_1510 = buffer.data(lk + 1510);
    const auto *lk_1511 = buffer.data(lk + 1511);
    const auto *lk_1512 = buffer.data(lk + 1512);
    const auto *lk_1515 = buffer.data(lk + 1515);
    const auto *lk_1517 = buffer.data(lk + 1517);
    const auto *lk_1518 = buffer.data(lk + 1518);
    const auto *lk_1521 = buffer.data(lk + 1521);
    const auto *lk_1522 = buffer.data(lk + 1522);
    const auto *lk_1524 = buffer.data(lk + 1524);
    const auto *lk_1526 = buffer.data(lk + 1526);
    const auto *lk_1527 = buffer.data(lk + 1527);
    const auto *lk_1529 = buffer.data(lk + 1529);
    const auto *lk_1530 = buffer.data(lk + 1530);
    const auto *lk_1532 = buffer.data(lk + 1532);
    const auto *lk_1533 = buffer.data(lk + 1533);
    const auto *lk_1535 = buffer.data(lk + 1535);
    const auto *lk_1536 = buffer.data(lk + 1536);
    const auto *lk_1537 = buffer.data(lk + 1537);
    const auto *lk_1539 = buffer.data(lk + 1539);
    const auto *lk_1540 = buffer.data(lk + 1540);
    const auto *lk_1541 = buffer.data(lk + 1541);
    const auto *lk_1542 = buffer.data(lk + 1542);
    const auto *lk_1543 = buffer.data(lk + 1543);
    const auto *lk_1544 = buffer.data(lk + 1544);
    const auto *lk_1545 = buffer.data(lk + 1545);
    const auto *lk_1546 = buffer.data(lk + 1546);
    const auto *lk_1547 = buffer.data(lk + 1547);

    const auto *ll_1575 = buffer.data(ll + 1575);
    const auto *ll_1791 = buffer.data(ll + 1791);
    const auto *ll_1792 = buffer.data(ll + 1792);
    const auto *ll_1793 = buffer.data(ll + 1793);
    const auto *ll_1794 = buffer.data(ll + 1794);
    const auto *ll_1795 = buffer.data(ll + 1795);
    const auto *ll_1796 = buffer.data(ll + 1796);
    const auto *ll_1797 = buffer.data(ll + 1797);
    const auto *ll_1798 = buffer.data(ll + 1798);
    const auto *ll_1799 = buffer.data(ll + 1799);
    const auto *ll_1800 = buffer.data(ll + 1800);
    const auto *ll_1803 = buffer.data(ll + 1803);
    const auto *ll_1805 = buffer.data(ll + 1805);
    const auto *ll_1806 = buffer.data(ll + 1806);
    const auto *ll_1809 = buffer.data(ll + 1809);
    const auto *ll_1810 = buffer.data(ll + 1810);
    const auto *ll_1812 = buffer.data(ll + 1812);
    const auto *ll_1814 = buffer.data(ll + 1814);
    const auto *ll_1815 = buffer.data(ll + 1815);
    const auto *ll_1817 = buffer.data(ll + 1817);
    const auto *ll_1818 = buffer.data(ll + 1818);
    const auto *ll_1820 = buffer.data(ll + 1820);
    const auto *ll_1821 = buffer.data(ll + 1821);
    const auto *ll_1823 = buffer.data(ll + 1823);
    const auto *ll_1824 = buffer.data(ll + 1824);
    const auto *ll_1825 = buffer.data(ll + 1825);
    const auto *ll_1827 = buffer.data(ll + 1827);
    const auto *ll_1836 = buffer.data(ll + 1836);
    const auto *ll_1837 = buffer.data(ll + 1837);
    const auto *ll_1838 = buffer.data(ll + 1838);
    const auto *ll_1839 = buffer.data(ll + 1839);
    const auto *ll_1840 = buffer.data(ll + 1840);
    const auto *ll_1841 = buffer.data(ll + 1841);
    const auto *ll_1842 = buffer.data(ll + 1842);
    const auto *ll_1843 = buffer.data(ll + 1843);
    const auto *ll_1844 = buffer.data(ll + 1844);
    const auto *ll_1845 = buffer.data(ll + 1845);
    const auto *ll_1848 = buffer.data(ll + 1848);
    const auto *ll_1850 = buffer.data(ll + 1850);
    const auto *ll_1851 = buffer.data(ll + 1851);
    const auto *ll_1854 = buffer.data(ll + 1854);
    const auto *ll_1855 = buffer.data(ll + 1855);
    const auto *ll_1857 = buffer.data(ll + 1857);
    const auto *ll_1859 = buffer.data(ll + 1859);
    const auto *ll_1860 = buffer.data(ll + 1860);
    const auto *ll_1862 = buffer.data(ll + 1862);
    const auto *ll_1863 = buffer.data(ll + 1863);
    const auto *ll_1865 = buffer.data(ll + 1865);
    const auto *ll_1866 = buffer.data(ll + 1866);
    const auto *ll_1868 = buffer.data(ll + 1868);
    const auto *ll_1869 = buffer.data(ll + 1869);
    const auto *ll_1870 = buffer.data(ll + 1870);
    const auto *ll_1872 = buffer.data(ll + 1872);
    const auto *ll_1881 = buffer.data(ll + 1881);
    const auto *ll_1882 = buffer.data(ll + 1882);
    const auto *ll_1883 = buffer.data(ll + 1883);
    const auto *ll_1884 = buffer.data(ll + 1884);
    const auto *ll_1885 = buffer.data(ll + 1885);
    const auto *ll_1886 = buffer.data(ll + 1886);
    const auto *ll_1887 = buffer.data(ll + 1887);
    const auto *ll_1888 = buffer.data(ll + 1888);
    const auto *ll_1889 = buffer.data(ll + 1889);
    const auto *ll_1890 = buffer.data(ll + 1890);
    const auto *ll_1893 = buffer.data(ll + 1893);
    const auto *ll_1895 = buffer.data(ll + 1895);
    const auto *ll_1896 = buffer.data(ll + 1896);
    const auto *ll_1899 = buffer.data(ll + 1899);
    const auto *ll_1900 = buffer.data(ll + 1900);
    const auto *ll_1902 = buffer.data(ll + 1902);
    const auto *ll_1904 = buffer.data(ll + 1904);
    const auto *ll_1905 = buffer.data(ll + 1905);
    const auto *ll_1907 = buffer.data(ll + 1907);
    const auto *ll_1908 = buffer.data(ll + 1908);
    const auto *ll_1910 = buffer.data(ll + 1910);
    const auto *ll_1911 = buffer.data(ll + 1911);
    const auto *ll_1913 = buffer.data(ll + 1913);
    const auto *ll_1914 = buffer.data(ll + 1914);
    const auto *ll_1915 = buffer.data(ll + 1915);
    const auto *ll_1917 = buffer.data(ll + 1917);
    const auto *ll_1926 = buffer.data(ll + 1926);
    const auto *ll_1927 = buffer.data(ll + 1927);
    const auto *ll_1928 = buffer.data(ll + 1928);
    const auto *ll_1929 = buffer.data(ll + 1929);
    const auto *ll_1930 = buffer.data(ll + 1930);
    const auto *ll_1931 = buffer.data(ll + 1931);
    const auto *ll_1932 = buffer.data(ll + 1932);
    const auto *ll_1933 = buffer.data(ll + 1933);
    const auto *ll_1934 = buffer.data(ll + 1934);

    const auto *mk_1438 = buffer.data(mk + 1438);
    const auto *mk_1439 = buffer.data(mk + 1439);
    const auto *mk_1440 = buffer.data(mk + 1440);
    const auto *mk_1442 = buffer.data(mk + 1442);
    const auto *mk_1443 = buffer.data(mk + 1443);
    const auto *mk_1445 = buffer.data(mk + 1445);
    const auto *mk_1446 = buffer.data(mk + 1446);
    const auto *mk_1449 = buffer.data(mk + 1449);
    const auto *mk_1450 = buffer.data(mk + 1450);
    const auto *mk_1454 = buffer.data(mk + 1454);
    const auto *mk_1455 = buffer.data(mk + 1455);
    const auto *mk_1460 = buffer.data(mk + 1460);
    const auto *mk_1468 = buffer.data(mk + 1468);
    const auto *mk_1469 = buffer.data(mk + 1469);
    const auto *mk_1470 = buffer.data(mk + 1470);
    const auto *mk_1471 = buffer.data(mk + 1471);
    const auto *mk_1472 = buffer.data(mk + 1472);
    const auto *mk_1473 = buffer.data(mk + 1473);
    const auto *mk_1474 = buffer.data(mk + 1474);
    const auto *mk_1475 = buffer.data(mk + 1475);
    const auto *mk_1476 = buffer.data(mk + 1476);
    const auto *mk_1478 = buffer.data(mk + 1478);
    const auto *mk_1479 = buffer.data(mk + 1479);
    const auto *mk_1481 = buffer.data(mk + 1481);
    const auto *mk_1482 = buffer.data(mk + 1482);
    const auto *mk_1485 = buffer.data(mk + 1485);
    const auto *mk_1486 = buffer.data(mk + 1486);
    const auto *mk_1490 = buffer.data(mk + 1490);
    const auto *mk_1491 = buffer.data(mk + 1491);
    const auto *mk_1496 = buffer.data(mk + 1496);
    const auto *mk_1504 = buffer.data(mk + 1504);
    const auto *mk_1505 = buffer.data(mk + 1505);
    const auto *mk_1506 = buffer.data(mk + 1506);
    const auto *mk_1507 = buffer.data(mk + 1507);
    const auto *mk_1508 = buffer.data(mk + 1508);
    const auto *mk_1509 = buffer.data(mk + 1509);
    const auto *mk_1510 = buffer.data(mk + 1510);
    const auto *mk_1511 = buffer.data(mk + 1511);
    const auto *mk_1512 = buffer.data(mk + 1512);
    const auto *mk_1514 = buffer.data(mk + 1514);
    const auto *mk_1515 = buffer.data(mk + 1515);
    const auto *mk_1517 = buffer.data(mk + 1517);
    const auto *mk_1518 = buffer.data(mk + 1518);
    const auto *mk_1521 = buffer.data(mk + 1521);
    const auto *mk_1522 = buffer.data(mk + 1522);
    const auto *mk_1526 = buffer.data(mk + 1526);
    const auto *mk_1527 = buffer.data(mk + 1527);
    const auto *mk_1532 = buffer.data(mk + 1532);
    const auto *mk_1540 = buffer.data(mk + 1540);
    const auto *mk_1541 = buffer.data(mk + 1541);
    const auto *mk_1542 = buffer.data(mk + 1542);
    const auto *mk_1543 = buffer.data(mk + 1543);
    const auto *mk_1544 = buffer.data(mk + 1544);
    const auto *mk_1545 = buffer.data(mk + 1545);
    const auto *mk_1546 = buffer.data(mk + 1546);
    const auto *mk_1547 = buffer.data(mk + 1547);

#pragma omp simd aligned(t_1789, t_1790, t_1791, t_1792, t_1793, t_1794, pa_x, pb_x, lk_1438, \
                         lk_1439, ll_1791, ll_1792, ll_1793, ll_1794, mk_1438, \
                         mk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1789[k] = f_13 * lk_1438[k]
                    + pb_x[k] * mk_1438[k];

        t_1790[k] = f_13 * lk_1439[k]
                    + pb_x[k] * mk_1439[k];

        t_1791[k] = pa_x[k] * ll_1791[k];

        t_1792[k] = pa_x[k] * ll_1792[k];

        t_1793[k] = pa_x[k] * ll_1793[k];

        t_1794[k] = pa_x[k] * ll_1794[k];
    }

#pragma omp simd aligned(t_1795, t_1796, t_1797, t_1798, t_1799, t_1800, pa_x, lk_1440, \
                         ll_1795, ll_1796, ll_1797, ll_1798, ll_1799, \
                         ll_1800 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1795[k] = pa_x[k] * ll_1795[k];

        t_1796[k] = pa_x[k] * ll_1796[k];

        t_1797[k] = pa_x[k] * ll_1797[k];

        t_1798[k] = pa_x[k] * ll_1798[k];

        t_1799[k] = pa_x[k] * ll_1799[k];

        t_1800[k] = f_19 * lk_1440[k]
                    + pa_x[k] * ll_1800[k];
    }

#pragma omp simd aligned(t_1801, t_1802, t_1803, t_1804, pa_x, pb_y, pb_z, lk_1116, lk_1152, \
                         lk_1154, lk_1443, ll_1803, mk_1440, mk_1442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1801[k] = f_16 * lk_1152[k]
                    + pb_y[k] * mk_1440[k];

        t_1802[k] = f_16 * lk_1116[k]
                    + pb_z[k] * mk_1440[k];

        t_1803[k] = f_18 * lk_1443[k]
                    + pa_x[k] * ll_1803[k];

        t_1804[k] = f_16 * lk_1154[k]
                    + pb_y[k] * mk_1442[k];
    }

#pragma omp simd aligned(t_1805, t_1806, t_1807, t_1808, pa_x, pb_y, pb_z, lk_1119, lk_1157, \
                         lk_1445, lk_1446, ll_1805, ll_1806, mk_1443, \
                         mk_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1805[k] = f_18 * lk_1445[k]
                    + pa_x[k] * ll_1805[k];

        t_1806[k] = f_17 * lk_1446[k]
                    + pa_x[k] * ll_1806[k];

        t_1807[k] = f_16 * lk_1119[k]
                    + pb_z[k] * mk_1443[k];

        t_1808[k] = f_16 * lk_1157[k]
                    + pb_y[k] * mk_1445[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, t_1812, pa_x, pb_z, lk_1122, lk_1449, \
                         lk_1450, lk_1452, ll_1809, ll_1810, ll_1812, \
                         mk_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = f_17 * lk_1449[k]
                    + pa_x[k] * ll_1809[k];

        t_1810[k] = f_16 * lk_1450[k]
                    + pa_x[k] * ll_1810[k];

        t_1811[k] = f_16 * lk_1122[k]
                    + pb_z[k] * mk_1446[k];

        t_1812[k] = f_16 * lk_1452[k]
                    + pa_x[k] * ll_1812[k];
    }

#pragma omp simd aligned(t_1813, t_1814, t_1815, t_1816, pa_x, pb_y, pb_z, lk_1126, lk_1161, \
                         lk_1454, lk_1455, ll_1814, ll_1815, mk_1449, \
                         mk_1450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1813[k] = f_16 * lk_1161[k]
                    + pb_y[k] * mk_1449[k];

        t_1814[k] = f_16 * lk_1454[k]
                    + pa_x[k] * ll_1814[k];

        t_1815[k] = f_15 * lk_1455[k]
                    + pa_x[k] * ll_1815[k];

        t_1816[k] = f_16 * lk_1126[k]
                    + pb_z[k] * mk_1450[k];
    }

#pragma omp simd aligned(t_1817, t_1818, t_1819, t_1820, pa_x, pb_y, lk_1166, lk_1457, \
                         lk_1458, lk_1460, ll_1817, ll_1818, ll_1820, \
                         mk_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1817[k] = f_15 * lk_1457[k]
                    + pa_x[k] * ll_1817[k];

        t_1818[k] = f_15 * lk_1458[k]
                    + pa_x[k] * ll_1818[k];

        t_1819[k] = f_16 * lk_1166[k]
                    + pb_y[k] * mk_1454[k];

        t_1820[k] = f_15 * lk_1460[k]
                    + pa_x[k] * ll_1820[k];
    }

#pragma omp simd aligned(t_1821, t_1822, t_1823, t_1824, pa_x, pb_z, lk_1131, lk_1461, \
                         lk_1463, lk_1464, ll_1821, ll_1823, ll_1824, \
                         mk_1455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1821[k] = f_14 * lk_1461[k]
                    + pa_x[k] * ll_1821[k];

        t_1822[k] = f_16 * lk_1131[k]
                    + pb_z[k] * mk_1455[k];

        t_1823[k] = f_14 * lk_1463[k]
                    + pa_x[k] * ll_1823[k];

        t_1824[k] = f_14 * lk_1464[k]
                    + pa_x[k] * ll_1824[k];
    }

#pragma omp simd aligned(t_1825, t_1826, t_1827, t_1828, pa_x, pb_x, pb_y, lk_1172, lk_1465, \
                         lk_1467, lk_1468, ll_1825, ll_1827, mk_1460, \
                         mk_1468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1825[k] = f_14 * lk_1465[k]
                    + pa_x[k] * ll_1825[k];

        t_1826[k] = f_16 * lk_1172[k]
                    + pb_y[k] * mk_1460[k];

        t_1827[k] = f_14 * lk_1467[k]
                    + pa_x[k] * ll_1827[k];

        t_1828[k] = f_13 * lk_1468[k]
                    + pb_x[k] * mk_1468[k];
    }

#pragma omp simd aligned(t_1829, t_1830, t_1831, t_1832, t_1833, pb_x, lk_1469, lk_1470, \
                         lk_1471, lk_1472, lk_1473, mk_1469, mk_1470, mk_1471, mk_1472, \
                         mk_1473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1829[k] = f_13 * lk_1469[k]
                    + pb_x[k] * mk_1469[k];

        t_1830[k] = f_13 * lk_1470[k]
                    + pb_x[k] * mk_1470[k];

        t_1831[k] = f_13 * lk_1471[k]
                    + pb_x[k] * mk_1471[k];

        t_1832[k] = f_13 * lk_1472[k]
                    + pb_x[k] * mk_1472[k];

        t_1833[k] = f_13 * lk_1473[k]
                    + pb_x[k] * mk_1473[k];
    }

#pragma omp simd aligned(t_1834, t_1835, t_1836, t_1837, t_1838, t_1839, pa_x, pb_x, lk_1474, \
                         lk_1475, ll_1836, ll_1837, ll_1838, ll_1839, mk_1474, \
                         mk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1834[k] = f_13 * lk_1474[k]
                    + pb_x[k] * mk_1474[k];

        t_1835[k] = f_13 * lk_1475[k]
                    + pb_x[k] * mk_1475[k];

        t_1836[k] = pa_x[k] * ll_1836[k];

        t_1837[k] = pa_x[k] * ll_1837[k];

        t_1838[k] = pa_x[k] * ll_1838[k];

        t_1839[k] = pa_x[k] * ll_1839[k];
    }

#pragma omp simd aligned(t_1840, t_1841, t_1842, t_1843, t_1844, t_1845, pa_x, lk_1476, \
                         ll_1840, ll_1841, ll_1842, ll_1843, ll_1844, \
                         ll_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1840[k] = pa_x[k] * ll_1840[k];

        t_1841[k] = pa_x[k] * ll_1841[k];

        t_1842[k] = pa_x[k] * ll_1842[k];

        t_1843[k] = pa_x[k] * ll_1843[k];

        t_1844[k] = pa_x[k] * ll_1844[k];

        t_1845[k] = f_19 * lk_1476[k]
                    + pa_x[k] * ll_1845[k];
    }

#pragma omp simd aligned(t_1846, t_1847, t_1848, t_1849, pa_x, pb_y, pb_z, lk_1152, lk_1188, \
                         lk_1190, lk_1479, ll_1848, mk_1476, mk_1478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1846[k] = f_15 * lk_1188[k]
                    + pb_y[k] * mk_1476[k];

        t_1847[k] = f_17 * lk_1152[k]
                    + pb_z[k] * mk_1476[k];

        t_1848[k] = f_18 * lk_1479[k]
                    + pa_x[k] * ll_1848[k];

        t_1849[k] = f_15 * lk_1190[k]
                    + pb_y[k] * mk_1478[k];
    }

#pragma omp simd aligned(t_1850, t_1851, t_1852, t_1853, pa_x, pb_y, pb_z, lk_1155, lk_1193, \
                         lk_1481, lk_1482, ll_1850, ll_1851, mk_1479, \
                         mk_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1850[k] = f_18 * lk_1481[k]
                    + pa_x[k] * ll_1850[k];

        t_1851[k] = f_17 * lk_1482[k]
                    + pa_x[k] * ll_1851[k];

        t_1852[k] = f_17 * lk_1155[k]
                    + pb_z[k] * mk_1479[k];

        t_1853[k] = f_15 * lk_1193[k]
                    + pb_y[k] * mk_1481[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, t_1857, pa_x, pb_z, lk_1158, lk_1485, \
                         lk_1486, lk_1488, ll_1854, ll_1855, ll_1857, \
                         mk_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_17 * lk_1485[k]
                    + pa_x[k] * ll_1854[k];

        t_1855[k] = f_16 * lk_1486[k]
                    + pa_x[k] * ll_1855[k];

        t_1856[k] = f_17 * lk_1158[k]
                    + pb_z[k] * mk_1482[k];

        t_1857[k] = f_16 * lk_1488[k]
                    + pa_x[k] * ll_1857[k];
    }

#pragma omp simd aligned(t_1858, t_1859, t_1860, t_1861, pa_x, pb_y, pb_z, lk_1162, lk_1197, \
                         lk_1490, lk_1491, ll_1859, ll_1860, mk_1485, \
                         mk_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1858[k] = f_15 * lk_1197[k]
                    + pb_y[k] * mk_1485[k];

        t_1859[k] = f_16 * lk_1490[k]
                    + pa_x[k] * ll_1859[k];

        t_1860[k] = f_15 * lk_1491[k]
                    + pa_x[k] * ll_1860[k];

        t_1861[k] = f_17 * lk_1162[k]
                    + pb_z[k] * mk_1486[k];
    }

#pragma omp simd aligned(t_1862, t_1863, t_1864, t_1865, pa_x, pb_y, lk_1202, lk_1493, \
                         lk_1494, lk_1496, ll_1862, ll_1863, ll_1865, \
                         mk_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1862[k] = f_15 * lk_1493[k]
                    + pa_x[k] * ll_1862[k];

        t_1863[k] = f_15 * lk_1494[k]
                    + pa_x[k] * ll_1863[k];

        t_1864[k] = f_15 * lk_1202[k]
                    + pb_y[k] * mk_1490[k];

        t_1865[k] = f_15 * lk_1496[k]
                    + pa_x[k] * ll_1865[k];
    }

#pragma omp simd aligned(t_1866, t_1867, t_1868, t_1869, pa_x, pb_z, lk_1167, lk_1497, \
                         lk_1499, lk_1500, ll_1866, ll_1868, ll_1869, \
                         mk_1491 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1866[k] = f_14 * lk_1497[k]
                    + pa_x[k] * ll_1866[k];

        t_1867[k] = f_17 * lk_1167[k]
                    + pb_z[k] * mk_1491[k];

        t_1868[k] = f_14 * lk_1499[k]
                    + pa_x[k] * ll_1868[k];

        t_1869[k] = f_14 * lk_1500[k]
                    + pa_x[k] * ll_1869[k];
    }

#pragma omp simd aligned(t_1870, t_1871, t_1872, t_1873, pa_x, pb_x, pb_y, lk_1208, lk_1501, \
                         lk_1503, lk_1504, ll_1870, ll_1872, mk_1496, \
                         mk_1504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1870[k] = f_14 * lk_1501[k]
                    + pa_x[k] * ll_1870[k];

        t_1871[k] = f_15 * lk_1208[k]
                    + pb_y[k] * mk_1496[k];

        t_1872[k] = f_14 * lk_1503[k]
                    + pa_x[k] * ll_1872[k];

        t_1873[k] = f_13 * lk_1504[k]
                    + pb_x[k] * mk_1504[k];
    }

#pragma omp simd aligned(t_1874, t_1875, t_1876, t_1877, t_1878, pb_x, lk_1505, lk_1506, \
                         lk_1507, lk_1508, lk_1509, mk_1505, mk_1506, mk_1507, mk_1508, \
                         mk_1509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1874[k] = f_13 * lk_1505[k]
                    + pb_x[k] * mk_1505[k];

        t_1875[k] = f_13 * lk_1506[k]
                    + pb_x[k] * mk_1506[k];

        t_1876[k] = f_13 * lk_1507[k]
                    + pb_x[k] * mk_1507[k];

        t_1877[k] = f_13 * lk_1508[k]
                    + pb_x[k] * mk_1508[k];

        t_1878[k] = f_13 * lk_1509[k]
                    + pb_x[k] * mk_1509[k];
    }

#pragma omp simd aligned(t_1879, t_1880, t_1881, t_1882, t_1883, t_1884, pa_x, pb_x, lk_1510, \
                         lk_1511, ll_1881, ll_1882, ll_1883, ll_1884, mk_1510, \
                         mk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1879[k] = f_13 * lk_1510[k]
                    + pb_x[k] * mk_1510[k];

        t_1880[k] = f_13 * lk_1511[k]
                    + pb_x[k] * mk_1511[k];

        t_1881[k] = pa_x[k] * ll_1881[k];

        t_1882[k] = pa_x[k] * ll_1882[k];

        t_1883[k] = pa_x[k] * ll_1883[k];

        t_1884[k] = pa_x[k] * ll_1884[k];
    }

#pragma omp simd aligned(t_1885, t_1886, t_1887, t_1888, t_1889, t_1890, pa_x, lk_1512, \
                         ll_1885, ll_1886, ll_1887, ll_1888, ll_1889, \
                         ll_1890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1885[k] = pa_x[k] * ll_1885[k];

        t_1886[k] = pa_x[k] * ll_1886[k];

        t_1887[k] = pa_x[k] * ll_1887[k];

        t_1888[k] = pa_x[k] * ll_1888[k];

        t_1889[k] = pa_x[k] * ll_1889[k];

        t_1890[k] = f_19 * lk_1512[k]
                    + pa_x[k] * ll_1890[k];
    }

#pragma omp simd aligned(t_1891, t_1892, t_1893, t_1894, pa_x, pb_y, pb_z, lk_1188, lk_1224, \
                         lk_1226, lk_1515, ll_1893, mk_1512, mk_1514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1891[k] = f_14 * lk_1224[k]
                    + pb_y[k] * mk_1512[k];

        t_1892[k] = f_18 * lk_1188[k]
                    + pb_z[k] * mk_1512[k];

        t_1893[k] = f_18 * lk_1515[k]
                    + pa_x[k] * ll_1893[k];

        t_1894[k] = f_14 * lk_1226[k]
                    + pb_y[k] * mk_1514[k];
    }

#pragma omp simd aligned(t_1895, t_1896, t_1897, t_1898, pa_x, pb_y, pb_z, lk_1191, lk_1229, \
                         lk_1517, lk_1518, ll_1895, ll_1896, mk_1515, \
                         mk_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1895[k] = f_18 * lk_1517[k]
                    + pa_x[k] * ll_1895[k];

        t_1896[k] = f_17 * lk_1518[k]
                    + pa_x[k] * ll_1896[k];

        t_1897[k] = f_18 * lk_1191[k]
                    + pb_z[k] * mk_1515[k];

        t_1898[k] = f_14 * lk_1229[k]
                    + pb_y[k] * mk_1517[k];
    }

#pragma omp simd aligned(t_1899, t_1900, t_1901, t_1902, pa_x, pb_z, lk_1194, lk_1521, \
                         lk_1522, lk_1524, ll_1899, ll_1900, ll_1902, \
                         mk_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1899[k] = f_17 * lk_1521[k]
                    + pa_x[k] * ll_1899[k];

        t_1900[k] = f_16 * lk_1522[k]
                    + pa_x[k] * ll_1900[k];

        t_1901[k] = f_18 * lk_1194[k]
                    + pb_z[k] * mk_1518[k];

        t_1902[k] = f_16 * lk_1524[k]
                    + pa_x[k] * ll_1902[k];
    }

#pragma omp simd aligned(t_1903, t_1904, t_1905, t_1906, pa_x, pb_y, pb_z, lk_1198, lk_1233, \
                         lk_1526, lk_1527, ll_1904, ll_1905, mk_1521, \
                         mk_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1903[k] = f_14 * lk_1233[k]
                    + pb_y[k] * mk_1521[k];

        t_1904[k] = f_16 * lk_1526[k]
                    + pa_x[k] * ll_1904[k];

        t_1905[k] = f_15 * lk_1527[k]
                    + pa_x[k] * ll_1905[k];

        t_1906[k] = f_18 * lk_1198[k]
                    + pb_z[k] * mk_1522[k];
    }

#pragma omp simd aligned(t_1907, t_1908, t_1909, t_1910, pa_x, pb_y, lk_1238, lk_1529, \
                         lk_1530, lk_1532, ll_1907, ll_1908, ll_1910, \
                         mk_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1907[k] = f_15 * lk_1529[k]
                    + pa_x[k] * ll_1907[k];

        t_1908[k] = f_15 * lk_1530[k]
                    + pa_x[k] * ll_1908[k];

        t_1909[k] = f_14 * lk_1238[k]
                    + pb_y[k] * mk_1526[k];

        t_1910[k] = f_15 * lk_1532[k]
                    + pa_x[k] * ll_1910[k];
    }

#pragma omp simd aligned(t_1911, t_1912, t_1913, t_1914, pa_x, pb_z, lk_1203, lk_1533, \
                         lk_1535, lk_1536, ll_1911, ll_1913, ll_1914, \
                         mk_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1911[k] = f_14 * lk_1533[k]
                    + pa_x[k] * ll_1911[k];

        t_1912[k] = f_18 * lk_1203[k]
                    + pb_z[k] * mk_1527[k];

        t_1913[k] = f_14 * lk_1535[k]
                    + pa_x[k] * ll_1913[k];

        t_1914[k] = f_14 * lk_1536[k]
                    + pa_x[k] * ll_1914[k];
    }

#pragma omp simd aligned(t_1915, t_1916, t_1917, t_1918, pa_x, pb_x, pb_y, lk_1244, lk_1537, \
                         lk_1539, lk_1540, ll_1915, ll_1917, mk_1532, \
                         mk_1540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1915[k] = f_14 * lk_1537[k]
                    + pa_x[k] * ll_1915[k];

        t_1916[k] = f_14 * lk_1244[k]
                    + pb_y[k] * mk_1532[k];

        t_1917[k] = f_14 * lk_1539[k]
                    + pa_x[k] * ll_1917[k];

        t_1918[k] = f_13 * lk_1540[k]
                    + pb_x[k] * mk_1540[k];
    }

#pragma omp simd aligned(t_1919, t_1920, t_1921, t_1922, t_1923, pb_x, lk_1541, lk_1542, \
                         lk_1543, lk_1544, lk_1545, mk_1541, mk_1542, mk_1543, mk_1544, \
                         mk_1545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1919[k] = f_13 * lk_1541[k]
                    + pb_x[k] * mk_1541[k];

        t_1920[k] = f_13 * lk_1542[k]
                    + pb_x[k] * mk_1542[k];

        t_1921[k] = f_13 * lk_1543[k]
                    + pb_x[k] * mk_1543[k];

        t_1922[k] = f_13 * lk_1544[k]
                    + pb_x[k] * mk_1544[k];

        t_1923[k] = f_13 * lk_1545[k]
                    + pb_x[k] * mk_1545[k];
    }

#pragma omp simd aligned(t_1924, t_1925, t_1926, t_1927, t_1928, t_1929, pa_x, pb_x, lk_1546, \
                         lk_1547, ll_1926, ll_1927, ll_1928, ll_1929, mk_1546, \
                         mk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1924[k] = f_13 * lk_1546[k]
                    + pb_x[k] * mk_1546[k];

        t_1925[k] = f_13 * lk_1547[k]
                    + pb_x[k] * mk_1547[k];

        t_1926[k] = pa_x[k] * ll_1926[k];

        t_1927[k] = pa_x[k] * ll_1927[k];

        t_1928[k] = pa_x[k] * ll_1928[k];

        t_1929[k] = pa_x[k] * ll_1929[k];
    }

#pragma omp simd aligned(t_1930, t_1931, t_1932, t_1933, t_1934, t_1935, pa_x, pa_y, ll_1575, \
                         ll_1930, ll_1931, ll_1932, ll_1933, ll_1934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1930[k] = pa_x[k] * ll_1930[k];

        t_1931[k] = pa_x[k] * ll_1931[k];

        t_1932[k] = pa_x[k] * ll_1932[k];

        t_1933[k] = pa_x[k] * ll_1933[k];

        t_1934[k] = pa_x[k] * ll_1934[k];

        t_1935[k] = pa_y[k] * ll_1575[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece15(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t lk,
                                             const size_t ll, const size_t mi0, const size_t mi1,
                                             const size_t mk, const size_t ncols,
                                             const double alpha, const double beta,
                                             const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_22 = 3.5 / p;

    auto *t_1936 = buffer.data(target + 1936);
    auto *t_1937 = buffer.data(target + 1937);
    auto *t_1938 = buffer.data(target + 1938);
    auto *t_1939 = buffer.data(target + 1939);
    auto *t_1940 = buffer.data(target + 1940);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lk_1227 = buffer.data(lk + 1227);
    const auto *lk_1230 = buffer.data(lk + 1230);
    const auto *lk_1234 = buffer.data(lk + 1234);
    const auto *lk_1239 = buffer.data(lk + 1239);
    const auto *lk_1260 = buffer.data(lk + 1260);
    const auto *lk_1262 = buffer.data(lk + 1262);
    const auto *lk_1263 = buffer.data(lk + 1263);
    const auto *lk_1265 = buffer.data(lk + 1265);
    const auto *lk_1266 = buffer.data(lk + 1266);
    const auto *lk_1269 = buffer.data(lk + 1269);
    const auto *lk_1270 = buffer.data(lk + 1270);
    const auto *lk_1274 = buffer.data(lk + 1274);
    const auto *lk_1275 = buffer.data(lk + 1275);
    const auto *lk_1280 = buffer.data(lk + 1280);
    const auto *lk_1296 = buffer.data(lk + 1296);
    const auto *lk_1298 = buffer.data(lk + 1298);
    const auto *lk_1299 = buffer.data(lk + 1299);
    const auto *lk_1301 = buffer.data(lk + 1301);
    const auto *lk_1305 = buffer.data(lk + 1305);
    const auto *lk_1310 = buffer.data(lk + 1310);
    const auto *lk_1316 = buffer.data(lk + 1316);
    const auto *lk_1324 = buffer.data(lk + 1324);
    const auto *lk_1331 = buffer.data(lk + 1331);
    const auto *lk_1334 = buffer.data(lk + 1334);
    const auto *lk_1337 = buffer.data(lk + 1337);
    const auto *lk_1551 = buffer.data(lk + 1551);
    const auto *lk_1554 = buffer.data(lk + 1554);
    const auto *lk_1558 = buffer.data(lk + 1558);
    const auto *lk_1560 = buffer.data(lk + 1560);
    const auto *lk_1563 = buffer.data(lk + 1563);
    const auto *lk_1565 = buffer.data(lk + 1565);
    const auto *lk_1566 = buffer.data(lk + 1566);
    const auto *lk_1569 = buffer.data(lk + 1569);
    const auto *lk_1571 = buffer.data(lk + 1571);
    const auto *lk_1572 = buffer.data(lk + 1572);
    const auto *lk_1573 = buffer.data(lk + 1573);
    const auto *lk_1576 = buffer.data(lk + 1576);
    const auto *lk_1577 = buffer.data(lk + 1577);
    const auto *lk_1578 = buffer.data(lk + 1578);
    const auto *lk_1579 = buffer.data(lk + 1579);
    const auto *lk_1580 = buffer.data(lk + 1580);
    const auto *lk_1581 = buffer.data(lk + 1581);
    const auto *lk_1582 = buffer.data(lk + 1582);
    const auto *lk_1584 = buffer.data(lk + 1584);
    const auto *lk_1587 = buffer.data(lk + 1587);
    const auto *lk_1589 = buffer.data(lk + 1589);
    const auto *lk_1590 = buffer.data(lk + 1590);
    const auto *lk_1593 = buffer.data(lk + 1593);
    const auto *lk_1594 = buffer.data(lk + 1594);
    const auto *lk_1596 = buffer.data(lk + 1596);
    const auto *lk_1598 = buffer.data(lk + 1598);
    const auto *lk_1599 = buffer.data(lk + 1599);
    const auto *lk_1601 = buffer.data(lk + 1601);
    const auto *lk_1602 = buffer.data(lk + 1602);
    const auto *lk_1604 = buffer.data(lk + 1604);
    const auto *lk_1605 = buffer.data(lk + 1605);
    const auto *lk_1607 = buffer.data(lk + 1607);
    const auto *lk_1608 = buffer.data(lk + 1608);
    const auto *lk_1609 = buffer.data(lk + 1609);
    const auto *lk_1611 = buffer.data(lk + 1611);
    const auto *lk_1612 = buffer.data(lk + 1612);
    const auto *lk_1613 = buffer.data(lk + 1613);
    const auto *lk_1614 = buffer.data(lk + 1614);
    const auto *lk_1615 = buffer.data(lk + 1615);
    const auto *lk_1616 = buffer.data(lk + 1616);
    const auto *lk_1617 = buffer.data(lk + 1617);
    const auto *lk_1619 = buffer.data(lk + 1619);

    const auto *ll_1577 = buffer.data(ll + 1577);
    const auto *ll_1580 = buffer.data(ll + 1580);
    const auto *ll_1584 = buffer.data(ll + 1584);
    const auto *ll_1589 = buffer.data(ll + 1589);
    const auto *ll_1595 = buffer.data(ll + 1595);
    const auto *ll_1602 = buffer.data(ll + 1602);
    const auto *ll_1610 = buffer.data(ll + 1610);
    const auto *ll_1620 = buffer.data(ll + 1620);
    const auto *ll_1621 = buffer.data(ll + 1621);
    const auto *ll_1623 = buffer.data(ll + 1623);
    const auto *ll_1625 = buffer.data(ll + 1625);
    const auto *ll_1626 = buffer.data(ll + 1626);
    const auto *ll_1938 = buffer.data(ll + 1938);
    const auto *ll_1941 = buffer.data(ll + 1941);
    const auto *ll_1945 = buffer.data(ll + 1945);
    const auto *ll_1947 = buffer.data(ll + 1947);
    const auto *ll_1950 = buffer.data(ll + 1950);
    const auto *ll_1952 = buffer.data(ll + 1952);
    const auto *ll_1953 = buffer.data(ll + 1953);
    const auto *ll_1956 = buffer.data(ll + 1956);
    const auto *ll_1958 = buffer.data(ll + 1958);
    const auto *ll_1959 = buffer.data(ll + 1959);
    const auto *ll_1960 = buffer.data(ll + 1960);
    const auto *ll_1971 = buffer.data(ll + 1971);
    const auto *ll_1972 = buffer.data(ll + 1972);
    const auto *ll_1973 = buffer.data(ll + 1973);
    const auto *ll_1974 = buffer.data(ll + 1974);
    const auto *ll_1975 = buffer.data(ll + 1975);
    const auto *ll_1976 = buffer.data(ll + 1976);
    const auto *ll_1977 = buffer.data(ll + 1977);
    const auto *ll_1978 = buffer.data(ll + 1978);
    const auto *ll_1979 = buffer.data(ll + 1979);
    const auto *ll_1980 = buffer.data(ll + 1980);
    const auto *ll_1983 = buffer.data(ll + 1983);
    const auto *ll_1985 = buffer.data(ll + 1985);
    const auto *ll_1986 = buffer.data(ll + 1986);
    const auto *ll_1989 = buffer.data(ll + 1989);
    const auto *ll_1990 = buffer.data(ll + 1990);
    const auto *ll_1992 = buffer.data(ll + 1992);
    const auto *ll_1994 = buffer.data(ll + 1994);
    const auto *ll_1995 = buffer.data(ll + 1995);
    const auto *ll_1997 = buffer.data(ll + 1997);
    const auto *ll_1998 = buffer.data(ll + 1998);
    const auto *ll_2000 = buffer.data(ll + 2000);
    const auto *ll_2001 = buffer.data(ll + 2001);
    const auto *ll_2003 = buffer.data(ll + 2003);
    const auto *ll_2004 = buffer.data(ll + 2004);
    const auto *ll_2005 = buffer.data(ll + 2005);
    const auto *ll_2007 = buffer.data(ll + 2007);
    const auto *ll_2016 = buffer.data(ll + 2016);
    const auto *ll_2017 = buffer.data(ll + 2017);
    const auto *ll_2018 = buffer.data(ll + 2018);
    const auto *ll_2019 = buffer.data(ll + 2019);
    const auto *ll_2020 = buffer.data(ll + 2020);
    const auto *ll_2021 = buffer.data(ll + 2021);
    const auto *ll_2022 = buffer.data(ll + 2022);
    const auto *ll_2024 = buffer.data(ll + 2024);

    const auto *mi0_1260 = buffer.data(mi0 + 1260);
    const auto *mi0_1263 = buffer.data(mi0 + 1263);
    const auto *mi0_1265 = buffer.data(mi0 + 1265);
    const auto *mi0_1266 = buffer.data(mi0 + 1266);
    const auto *mi0_1269 = buffer.data(mi0 + 1269);
    const auto *mi0_1270 = buffer.data(mi0 + 1270);
    const auto *mi0_1272 = buffer.data(mi0 + 1272);
    const auto *mi0_1274 = buffer.data(mi0 + 1274);
    const auto *mi0_1275 = buffer.data(mi0 + 1275);
    const auto *mi0_1277 = buffer.data(mi0 + 1277);
    const auto *mi0_1278 = buffer.data(mi0 + 1278);
    const auto *mi0_1280 = buffer.data(mi0 + 1280);
    const auto *mi0_1281 = buffer.data(mi0 + 1281);
    const auto *mi0_1282 = buffer.data(mi0 + 1282);
    const auto *mi0_1283 = buffer.data(mi0 + 1283);
    const auto *mi0_1284 = buffer.data(mi0 + 1284);
    const auto *mi0_1285 = buffer.data(mi0 + 1285);
    const auto *mi0_1287 = buffer.data(mi0 + 1287);

    const auto *mi1_1260 = buffer.data(mi1 + 1260);
    const auto *mi1_1263 = buffer.data(mi1 + 1263);
    const auto *mi1_1265 = buffer.data(mi1 + 1265);
    const auto *mi1_1266 = buffer.data(mi1 + 1266);
    const auto *mi1_1269 = buffer.data(mi1 + 1269);
    const auto *mi1_1270 = buffer.data(mi1 + 1270);
    const auto *mi1_1272 = buffer.data(mi1 + 1272);
    const auto *mi1_1274 = buffer.data(mi1 + 1274);
    const auto *mi1_1275 = buffer.data(mi1 + 1275);
    const auto *mi1_1277 = buffer.data(mi1 + 1277);
    const auto *mi1_1278 = buffer.data(mi1 + 1278);
    const auto *mi1_1280 = buffer.data(mi1 + 1280);
    const auto *mi1_1281 = buffer.data(mi1 + 1281);
    const auto *mi1_1282 = buffer.data(mi1 + 1282);
    const auto *mi1_1283 = buffer.data(mi1 + 1283);
    const auto *mi1_1284 = buffer.data(mi1 + 1284);
    const auto *mi1_1285 = buffer.data(mi1 + 1285);
    const auto *mi1_1287 = buffer.data(mi1 + 1287);

    const auto *mk_1548 = buffer.data(mk + 1548);
    const auto *mk_1550 = buffer.data(mk + 1550);
    const auto *mk_1551 = buffer.data(mk + 1551);
    const auto *mk_1553 = buffer.data(mk + 1553);
    const auto *mk_1554 = buffer.data(mk + 1554);
    const auto *mk_1557 = buffer.data(mk + 1557);
    const auto *mk_1558 = buffer.data(mk + 1558);
    const auto *mk_1562 = buffer.data(mk + 1562);
    const auto *mk_1563 = buffer.data(mk + 1563);
    const auto *mk_1568 = buffer.data(mk + 1568);
    const auto *mk_1576 = buffer.data(mk + 1576);
    const auto *mk_1577 = buffer.data(mk + 1577);
    const auto *mk_1578 = buffer.data(mk + 1578);
    const auto *mk_1579 = buffer.data(mk + 1579);
    const auto *mk_1580 = buffer.data(mk + 1580);
    const auto *mk_1581 = buffer.data(mk + 1581);
    const auto *mk_1582 = buffer.data(mk + 1582);
    const auto *mk_1584 = buffer.data(mk + 1584);
    const auto *mk_1586 = buffer.data(mk + 1586);
    const auto *mk_1587 = buffer.data(mk + 1587);
    const auto *mk_1589 = buffer.data(mk + 1589);
    const auto *mk_1590 = buffer.data(mk + 1590);
    const auto *mk_1593 = buffer.data(mk + 1593);
    const auto *mk_1594 = buffer.data(mk + 1594);
    const auto *mk_1598 = buffer.data(mk + 1598);
    const auto *mk_1599 = buffer.data(mk + 1599);
    const auto *mk_1604 = buffer.data(mk + 1604);
    const auto *mk_1611 = buffer.data(mk + 1611);
    const auto *mk_1612 = buffer.data(mk + 1612);
    const auto *mk_1613 = buffer.data(mk + 1613);
    const auto *mk_1614 = buffer.data(mk + 1614);
    const auto *mk_1615 = buffer.data(mk + 1615);
    const auto *mk_1616 = buffer.data(mk + 1616);
    const auto *mk_1617 = buffer.data(mk + 1617);
    const auto *mk_1619 = buffer.data(mk + 1619);
    const auto *mk_1620 = buffer.data(mk + 1620);
    const auto *mk_1621 = buffer.data(mk + 1621);
    const auto *mk_1623 = buffer.data(mk + 1623);
    const auto *mk_1625 = buffer.data(mk + 1625);
    const auto *mk_1626 = buffer.data(mk + 1626);
    const auto *mk_1629 = buffer.data(mk + 1629);
    const auto *mk_1630 = buffer.data(mk + 1630);
    const auto *mk_1632 = buffer.data(mk + 1632);
    const auto *mk_1634 = buffer.data(mk + 1634);
    const auto *mk_1635 = buffer.data(mk + 1635);
    const auto *mk_1637 = buffer.data(mk + 1637);
    const auto *mk_1638 = buffer.data(mk + 1638);
    const auto *mk_1640 = buffer.data(mk + 1640);
    const auto *mk_1641 = buffer.data(mk + 1641);
    const auto *mk_1643 = buffer.data(mk + 1643);
    const auto *mk_1644 = buffer.data(mk + 1644);
    const auto *mk_1645 = buffer.data(mk + 1645);
    const auto *mk_1647 = buffer.data(mk + 1647);
    const auto *mk_1648 = buffer.data(mk + 1648);
    const auto *mk_1649 = buffer.data(mk + 1649);
    const auto *mk_1650 = buffer.data(mk + 1650);
    const auto *mk_1651 = buffer.data(mk + 1651);
    const auto *mk_1652 = buffer.data(mk + 1652);
    const auto *mk_1653 = buffer.data(mk + 1653);
    const auto *mk_1654 = buffer.data(mk + 1654);
    const auto *mk_1655 = buffer.data(mk + 1655);
    const auto *mk_1656 = buffer.data(mk + 1656);
    const auto *mk_1658 = buffer.data(mk + 1658);
    const auto *mk_1659 = buffer.data(mk + 1659);
    const auto *mk_1661 = buffer.data(mk + 1661);

#pragma omp simd aligned(t_1936, t_1937, t_1938, t_1939, t_1940, pa_x, pa_y, pb_y, lk_1260, \
                         lk_1262, lk_1551, ll_1577, ll_1580, ll_1938, mk_1548, \
                         mk_1550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1936[k] = f_13 * lk_1260[k]
                    + pb_y[k] * mk_1548[k];

        t_1937[k] = pa_y[k] * ll_1577[k];

        t_1938[k] = f_18 * lk_1551[k]
                    + pa_x[k] * ll_1938[k];

        t_1939[k] = f_13 * lk_1262[k]
                    + pb_y[k] * mk_1550[k];

        t_1940[k] = pa_y[k] * ll_1580[k];
    }

#pragma omp simd aligned(t_1941, t_1942, t_1943, t_1944, pa_x, pa_y, pb_y, pb_z, lk_1227, \
                         lk_1265, lk_1554, ll_1584, ll_1941, mk_1551, \
                         mk_1553 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1941[k] = f_17 * lk_1554[k]
                    + pa_x[k] * ll_1941[k];

        t_1942[k] = f_22 * lk_1227[k]
                    + pb_z[k] * mk_1551[k];

        t_1943[k] = f_13 * lk_1265[k]
                    + pb_y[k] * mk_1553[k];

        t_1944[k] = pa_y[k] * ll_1584[k];
    }

#pragma omp simd aligned(t_1945, t_1946, t_1947, t_1948, pa_x, pb_y, pb_z, lk_1230, lk_1269, \
                         lk_1558, lk_1560, ll_1945, ll_1947, mk_1554, \
                         mk_1557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1945[k] = f_16 * lk_1558[k]
                    + pa_x[k] * ll_1945[k];

        t_1946[k] = f_22 * lk_1230[k]
                    + pb_z[k] * mk_1554[k];

        t_1947[k] = f_16 * lk_1560[k]
                    + pa_x[k] * ll_1947[k];

        t_1948[k] = f_13 * lk_1269[k]
                    + pb_y[k] * mk_1557[k];
    }

#pragma omp simd aligned(t_1949, t_1950, t_1951, t_1952, pa_x, pa_y, pb_z, lk_1234, lk_1563, \
                         lk_1565, ll_1589, ll_1950, ll_1952, mk_1558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1949[k] = pa_y[k] * ll_1589[k];

        t_1950[k] = f_15 * lk_1563[k]
                    + pa_x[k] * ll_1950[k];

        t_1951[k] = f_22 * lk_1234[k]
                    + pb_z[k] * mk_1558[k];

        t_1952[k] = f_15 * lk_1565[k]
                    + pa_x[k] * ll_1952[k];
    }

#pragma omp simd aligned(t_1953, t_1954, t_1955, t_1956, pa_x, pa_y, pb_y, lk_1274, lk_1566, \
                         lk_1569, ll_1595, ll_1953, ll_1956, mk_1562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1953[k] = f_15 * lk_1566[k]
                    + pa_x[k] * ll_1953[k];

        t_1954[k] = f_13 * lk_1274[k]
                    + pb_y[k] * mk_1562[k];

        t_1955[k] = pa_y[k] * ll_1595[k];

        t_1956[k] = f_14 * lk_1569[k]
                    + pa_x[k] * ll_1956[k];
    }

#pragma omp simd aligned(t_1957, t_1958, t_1959, t_1960, pa_x, pb_z, lk_1239, lk_1571, \
                         lk_1572, lk_1573, ll_1958, ll_1959, ll_1960, \
                         mk_1563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1957[k] = f_22 * lk_1239[k]
                    + pb_z[k] * mk_1563[k];

        t_1958[k] = f_14 * lk_1571[k]
                    + pa_x[k] * ll_1958[k];

        t_1959[k] = f_14 * lk_1572[k]
                    + pa_x[k] * ll_1959[k];

        t_1960[k] = f_14 * lk_1573[k]
                    + pa_x[k] * ll_1960[k];
    }

#pragma omp simd aligned(t_1961, t_1962, t_1963, t_1964, pa_y, pb_x, pb_y, lk_1280, lk_1576, \
                         lk_1577, ll_1602, mk_1568, mk_1576, mk_1577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1961[k] = f_13 * lk_1280[k]
                    + pb_y[k] * mk_1568[k];

        t_1962[k] = pa_y[k] * ll_1602[k];

        t_1963[k] = f_13 * lk_1576[k]
                    + pb_x[k] * mk_1576[k];

        t_1964[k] = f_13 * lk_1577[k]
                    + pb_x[k] * mk_1577[k];
    }

#pragma omp simd aligned(t_1965, t_1966, t_1967, t_1968, t_1969, pb_x, lk_1578, lk_1579, \
                         lk_1580, lk_1581, lk_1582, mk_1578, mk_1579, mk_1580, mk_1581, \
                         mk_1582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1965[k] = f_13 * lk_1578[k]
                    + pb_x[k] * mk_1578[k];

        t_1966[k] = f_13 * lk_1579[k]
                    + pb_x[k] * mk_1579[k];

        t_1967[k] = f_13 * lk_1580[k]
                    + pb_x[k] * mk_1580[k];

        t_1968[k] = f_13 * lk_1581[k]
                    + pb_x[k] * mk_1581[k];

        t_1969[k] = f_13 * lk_1582[k]
                    + pb_x[k] * mk_1582[k];
    }

#pragma omp simd aligned(t_1970, t_1971, t_1972, t_1973, t_1974, t_1975, t_1976, pa_x, pa_y, \
                         ll_1610, ll_1971, ll_1972, ll_1973, ll_1974, ll_1975, \
                         ll_1976 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1970[k] = pa_y[k] * ll_1610[k];

        t_1971[k] = pa_x[k] * ll_1971[k];

        t_1972[k] = pa_x[k] * ll_1972[k];

        t_1973[k] = pa_x[k] * ll_1973[k];

        t_1974[k] = pa_x[k] * ll_1974[k];

        t_1975[k] = pa_x[k] * ll_1975[k];

        t_1976[k] = pa_x[k] * ll_1976[k];
    }

#pragma omp simd aligned(t_1977, t_1978, t_1979, t_1980, t_1981, t_1982, pa_x, pb_y, pb_z, \
                         lk_1260, lk_1584, ll_1977, ll_1978, ll_1979, ll_1980, \
                         mk_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1977[k] = pa_x[k] * ll_1977[k];

        t_1978[k] = pa_x[k] * ll_1978[k];

        t_1979[k] = pa_x[k] * ll_1979[k];

        t_1980[k] = f_19 * lk_1584[k]
                    + pa_x[k] * ll_1980[k];

        t_1981[k] = pb_y[k] * mk_1584[k];

        t_1982[k] = f_19 * lk_1260[k]
                    + pb_z[k] * mk_1584[k];
    }

#pragma omp simd aligned(t_1983, t_1984, t_1985, t_1986, pa_x, pb_y, lk_1587, lk_1589, \
                         lk_1590, ll_1983, ll_1985, ll_1986, mk_1586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1983[k] = f_18 * lk_1587[k]
                    + pa_x[k] * ll_1983[k];

        t_1984[k] = pb_y[k] * mk_1586[k];

        t_1985[k] = f_18 * lk_1589[k]
                    + pa_x[k] * ll_1985[k];

        t_1986[k] = f_17 * lk_1590[k]
                    + pa_x[k] * ll_1986[k];
    }

#pragma omp simd aligned(t_1987, t_1988, t_1989, t_1990, pa_x, pb_y, pb_z, lk_1263, lk_1593, \
                         lk_1594, ll_1989, ll_1990, mk_1587, mk_1589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1987[k] = f_19 * lk_1263[k]
                    + pb_z[k] * mk_1587[k];

        t_1988[k] = pb_y[k] * mk_1589[k];

        t_1989[k] = f_17 * lk_1593[k]
                    + pa_x[k] * ll_1989[k];

        t_1990[k] = f_16 * lk_1594[k]
                    + pa_x[k] * ll_1990[k];
    }

#pragma omp simd aligned(t_1991, t_1992, t_1993, t_1994, pa_x, pb_y, pb_z, lk_1266, lk_1596, \
                         lk_1598, ll_1992, ll_1994, mk_1590, mk_1593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1991[k] = f_19 * lk_1266[k]
                    + pb_z[k] * mk_1590[k];

        t_1992[k] = f_16 * lk_1596[k]
                    + pa_x[k] * ll_1992[k];

        t_1993[k] = pb_y[k] * mk_1593[k];

        t_1994[k] = f_16 * lk_1598[k]
                    + pa_x[k] * ll_1994[k];
    }

#pragma omp simd aligned(t_1995, t_1996, t_1997, t_1998, pa_x, pb_z, lk_1270, lk_1599, \
                         lk_1601, lk_1602, ll_1995, ll_1997, ll_1998, \
                         mk_1594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1995[k] = f_15 * lk_1599[k]
                    + pa_x[k] * ll_1995[k];

        t_1996[k] = f_19 * lk_1270[k]
                    + pb_z[k] * mk_1594[k];

        t_1997[k] = f_15 * lk_1601[k]
                    + pa_x[k] * ll_1997[k];

        t_1998[k] = f_15 * lk_1602[k]
                    + pa_x[k] * ll_1998[k];
    }

#pragma omp simd aligned(t_1999, t_2000, t_2001, t_2002, pa_x, pb_y, pb_z, lk_1275, lk_1604, \
                         lk_1605, ll_2000, ll_2001, mk_1598, mk_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1999[k] = pb_y[k] * mk_1598[k];

        t_2000[k] = f_15 * lk_1604[k]
                    + pa_x[k] * ll_2000[k];

        t_2001[k] = f_14 * lk_1605[k]
                    + pa_x[k] * ll_2001[k];

        t_2002[k] = f_19 * lk_1275[k]
                    + pb_z[k] * mk_1599[k];
    }

#pragma omp simd aligned(t_2003, t_2004, t_2005, t_2006, t_2007, pa_x, pb_y, lk_1607, lk_1608, \
                         lk_1609, lk_1611, ll_2003, ll_2004, ll_2005, ll_2007, \
                         mk_1604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2003[k] = f_14 * lk_1607[k]
                    + pa_x[k] * ll_2003[k];

        t_2004[k] = f_14 * lk_1608[k]
                    + pa_x[k] * ll_2004[k];

        t_2005[k] = f_14 * lk_1609[k]
                    + pa_x[k] * ll_2005[k];

        t_2006[k] = pb_y[k] * mk_1604[k];

        t_2007[k] = f_14 * lk_1611[k]
                    + pa_x[k] * ll_2007[k];
    }

#pragma omp simd aligned(t_2008, t_2009, t_2010, t_2011, t_2012, pb_x, lk_1612, lk_1613, \
                         lk_1614, lk_1615, lk_1616, mk_1612, mk_1613, mk_1614, mk_1615, \
                         mk_1616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2008[k] = f_13 * lk_1612[k]
                    + pb_x[k] * mk_1612[k];

        t_2009[k] = f_13 * lk_1613[k]
                    + pb_x[k] * mk_1613[k];

        t_2010[k] = f_13 * lk_1614[k]
                    + pb_x[k] * mk_1614[k];

        t_2011[k] = f_13 * lk_1615[k]
                    + pb_x[k] * mk_1615[k];

        t_2012[k] = f_13 * lk_1616[k]
                    + pb_x[k] * mk_1616[k];
    }

#pragma omp simd aligned(t_2013, t_2014, t_2015, t_2016, t_2017, pa_x, pb_x, pb_y, lk_1617, \
                         lk_1619, ll_2016, ll_2017, mk_1611, mk_1617, \
                         mk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2013[k] = f_13 * lk_1617[k]
                    + pb_x[k] * mk_1617[k];

        t_2014[k] = pb_y[k] * mk_1611[k];

        t_2015[k] = f_13 * lk_1619[k]
                    + pb_x[k] * mk_1619[k];

        t_2016[k] = pa_x[k] * ll_2016[k];

        t_2017[k] = pa_x[k] * ll_2017[k];
    }

#pragma omp simd aligned(t_2018, t_2019, t_2020, t_2021, t_2022, t_2023, t_2024, pa_x, pb_y, \
                         ll_2018, ll_2019, ll_2020, ll_2021, ll_2022, ll_2024, \
                         mk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2018[k] = pa_x[k] * ll_2018[k];

        t_2019[k] = pa_x[k] * ll_2019[k];

        t_2020[k] = pa_x[k] * ll_2020[k];

        t_2021[k] = pa_x[k] * ll_2021[k];

        t_2022[k] = pa_x[k] * ll_2022[k];

        t_2023[k] = pb_y[k] * mk_1619[k];

        t_2024[k] = pa_x[k] * ll_2024[k];
    }

#pragma omp simd aligned(t_2025, t_2026, t_2027, t_2028, t_2029, pb_x, pb_y, pb_z, lk_1296, \
                         mi0_1260, mi0_1263, mi1_1260, mi1_1263, mk_1620, mk_1621, \
                         mk_1623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2025[k] = f_1 * mi0_1260[k]
                    - f_2 * mi1_1260[k]
                    + pb_x[k] * mk_1620[k];

        t_2026[k] = f_0 * lk_1296[k]
                    + pb_y[k] * mk_1620[k];

        t_2027[k] = pb_z[k] * mk_1620[k];

        t_2028[k] = f_11 * mi0_1263[k]
                    - f_12 * mi1_1263[k]
                    + pb_x[k] * mk_1623[k];

        t_2029[k] = pb_z[k] * mk_1621[k];
    }

#pragma omp simd aligned(t_2030, t_2031, t_2032, t_2033, pb_x, pb_y, pb_z, lk_1301, mi0_1265, \
                         mi0_1266, mi1_1265, mi1_1266, mk_1623, mk_1625, \
                         mk_1626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2030[k] = f_11 * mi0_1265[k]
                    - f_12 * mi1_1265[k]
                    + pb_x[k] * mk_1625[k];

        t_2031[k] = f_9 * mi0_1266[k]
                    - f_10 * mi1_1266[k]
                    + pb_x[k] * mk_1626[k];

        t_2032[k] = pb_z[k] * mk_1623[k];

        t_2033[k] = f_0 * lk_1301[k]
                    + pb_y[k] * mk_1625[k];
    }

#pragma omp simd aligned(t_2034, t_2035, t_2036, t_2037, pb_x, pb_z, mi0_1269, mi0_1270, \
                         mi0_1272, mi1_1269, mi1_1270, mi1_1272, mk_1626, mk_1629, mk_1630, \
                         mk_1632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2034[k] = f_9 * mi0_1269[k]
                    - f_10 * mi1_1269[k]
                    + pb_x[k] * mk_1629[k];

        t_2035[k] = f_7 * mi0_1270[k]
                    - f_8 * mi1_1270[k]
                    + pb_x[k] * mk_1630[k];

        t_2036[k] = pb_z[k] * mk_1626[k];

        t_2037[k] = f_7 * mi0_1272[k]
                    - f_8 * mi1_1272[k]
                    + pb_x[k] * mk_1632[k];
    }

#pragma omp simd aligned(t_2038, t_2039, t_2040, t_2041, pb_x, pb_y, pb_z, lk_1305, mi0_1274, \
                         mi0_1275, mi1_1274, mi1_1275, mk_1629, mk_1630, mk_1634, \
                         mk_1635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2038[k] = f_0 * lk_1305[k]
                    + pb_y[k] * mk_1629[k];

        t_2039[k] = f_7 * mi0_1274[k]
                    - f_8 * mi1_1274[k]
                    + pb_x[k] * mk_1634[k];

        t_2040[k] = f_5 * mi0_1275[k]
                    - f_6 * mi1_1275[k]
                    + pb_x[k] * mk_1635[k];

        t_2041[k] = pb_z[k] * mk_1630[k];
    }

#pragma omp simd aligned(t_2042, t_2043, t_2044, pb_x, pb_y, lk_1310, mi0_1277, mi0_1278, \
                         mi1_1277, mi1_1278, mk_1634, mk_1637, \
                         mk_1638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2042[k] = f_5 * mi0_1277[k]
                    - f_6 * mi1_1277[k]
                    + pb_x[k] * mk_1637[k];

        t_2043[k] = f_5 * mi0_1278[k]
                    - f_6 * mi1_1278[k]
                    + pb_x[k] * mk_1638[k];

        t_2044[k] = f_0 * lk_1310[k]
                    + pb_y[k] * mk_1634[k];
    }

#pragma omp simd aligned(t_2045, t_2046, t_2047, t_2048, pb_x, pb_z, mi0_1280, mi0_1281, \
                         mi0_1283, mi1_1280, mi1_1281, mi1_1283, mk_1635, mk_1640, mk_1641, \
                         mk_1643 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2045[k] = f_5 * mi0_1280[k]
                    - f_6 * mi1_1280[k]
                    + pb_x[k] * mk_1640[k];

        t_2046[k] = f_3 * mi0_1281[k]
                    - f_4 * mi1_1281[k]
                    + pb_x[k] * mk_1641[k];

        t_2047[k] = pb_z[k] * mk_1635[k];

        t_2048[k] = f_3 * mi0_1283[k]
                    - f_4 * mi1_1283[k]
                    + pb_x[k] * mk_1643[k];
    }

#pragma omp simd aligned(t_2049, t_2050, t_2051, pb_x, pb_y, lk_1316, mi0_1284, mi0_1285, \
                         mi1_1284, mi1_1285, mk_1640, mk_1644, \
                         mk_1645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2049[k] = f_3 * mi0_1284[k]
                    - f_4 * mi1_1284[k]
                    + pb_x[k] * mk_1644[k];

        t_2050[k] = f_3 * mi0_1285[k]
                    - f_4 * mi1_1285[k]
                    + pb_x[k] * mk_1645[k];

        t_2051[k] = f_0 * lk_1316[k]
                    + pb_y[k] * mk_1640[k];
    }

#pragma omp simd aligned(t_2052, t_2053, t_2054, t_2055, t_2056, t_2057, pb_x, mi0_1287, \
                         mi1_1287, mk_1647, mk_1648, mk_1649, mk_1650, mk_1651, \
                         mk_1652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2052[k] = f_3 * mi0_1287[k]
                    - f_4 * mi1_1287[k]
                    + pb_x[k] * mk_1647[k];

        t_2053[k] = pb_x[k] * mk_1648[k];

        t_2054[k] = pb_x[k] * mk_1649[k];

        t_2055[k] = pb_x[k] * mk_1650[k];

        t_2056[k] = pb_x[k] * mk_1651[k];

        t_2057[k] = pb_x[k] * mk_1652[k];
    }

#pragma omp simd aligned(t_2058, t_2059, t_2060, t_2061, t_2062, pb_x, pb_y, pb_z, lk_1324, \
                         mi0_1281, mi1_1281, mk_1648, mk_1653, mk_1654, \
                         mk_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2058[k] = pb_x[k] * mk_1653[k];

        t_2059[k] = pb_x[k] * mk_1654[k];

        t_2060[k] = pb_x[k] * mk_1655[k];

        t_2061[k] = f_0 * lk_1324[k]
                    + f_1 * mi0_1281[k]
                    - f_2 * mi1_1281[k]
                    + pb_y[k] * mk_1648[k];

        t_2062[k] = pb_z[k] * mk_1648[k];
    }

#pragma omp simd aligned(t_2063, t_2064, t_2065, pb_z, mi0_1281, mi0_1282, mi0_1283, mi1_1281, \
                         mi1_1282, mi1_1283, mk_1649, mk_1650, \
                         mk_1651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2063[k] = f_3 * mi0_1281[k]
                    - f_4 * mi1_1281[k]
                    + pb_z[k] * mk_1649[k];

        t_2064[k] = f_5 * mi0_1282[k]
                    - f_6 * mi1_1282[k]
                    + pb_z[k] * mk_1650[k];

        t_2065[k] = f_7 * mi0_1283[k]
                    - f_8 * mi1_1283[k]
                    + pb_z[k] * mk_1651[k];
    }

#pragma omp simd aligned(t_2066, t_2067, t_2068, t_2069, pb_y, pb_z, lk_1331, mi0_1284, \
                         mi0_1285, mi0_1287, mi1_1284, mi1_1285, mi1_1287, mk_1652, mk_1653, \
                         mk_1655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2066[k] = f_9 * mi0_1284[k]
                    - f_10 * mi1_1284[k]
                    + pb_z[k] * mk_1652[k];

        t_2067[k] = f_11 * mi0_1285[k]
                    - f_12 * mi1_1285[k]
                    + pb_z[k] * mk_1653[k];

        t_2068[k] = f_0 * lk_1331[k]
                    + pb_y[k] * mk_1655[k];

        t_2069[k] = f_1 * mi0_1287[k]
                    - f_2 * mi1_1287[k]
                    + pb_z[k] * mk_1655[k];
    }

#pragma omp simd aligned(t_2070, t_2071, t_2072, t_2073, t_2074, pa_z, pb_y, pb_z, lk_1296, \
                         lk_1334, ll_1620, ll_1621, ll_1623, mk_1656, \
                         mk_1658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2070[k] = pa_z[k] * ll_1620[k];

        t_2071[k] = pa_z[k] * ll_1621[k];

        t_2072[k] = f_13 * lk_1296[k]
                    + pb_z[k] * mk_1656[k];

        t_2073[k] = pa_z[k] * ll_1623[k];

        t_2074[k] = f_19 * lk_1334[k]
                    + pb_y[k] * mk_1658[k];
    }

#pragma omp simd aligned(t_2075, t_2076, t_2077, t_2078, pa_z, pb_y, pb_z, lk_1298, lk_1299, \
                         lk_1337, ll_1625, ll_1626, mk_1659, mk_1661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2075[k] = f_14 * lk_1298[k]
                    + pa_z[k] * ll_1625[k];

        t_2076[k] = pa_z[k] * ll_1626[k];

        t_2077[k] = f_13 * lk_1299[k]
                    + pb_z[k] * mk_1659[k];

        t_2078[k] = f_19 * lk_1337[k]
                    + pb_y[k] * mk_1661[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece16(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kl0,
                                             const size_t kl1, const size_t lk, const size_t ll,
                                             const size_t mi0, const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;
    const auto f_23 = 3.0 / alpha;
    const auto f_24 = 3.0 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_27 = 2.5 / alpha;
    const auto f_28 = 2.5 * beta / (alpha * p);

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
    auto *t_2184 = buffer.data(target + 2184);
    auto *t_2185 = buffer.data(target + 2185);
    auto *t_2186 = buffer.data(target + 2186);
    auto *t_2187 = buffer.data(target + 2187);
    auto *t_2188 = buffer.data(target + 2188);
    auto *t_2189 = buffer.data(target + 2189);
    auto *t_2190 = buffer.data(target + 2190);
    auto *t_2191 = buffer.data(target + 2191);
    auto *t_2192 = buffer.data(target + 2192);
    auto *t_2193 = buffer.data(target + 2193);
    auto *t_2194 = buffer.data(target + 2194);
    auto *t_2195 = buffer.data(target + 2195);
    auto *t_2196 = buffer.data(target + 2196);
    auto *t_2197 = buffer.data(target + 2197);
    auto *t_2198 = buffer.data(target + 2198);
    auto *t_2199 = buffer.data(target + 2199);
    auto *t_2200 = buffer.data(target + 2200);
    auto *t_2201 = buffer.data(target + 2201);
    auto *t_2202 = buffer.data(target + 2202);
    auto *t_2203 = buffer.data(target + 2203);
    auto *t_2204 = buffer.data(target + 2204);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_1296 = buffer.data(kl0 + 1296);
    const auto *kl0_1341 = buffer.data(kl0 + 1341);
    const auto *kl0_1394 = buffer.data(kl0 + 1394);
    const auto *kl0_1439 = buffer.data(kl0 + 1439);

    const auto *kl1_1296 = buffer.data(kl1 + 1296);
    const auto *kl1_1341 = buffer.data(kl1 + 1341);
    const auto *kl1_1394 = buffer.data(kl1 + 1394);
    const auto *kl1_1439 = buffer.data(kl1 + 1439);

    const auto *lk_1301 = buffer.data(lk + 1301);
    const auto *lk_1302 = buffer.data(lk + 1302);
    const auto *lk_1303 = buffer.data(lk + 1303);
    const auto *lk_1305 = buffer.data(lk + 1305);
    const auto *lk_1306 = buffer.data(lk + 1306);
    const auto *lk_1307 = buffer.data(lk + 1307);
    const auto *lk_1308 = buffer.data(lk + 1308);
    const auto *lk_1310 = buffer.data(lk + 1310);
    const auto *lk_1311 = buffer.data(lk + 1311);
    const auto *lk_1312 = buffer.data(lk + 1312);
    const auto *lk_1313 = buffer.data(lk + 1313);
    const auto *lk_1314 = buffer.data(lk + 1314);
    const auto *lk_1316 = buffer.data(lk + 1316);
    const auto *lk_1324 = buffer.data(lk + 1324);
    const auto *lk_1325 = buffer.data(lk + 1325);
    const auto *lk_1326 = buffer.data(lk + 1326);
    const auto *lk_1327 = buffer.data(lk + 1327);
    const auto *lk_1328 = buffer.data(lk + 1328);
    const auto *lk_1329 = buffer.data(lk + 1329);
    const auto *lk_1331 = buffer.data(lk + 1331);
    const auto *lk_1332 = buffer.data(lk + 1332);
    const auto *lk_1335 = buffer.data(lk + 1335);
    const auto *lk_1338 = buffer.data(lk + 1338);
    const auto *lk_1341 = buffer.data(lk + 1341);
    const auto *lk_1342 = buffer.data(lk + 1342);
    const auto *lk_1346 = buffer.data(lk + 1346);
    const auto *lk_1347 = buffer.data(lk + 1347);
    const auto *lk_1352 = buffer.data(lk + 1352);
    const auto *lk_1360 = buffer.data(lk + 1360);
    const auto *lk_1367 = buffer.data(lk + 1367);
    const auto *lk_1368 = buffer.data(lk + 1368);
    const auto *lk_1370 = buffer.data(lk + 1370);
    const auto *lk_1371 = buffer.data(lk + 1371);
    const auto *lk_1373 = buffer.data(lk + 1373);
    const auto *lk_1374 = buffer.data(lk + 1374);
    const auto *lk_1377 = buffer.data(lk + 1377);
    const auto *lk_1378 = buffer.data(lk + 1378);
    const auto *lk_1382 = buffer.data(lk + 1382);
    const auto *lk_1383 = buffer.data(lk + 1383);
    const auto *lk_1388 = buffer.data(lk + 1388);
    const auto *lk_1396 = buffer.data(lk + 1396);
    const auto *lk_1398 = buffer.data(lk + 1398);
    const auto *lk_1399 = buffer.data(lk + 1399);
    const auto *lk_1400 = buffer.data(lk + 1400);
    const auto *lk_1401 = buffer.data(lk + 1401);
    const auto *lk_1402 = buffer.data(lk + 1402);
    const auto *lk_1403 = buffer.data(lk + 1403);
    const auto *lk_1404 = buffer.data(lk + 1404);
    const auto *lk_1406 = buffer.data(lk + 1406);
    const auto *lk_1409 = buffer.data(lk + 1409);
    const auto *lk_1413 = buffer.data(lk + 1413);
    const auto *lk_1418 = buffer.data(lk + 1418);
    const auto *lk_1424 = buffer.data(lk + 1424);
    const auto *lk_1434 = buffer.data(lk + 1434);
    const auto *lk_1435 = buffer.data(lk + 1435);
    const auto *lk_1436 = buffer.data(lk + 1436);
    const auto *lk_1437 = buffer.data(lk + 1437);
    const auto *lk_1438 = buffer.data(lk + 1438);
    const auto *lk_1439 = buffer.data(lk + 1439);

    const auto *ll_1629 = buffer.data(ll + 1629);
    const auto *ll_1630 = buffer.data(ll + 1630);
    const auto *ll_1632 = buffer.data(ll + 1632);
    const auto *ll_1634 = buffer.data(ll + 1634);
    const auto *ll_1635 = buffer.data(ll + 1635);
    const auto *ll_1637 = buffer.data(ll + 1637);
    const auto *ll_1638 = buffer.data(ll + 1638);
    const auto *ll_1640 = buffer.data(ll + 1640);
    const auto *ll_1641 = buffer.data(ll + 1641);
    const auto *ll_1643 = buffer.data(ll + 1643);
    const auto *ll_1644 = buffer.data(ll + 1644);
    const auto *ll_1645 = buffer.data(ll + 1645);
    const auto *ll_1647 = buffer.data(ll + 1647);
    const auto *ll_1656 = buffer.data(ll + 1656);
    const auto *ll_1658 = buffer.data(ll + 1658);
    const auto *ll_1659 = buffer.data(ll + 1659);
    const auto *ll_1660 = buffer.data(ll + 1660);
    const auto *ll_1661 = buffer.data(ll + 1661);
    const auto *ll_1662 = buffer.data(ll + 1662);
    const auto *ll_1664 = buffer.data(ll + 1664);
    const auto *ll_1701 = buffer.data(ll + 1701);
    const auto *ll_1746 = buffer.data(ll + 1746);
    const auto *ll_1754 = buffer.data(ll + 1754);
    const auto *ll_1799 = buffer.data(ll + 1799);

    const auto *mi0_1316 = buffer.data(mi0 + 1316);
    const auto *mi0_1319 = buffer.data(mi0 + 1319);
    const auto *mi0_1321 = buffer.data(mi0 + 1321);
    const auto *mi0_1322 = buffer.data(mi0 + 1322);
    const auto *mi0_1325 = buffer.data(mi0 + 1325);
    const auto *mi0_1326 = buffer.data(mi0 + 1326);
    const auto *mi0_1328 = buffer.data(mi0 + 1328);
    const auto *mi0_1330 = buffer.data(mi0 + 1330);
    const auto *mi0_1331 = buffer.data(mi0 + 1331);
    const auto *mi0_1333 = buffer.data(mi0 + 1333);
    const auto *mi0_1334 = buffer.data(mi0 + 1334);
    const auto *mi0_1336 = buffer.data(mi0 + 1336);
    const auto *mi0_1337 = buffer.data(mi0 + 1337);
    const auto *mi0_1339 = buffer.data(mi0 + 1339);
    const auto *mi0_1340 = buffer.data(mi0 + 1340);
    const auto *mi0_1341 = buffer.data(mi0 + 1341);
    const auto *mi0_1342 = buffer.data(mi0 + 1342);
    const auto *mi0_1343 = buffer.data(mi0 + 1343);
    const auto *mi0_1344 = buffer.data(mi0 + 1344);
    const auto *mi0_1347 = buffer.data(mi0 + 1347);
    const auto *mi0_1349 = buffer.data(mi0 + 1349);
    const auto *mi0_1350 = buffer.data(mi0 + 1350);
    const auto *mi0_1353 = buffer.data(mi0 + 1353);
    const auto *mi0_1354 = buffer.data(mi0 + 1354);
    const auto *mi0_1356 = buffer.data(mi0 + 1356);
    const auto *mi0_1358 = buffer.data(mi0 + 1358);
    const auto *mi0_1359 = buffer.data(mi0 + 1359);
    const auto *mi0_1361 = buffer.data(mi0 + 1361);
    const auto *mi0_1362 = buffer.data(mi0 + 1362);
    const auto *mi0_1364 = buffer.data(mi0 + 1364);
    const auto *mi0_1365 = buffer.data(mi0 + 1365);
    const auto *mi0_1367 = buffer.data(mi0 + 1367);
    const auto *mi0_1368 = buffer.data(mi0 + 1368);
    const auto *mi0_1369 = buffer.data(mi0 + 1369);
    const auto *mi0_1370 = buffer.data(mi0 + 1370);
    const auto *mi0_1371 = buffer.data(mi0 + 1371);

    const auto *mi1_1316 = buffer.data(mi1 + 1316);
    const auto *mi1_1319 = buffer.data(mi1 + 1319);
    const auto *mi1_1321 = buffer.data(mi1 + 1321);
    const auto *mi1_1322 = buffer.data(mi1 + 1322);
    const auto *mi1_1325 = buffer.data(mi1 + 1325);
    const auto *mi1_1326 = buffer.data(mi1 + 1326);
    const auto *mi1_1328 = buffer.data(mi1 + 1328);
    const auto *mi1_1330 = buffer.data(mi1 + 1330);
    const auto *mi1_1331 = buffer.data(mi1 + 1331);
    const auto *mi1_1333 = buffer.data(mi1 + 1333);
    const auto *mi1_1334 = buffer.data(mi1 + 1334);
    const auto *mi1_1336 = buffer.data(mi1 + 1336);
    const auto *mi1_1337 = buffer.data(mi1 + 1337);
    const auto *mi1_1339 = buffer.data(mi1 + 1339);
    const auto *mi1_1340 = buffer.data(mi1 + 1340);
    const auto *mi1_1341 = buffer.data(mi1 + 1341);
    const auto *mi1_1342 = buffer.data(mi1 + 1342);
    const auto *mi1_1343 = buffer.data(mi1 + 1343);
    const auto *mi1_1344 = buffer.data(mi1 + 1344);
    const auto *mi1_1347 = buffer.data(mi1 + 1347);
    const auto *mi1_1349 = buffer.data(mi1 + 1349);
    const auto *mi1_1350 = buffer.data(mi1 + 1350);
    const auto *mi1_1353 = buffer.data(mi1 + 1353);
    const auto *mi1_1354 = buffer.data(mi1 + 1354);
    const auto *mi1_1356 = buffer.data(mi1 + 1356);
    const auto *mi1_1358 = buffer.data(mi1 + 1358);
    const auto *mi1_1359 = buffer.data(mi1 + 1359);
    const auto *mi1_1361 = buffer.data(mi1 + 1361);
    const auto *mi1_1362 = buffer.data(mi1 + 1362);
    const auto *mi1_1364 = buffer.data(mi1 + 1364);
    const auto *mi1_1365 = buffer.data(mi1 + 1365);
    const auto *mi1_1367 = buffer.data(mi1 + 1367);
    const auto *mi1_1368 = buffer.data(mi1 + 1368);
    const auto *mi1_1369 = buffer.data(mi1 + 1369);
    const auto *mi1_1370 = buffer.data(mi1 + 1370);
    const auto *mi1_1371 = buffer.data(mi1 + 1371);

    const auto *mk_1662 = buffer.data(mk + 1662);
    const auto *mk_1665 = buffer.data(mk + 1665);
    const auto *mk_1666 = buffer.data(mk + 1666);
    const auto *mk_1670 = buffer.data(mk + 1670);
    const auto *mk_1671 = buffer.data(mk + 1671);
    const auto *mk_1676 = buffer.data(mk + 1676);
    const auto *mk_1684 = buffer.data(mk + 1684);
    const auto *mk_1685 = buffer.data(mk + 1685);
    const auto *mk_1686 = buffer.data(mk + 1686);
    const auto *mk_1687 = buffer.data(mk + 1687);
    const auto *mk_1688 = buffer.data(mk + 1688);
    const auto *mk_1689 = buffer.data(mk + 1689);
    const auto *mk_1690 = buffer.data(mk + 1690);
    const auto *mk_1691 = buffer.data(mk + 1691);
    const auto *mk_1692 = buffer.data(mk + 1692);
    const auto *mk_1694 = buffer.data(mk + 1694);
    const auto *mk_1695 = buffer.data(mk + 1695);
    const auto *mk_1697 = buffer.data(mk + 1697);
    const auto *mk_1698 = buffer.data(mk + 1698);
    const auto *mk_1701 = buffer.data(mk + 1701);
    const auto *mk_1702 = buffer.data(mk + 1702);
    const auto *mk_1704 = buffer.data(mk + 1704);
    const auto *mk_1706 = buffer.data(mk + 1706);
    const auto *mk_1707 = buffer.data(mk + 1707);
    const auto *mk_1709 = buffer.data(mk + 1709);
    const auto *mk_1710 = buffer.data(mk + 1710);
    const auto *mk_1712 = buffer.data(mk + 1712);
    const auto *mk_1713 = buffer.data(mk + 1713);
    const auto *mk_1715 = buffer.data(mk + 1715);
    const auto *mk_1716 = buffer.data(mk + 1716);
    const auto *mk_1717 = buffer.data(mk + 1717);
    const auto *mk_1719 = buffer.data(mk + 1719);
    const auto *mk_1720 = buffer.data(mk + 1720);
    const auto *mk_1721 = buffer.data(mk + 1721);
    const auto *mk_1722 = buffer.data(mk + 1722);
    const auto *mk_1723 = buffer.data(mk + 1723);
    const auto *mk_1724 = buffer.data(mk + 1724);
    const auto *mk_1725 = buffer.data(mk + 1725);
    const auto *mk_1726 = buffer.data(mk + 1726);
    const auto *mk_1727 = buffer.data(mk + 1727);
    const auto *mk_1728 = buffer.data(mk + 1728);
    const auto *mk_1730 = buffer.data(mk + 1730);
    const auto *mk_1731 = buffer.data(mk + 1731);
    const auto *mk_1733 = buffer.data(mk + 1733);
    const auto *mk_1734 = buffer.data(mk + 1734);
    const auto *mk_1737 = buffer.data(mk + 1737);
    const auto *mk_1738 = buffer.data(mk + 1738);
    const auto *mk_1740 = buffer.data(mk + 1740);
    const auto *mk_1742 = buffer.data(mk + 1742);
    const auto *mk_1743 = buffer.data(mk + 1743);
    const auto *mk_1745 = buffer.data(mk + 1745);
    const auto *mk_1746 = buffer.data(mk + 1746);
    const auto *mk_1748 = buffer.data(mk + 1748);
    const auto *mk_1749 = buffer.data(mk + 1749);
    const auto *mk_1751 = buffer.data(mk + 1751);
    const auto *mk_1752 = buffer.data(mk + 1752);
    const auto *mk_1753 = buffer.data(mk + 1753);
    const auto *mk_1755 = buffer.data(mk + 1755);
    const auto *mk_1756 = buffer.data(mk + 1756);
    const auto *mk_1757 = buffer.data(mk + 1757);
    const auto *mk_1758 = buffer.data(mk + 1758);
    const auto *mk_1759 = buffer.data(mk + 1759);
    const auto *mk_1760 = buffer.data(mk + 1760);
    const auto *mk_1761 = buffer.data(mk + 1761);
    const auto *mk_1762 = buffer.data(mk + 1762);
    const auto *mk_1763 = buffer.data(mk + 1763);

#pragma omp simd aligned(t_2079, t_2080, t_2081, t_2082, pa_z, pb_z, lk_1301, lk_1302, \
                         lk_1303, ll_1629, ll_1630, ll_1632, mk_1662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2079[k] = f_15 * lk_1301[k]
                    + pa_z[k] * ll_1629[k];

        t_2080[k] = pa_z[k] * ll_1630[k];

        t_2081[k] = f_13 * lk_1302[k]
                    + pb_z[k] * mk_1662[k];

        t_2082[k] = f_14 * lk_1303[k]
                    + pa_z[k] * ll_1632[k];
    }

#pragma omp simd aligned(t_2083, t_2084, t_2085, t_2086, pa_z, pb_y, pb_z, lk_1305, lk_1306, \
                         lk_1341, ll_1634, ll_1635, mk_1665, mk_1666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2083[k] = f_19 * lk_1341[k]
                    + pb_y[k] * mk_1665[k];

        t_2084[k] = f_16 * lk_1305[k]
                    + pa_z[k] * ll_1634[k];

        t_2085[k] = pa_z[k] * ll_1635[k];

        t_2086[k] = f_13 * lk_1306[k]
                    + pb_z[k] * mk_1666[k];
    }

#pragma omp simd aligned(t_2087, t_2088, t_2089, t_2090, t_2091, pa_z, pb_y, lk_1307, lk_1308, \
                         lk_1310, lk_1346, ll_1637, ll_1638, ll_1640, ll_1641, \
                         mk_1670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2087[k] = f_14 * lk_1307[k]
                    + pa_z[k] * ll_1637[k];

        t_2088[k] = f_15 * lk_1308[k]
                    + pa_z[k] * ll_1638[k];

        t_2089[k] = f_19 * lk_1346[k]
                    + pb_y[k] * mk_1670[k];

        t_2090[k] = f_17 * lk_1310[k]
                    + pa_z[k] * ll_1640[k];

        t_2091[k] = pa_z[k] * ll_1641[k];
    }

#pragma omp simd aligned(t_2092, t_2093, t_2094, t_2095, pa_z, pb_z, lk_1311, lk_1312, \
                         lk_1313, lk_1314, ll_1643, ll_1644, ll_1645, \
                         mk_1671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2092[k] = f_13 * lk_1311[k]
                    + pb_z[k] * mk_1671[k];

        t_2093[k] = f_14 * lk_1312[k]
                    + pa_z[k] * ll_1643[k];

        t_2094[k] = f_15 * lk_1313[k]
                    + pa_z[k] * ll_1644[k];

        t_2095[k] = f_16 * lk_1314[k]
                    + pa_z[k] * ll_1645[k];
    }

#pragma omp simd aligned(t_2096, t_2097, t_2098, t_2099, t_2100, pa_z, pb_x, pb_y, lk_1316, \
                         lk_1352, ll_1647, mk_1676, mk_1684, mk_1685, \
                         mk_1686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2096[k] = f_19 * lk_1352[k]
                    + pb_y[k] * mk_1676[k];

        t_2097[k] = f_18 * lk_1316[k]
                    + pa_z[k] * ll_1647[k];

        t_2098[k] = pb_x[k] * mk_1684[k];

        t_2099[k] = pb_x[k] * mk_1685[k];

        t_2100[k] = pb_x[k] * mk_1686[k];
    }

#pragma omp simd aligned(t_2101, t_2102, t_2103, t_2104, t_2105, t_2106, pa_z, pb_x, ll_1656, \
                         mk_1687, mk_1688, mk_1689, mk_1690, mk_1691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2101[k] = pb_x[k] * mk_1687[k];

        t_2102[k] = pb_x[k] * mk_1688[k];

        t_2103[k] = pb_x[k] * mk_1689[k];

        t_2104[k] = pb_x[k] * mk_1690[k];

        t_2105[k] = pb_x[k] * mk_1691[k];

        t_2106[k] = pa_z[k] * ll_1656[k];
    }

#pragma omp simd aligned(t_2107, t_2108, t_2109, t_2110, pa_z, pb_z, lk_1324, lk_1325, \
                         lk_1326, lk_1327, ll_1658, ll_1659, ll_1660, \
                         mk_1684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2107[k] = f_13 * lk_1324[k]
                    + pb_z[k] * mk_1684[k];

        t_2108[k] = f_14 * lk_1325[k]
                    + pa_z[k] * ll_1658[k];

        t_2109[k] = f_15 * lk_1326[k]
                    + pa_z[k] * ll_1659[k];

        t_2110[k] = f_16 * lk_1327[k]
                    + pa_z[k] * ll_1660[k];
    }

#pragma omp simd aligned(t_2111, t_2112, t_2113, t_2114, pa_z, pb_y, lk_1328, lk_1329, \
                         lk_1331, lk_1367, ll_1661, ll_1662, ll_1664, \
                         mk_1691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2111[k] = f_17 * lk_1328[k]
                    + pa_z[k] * ll_1661[k];

        t_2112[k] = f_18 * lk_1329[k]
                    + pa_z[k] * ll_1662[k];

        t_2113[k] = f_19 * lk_1367[k]
                    + pb_y[k] * mk_1691[k];

        t_2114[k] = f_19 * lk_1331[k]
                    + pa_z[k] * ll_1664[k];
    }

#pragma omp simd aligned(t_2115, t_2116, t_2117, t_2118, pb_x, pb_y, pb_z, lk_1332, lk_1368, \
                         mi0_1316, mi0_1319, mi1_1316, mi1_1319, mk_1692, \
                         mk_1695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2115[k] = f_1 * mi0_1316[k]
                    - f_2 * mi1_1316[k]
                    + pb_x[k] * mk_1692[k];

        t_2116[k] = f_22 * lk_1368[k]
                    + pb_y[k] * mk_1692[k];

        t_2117[k] = f_14 * lk_1332[k]
                    + pb_z[k] * mk_1692[k];

        t_2118[k] = f_11 * mi0_1319[k]
                    - f_12 * mi1_1319[k]
                    + pb_x[k] * mk_1695[k];
    }

#pragma omp simd aligned(t_2119, t_2120, t_2121, pb_x, pb_y, lk_1370, mi0_1321, mi0_1322, \
                         mi1_1321, mi1_1322, mk_1694, mk_1697, \
                         mk_1698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2119[k] = f_22 * lk_1370[k]
                    + pb_y[k] * mk_1694[k];

        t_2120[k] = f_11 * mi0_1321[k]
                    - f_12 * mi1_1321[k]
                    + pb_x[k] * mk_1697[k];

        t_2121[k] = f_9 * mi0_1322[k]
                    - f_10 * mi1_1322[k]
                    + pb_x[k] * mk_1698[k];
    }

#pragma omp simd aligned(t_2122, t_2123, t_2124, pb_x, pb_y, pb_z, lk_1335, lk_1373, mi0_1325, \
                         mi1_1325, mk_1695, mk_1697, mk_1701 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2122[k] = f_14 * lk_1335[k]
                    + pb_z[k] * mk_1695[k];

        t_2123[k] = f_22 * lk_1373[k]
                    + pb_y[k] * mk_1697[k];

        t_2124[k] = f_9 * mi0_1325[k]
                    - f_10 * mi1_1325[k]
                    + pb_x[k] * mk_1701[k];
    }

#pragma omp simd aligned(t_2125, t_2126, t_2127, pb_x, pb_z, lk_1338, mi0_1326, mi0_1328, \
                         mi1_1326, mi1_1328, mk_1698, mk_1702, \
                         mk_1704 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2125[k] = f_7 * mi0_1326[k]
                    - f_8 * mi1_1326[k]
                    + pb_x[k] * mk_1702[k];

        t_2126[k] = f_14 * lk_1338[k]
                    + pb_z[k] * mk_1698[k];

        t_2127[k] = f_7 * mi0_1328[k]
                    - f_8 * mi1_1328[k]
                    + pb_x[k] * mk_1704[k];
    }

#pragma omp simd aligned(t_2128, t_2129, t_2130, pb_x, pb_y, lk_1377, mi0_1330, mi0_1331, \
                         mi1_1330, mi1_1331, mk_1701, mk_1706, \
                         mk_1707 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2128[k] = f_22 * lk_1377[k]
                    + pb_y[k] * mk_1701[k];

        t_2129[k] = f_7 * mi0_1330[k]
                    - f_8 * mi1_1330[k]
                    + pb_x[k] * mk_1706[k];

        t_2130[k] = f_5 * mi0_1331[k]
                    - f_6 * mi1_1331[k]
                    + pb_x[k] * mk_1707[k];
    }

#pragma omp simd aligned(t_2131, t_2132, t_2133, pb_x, pb_z, lk_1342, mi0_1333, mi0_1334, \
                         mi1_1333, mi1_1334, mk_1702, mk_1709, \
                         mk_1710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2131[k] = f_14 * lk_1342[k]
                    + pb_z[k] * mk_1702[k];

        t_2132[k] = f_5 * mi0_1333[k]
                    - f_6 * mi1_1333[k]
                    + pb_x[k] * mk_1709[k];

        t_2133[k] = f_5 * mi0_1334[k]
                    - f_6 * mi1_1334[k]
                    + pb_x[k] * mk_1710[k];
    }

#pragma omp simd aligned(t_2134, t_2135, t_2136, pb_x, pb_y, lk_1382, mi0_1336, mi0_1337, \
                         mi1_1336, mi1_1337, mk_1706, mk_1712, \
                         mk_1713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2134[k] = f_22 * lk_1382[k]
                    + pb_y[k] * mk_1706[k];

        t_2135[k] = f_5 * mi0_1336[k]
                    - f_6 * mi1_1336[k]
                    + pb_x[k] * mk_1712[k];

        t_2136[k] = f_3 * mi0_1337[k]
                    - f_4 * mi1_1337[k]
                    + pb_x[k] * mk_1713[k];
    }

#pragma omp simd aligned(t_2137, t_2138, t_2139, pb_x, pb_z, lk_1347, mi0_1339, mi0_1340, \
                         mi1_1339, mi1_1340, mk_1707, mk_1715, \
                         mk_1716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2137[k] = f_14 * lk_1347[k]
                    + pb_z[k] * mk_1707[k];

        t_2138[k] = f_3 * mi0_1339[k]
                    - f_4 * mi1_1339[k]
                    + pb_x[k] * mk_1715[k];

        t_2139[k] = f_3 * mi0_1340[k]
                    - f_4 * mi1_1340[k]
                    + pb_x[k] * mk_1716[k];
    }

#pragma omp simd aligned(t_2140, t_2141, t_2142, t_2143, pb_x, pb_y, lk_1388, mi0_1341, \
                         mi0_1343, mi1_1341, mi1_1343, mk_1712, mk_1717, mk_1719, \
                         mk_1720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2140[k] = f_3 * mi0_1341[k]
                    - f_4 * mi1_1341[k]
                    + pb_x[k] * mk_1717[k];

        t_2141[k] = f_22 * lk_1388[k]
                    + pb_y[k] * mk_1712[k];

        t_2142[k] = f_3 * mi0_1343[k]
                    - f_4 * mi1_1343[k]
                    + pb_x[k] * mk_1719[k];

        t_2143[k] = pb_x[k] * mk_1720[k];
    }

#pragma omp simd aligned(t_2144, t_2145, t_2146, t_2147, t_2148, t_2149, t_2150, pb_x, \
                         mk_1721, mk_1722, mk_1723, mk_1724, mk_1725, mk_1726, \
                         mk_1727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2144[k] = pb_x[k] * mk_1721[k];

        t_2145[k] = pb_x[k] * mk_1722[k];

        t_2146[k] = pb_x[k] * mk_1723[k];

        t_2147[k] = pb_x[k] * mk_1724[k];

        t_2148[k] = pb_x[k] * mk_1725[k];

        t_2149[k] = pb_x[k] * mk_1726[k];

        t_2150[k] = pb_x[k] * mk_1727[k];
    }

#pragma omp simd aligned(t_2151, t_2152, t_2153, pa_z, pb_y, pb_z, kl0_1296, kl1_1296, \
                         lk_1360, lk_1398, ll_1701, mi0_1339, mi1_1339, mk_1720, \
                         mk_1722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2151[k] = f_20 * kl0_1296[k]
                    - f_21 * kl1_1296[k]
                    + pa_z[k] * ll_1701[k];

        t_2152[k] = f_14 * lk_1360[k]
                    + pb_z[k] * mk_1720[k];

        t_2153[k] = f_22 * lk_1398[k]
                    + f_11 * mi0_1339[k]
                    - f_12 * mi1_1339[k]
                    + pb_y[k] * mk_1722[k];
    }

#pragma omp simd aligned(t_2154, t_2155, t_2156, pb_y, lk_1399, lk_1400, lk_1401, mi0_1340, \
                         mi0_1341, mi0_1342, mi1_1340, mi1_1341, mi1_1342, mk_1723, mk_1724, \
                         mk_1725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2154[k] = f_22 * lk_1399[k]
                    + f_9 * mi0_1340[k]
                    - f_10 * mi1_1340[k]
                    + pb_y[k] * mk_1723[k];

        t_2155[k] = f_22 * lk_1400[k]
                    + f_7 * mi0_1341[k]
                    - f_8 * mi1_1341[k]
                    + pb_y[k] * mk_1724[k];

        t_2156[k] = f_22 * lk_1401[k]
                    + f_5 * mi0_1342[k]
                    - f_6 * mi1_1342[k]
                    + pb_y[k] * mk_1725[k];
    }

#pragma omp simd aligned(t_2157, t_2158, t_2159, pa_y, pb_y, kl0_1394, kl1_1394, lk_1402, \
                         lk_1403, ll_1754, mi0_1343, mi1_1343, mk_1726, \
                         mk_1727 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2157[k] = f_22 * lk_1402[k]
                    + f_3 * mi0_1343[k]
                    - f_4 * mi1_1343[k]
                    + pb_y[k] * mk_1726[k];

        t_2158[k] = f_22 * lk_1403[k]
                    + pb_y[k] * mk_1727[k];

        t_2159[k] = f_23 * kl0_1394[k]
                    - f_24 * kl1_1394[k]
                    + pa_y[k] * ll_1754[k];
    }

#pragma omp simd aligned(t_2160, t_2161, t_2162, t_2163, pb_x, pb_y, pb_z, lk_1368, lk_1404, \
                         mi0_1344, mi0_1347, mi1_1344, mi1_1347, mk_1728, \
                         mk_1731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2160[k] = f_1 * mi0_1344[k]
                    - f_2 * mi1_1344[k]
                    + pb_x[k] * mk_1728[k];

        t_2161[k] = f_18 * lk_1404[k]
                    + pb_y[k] * mk_1728[k];

        t_2162[k] = f_15 * lk_1368[k]
                    + pb_z[k] * mk_1728[k];

        t_2163[k] = f_11 * mi0_1347[k]
                    - f_12 * mi1_1347[k]
                    + pb_x[k] * mk_1731[k];
    }

#pragma omp simd aligned(t_2164, t_2165, t_2166, pb_x, pb_y, lk_1406, mi0_1349, mi0_1350, \
                         mi1_1349, mi1_1350, mk_1730, mk_1733, \
                         mk_1734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2164[k] = f_18 * lk_1406[k]
                    + pb_y[k] * mk_1730[k];

        t_2165[k] = f_11 * mi0_1349[k]
                    - f_12 * mi1_1349[k]
                    + pb_x[k] * mk_1733[k];

        t_2166[k] = f_9 * mi0_1350[k]
                    - f_10 * mi1_1350[k]
                    + pb_x[k] * mk_1734[k];
    }

#pragma omp simd aligned(t_2167, t_2168, t_2169, pb_x, pb_y, pb_z, lk_1371, lk_1409, mi0_1353, \
                         mi1_1353, mk_1731, mk_1733, mk_1737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2167[k] = f_15 * lk_1371[k]
                    + pb_z[k] * mk_1731[k];

        t_2168[k] = f_18 * lk_1409[k]
                    + pb_y[k] * mk_1733[k];

        t_2169[k] = f_9 * mi0_1353[k]
                    - f_10 * mi1_1353[k]
                    + pb_x[k] * mk_1737[k];
    }

#pragma omp simd aligned(t_2170, t_2171, t_2172, pb_x, pb_z, lk_1374, mi0_1354, mi0_1356, \
                         mi1_1354, mi1_1356, mk_1734, mk_1738, \
                         mk_1740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2170[k] = f_7 * mi0_1354[k]
                    - f_8 * mi1_1354[k]
                    + pb_x[k] * mk_1738[k];

        t_2171[k] = f_15 * lk_1374[k]
                    + pb_z[k] * mk_1734[k];

        t_2172[k] = f_7 * mi0_1356[k]
                    - f_8 * mi1_1356[k]
                    + pb_x[k] * mk_1740[k];
    }

#pragma omp simd aligned(t_2173, t_2174, t_2175, pb_x, pb_y, lk_1413, mi0_1358, mi0_1359, \
                         mi1_1358, mi1_1359, mk_1737, mk_1742, \
                         mk_1743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2173[k] = f_18 * lk_1413[k]
                    + pb_y[k] * mk_1737[k];

        t_2174[k] = f_7 * mi0_1358[k]
                    - f_8 * mi1_1358[k]
                    + pb_x[k] * mk_1742[k];

        t_2175[k] = f_5 * mi0_1359[k]
                    - f_6 * mi1_1359[k]
                    + pb_x[k] * mk_1743[k];
    }

#pragma omp simd aligned(t_2176, t_2177, t_2178, pb_x, pb_z, lk_1378, mi0_1361, mi0_1362, \
                         mi1_1361, mi1_1362, mk_1738, mk_1745, \
                         mk_1746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2176[k] = f_15 * lk_1378[k]
                    + pb_z[k] * mk_1738[k];

        t_2177[k] = f_5 * mi0_1361[k]
                    - f_6 * mi1_1361[k]
                    + pb_x[k] * mk_1745[k];

        t_2178[k] = f_5 * mi0_1362[k]
                    - f_6 * mi1_1362[k]
                    + pb_x[k] * mk_1746[k];
    }

#pragma omp simd aligned(t_2179, t_2180, t_2181, pb_x, pb_y, lk_1418, mi0_1364, mi0_1365, \
                         mi1_1364, mi1_1365, mk_1742, mk_1748, \
                         mk_1749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2179[k] = f_18 * lk_1418[k]
                    + pb_y[k] * mk_1742[k];

        t_2180[k] = f_5 * mi0_1364[k]
                    - f_6 * mi1_1364[k]
                    + pb_x[k] * mk_1748[k];

        t_2181[k] = f_3 * mi0_1365[k]
                    - f_4 * mi1_1365[k]
                    + pb_x[k] * mk_1749[k];
    }

#pragma omp simd aligned(t_2182, t_2183, t_2184, pb_x, pb_z, lk_1383, mi0_1367, mi0_1368, \
                         mi1_1367, mi1_1368, mk_1743, mk_1751, \
                         mk_1752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2182[k] = f_15 * lk_1383[k]
                    + pb_z[k] * mk_1743[k];

        t_2183[k] = f_3 * mi0_1367[k]
                    - f_4 * mi1_1367[k]
                    + pb_x[k] * mk_1751[k];

        t_2184[k] = f_3 * mi0_1368[k]
                    - f_4 * mi1_1368[k]
                    + pb_x[k] * mk_1752[k];
    }

#pragma omp simd aligned(t_2185, t_2186, t_2187, t_2188, pb_x, pb_y, lk_1424, mi0_1369, \
                         mi0_1371, mi1_1369, mi1_1371, mk_1748, mk_1753, mk_1755, \
                         mk_1756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2185[k] = f_3 * mi0_1369[k]
                    - f_4 * mi1_1369[k]
                    + pb_x[k] * mk_1753[k];

        t_2186[k] = f_18 * lk_1424[k]
                    + pb_y[k] * mk_1748[k];

        t_2187[k] = f_3 * mi0_1371[k]
                    - f_4 * mi1_1371[k]
                    + pb_x[k] * mk_1755[k];

        t_2188[k] = pb_x[k] * mk_1756[k];
    }

#pragma omp simd aligned(t_2189, t_2190, t_2191, t_2192, t_2193, t_2194, t_2195, pb_x, \
                         mk_1757, mk_1758, mk_1759, mk_1760, mk_1761, mk_1762, \
                         mk_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2189[k] = pb_x[k] * mk_1757[k];

        t_2190[k] = pb_x[k] * mk_1758[k];

        t_2191[k] = pb_x[k] * mk_1759[k];

        t_2192[k] = pb_x[k] * mk_1760[k];

        t_2193[k] = pb_x[k] * mk_1761[k];

        t_2194[k] = pb_x[k] * mk_1762[k];

        t_2195[k] = pb_x[k] * mk_1763[k];
    }

#pragma omp simd aligned(t_2196, t_2197, t_2198, pa_z, pb_y, pb_z, kl0_1341, kl1_1341, \
                         lk_1396, lk_1434, ll_1746, mi0_1367, mi1_1367, mk_1756, \
                         mk_1758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2196[k] = f_25 * kl0_1341[k]
                    - f_26 * kl1_1341[k]
                    + pa_z[k] * ll_1746[k];

        t_2197[k] = f_15 * lk_1396[k]
                    + pb_z[k] * mk_1756[k];

        t_2198[k] = f_18 * lk_1434[k]
                    + f_11 * mi0_1367[k]
                    - f_12 * mi1_1367[k]
                    + pb_y[k] * mk_1758[k];
    }

#pragma omp simd aligned(t_2199, t_2200, t_2201, pb_y, lk_1435, lk_1436, lk_1437, mi0_1368, \
                         mi0_1369, mi0_1370, mi1_1368, mi1_1369, mi1_1370, mk_1759, mk_1760, \
                         mk_1761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2199[k] = f_18 * lk_1435[k]
                    + f_9 * mi0_1368[k]
                    - f_10 * mi1_1368[k]
                    + pb_y[k] * mk_1759[k];

        t_2200[k] = f_18 * lk_1436[k]
                    + f_7 * mi0_1369[k]
                    - f_8 * mi1_1369[k]
                    + pb_y[k] * mk_1760[k];

        t_2201[k] = f_18 * lk_1437[k]
                    + f_5 * mi0_1370[k]
                    - f_6 * mi1_1370[k]
                    + pb_y[k] * mk_1761[k];
    }

#pragma omp simd aligned(t_2202, t_2203, t_2204, pa_y, pb_y, kl0_1439, kl1_1439, lk_1438, \
                         lk_1439, ll_1799, mi0_1371, mi1_1371, mk_1762, \
                         mk_1763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2202[k] = f_18 * lk_1438[k]
                    + f_3 * mi0_1371[k]
                    - f_4 * mi1_1371[k]
                    + pb_y[k] * mk_1762[k];

        t_2203[k] = f_18 * lk_1439[k]
                    + pb_y[k] * mk_1763[k];

        t_2204[k] = f_27 * kl0_1439[k]
                    - f_28 * kl1_1439[k]
                    + pa_y[k] * ll_1799[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece17(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kl0,
                                             const size_t kl1, const size_t lk, const size_t ll,
                                             const size_t mi0, const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_15 = 1.5 / p;
    const auto f_16 = 2.0 / p;
    const auto f_17 = 2.5 / p;
    const auto f_18 = 3.0 / p;
    const auto f_29 = 1.5 / alpha;
    const auto f_30 = 1.5 * beta / (alpha * p);
    const auto f_31 = 2.0 / alpha;
    const auto f_32 = 2.0 * beta / (alpha * p);

    auto *t_2205 = buffer.data(target + 2205);
    auto *t_2206 = buffer.data(target + 2206);
    auto *t_2207 = buffer.data(target + 2207);
    auto *t_2208 = buffer.data(target + 2208);
    auto *t_2209 = buffer.data(target + 2209);
    auto *t_2210 = buffer.data(target + 2210);
    auto *t_2211 = buffer.data(target + 2211);
    auto *t_2212 = buffer.data(target + 2212);
    auto *t_2213 = buffer.data(target + 2213);
    auto *t_2214 = buffer.data(target + 2214);
    auto *t_2215 = buffer.data(target + 2215);
    auto *t_2216 = buffer.data(target + 2216);
    auto *t_2217 = buffer.data(target + 2217);
    auto *t_2218 = buffer.data(target + 2218);
    auto *t_2219 = buffer.data(target + 2219);
    auto *t_2220 = buffer.data(target + 2220);
    auto *t_2221 = buffer.data(target + 2221);
    auto *t_2222 = buffer.data(target + 2222);
    auto *t_2223 = buffer.data(target + 2223);
    auto *t_2224 = buffer.data(target + 2224);
    auto *t_2225 = buffer.data(target + 2225);
    auto *t_2226 = buffer.data(target + 2226);
    auto *t_2227 = buffer.data(target + 2227);
    auto *t_2228 = buffer.data(target + 2228);
    auto *t_2229 = buffer.data(target + 2229);
    auto *t_2230 = buffer.data(target + 2230);
    auto *t_2231 = buffer.data(target + 2231);
    auto *t_2232 = buffer.data(target + 2232);
    auto *t_2233 = buffer.data(target + 2233);
    auto *t_2234 = buffer.data(target + 2234);
    auto *t_2235 = buffer.data(target + 2235);
    auto *t_2236 = buffer.data(target + 2236);
    auto *t_2237 = buffer.data(target + 2237);
    auto *t_2238 = buffer.data(target + 2238);
    auto *t_2239 = buffer.data(target + 2239);
    auto *t_2240 = buffer.data(target + 2240);
    auto *t_2241 = buffer.data(target + 2241);
    auto *t_2242 = buffer.data(target + 2242);
    auto *t_2243 = buffer.data(target + 2243);
    auto *t_2244 = buffer.data(target + 2244);
    auto *t_2245 = buffer.data(target + 2245);
    auto *t_2246 = buffer.data(target + 2246);
    auto *t_2247 = buffer.data(target + 2247);
    auto *t_2248 = buffer.data(target + 2248);
    auto *t_2249 = buffer.data(target + 2249);
    auto *t_2250 = buffer.data(target + 2250);
    auto *t_2251 = buffer.data(target + 2251);
    auto *t_2252 = buffer.data(target + 2252);
    auto *t_2253 = buffer.data(target + 2253);
    auto *t_2254 = buffer.data(target + 2254);
    auto *t_2255 = buffer.data(target + 2255);
    auto *t_2256 = buffer.data(target + 2256);
    auto *t_2257 = buffer.data(target + 2257);
    auto *t_2258 = buffer.data(target + 2258);
    auto *t_2259 = buffer.data(target + 2259);
    auto *t_2260 = buffer.data(target + 2260);
    auto *t_2261 = buffer.data(target + 2261);
    auto *t_2262 = buffer.data(target + 2262);
    auto *t_2263 = buffer.data(target + 2263);
    auto *t_2264 = buffer.data(target + 2264);
    auto *t_2265 = buffer.data(target + 2265);
    auto *t_2266 = buffer.data(target + 2266);
    auto *t_2267 = buffer.data(target + 2267);
    auto *t_2268 = buffer.data(target + 2268);
    auto *t_2269 = buffer.data(target + 2269);
    auto *t_2270 = buffer.data(target + 2270);
    auto *t_2271 = buffer.data(target + 2271);
    auto *t_2272 = buffer.data(target + 2272);
    auto *t_2273 = buffer.data(target + 2273);
    auto *t_2274 = buffer.data(target + 2274);
    auto *t_2275 = buffer.data(target + 2275);
    auto *t_2276 = buffer.data(target + 2276);
    auto *t_2277 = buffer.data(target + 2277);
    auto *t_2278 = buffer.data(target + 2278);
    auto *t_2279 = buffer.data(target + 2279);
    auto *t_2280 = buffer.data(target + 2280);
    auto *t_2281 = buffer.data(target + 2281);
    auto *t_2282 = buffer.data(target + 2282);
    auto *t_2283 = buffer.data(target + 2283);
    auto *t_2284 = buffer.data(target + 2284);
    auto *t_2285 = buffer.data(target + 2285);
    auto *t_2286 = buffer.data(target + 2286);
    auto *t_2287 = buffer.data(target + 2287);
    auto *t_2288 = buffer.data(target + 2288);
    auto *t_2289 = buffer.data(target + 2289);
    auto *t_2290 = buffer.data(target + 2290);
    auto *t_2291 = buffer.data(target + 2291);
    auto *t_2292 = buffer.data(target + 2292);
    auto *t_2293 = buffer.data(target + 2293);
    auto *t_2294 = buffer.data(target + 2294);
    auto *t_2295 = buffer.data(target + 2295);
    auto *t_2296 = buffer.data(target + 2296);
    auto *t_2297 = buffer.data(target + 2297);
    auto *t_2298 = buffer.data(target + 2298);
    auto *t_2299 = buffer.data(target + 2299);
    auto *t_2300 = buffer.data(target + 2300);
    auto *t_2301 = buffer.data(target + 2301);
    auto *t_2302 = buffer.data(target + 2302);
    auto *t_2303 = buffer.data(target + 2303);
    auto *t_2304 = buffer.data(target + 2304);
    auto *t_2305 = buffer.data(target + 2305);
    auto *t_2306 = buffer.data(target + 2306);
    auto *t_2307 = buffer.data(target + 2307);
    auto *t_2308 = buffer.data(target + 2308);
    auto *t_2309 = buffer.data(target + 2309);
    auto *t_2310 = buffer.data(target + 2310);
    auto *t_2311 = buffer.data(target + 2311);
    auto *t_2312 = buffer.data(target + 2312);
    auto *t_2313 = buffer.data(target + 2313);
    auto *t_2314 = buffer.data(target + 2314);
    auto *t_2315 = buffer.data(target + 2315);
    auto *t_2316 = buffer.data(target + 2316);
    auto *t_2317 = buffer.data(target + 2317);
    auto *t_2318 = buffer.data(target + 2318);
    auto *t_2319 = buffer.data(target + 2319);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_1386 = buffer.data(kl0 + 1386);
    const auto *kl0_1431 = buffer.data(kl0 + 1431);
    const auto *kl0_1484 = buffer.data(kl0 + 1484);
    const auto *kl0_1529 = buffer.data(kl0 + 1529);

    const auto *kl1_1386 = buffer.data(kl1 + 1386);
    const auto *kl1_1431 = buffer.data(kl1 + 1431);
    const auto *kl1_1484 = buffer.data(kl1 + 1484);
    const auto *kl1_1529 = buffer.data(kl1 + 1529);

    const auto *lk_1404 = buffer.data(lk + 1404);
    const auto *lk_1407 = buffer.data(lk + 1407);
    const auto *lk_1410 = buffer.data(lk + 1410);
    const auto *lk_1414 = buffer.data(lk + 1414);
    const auto *lk_1419 = buffer.data(lk + 1419);
    const auto *lk_1432 = buffer.data(lk + 1432);
    const auto *lk_1440 = buffer.data(lk + 1440);
    const auto *lk_1442 = buffer.data(lk + 1442);
    const auto *lk_1443 = buffer.data(lk + 1443);
    const auto *lk_1445 = buffer.data(lk + 1445);
    const auto *lk_1446 = buffer.data(lk + 1446);
    const auto *lk_1449 = buffer.data(lk + 1449);
    const auto *lk_1450 = buffer.data(lk + 1450);
    const auto *lk_1454 = buffer.data(lk + 1454);
    const auto *lk_1455 = buffer.data(lk + 1455);
    const auto *lk_1460 = buffer.data(lk + 1460);
    const auto *lk_1468 = buffer.data(lk + 1468);
    const auto *lk_1470 = buffer.data(lk + 1470);
    const auto *lk_1471 = buffer.data(lk + 1471);
    const auto *lk_1472 = buffer.data(lk + 1472);
    const auto *lk_1473 = buffer.data(lk + 1473);
    const auto *lk_1474 = buffer.data(lk + 1474);
    const auto *lk_1475 = buffer.data(lk + 1475);
    const auto *lk_1476 = buffer.data(lk + 1476);
    const auto *lk_1478 = buffer.data(lk + 1478);
    const auto *lk_1479 = buffer.data(lk + 1479);
    const auto *lk_1481 = buffer.data(lk + 1481);
    const auto *lk_1482 = buffer.data(lk + 1482);
    const auto *lk_1485 = buffer.data(lk + 1485);
    const auto *lk_1486 = buffer.data(lk + 1486);
    const auto *lk_1490 = buffer.data(lk + 1490);
    const auto *lk_1491 = buffer.data(lk + 1491);
    const auto *lk_1496 = buffer.data(lk + 1496);
    const auto *lk_1506 = buffer.data(lk + 1506);
    const auto *lk_1507 = buffer.data(lk + 1507);
    const auto *lk_1508 = buffer.data(lk + 1508);
    const auto *lk_1509 = buffer.data(lk + 1509);
    const auto *lk_1510 = buffer.data(lk + 1510);
    const auto *lk_1511 = buffer.data(lk + 1511);
    const auto *lk_1512 = buffer.data(lk + 1512);
    const auto *lk_1514 = buffer.data(lk + 1514);
    const auto *lk_1517 = buffer.data(lk + 1517);
    const auto *lk_1521 = buffer.data(lk + 1521);
    const auto *lk_1526 = buffer.data(lk + 1526);

    const auto *ll_1791 = buffer.data(ll + 1791);
    const auto *ll_1836 = buffer.data(ll + 1836);
    const auto *ll_1844 = buffer.data(ll + 1844);
    const auto *ll_1889 = buffer.data(ll + 1889);

    const auto *mi0_1372 = buffer.data(mi0 + 1372);
    const auto *mi0_1375 = buffer.data(mi0 + 1375);
    const auto *mi0_1377 = buffer.data(mi0 + 1377);
    const auto *mi0_1378 = buffer.data(mi0 + 1378);
    const auto *mi0_1381 = buffer.data(mi0 + 1381);
    const auto *mi0_1382 = buffer.data(mi0 + 1382);
    const auto *mi0_1384 = buffer.data(mi0 + 1384);
    const auto *mi0_1386 = buffer.data(mi0 + 1386);
    const auto *mi0_1387 = buffer.data(mi0 + 1387);
    const auto *mi0_1389 = buffer.data(mi0 + 1389);
    const auto *mi0_1390 = buffer.data(mi0 + 1390);
    const auto *mi0_1392 = buffer.data(mi0 + 1392);
    const auto *mi0_1393 = buffer.data(mi0 + 1393);
    const auto *mi0_1395 = buffer.data(mi0 + 1395);
    const auto *mi0_1396 = buffer.data(mi0 + 1396);
    const auto *mi0_1397 = buffer.data(mi0 + 1397);
    const auto *mi0_1398 = buffer.data(mi0 + 1398);
    const auto *mi0_1399 = buffer.data(mi0 + 1399);
    const auto *mi0_1400 = buffer.data(mi0 + 1400);
    const auto *mi0_1403 = buffer.data(mi0 + 1403);
    const auto *mi0_1405 = buffer.data(mi0 + 1405);
    const auto *mi0_1406 = buffer.data(mi0 + 1406);
    const auto *mi0_1409 = buffer.data(mi0 + 1409);
    const auto *mi0_1410 = buffer.data(mi0 + 1410);
    const auto *mi0_1412 = buffer.data(mi0 + 1412);
    const auto *mi0_1414 = buffer.data(mi0 + 1414);
    const auto *mi0_1415 = buffer.data(mi0 + 1415);
    const auto *mi0_1417 = buffer.data(mi0 + 1417);
    const auto *mi0_1418 = buffer.data(mi0 + 1418);
    const auto *mi0_1420 = buffer.data(mi0 + 1420);
    const auto *mi0_1421 = buffer.data(mi0 + 1421);
    const auto *mi0_1423 = buffer.data(mi0 + 1423);
    const auto *mi0_1424 = buffer.data(mi0 + 1424);
    const auto *mi0_1425 = buffer.data(mi0 + 1425);
    const auto *mi0_1426 = buffer.data(mi0 + 1426);
    const auto *mi0_1427 = buffer.data(mi0 + 1427);
    const auto *mi0_1428 = buffer.data(mi0 + 1428);
    const auto *mi0_1431 = buffer.data(mi0 + 1431);
    const auto *mi0_1433 = buffer.data(mi0 + 1433);
    const auto *mi0_1434 = buffer.data(mi0 + 1434);
    const auto *mi0_1437 = buffer.data(mi0 + 1437);
    const auto *mi0_1438 = buffer.data(mi0 + 1438);
    const auto *mi0_1440 = buffer.data(mi0 + 1440);
    const auto *mi0_1442 = buffer.data(mi0 + 1442);
    const auto *mi0_1443 = buffer.data(mi0 + 1443);
    const auto *mi0_1445 = buffer.data(mi0 + 1445);
    const auto *mi0_1446 = buffer.data(mi0 + 1446);
    const auto *mi0_1448 = buffer.data(mi0 + 1448);
    const auto *mi0_1449 = buffer.data(mi0 + 1449);
    const auto *mi0_1451 = buffer.data(mi0 + 1451);
    const auto *mi0_1452 = buffer.data(mi0 + 1452);

    const auto *mi1_1372 = buffer.data(mi1 + 1372);
    const auto *mi1_1375 = buffer.data(mi1 + 1375);
    const auto *mi1_1377 = buffer.data(mi1 + 1377);
    const auto *mi1_1378 = buffer.data(mi1 + 1378);
    const auto *mi1_1381 = buffer.data(mi1 + 1381);
    const auto *mi1_1382 = buffer.data(mi1 + 1382);
    const auto *mi1_1384 = buffer.data(mi1 + 1384);
    const auto *mi1_1386 = buffer.data(mi1 + 1386);
    const auto *mi1_1387 = buffer.data(mi1 + 1387);
    const auto *mi1_1389 = buffer.data(mi1 + 1389);
    const auto *mi1_1390 = buffer.data(mi1 + 1390);
    const auto *mi1_1392 = buffer.data(mi1 + 1392);
    const auto *mi1_1393 = buffer.data(mi1 + 1393);
    const auto *mi1_1395 = buffer.data(mi1 + 1395);
    const auto *mi1_1396 = buffer.data(mi1 + 1396);
    const auto *mi1_1397 = buffer.data(mi1 + 1397);
    const auto *mi1_1398 = buffer.data(mi1 + 1398);
    const auto *mi1_1399 = buffer.data(mi1 + 1399);
    const auto *mi1_1400 = buffer.data(mi1 + 1400);
    const auto *mi1_1403 = buffer.data(mi1 + 1403);
    const auto *mi1_1405 = buffer.data(mi1 + 1405);
    const auto *mi1_1406 = buffer.data(mi1 + 1406);
    const auto *mi1_1409 = buffer.data(mi1 + 1409);
    const auto *mi1_1410 = buffer.data(mi1 + 1410);
    const auto *mi1_1412 = buffer.data(mi1 + 1412);
    const auto *mi1_1414 = buffer.data(mi1 + 1414);
    const auto *mi1_1415 = buffer.data(mi1 + 1415);
    const auto *mi1_1417 = buffer.data(mi1 + 1417);
    const auto *mi1_1418 = buffer.data(mi1 + 1418);
    const auto *mi1_1420 = buffer.data(mi1 + 1420);
    const auto *mi1_1421 = buffer.data(mi1 + 1421);
    const auto *mi1_1423 = buffer.data(mi1 + 1423);
    const auto *mi1_1424 = buffer.data(mi1 + 1424);
    const auto *mi1_1425 = buffer.data(mi1 + 1425);
    const auto *mi1_1426 = buffer.data(mi1 + 1426);
    const auto *mi1_1427 = buffer.data(mi1 + 1427);
    const auto *mi1_1428 = buffer.data(mi1 + 1428);
    const auto *mi1_1431 = buffer.data(mi1 + 1431);
    const auto *mi1_1433 = buffer.data(mi1 + 1433);
    const auto *mi1_1434 = buffer.data(mi1 + 1434);
    const auto *mi1_1437 = buffer.data(mi1 + 1437);
    const auto *mi1_1438 = buffer.data(mi1 + 1438);
    const auto *mi1_1440 = buffer.data(mi1 + 1440);
    const auto *mi1_1442 = buffer.data(mi1 + 1442);
    const auto *mi1_1443 = buffer.data(mi1 + 1443);
    const auto *mi1_1445 = buffer.data(mi1 + 1445);
    const auto *mi1_1446 = buffer.data(mi1 + 1446);
    const auto *mi1_1448 = buffer.data(mi1 + 1448);
    const auto *mi1_1449 = buffer.data(mi1 + 1449);
    const auto *mi1_1451 = buffer.data(mi1 + 1451);
    const auto *mi1_1452 = buffer.data(mi1 + 1452);

    const auto *mk_1764 = buffer.data(mk + 1764);
    const auto *mk_1766 = buffer.data(mk + 1766);
    const auto *mk_1767 = buffer.data(mk + 1767);
    const auto *mk_1769 = buffer.data(mk + 1769);
    const auto *mk_1770 = buffer.data(mk + 1770);
    const auto *mk_1773 = buffer.data(mk + 1773);
    const auto *mk_1774 = buffer.data(mk + 1774);
    const auto *mk_1776 = buffer.data(mk + 1776);
    const auto *mk_1778 = buffer.data(mk + 1778);
    const auto *mk_1779 = buffer.data(mk + 1779);
    const auto *mk_1781 = buffer.data(mk + 1781);
    const auto *mk_1782 = buffer.data(mk + 1782);
    const auto *mk_1784 = buffer.data(mk + 1784);
    const auto *mk_1785 = buffer.data(mk + 1785);
    const auto *mk_1787 = buffer.data(mk + 1787);
    const auto *mk_1788 = buffer.data(mk + 1788);
    const auto *mk_1789 = buffer.data(mk + 1789);
    const auto *mk_1791 = buffer.data(mk + 1791);
    const auto *mk_1792 = buffer.data(mk + 1792);
    const auto *mk_1793 = buffer.data(mk + 1793);
    const auto *mk_1794 = buffer.data(mk + 1794);
    const auto *mk_1795 = buffer.data(mk + 1795);
    const auto *mk_1796 = buffer.data(mk + 1796);
    const auto *mk_1797 = buffer.data(mk + 1797);
    const auto *mk_1798 = buffer.data(mk + 1798);
    const auto *mk_1799 = buffer.data(mk + 1799);
    const auto *mk_1800 = buffer.data(mk + 1800);
    const auto *mk_1802 = buffer.data(mk + 1802);
    const auto *mk_1803 = buffer.data(mk + 1803);
    const auto *mk_1805 = buffer.data(mk + 1805);
    const auto *mk_1806 = buffer.data(mk + 1806);
    const auto *mk_1809 = buffer.data(mk + 1809);
    const auto *mk_1810 = buffer.data(mk + 1810);
    const auto *mk_1812 = buffer.data(mk + 1812);
    const auto *mk_1814 = buffer.data(mk + 1814);
    const auto *mk_1815 = buffer.data(mk + 1815);
    const auto *mk_1817 = buffer.data(mk + 1817);
    const auto *mk_1818 = buffer.data(mk + 1818);
    const auto *mk_1820 = buffer.data(mk + 1820);
    const auto *mk_1821 = buffer.data(mk + 1821);
    const auto *mk_1823 = buffer.data(mk + 1823);
    const auto *mk_1824 = buffer.data(mk + 1824);
    const auto *mk_1825 = buffer.data(mk + 1825);
    const auto *mk_1827 = buffer.data(mk + 1827);
    const auto *mk_1828 = buffer.data(mk + 1828);
    const auto *mk_1829 = buffer.data(mk + 1829);
    const auto *mk_1830 = buffer.data(mk + 1830);
    const auto *mk_1831 = buffer.data(mk + 1831);
    const auto *mk_1832 = buffer.data(mk + 1832);
    const auto *mk_1833 = buffer.data(mk + 1833);
    const auto *mk_1834 = buffer.data(mk + 1834);
    const auto *mk_1835 = buffer.data(mk + 1835);
    const auto *mk_1836 = buffer.data(mk + 1836);
    const auto *mk_1838 = buffer.data(mk + 1838);
    const auto *mk_1839 = buffer.data(mk + 1839);
    const auto *mk_1841 = buffer.data(mk + 1841);
    const auto *mk_1842 = buffer.data(mk + 1842);
    const auto *mk_1845 = buffer.data(mk + 1845);
    const auto *mk_1846 = buffer.data(mk + 1846);
    const auto *mk_1848 = buffer.data(mk + 1848);
    const auto *mk_1850 = buffer.data(mk + 1850);
    const auto *mk_1851 = buffer.data(mk + 1851);
    const auto *mk_1853 = buffer.data(mk + 1853);
    const auto *mk_1854 = buffer.data(mk + 1854);
    const auto *mk_1856 = buffer.data(mk + 1856);
    const auto *mk_1857 = buffer.data(mk + 1857);
    const auto *mk_1859 = buffer.data(mk + 1859);
    const auto *mk_1860 = buffer.data(mk + 1860);

#pragma omp simd aligned(t_2205, t_2206, t_2207, t_2208, pb_x, pb_y, pb_z, lk_1404, lk_1440, \
                         mi0_1372, mi0_1375, mi1_1372, mi1_1375, mk_1764, \
                         mk_1767 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2205[k] = f_1 * mi0_1372[k]
                    - f_2 * mi1_1372[k]
                    + pb_x[k] * mk_1764[k];

        t_2206[k] = f_17 * lk_1440[k]
                    + pb_y[k] * mk_1764[k];

        t_2207[k] = f_16 * lk_1404[k]
                    + pb_z[k] * mk_1764[k];

        t_2208[k] = f_11 * mi0_1375[k]
                    - f_12 * mi1_1375[k]
                    + pb_x[k] * mk_1767[k];
    }

#pragma omp simd aligned(t_2209, t_2210, t_2211, pb_x, pb_y, lk_1442, mi0_1377, mi0_1378, \
                         mi1_1377, mi1_1378, mk_1766, mk_1769, \
                         mk_1770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2209[k] = f_17 * lk_1442[k]
                    + pb_y[k] * mk_1766[k];

        t_2210[k] = f_11 * mi0_1377[k]
                    - f_12 * mi1_1377[k]
                    + pb_x[k] * mk_1769[k];

        t_2211[k] = f_9 * mi0_1378[k]
                    - f_10 * mi1_1378[k]
                    + pb_x[k] * mk_1770[k];
    }

#pragma omp simd aligned(t_2212, t_2213, t_2214, pb_x, pb_y, pb_z, lk_1407, lk_1445, mi0_1381, \
                         mi1_1381, mk_1767, mk_1769, mk_1773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2212[k] = f_16 * lk_1407[k]
                    + pb_z[k] * mk_1767[k];

        t_2213[k] = f_17 * lk_1445[k]
                    + pb_y[k] * mk_1769[k];

        t_2214[k] = f_9 * mi0_1381[k]
                    - f_10 * mi1_1381[k]
                    + pb_x[k] * mk_1773[k];
    }

#pragma omp simd aligned(t_2215, t_2216, t_2217, pb_x, pb_z, lk_1410, mi0_1382, mi0_1384, \
                         mi1_1382, mi1_1384, mk_1770, mk_1774, \
                         mk_1776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2215[k] = f_7 * mi0_1382[k]
                    - f_8 * mi1_1382[k]
                    + pb_x[k] * mk_1774[k];

        t_2216[k] = f_16 * lk_1410[k]
                    + pb_z[k] * mk_1770[k];

        t_2217[k] = f_7 * mi0_1384[k]
                    - f_8 * mi1_1384[k]
                    + pb_x[k] * mk_1776[k];
    }

#pragma omp simd aligned(t_2218, t_2219, t_2220, pb_x, pb_y, lk_1449, mi0_1386, mi0_1387, \
                         mi1_1386, mi1_1387, mk_1773, mk_1778, \
                         mk_1779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2218[k] = f_17 * lk_1449[k]
                    + pb_y[k] * mk_1773[k];

        t_2219[k] = f_7 * mi0_1386[k]
                    - f_8 * mi1_1386[k]
                    + pb_x[k] * mk_1778[k];

        t_2220[k] = f_5 * mi0_1387[k]
                    - f_6 * mi1_1387[k]
                    + pb_x[k] * mk_1779[k];
    }

#pragma omp simd aligned(t_2221, t_2222, t_2223, pb_x, pb_z, lk_1414, mi0_1389, mi0_1390, \
                         mi1_1389, mi1_1390, mk_1774, mk_1781, \
                         mk_1782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2221[k] = f_16 * lk_1414[k]
                    + pb_z[k] * mk_1774[k];

        t_2222[k] = f_5 * mi0_1389[k]
                    - f_6 * mi1_1389[k]
                    + pb_x[k] * mk_1781[k];

        t_2223[k] = f_5 * mi0_1390[k]
                    - f_6 * mi1_1390[k]
                    + pb_x[k] * mk_1782[k];
    }

#pragma omp simd aligned(t_2224, t_2225, t_2226, pb_x, pb_y, lk_1454, mi0_1392, mi0_1393, \
                         mi1_1392, mi1_1393, mk_1778, mk_1784, \
                         mk_1785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2224[k] = f_17 * lk_1454[k]
                    + pb_y[k] * mk_1778[k];

        t_2225[k] = f_5 * mi0_1392[k]
                    - f_6 * mi1_1392[k]
                    + pb_x[k] * mk_1784[k];

        t_2226[k] = f_3 * mi0_1393[k]
                    - f_4 * mi1_1393[k]
                    + pb_x[k] * mk_1785[k];
    }

#pragma omp simd aligned(t_2227, t_2228, t_2229, pb_x, pb_z, lk_1419, mi0_1395, mi0_1396, \
                         mi1_1395, mi1_1396, mk_1779, mk_1787, \
                         mk_1788 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2227[k] = f_16 * lk_1419[k]
                    + pb_z[k] * mk_1779[k];

        t_2228[k] = f_3 * mi0_1395[k]
                    - f_4 * mi1_1395[k]
                    + pb_x[k] * mk_1787[k];

        t_2229[k] = f_3 * mi0_1396[k]
                    - f_4 * mi1_1396[k]
                    + pb_x[k] * mk_1788[k];
    }

#pragma omp simd aligned(t_2230, t_2231, t_2232, t_2233, pb_x, pb_y, lk_1460, mi0_1397, \
                         mi0_1399, mi1_1397, mi1_1399, mk_1784, mk_1789, mk_1791, \
                         mk_1792 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2230[k] = f_3 * mi0_1397[k]
                    - f_4 * mi1_1397[k]
                    + pb_x[k] * mk_1789[k];

        t_2231[k] = f_17 * lk_1460[k]
                    + pb_y[k] * mk_1784[k];

        t_2232[k] = f_3 * mi0_1399[k]
                    - f_4 * mi1_1399[k]
                    + pb_x[k] * mk_1791[k];

        t_2233[k] = pb_x[k] * mk_1792[k];
    }

#pragma omp simd aligned(t_2234, t_2235, t_2236, t_2237, t_2238, t_2239, t_2240, pb_x, \
                         mk_1793, mk_1794, mk_1795, mk_1796, mk_1797, mk_1798, \
                         mk_1799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2234[k] = pb_x[k] * mk_1793[k];

        t_2235[k] = pb_x[k] * mk_1794[k];

        t_2236[k] = pb_x[k] * mk_1795[k];

        t_2237[k] = pb_x[k] * mk_1796[k];

        t_2238[k] = pb_x[k] * mk_1797[k];

        t_2239[k] = pb_x[k] * mk_1798[k];

        t_2240[k] = pb_x[k] * mk_1799[k];
    }

#pragma omp simd aligned(t_2241, t_2242, t_2243, pa_z, pb_y, pb_z, kl0_1386, kl1_1386, \
                         lk_1432, lk_1470, ll_1791, mi0_1395, mi1_1395, mk_1792, \
                         mk_1794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2241[k] = f_29 * kl0_1386[k]
                    - f_30 * kl1_1386[k]
                    + pa_z[k] * ll_1791[k];

        t_2242[k] = f_16 * lk_1432[k]
                    + pb_z[k] * mk_1792[k];

        t_2243[k] = f_17 * lk_1470[k]
                    + f_11 * mi0_1395[k]
                    - f_12 * mi1_1395[k]
                    + pb_y[k] * mk_1794[k];
    }

#pragma omp simd aligned(t_2244, t_2245, t_2246, pb_y, lk_1471, lk_1472, lk_1473, mi0_1396, \
                         mi0_1397, mi0_1398, mi1_1396, mi1_1397, mi1_1398, mk_1795, mk_1796, \
                         mk_1797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2244[k] = f_17 * lk_1471[k]
                    + f_9 * mi0_1396[k]
                    - f_10 * mi1_1396[k]
                    + pb_y[k] * mk_1795[k];

        t_2245[k] = f_17 * lk_1472[k]
                    + f_7 * mi0_1397[k]
                    - f_8 * mi1_1397[k]
                    + pb_y[k] * mk_1796[k];

        t_2246[k] = f_17 * lk_1473[k]
                    + f_5 * mi0_1398[k]
                    - f_6 * mi1_1398[k]
                    + pb_y[k] * mk_1797[k];
    }

#pragma omp simd aligned(t_2247, t_2248, t_2249, pa_y, pb_y, kl0_1484, kl1_1484, lk_1474, \
                         lk_1475, ll_1844, mi0_1399, mi1_1399, mk_1798, \
                         mk_1799 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2247[k] = f_17 * lk_1474[k]
                    + f_3 * mi0_1399[k]
                    - f_4 * mi1_1399[k]
                    + pb_y[k] * mk_1798[k];

        t_2248[k] = f_17 * lk_1475[k]
                    + pb_y[k] * mk_1799[k];

        t_2249[k] = f_31 * kl0_1484[k]
                    - f_32 * kl1_1484[k]
                    + pa_y[k] * ll_1844[k];
    }

#pragma omp simd aligned(t_2250, t_2251, t_2252, t_2253, pb_x, pb_y, pb_z, lk_1440, lk_1476, \
                         mi0_1400, mi0_1403, mi1_1400, mi1_1403, mk_1800, \
                         mk_1803 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2250[k] = f_1 * mi0_1400[k]
                    - f_2 * mi1_1400[k]
                    + pb_x[k] * mk_1800[k];

        t_2251[k] = f_16 * lk_1476[k]
                    + pb_y[k] * mk_1800[k];

        t_2252[k] = f_17 * lk_1440[k]
                    + pb_z[k] * mk_1800[k];

        t_2253[k] = f_11 * mi0_1403[k]
                    - f_12 * mi1_1403[k]
                    + pb_x[k] * mk_1803[k];
    }

#pragma omp simd aligned(t_2254, t_2255, t_2256, pb_x, pb_y, lk_1478, mi0_1405, mi0_1406, \
                         mi1_1405, mi1_1406, mk_1802, mk_1805, \
                         mk_1806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2254[k] = f_16 * lk_1478[k]
                    + pb_y[k] * mk_1802[k];

        t_2255[k] = f_11 * mi0_1405[k]
                    - f_12 * mi1_1405[k]
                    + pb_x[k] * mk_1805[k];

        t_2256[k] = f_9 * mi0_1406[k]
                    - f_10 * mi1_1406[k]
                    + pb_x[k] * mk_1806[k];
    }

#pragma omp simd aligned(t_2257, t_2258, t_2259, pb_x, pb_y, pb_z, lk_1443, lk_1481, mi0_1409, \
                         mi1_1409, mk_1803, mk_1805, mk_1809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2257[k] = f_17 * lk_1443[k]
                    + pb_z[k] * mk_1803[k];

        t_2258[k] = f_16 * lk_1481[k]
                    + pb_y[k] * mk_1805[k];

        t_2259[k] = f_9 * mi0_1409[k]
                    - f_10 * mi1_1409[k]
                    + pb_x[k] * mk_1809[k];
    }

#pragma omp simd aligned(t_2260, t_2261, t_2262, pb_x, pb_z, lk_1446, mi0_1410, mi0_1412, \
                         mi1_1410, mi1_1412, mk_1806, mk_1810, \
                         mk_1812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2260[k] = f_7 * mi0_1410[k]
                    - f_8 * mi1_1410[k]
                    + pb_x[k] * mk_1810[k];

        t_2261[k] = f_17 * lk_1446[k]
                    + pb_z[k] * mk_1806[k];

        t_2262[k] = f_7 * mi0_1412[k]
                    - f_8 * mi1_1412[k]
                    + pb_x[k] * mk_1812[k];
    }

#pragma omp simd aligned(t_2263, t_2264, t_2265, pb_x, pb_y, lk_1485, mi0_1414, mi0_1415, \
                         mi1_1414, mi1_1415, mk_1809, mk_1814, \
                         mk_1815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2263[k] = f_16 * lk_1485[k]
                    + pb_y[k] * mk_1809[k];

        t_2264[k] = f_7 * mi0_1414[k]
                    - f_8 * mi1_1414[k]
                    + pb_x[k] * mk_1814[k];

        t_2265[k] = f_5 * mi0_1415[k]
                    - f_6 * mi1_1415[k]
                    + pb_x[k] * mk_1815[k];
    }

#pragma omp simd aligned(t_2266, t_2267, t_2268, pb_x, pb_z, lk_1450, mi0_1417, mi0_1418, \
                         mi1_1417, mi1_1418, mk_1810, mk_1817, \
                         mk_1818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2266[k] = f_17 * lk_1450[k]
                    + pb_z[k] * mk_1810[k];

        t_2267[k] = f_5 * mi0_1417[k]
                    - f_6 * mi1_1417[k]
                    + pb_x[k] * mk_1817[k];

        t_2268[k] = f_5 * mi0_1418[k]
                    - f_6 * mi1_1418[k]
                    + pb_x[k] * mk_1818[k];
    }

#pragma omp simd aligned(t_2269, t_2270, t_2271, pb_x, pb_y, lk_1490, mi0_1420, mi0_1421, \
                         mi1_1420, mi1_1421, mk_1814, mk_1820, \
                         mk_1821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2269[k] = f_16 * lk_1490[k]
                    + pb_y[k] * mk_1814[k];

        t_2270[k] = f_5 * mi0_1420[k]
                    - f_6 * mi1_1420[k]
                    + pb_x[k] * mk_1820[k];

        t_2271[k] = f_3 * mi0_1421[k]
                    - f_4 * mi1_1421[k]
                    + pb_x[k] * mk_1821[k];
    }

#pragma omp simd aligned(t_2272, t_2273, t_2274, pb_x, pb_z, lk_1455, mi0_1423, mi0_1424, \
                         mi1_1423, mi1_1424, mk_1815, mk_1823, \
                         mk_1824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2272[k] = f_17 * lk_1455[k]
                    + pb_z[k] * mk_1815[k];

        t_2273[k] = f_3 * mi0_1423[k]
                    - f_4 * mi1_1423[k]
                    + pb_x[k] * mk_1823[k];

        t_2274[k] = f_3 * mi0_1424[k]
                    - f_4 * mi1_1424[k]
                    + pb_x[k] * mk_1824[k];
    }

#pragma omp simd aligned(t_2275, t_2276, t_2277, t_2278, pb_x, pb_y, lk_1496, mi0_1425, \
                         mi0_1427, mi1_1425, mi1_1427, mk_1820, mk_1825, mk_1827, \
                         mk_1828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2275[k] = f_3 * mi0_1425[k]
                    - f_4 * mi1_1425[k]
                    + pb_x[k] * mk_1825[k];

        t_2276[k] = f_16 * lk_1496[k]
                    + pb_y[k] * mk_1820[k];

        t_2277[k] = f_3 * mi0_1427[k]
                    - f_4 * mi1_1427[k]
                    + pb_x[k] * mk_1827[k];

        t_2278[k] = pb_x[k] * mk_1828[k];
    }

#pragma omp simd aligned(t_2279, t_2280, t_2281, t_2282, t_2283, t_2284, t_2285, pb_x, \
                         mk_1829, mk_1830, mk_1831, mk_1832, mk_1833, mk_1834, \
                         mk_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2279[k] = pb_x[k] * mk_1829[k];

        t_2280[k] = pb_x[k] * mk_1830[k];

        t_2281[k] = pb_x[k] * mk_1831[k];

        t_2282[k] = pb_x[k] * mk_1832[k];

        t_2283[k] = pb_x[k] * mk_1833[k];

        t_2284[k] = pb_x[k] * mk_1834[k];

        t_2285[k] = pb_x[k] * mk_1835[k];
    }

#pragma omp simd aligned(t_2286, t_2287, t_2288, pa_z, pb_y, pb_z, kl0_1431, kl1_1431, \
                         lk_1468, lk_1506, ll_1836, mi0_1423, mi1_1423, mk_1828, \
                         mk_1830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2286[k] = f_31 * kl0_1431[k]
                    - f_32 * kl1_1431[k]
                    + pa_z[k] * ll_1836[k];

        t_2287[k] = f_17 * lk_1468[k]
                    + pb_z[k] * mk_1828[k];

        t_2288[k] = f_16 * lk_1506[k]
                    + f_11 * mi0_1423[k]
                    - f_12 * mi1_1423[k]
                    + pb_y[k] * mk_1830[k];
    }

#pragma omp simd aligned(t_2289, t_2290, t_2291, pb_y, lk_1507, lk_1508, lk_1509, mi0_1424, \
                         mi0_1425, mi0_1426, mi1_1424, mi1_1425, mi1_1426, mk_1831, mk_1832, \
                         mk_1833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2289[k] = f_16 * lk_1507[k]
                    + f_9 * mi0_1424[k]
                    - f_10 * mi1_1424[k]
                    + pb_y[k] * mk_1831[k];

        t_2290[k] = f_16 * lk_1508[k]
                    + f_7 * mi0_1425[k]
                    - f_8 * mi1_1425[k]
                    + pb_y[k] * mk_1832[k];

        t_2291[k] = f_16 * lk_1509[k]
                    + f_5 * mi0_1426[k]
                    - f_6 * mi1_1426[k]
                    + pb_y[k] * mk_1833[k];
    }

#pragma omp simd aligned(t_2292, t_2293, t_2294, pa_y, pb_y, kl0_1529, kl1_1529, lk_1510, \
                         lk_1511, ll_1889, mi0_1427, mi1_1427, mk_1834, \
                         mk_1835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2292[k] = f_16 * lk_1510[k]
                    + f_3 * mi0_1427[k]
                    - f_4 * mi1_1427[k]
                    + pb_y[k] * mk_1834[k];

        t_2293[k] = f_16 * lk_1511[k]
                    + pb_y[k] * mk_1835[k];

        t_2294[k] = f_29 * kl0_1529[k]
                    - f_30 * kl1_1529[k]
                    + pa_y[k] * ll_1889[k];
    }

#pragma omp simd aligned(t_2295, t_2296, t_2297, t_2298, pb_x, pb_y, pb_z, lk_1476, lk_1512, \
                         mi0_1428, mi0_1431, mi1_1428, mi1_1431, mk_1836, \
                         mk_1839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2295[k] = f_1 * mi0_1428[k]
                    - f_2 * mi1_1428[k]
                    + pb_x[k] * mk_1836[k];

        t_2296[k] = f_15 * lk_1512[k]
                    + pb_y[k] * mk_1836[k];

        t_2297[k] = f_18 * lk_1476[k]
                    + pb_z[k] * mk_1836[k];

        t_2298[k] = f_11 * mi0_1431[k]
                    - f_12 * mi1_1431[k]
                    + pb_x[k] * mk_1839[k];
    }

#pragma omp simd aligned(t_2299, t_2300, t_2301, pb_x, pb_y, lk_1514, mi0_1433, mi0_1434, \
                         mi1_1433, mi1_1434, mk_1838, mk_1841, \
                         mk_1842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2299[k] = f_15 * lk_1514[k]
                    + pb_y[k] * mk_1838[k];

        t_2300[k] = f_11 * mi0_1433[k]
                    - f_12 * mi1_1433[k]
                    + pb_x[k] * mk_1841[k];

        t_2301[k] = f_9 * mi0_1434[k]
                    - f_10 * mi1_1434[k]
                    + pb_x[k] * mk_1842[k];
    }

#pragma omp simd aligned(t_2302, t_2303, t_2304, pb_x, pb_y, pb_z, lk_1479, lk_1517, mi0_1437, \
                         mi1_1437, mk_1839, mk_1841, mk_1845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2302[k] = f_18 * lk_1479[k]
                    + pb_z[k] * mk_1839[k];

        t_2303[k] = f_15 * lk_1517[k]
                    + pb_y[k] * mk_1841[k];

        t_2304[k] = f_9 * mi0_1437[k]
                    - f_10 * mi1_1437[k]
                    + pb_x[k] * mk_1845[k];
    }

#pragma omp simd aligned(t_2305, t_2306, t_2307, pb_x, pb_z, lk_1482, mi0_1438, mi0_1440, \
                         mi1_1438, mi1_1440, mk_1842, mk_1846, \
                         mk_1848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2305[k] = f_7 * mi0_1438[k]
                    - f_8 * mi1_1438[k]
                    + pb_x[k] * mk_1846[k];

        t_2306[k] = f_18 * lk_1482[k]
                    + pb_z[k] * mk_1842[k];

        t_2307[k] = f_7 * mi0_1440[k]
                    - f_8 * mi1_1440[k]
                    + pb_x[k] * mk_1848[k];
    }

#pragma omp simd aligned(t_2308, t_2309, t_2310, pb_x, pb_y, lk_1521, mi0_1442, mi0_1443, \
                         mi1_1442, mi1_1443, mk_1845, mk_1850, \
                         mk_1851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2308[k] = f_15 * lk_1521[k]
                    + pb_y[k] * mk_1845[k];

        t_2309[k] = f_7 * mi0_1442[k]
                    - f_8 * mi1_1442[k]
                    + pb_x[k] * mk_1850[k];

        t_2310[k] = f_5 * mi0_1443[k]
                    - f_6 * mi1_1443[k]
                    + pb_x[k] * mk_1851[k];
    }

#pragma omp simd aligned(t_2311, t_2312, t_2313, pb_x, pb_z, lk_1486, mi0_1445, mi0_1446, \
                         mi1_1445, mi1_1446, mk_1846, mk_1853, \
                         mk_1854 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2311[k] = f_18 * lk_1486[k]
                    + pb_z[k] * mk_1846[k];

        t_2312[k] = f_5 * mi0_1445[k]
                    - f_6 * mi1_1445[k]
                    + pb_x[k] * mk_1853[k];

        t_2313[k] = f_5 * mi0_1446[k]
                    - f_6 * mi1_1446[k]
                    + pb_x[k] * mk_1854[k];
    }

#pragma omp simd aligned(t_2314, t_2315, t_2316, pb_x, pb_y, lk_1526, mi0_1448, mi0_1449, \
                         mi1_1448, mi1_1449, mk_1850, mk_1856, \
                         mk_1857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2314[k] = f_15 * lk_1526[k]
                    + pb_y[k] * mk_1850[k];

        t_2315[k] = f_5 * mi0_1448[k]
                    - f_6 * mi1_1448[k]
                    + pb_x[k] * mk_1856[k];

        t_2316[k] = f_3 * mi0_1449[k]
                    - f_4 * mi1_1449[k]
                    + pb_x[k] * mk_1857[k];
    }

#pragma omp simd aligned(t_2317, t_2318, t_2319, pb_x, pb_z, lk_1491, mi0_1451, mi0_1452, \
                         mi1_1451, mi1_1452, mk_1851, mk_1859, \
                         mk_1860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2317[k] = f_18 * lk_1491[k]
                    + pb_z[k] * mk_1851[k];

        t_2318[k] = f_3 * mi0_1451[k]
                    - f_4 * mi1_1451[k]
                    + pb_x[k] * mk_1859[k];

        t_2319[k] = f_3 * mi0_1452[k]
                    - f_4 * mi1_1452[k]
                    + pb_x[k] * mk_1860[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece18(CSimdMatrix &buffer, const size_t target,
                                             const size_t pa, const size_t pb, const size_t kl0,
                                             const size_t kl1, const size_t lk, const size_t ll,
                                             const size_t mi0, const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
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
    const auto f_18 = 3.0 / p;
    const auto f_19 = 4.0 / p;
    const auto f_20 = 0.5 / alpha;
    const auto f_21 = 0.5 * beta / (alpha * p);
    const auto f_22 = 3.5 / p;
    const auto f_23 = 3.0 / alpha;
    const auto f_24 = 3.0 * beta / (alpha * p);
    const auto f_25 = 1.0 / alpha;
    const auto f_26 = beta / (alpha * p);
    const auto f_27 = 2.5 / alpha;
    const auto f_28 = 2.5 * beta / (alpha * p);

    auto *t_2320 = buffer.data(target + 2320);
    auto *t_2321 = buffer.data(target + 2321);
    auto *t_2322 = buffer.data(target + 2322);
    auto *t_2323 = buffer.data(target + 2323);
    auto *t_2324 = buffer.data(target + 2324);
    auto *t_2325 = buffer.data(target + 2325);
    auto *t_2326 = buffer.data(target + 2326);
    auto *t_2327 = buffer.data(target + 2327);
    auto *t_2328 = buffer.data(target + 2328);
    auto *t_2329 = buffer.data(target + 2329);
    auto *t_2330 = buffer.data(target + 2330);
    auto *t_2331 = buffer.data(target + 2331);
    auto *t_2332 = buffer.data(target + 2332);
    auto *t_2333 = buffer.data(target + 2333);
    auto *t_2334 = buffer.data(target + 2334);
    auto *t_2335 = buffer.data(target + 2335);
    auto *t_2336 = buffer.data(target + 2336);
    auto *t_2337 = buffer.data(target + 2337);
    auto *t_2338 = buffer.data(target + 2338);
    auto *t_2339 = buffer.data(target + 2339);
    auto *t_2340 = buffer.data(target + 2340);
    auto *t_2341 = buffer.data(target + 2341);
    auto *t_2342 = buffer.data(target + 2342);
    auto *t_2343 = buffer.data(target + 2343);
    auto *t_2344 = buffer.data(target + 2344);
    auto *t_2345 = buffer.data(target + 2345);
    auto *t_2346 = buffer.data(target + 2346);
    auto *t_2347 = buffer.data(target + 2347);
    auto *t_2348 = buffer.data(target + 2348);
    auto *t_2349 = buffer.data(target + 2349);
    auto *t_2350 = buffer.data(target + 2350);
    auto *t_2351 = buffer.data(target + 2351);
    auto *t_2352 = buffer.data(target + 2352);
    auto *t_2353 = buffer.data(target + 2353);
    auto *t_2354 = buffer.data(target + 2354);
    auto *t_2355 = buffer.data(target + 2355);
    auto *t_2356 = buffer.data(target + 2356);
    auto *t_2357 = buffer.data(target + 2357);
    auto *t_2358 = buffer.data(target + 2358);
    auto *t_2359 = buffer.data(target + 2359);
    auto *t_2360 = buffer.data(target + 2360);
    auto *t_2361 = buffer.data(target + 2361);
    auto *t_2362 = buffer.data(target + 2362);
    auto *t_2363 = buffer.data(target + 2363);
    auto *t_2364 = buffer.data(target + 2364);
    auto *t_2365 = buffer.data(target + 2365);
    auto *t_2366 = buffer.data(target + 2366);
    auto *t_2367 = buffer.data(target + 2367);
    auto *t_2368 = buffer.data(target + 2368);
    auto *t_2369 = buffer.data(target + 2369);
    auto *t_2370 = buffer.data(target + 2370);
    auto *t_2371 = buffer.data(target + 2371);
    auto *t_2372 = buffer.data(target + 2372);
    auto *t_2373 = buffer.data(target + 2373);
    auto *t_2374 = buffer.data(target + 2374);
    auto *t_2375 = buffer.data(target + 2375);
    auto *t_2376 = buffer.data(target + 2376);
    auto *t_2377 = buffer.data(target + 2377);
    auto *t_2378 = buffer.data(target + 2378);
    auto *t_2379 = buffer.data(target + 2379);
    auto *t_2380 = buffer.data(target + 2380);
    auto *t_2381 = buffer.data(target + 2381);
    auto *t_2382 = buffer.data(target + 2382);
    auto *t_2383 = buffer.data(target + 2383);
    auto *t_2384 = buffer.data(target + 2384);
    auto *t_2385 = buffer.data(target + 2385);
    auto *t_2386 = buffer.data(target + 2386);
    auto *t_2387 = buffer.data(target + 2387);
    auto *t_2388 = buffer.data(target + 2388);
    auto *t_2389 = buffer.data(target + 2389);
    auto *t_2390 = buffer.data(target + 2390);
    auto *t_2391 = buffer.data(target + 2391);
    auto *t_2392 = buffer.data(target + 2392);
    auto *t_2393 = buffer.data(target + 2393);
    auto *t_2394 = buffer.data(target + 2394);
    auto *t_2395 = buffer.data(target + 2395);
    auto *t_2396 = buffer.data(target + 2396);
    auto *t_2397 = buffer.data(target + 2397);
    auto *t_2398 = buffer.data(target + 2398);
    auto *t_2399 = buffer.data(target + 2399);
    auto *t_2400 = buffer.data(target + 2400);
    auto *t_2401 = buffer.data(target + 2401);
    auto *t_2402 = buffer.data(target + 2402);
    auto *t_2403 = buffer.data(target + 2403);
    auto *t_2404 = buffer.data(target + 2404);
    auto *t_2405 = buffer.data(target + 2405);
    auto *t_2406 = buffer.data(target + 2406);
    auto *t_2407 = buffer.data(target + 2407);
    auto *t_2408 = buffer.data(target + 2408);
    auto *t_2409 = buffer.data(target + 2409);
    auto *t_2410 = buffer.data(target + 2410);
    auto *t_2411 = buffer.data(target + 2411);
    auto *t_2412 = buffer.data(target + 2412);
    auto *t_2413 = buffer.data(target + 2413);
    auto *t_2414 = buffer.data(target + 2414);
    auto *t_2415 = buffer.data(target + 2415);
    auto *t_2416 = buffer.data(target + 2416);
    auto *t_2417 = buffer.data(target + 2417);
    auto *t_2418 = buffer.data(target + 2418);
    auto *t_2419 = buffer.data(target + 2419);
    auto *t_2420 = buffer.data(target + 2420);
    auto *t_2421 = buffer.data(target + 2421);
    auto *t_2422 = buffer.data(target + 2422);
    auto *t_2423 = buffer.data(target + 2423);
    auto *t_2424 = buffer.data(target + 2424);
    auto *t_2425 = buffer.data(target + 2425);
    auto *t_2426 = buffer.data(target + 2426);
    auto *t_2427 = buffer.data(target + 2427);
    auto *t_2428 = buffer.data(target + 2428);
    auto *t_2429 = buffer.data(target + 2429);
    auto *t_2430 = buffer.data(target + 2430);
    auto *t_2431 = buffer.data(target + 2431);
    auto *t_2432 = buffer.data(target + 2432);
    auto *t_2433 = buffer.data(target + 2433);
    auto *t_2434 = buffer.data(target + 2434);
    auto *t_2435 = buffer.data(target + 2435);
    auto *t_2436 = buffer.data(target + 2436);
    auto *t_2437 = buffer.data(target + 2437);
    auto *t_2438 = buffer.data(target + 2438);
    auto *t_2439 = buffer.data(target + 2439);
    auto *t_2440 = buffer.data(target + 2440);
    auto *t_2441 = buffer.data(target + 2441);
    auto *t_2442 = buffer.data(target + 2442);
    auto *t_2443 = buffer.data(target + 2443);
    auto *t_2444 = buffer.data(target + 2444);
    auto *t_2445 = buffer.data(target + 2445);
    auto *t_2446 = buffer.data(target + 2446);
    auto *t_2447 = buffer.data(target + 2447);
    auto *t_2448 = buffer.data(target + 2448);
    auto *t_2449 = buffer.data(target + 2449);
    auto *t_2450 = buffer.data(target + 2450);
    auto *t_2451 = buffer.data(target + 2451);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kl0_1476 = buffer.data(kl0 + 1476);
    const auto *kl0_1521 = buffer.data(kl0 + 1521);
    const auto *kl0_1574 = buffer.data(kl0 + 1574);
    const auto *kl0_1619 = buffer.data(kl0 + 1619);

    const auto *kl1_1476 = buffer.data(kl1 + 1476);
    const auto *kl1_1521 = buffer.data(kl1 + 1521);
    const auto *kl1_1574 = buffer.data(kl1 + 1574);
    const auto *kl1_1619 = buffer.data(kl1 + 1619);

    const auto *lk_1504 = buffer.data(lk + 1504);
    const auto *lk_1512 = buffer.data(lk + 1512);
    const auto *lk_1515 = buffer.data(lk + 1515);
    const auto *lk_1518 = buffer.data(lk + 1518);
    const auto *lk_1522 = buffer.data(lk + 1522);
    const auto *lk_1527 = buffer.data(lk + 1527);
    const auto *lk_1532 = buffer.data(lk + 1532);
    const auto *lk_1540 = buffer.data(lk + 1540);
    const auto *lk_1542 = buffer.data(lk + 1542);
    const auto *lk_1543 = buffer.data(lk + 1543);
    const auto *lk_1544 = buffer.data(lk + 1544);
    const auto *lk_1545 = buffer.data(lk + 1545);
    const auto *lk_1546 = buffer.data(lk + 1546);
    const auto *lk_1547 = buffer.data(lk + 1547);
    const auto *lk_1548 = buffer.data(lk + 1548);
    const auto *lk_1550 = buffer.data(lk + 1550);
    const auto *lk_1551 = buffer.data(lk + 1551);
    const auto *lk_1553 = buffer.data(lk + 1553);
    const auto *lk_1554 = buffer.data(lk + 1554);
    const auto *lk_1557 = buffer.data(lk + 1557);
    const auto *lk_1558 = buffer.data(lk + 1558);
    const auto *lk_1562 = buffer.data(lk + 1562);
    const auto *lk_1563 = buffer.data(lk + 1563);
    const auto *lk_1568 = buffer.data(lk + 1568);
    const auto *lk_1576 = buffer.data(lk + 1576);
    const auto *lk_1578 = buffer.data(lk + 1578);
    const auto *lk_1579 = buffer.data(lk + 1579);
    const auto *lk_1580 = buffer.data(lk + 1580);
    const auto *lk_1581 = buffer.data(lk + 1581);
    const auto *lk_1582 = buffer.data(lk + 1582);
    const auto *lk_1583 = buffer.data(lk + 1583);
    const auto *lk_1584 = buffer.data(lk + 1584);
    const auto *lk_1585 = buffer.data(lk + 1585);
    const auto *lk_1586 = buffer.data(lk + 1586);
    const auto *lk_1587 = buffer.data(lk + 1587);
    const auto *lk_1589 = buffer.data(lk + 1589);
    const auto *lk_1590 = buffer.data(lk + 1590);
    const auto *lk_1592 = buffer.data(lk + 1592);
    const auto *lk_1593 = buffer.data(lk + 1593);
    const auto *lk_1594 = buffer.data(lk + 1594);
    const auto *lk_1596 = buffer.data(lk + 1596);
    const auto *lk_1597 = buffer.data(lk + 1597);
    const auto *lk_1598 = buffer.data(lk + 1598);
    const auto *lk_1599 = buffer.data(lk + 1599);
    const auto *lk_1601 = buffer.data(lk + 1601);
    const auto *lk_1602 = buffer.data(lk + 1602);
    const auto *lk_1603 = buffer.data(lk + 1603);
    const auto *lk_1604 = buffer.data(lk + 1604);
    const auto *lk_1612 = buffer.data(lk + 1612);
    const auto *lk_1614 = buffer.data(lk + 1614);
    const auto *lk_1615 = buffer.data(lk + 1615);
    const auto *lk_1616 = buffer.data(lk + 1616);
    const auto *lk_1617 = buffer.data(lk + 1617);
    const auto *lk_1618 = buffer.data(lk + 1618);
    const auto *lk_1619 = buffer.data(lk + 1619);

    const auto *ll_1881 = buffer.data(ll + 1881);
    const auto *ll_1926 = buffer.data(ll + 1926);
    const auto *ll_1934 = buffer.data(ll + 1934);
    const auto *ll_1979 = buffer.data(ll + 1979);
    const auto *ll_1980 = buffer.data(ll + 1980);
    const auto *ll_1982 = buffer.data(ll + 1982);
    const auto *ll_1983 = buffer.data(ll + 1983);
    const auto *ll_1985 = buffer.data(ll + 1985);
    const auto *ll_1986 = buffer.data(ll + 1986);
    const auto *ll_1989 = buffer.data(ll + 1989);
    const auto *ll_1990 = buffer.data(ll + 1990);
    const auto *ll_1992 = buffer.data(ll + 1992);
    const auto *ll_1994 = buffer.data(ll + 1994);
    const auto *ll_1995 = buffer.data(ll + 1995);
    const auto *ll_1997 = buffer.data(ll + 1997);
    const auto *ll_1998 = buffer.data(ll + 1998);
    const auto *ll_2000 = buffer.data(ll + 2000);
    const auto *ll_2001 = buffer.data(ll + 2001);
    const auto *ll_2003 = buffer.data(ll + 2003);
    const auto *ll_2004 = buffer.data(ll + 2004);
    const auto *ll_2005 = buffer.data(ll + 2005);
    const auto *ll_2007 = buffer.data(ll + 2007);
    const auto *ll_2016 = buffer.data(ll + 2016);
    const auto *ll_2018 = buffer.data(ll + 2018);
    const auto *ll_2019 = buffer.data(ll + 2019);
    const auto *ll_2020 = buffer.data(ll + 2020);
    const auto *ll_2021 = buffer.data(ll + 2021);
    const auto *ll_2022 = buffer.data(ll + 2022);
    const auto *ll_2024 = buffer.data(ll + 2024);

    const auto *mi0_1451 = buffer.data(mi0 + 1451);
    const auto *mi0_1452 = buffer.data(mi0 + 1452);
    const auto *mi0_1453 = buffer.data(mi0 + 1453);
    const auto *mi0_1454 = buffer.data(mi0 + 1454);
    const auto *mi0_1455 = buffer.data(mi0 + 1455);
    const auto *mi0_1456 = buffer.data(mi0 + 1456);
    const auto *mi0_1459 = buffer.data(mi0 + 1459);
    const auto *mi0_1461 = buffer.data(mi0 + 1461);
    const auto *mi0_1462 = buffer.data(mi0 + 1462);
    const auto *mi0_1465 = buffer.data(mi0 + 1465);
    const auto *mi0_1466 = buffer.data(mi0 + 1466);
    const auto *mi0_1468 = buffer.data(mi0 + 1468);
    const auto *mi0_1470 = buffer.data(mi0 + 1470);
    const auto *mi0_1471 = buffer.data(mi0 + 1471);
    const auto *mi0_1473 = buffer.data(mi0 + 1473);
    const auto *mi0_1474 = buffer.data(mi0 + 1474);
    const auto *mi0_1476 = buffer.data(mi0 + 1476);
    const auto *mi0_1477 = buffer.data(mi0 + 1477);
    const auto *mi0_1479 = buffer.data(mi0 + 1479);
    const auto *mi0_1480 = buffer.data(mi0 + 1480);
    const auto *mi0_1481 = buffer.data(mi0 + 1481);
    const auto *mi0_1482 = buffer.data(mi0 + 1482);
    const auto *mi0_1483 = buffer.data(mi0 + 1483);
    const auto *mi0_1512 = buffer.data(mi0 + 1512);
    const auto *mi0_1515 = buffer.data(mi0 + 1515);
    const auto *mi0_1517 = buffer.data(mi0 + 1517);
    const auto *mi0_1518 = buffer.data(mi0 + 1518);
    const auto *mi0_1521 = buffer.data(mi0 + 1521);
    const auto *mi0_1522 = buffer.data(mi0 + 1522);
    const auto *mi0_1524 = buffer.data(mi0 + 1524);
    const auto *mi0_1526 = buffer.data(mi0 + 1526);
    const auto *mi0_1527 = buffer.data(mi0 + 1527);
    const auto *mi0_1529 = buffer.data(mi0 + 1529);
    const auto *mi0_1530 = buffer.data(mi0 + 1530);
    const auto *mi0_1532 = buffer.data(mi0 + 1532);
    const auto *mi0_1533 = buffer.data(mi0 + 1533);

    const auto *mi1_1451 = buffer.data(mi1 + 1451);
    const auto *mi1_1452 = buffer.data(mi1 + 1452);
    const auto *mi1_1453 = buffer.data(mi1 + 1453);
    const auto *mi1_1454 = buffer.data(mi1 + 1454);
    const auto *mi1_1455 = buffer.data(mi1 + 1455);
    const auto *mi1_1456 = buffer.data(mi1 + 1456);
    const auto *mi1_1459 = buffer.data(mi1 + 1459);
    const auto *mi1_1461 = buffer.data(mi1 + 1461);
    const auto *mi1_1462 = buffer.data(mi1 + 1462);
    const auto *mi1_1465 = buffer.data(mi1 + 1465);
    const auto *mi1_1466 = buffer.data(mi1 + 1466);
    const auto *mi1_1468 = buffer.data(mi1 + 1468);
    const auto *mi1_1470 = buffer.data(mi1 + 1470);
    const auto *mi1_1471 = buffer.data(mi1 + 1471);
    const auto *mi1_1473 = buffer.data(mi1 + 1473);
    const auto *mi1_1474 = buffer.data(mi1 + 1474);
    const auto *mi1_1476 = buffer.data(mi1 + 1476);
    const auto *mi1_1477 = buffer.data(mi1 + 1477);
    const auto *mi1_1479 = buffer.data(mi1 + 1479);
    const auto *mi1_1480 = buffer.data(mi1 + 1480);
    const auto *mi1_1481 = buffer.data(mi1 + 1481);
    const auto *mi1_1482 = buffer.data(mi1 + 1482);
    const auto *mi1_1483 = buffer.data(mi1 + 1483);
    const auto *mi1_1512 = buffer.data(mi1 + 1512);
    const auto *mi1_1515 = buffer.data(mi1 + 1515);
    const auto *mi1_1517 = buffer.data(mi1 + 1517);
    const auto *mi1_1518 = buffer.data(mi1 + 1518);
    const auto *mi1_1521 = buffer.data(mi1 + 1521);
    const auto *mi1_1522 = buffer.data(mi1 + 1522);
    const auto *mi1_1524 = buffer.data(mi1 + 1524);
    const auto *mi1_1526 = buffer.data(mi1 + 1526);
    const auto *mi1_1527 = buffer.data(mi1 + 1527);
    const auto *mi1_1529 = buffer.data(mi1 + 1529);
    const auto *mi1_1530 = buffer.data(mi1 + 1530);
    const auto *mi1_1532 = buffer.data(mi1 + 1532);
    const auto *mi1_1533 = buffer.data(mi1 + 1533);

    const auto *mk_1856 = buffer.data(mk + 1856);
    const auto *mk_1861 = buffer.data(mk + 1861);
    const auto *mk_1863 = buffer.data(mk + 1863);
    const auto *mk_1864 = buffer.data(mk + 1864);
    const auto *mk_1865 = buffer.data(mk + 1865);
    const auto *mk_1866 = buffer.data(mk + 1866);
    const auto *mk_1867 = buffer.data(mk + 1867);
    const auto *mk_1868 = buffer.data(mk + 1868);
    const auto *mk_1869 = buffer.data(mk + 1869);
    const auto *mk_1870 = buffer.data(mk + 1870);
    const auto *mk_1871 = buffer.data(mk + 1871);
    const auto *mk_1872 = buffer.data(mk + 1872);
    const auto *mk_1874 = buffer.data(mk + 1874);
    const auto *mk_1875 = buffer.data(mk + 1875);
    const auto *mk_1877 = buffer.data(mk + 1877);
    const auto *mk_1878 = buffer.data(mk + 1878);
    const auto *mk_1881 = buffer.data(mk + 1881);
    const auto *mk_1882 = buffer.data(mk + 1882);
    const auto *mk_1884 = buffer.data(mk + 1884);
    const auto *mk_1886 = buffer.data(mk + 1886);
    const auto *mk_1887 = buffer.data(mk + 1887);
    const auto *mk_1889 = buffer.data(mk + 1889);
    const auto *mk_1890 = buffer.data(mk + 1890);
    const auto *mk_1892 = buffer.data(mk + 1892);
    const auto *mk_1893 = buffer.data(mk + 1893);
    const auto *mk_1895 = buffer.data(mk + 1895);
    const auto *mk_1896 = buffer.data(mk + 1896);
    const auto *mk_1897 = buffer.data(mk + 1897);
    const auto *mk_1899 = buffer.data(mk + 1899);
    const auto *mk_1900 = buffer.data(mk + 1900);
    const auto *mk_1901 = buffer.data(mk + 1901);
    const auto *mk_1902 = buffer.data(mk + 1902);
    const auto *mk_1903 = buffer.data(mk + 1903);
    const auto *mk_1904 = buffer.data(mk + 1904);
    const auto *mk_1905 = buffer.data(mk + 1905);
    const auto *mk_1906 = buffer.data(mk + 1906);
    const auto *mk_1907 = buffer.data(mk + 1907);
    const auto *mk_1908 = buffer.data(mk + 1908);
    const auto *mk_1910 = buffer.data(mk + 1910);
    const auto *mk_1911 = buffer.data(mk + 1911);
    const auto *mk_1913 = buffer.data(mk + 1913);
    const auto *mk_1914 = buffer.data(mk + 1914);
    const auto *mk_1917 = buffer.data(mk + 1917);
    const auto *mk_1918 = buffer.data(mk + 1918);
    const auto *mk_1922 = buffer.data(mk + 1922);
    const auto *mk_1923 = buffer.data(mk + 1923);
    const auto *mk_1928 = buffer.data(mk + 1928);
    const auto *mk_1936 = buffer.data(mk + 1936);
    const auto *mk_1937 = buffer.data(mk + 1937);
    const auto *mk_1938 = buffer.data(mk + 1938);
    const auto *mk_1939 = buffer.data(mk + 1939);
    const auto *mk_1940 = buffer.data(mk + 1940);
    const auto *mk_1941 = buffer.data(mk + 1941);
    const auto *mk_1942 = buffer.data(mk + 1942);
    const auto *mk_1943 = buffer.data(mk + 1943);
    const auto *mk_1944 = buffer.data(mk + 1944);
    const auto *mk_1946 = buffer.data(mk + 1946);
    const auto *mk_1947 = buffer.data(mk + 1947);
    const auto *mk_1949 = buffer.data(mk + 1949);
    const auto *mk_1950 = buffer.data(mk + 1950);
    const auto *mk_1953 = buffer.data(mk + 1953);
    const auto *mk_1954 = buffer.data(mk + 1954);
    const auto *mk_1956 = buffer.data(mk + 1956);
    const auto *mk_1958 = buffer.data(mk + 1958);
    const auto *mk_1959 = buffer.data(mk + 1959);
    const auto *mk_1961 = buffer.data(mk + 1961);
    const auto *mk_1962 = buffer.data(mk + 1962);
    const auto *mk_1964 = buffer.data(mk + 1964);
    const auto *mk_1965 = buffer.data(mk + 1965);

#pragma omp simd aligned(t_2320, t_2321, t_2322, t_2323, pb_x, pb_y, lk_1532, mi0_1453, \
                         mi0_1455, mi1_1453, mi1_1455, mk_1856, mk_1861, mk_1863, \
                         mk_1864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2320[k] = f_3 * mi0_1453[k]
                    - f_4 * mi1_1453[k]
                    + pb_x[k] * mk_1861[k];

        t_2321[k] = f_15 * lk_1532[k]
                    + pb_y[k] * mk_1856[k];

        t_2322[k] = f_3 * mi0_1455[k]
                    - f_4 * mi1_1455[k]
                    + pb_x[k] * mk_1863[k];

        t_2323[k] = pb_x[k] * mk_1864[k];
    }

#pragma omp simd aligned(t_2324, t_2325, t_2326, t_2327, t_2328, t_2329, t_2330, pb_x, \
                         mk_1865, mk_1866, mk_1867, mk_1868, mk_1869, mk_1870, \
                         mk_1871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2324[k] = pb_x[k] * mk_1865[k];

        t_2325[k] = pb_x[k] * mk_1866[k];

        t_2326[k] = pb_x[k] * mk_1867[k];

        t_2327[k] = pb_x[k] * mk_1868[k];

        t_2328[k] = pb_x[k] * mk_1869[k];

        t_2329[k] = pb_x[k] * mk_1870[k];

        t_2330[k] = pb_x[k] * mk_1871[k];
    }

#pragma omp simd aligned(t_2331, t_2332, t_2333, pa_z, pb_y, pb_z, kl0_1476, kl1_1476, \
                         lk_1504, lk_1542, ll_1881, mi0_1451, mi1_1451, mk_1864, \
                         mk_1866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2331[k] = f_27 * kl0_1476[k]
                    - f_28 * kl1_1476[k]
                    + pa_z[k] * ll_1881[k];

        t_2332[k] = f_18 * lk_1504[k]
                    + pb_z[k] * mk_1864[k];

        t_2333[k] = f_15 * lk_1542[k]
                    + f_11 * mi0_1451[k]
                    - f_12 * mi1_1451[k]
                    + pb_y[k] * mk_1866[k];
    }

#pragma omp simd aligned(t_2334, t_2335, t_2336, pb_y, lk_1543, lk_1544, lk_1545, mi0_1452, \
                         mi0_1453, mi0_1454, mi1_1452, mi1_1453, mi1_1454, mk_1867, mk_1868, \
                         mk_1869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2334[k] = f_15 * lk_1543[k]
                    + f_9 * mi0_1452[k]
                    - f_10 * mi1_1452[k]
                    + pb_y[k] * mk_1867[k];

        t_2335[k] = f_15 * lk_1544[k]
                    + f_7 * mi0_1453[k]
                    - f_8 * mi1_1453[k]
                    + pb_y[k] * mk_1868[k];

        t_2336[k] = f_15 * lk_1545[k]
                    + f_5 * mi0_1454[k]
                    - f_6 * mi1_1454[k]
                    + pb_y[k] * mk_1869[k];
    }

#pragma omp simd aligned(t_2337, t_2338, t_2339, pa_y, pb_y, kl0_1574, kl1_1574, lk_1546, \
                         lk_1547, ll_1934, mi0_1455, mi1_1455, mk_1870, \
                         mk_1871 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2337[k] = f_15 * lk_1546[k]
                    + f_3 * mi0_1455[k]
                    - f_4 * mi1_1455[k]
                    + pb_y[k] * mk_1870[k];

        t_2338[k] = f_15 * lk_1547[k]
                    + pb_y[k] * mk_1871[k];

        t_2339[k] = f_25 * kl0_1574[k]
                    - f_26 * kl1_1574[k]
                    + pa_y[k] * ll_1934[k];
    }

#pragma omp simd aligned(t_2340, t_2341, t_2342, t_2343, pb_x, pb_y, pb_z, lk_1512, lk_1548, \
                         mi0_1456, mi0_1459, mi1_1456, mi1_1459, mk_1872, \
                         mk_1875 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2340[k] = f_1 * mi0_1456[k]
                    - f_2 * mi1_1456[k]
                    + pb_x[k] * mk_1872[k];

        t_2341[k] = f_14 * lk_1548[k]
                    + pb_y[k] * mk_1872[k];

        t_2342[k] = f_22 * lk_1512[k]
                    + pb_z[k] * mk_1872[k];

        t_2343[k] = f_11 * mi0_1459[k]
                    - f_12 * mi1_1459[k]
                    + pb_x[k] * mk_1875[k];
    }

#pragma omp simd aligned(t_2344, t_2345, t_2346, pb_x, pb_y, lk_1550, mi0_1461, mi0_1462, \
                         mi1_1461, mi1_1462, mk_1874, mk_1877, \
                         mk_1878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2344[k] = f_14 * lk_1550[k]
                    + pb_y[k] * mk_1874[k];

        t_2345[k] = f_11 * mi0_1461[k]
                    - f_12 * mi1_1461[k]
                    + pb_x[k] * mk_1877[k];

        t_2346[k] = f_9 * mi0_1462[k]
                    - f_10 * mi1_1462[k]
                    + pb_x[k] * mk_1878[k];
    }

#pragma omp simd aligned(t_2347, t_2348, t_2349, pb_x, pb_y, pb_z, lk_1515, lk_1553, mi0_1465, \
                         mi1_1465, mk_1875, mk_1877, mk_1881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2347[k] = f_22 * lk_1515[k]
                    + pb_z[k] * mk_1875[k];

        t_2348[k] = f_14 * lk_1553[k]
                    + pb_y[k] * mk_1877[k];

        t_2349[k] = f_9 * mi0_1465[k]
                    - f_10 * mi1_1465[k]
                    + pb_x[k] * mk_1881[k];
    }

#pragma omp simd aligned(t_2350, t_2351, t_2352, pb_x, pb_z, lk_1518, mi0_1466, mi0_1468, \
                         mi1_1466, mi1_1468, mk_1878, mk_1882, \
                         mk_1884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2350[k] = f_7 * mi0_1466[k]
                    - f_8 * mi1_1466[k]
                    + pb_x[k] * mk_1882[k];

        t_2351[k] = f_22 * lk_1518[k]
                    + pb_z[k] * mk_1878[k];

        t_2352[k] = f_7 * mi0_1468[k]
                    - f_8 * mi1_1468[k]
                    + pb_x[k] * mk_1884[k];
    }

#pragma omp simd aligned(t_2353, t_2354, t_2355, pb_x, pb_y, lk_1557, mi0_1470, mi0_1471, \
                         mi1_1470, mi1_1471, mk_1881, mk_1886, \
                         mk_1887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2353[k] = f_14 * lk_1557[k]
                    + pb_y[k] * mk_1881[k];

        t_2354[k] = f_7 * mi0_1470[k]
                    - f_8 * mi1_1470[k]
                    + pb_x[k] * mk_1886[k];

        t_2355[k] = f_5 * mi0_1471[k]
                    - f_6 * mi1_1471[k]
                    + pb_x[k] * mk_1887[k];
    }

#pragma omp simd aligned(t_2356, t_2357, t_2358, pb_x, pb_z, lk_1522, mi0_1473, mi0_1474, \
                         mi1_1473, mi1_1474, mk_1882, mk_1889, \
                         mk_1890 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2356[k] = f_22 * lk_1522[k]
                    + pb_z[k] * mk_1882[k];

        t_2357[k] = f_5 * mi0_1473[k]
                    - f_6 * mi1_1473[k]
                    + pb_x[k] * mk_1889[k];

        t_2358[k] = f_5 * mi0_1474[k]
                    - f_6 * mi1_1474[k]
                    + pb_x[k] * mk_1890[k];
    }

#pragma omp simd aligned(t_2359, t_2360, t_2361, pb_x, pb_y, lk_1562, mi0_1476, mi0_1477, \
                         mi1_1476, mi1_1477, mk_1886, mk_1892, \
                         mk_1893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2359[k] = f_14 * lk_1562[k]
                    + pb_y[k] * mk_1886[k];

        t_2360[k] = f_5 * mi0_1476[k]
                    - f_6 * mi1_1476[k]
                    + pb_x[k] * mk_1892[k];

        t_2361[k] = f_3 * mi0_1477[k]
                    - f_4 * mi1_1477[k]
                    + pb_x[k] * mk_1893[k];
    }

#pragma omp simd aligned(t_2362, t_2363, t_2364, pb_x, pb_z, lk_1527, mi0_1479, mi0_1480, \
                         mi1_1479, mi1_1480, mk_1887, mk_1895, \
                         mk_1896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2362[k] = f_22 * lk_1527[k]
                    + pb_z[k] * mk_1887[k];

        t_2363[k] = f_3 * mi0_1479[k]
                    - f_4 * mi1_1479[k]
                    + pb_x[k] * mk_1895[k];

        t_2364[k] = f_3 * mi0_1480[k]
                    - f_4 * mi1_1480[k]
                    + pb_x[k] * mk_1896[k];
    }

#pragma omp simd aligned(t_2365, t_2366, t_2367, t_2368, pb_x, pb_y, lk_1568, mi0_1481, \
                         mi0_1483, mi1_1481, mi1_1483, mk_1892, mk_1897, mk_1899, \
                         mk_1900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2365[k] = f_3 * mi0_1481[k]
                    - f_4 * mi1_1481[k]
                    + pb_x[k] * mk_1897[k];

        t_2366[k] = f_14 * lk_1568[k]
                    + pb_y[k] * mk_1892[k];

        t_2367[k] = f_3 * mi0_1483[k]
                    - f_4 * mi1_1483[k]
                    + pb_x[k] * mk_1899[k];

        t_2368[k] = pb_x[k] * mk_1900[k];
    }

#pragma omp simd aligned(t_2369, t_2370, t_2371, t_2372, t_2373, t_2374, t_2375, pb_x, \
                         mk_1901, mk_1902, mk_1903, mk_1904, mk_1905, mk_1906, \
                         mk_1907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2369[k] = pb_x[k] * mk_1901[k];

        t_2370[k] = pb_x[k] * mk_1902[k];

        t_2371[k] = pb_x[k] * mk_1903[k];

        t_2372[k] = pb_x[k] * mk_1904[k];

        t_2373[k] = pb_x[k] * mk_1905[k];

        t_2374[k] = pb_x[k] * mk_1906[k];

        t_2375[k] = pb_x[k] * mk_1907[k];
    }

#pragma omp simd aligned(t_2376, t_2377, t_2378, pa_z, pb_y, pb_z, kl0_1521, kl1_1521, \
                         lk_1540, lk_1578, ll_1926, mi0_1479, mi1_1479, mk_1900, \
                         mk_1902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2376[k] = f_23 * kl0_1521[k]
                    - f_24 * kl1_1521[k]
                    + pa_z[k] * ll_1926[k];

        t_2377[k] = f_22 * lk_1540[k]
                    + pb_z[k] * mk_1900[k];

        t_2378[k] = f_14 * lk_1578[k]
                    + f_11 * mi0_1479[k]
                    - f_12 * mi1_1479[k]
                    + pb_y[k] * mk_1902[k];
    }

#pragma omp simd aligned(t_2379, t_2380, t_2381, pb_y, lk_1579, lk_1580, lk_1581, mi0_1480, \
                         mi0_1481, mi0_1482, mi1_1480, mi1_1481, mi1_1482, mk_1903, mk_1904, \
                         mk_1905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2379[k] = f_14 * lk_1579[k]
                    + f_9 * mi0_1480[k]
                    - f_10 * mi1_1480[k]
                    + pb_y[k] * mk_1903[k];

        t_2380[k] = f_14 * lk_1580[k]
                    + f_7 * mi0_1481[k]
                    - f_8 * mi1_1481[k]
                    + pb_y[k] * mk_1904[k];

        t_2381[k] = f_14 * lk_1581[k]
                    + f_5 * mi0_1482[k]
                    - f_6 * mi1_1482[k]
                    + pb_y[k] * mk_1905[k];
    }

#pragma omp simd aligned(t_2382, t_2383, t_2384, t_2385, pa_y, pb_y, kl0_1619, kl1_1619, \
                         lk_1582, lk_1583, ll_1979, ll_1980, mi0_1483, mi1_1483, mk_1906, \
                         mk_1907 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2382[k] = f_14 * lk_1582[k]
                    + f_3 * mi0_1483[k]
                    - f_4 * mi1_1483[k]
                    + pb_y[k] * mk_1906[k];

        t_2383[k] = f_14 * lk_1583[k]
                    + pb_y[k] * mk_1907[k];

        t_2384[k] = f_20 * kl0_1619[k]
                    - f_21 * kl1_1619[k]
                    + pa_y[k] * ll_1979[k];

        t_2385[k] = pa_y[k] * ll_1980[k];
    }

#pragma omp simd aligned(t_2386, t_2387, t_2388, t_2389, t_2390, pa_y, pb_y, lk_1584, lk_1585, \
                         lk_1586, ll_1982, ll_1983, ll_1985, mk_1908, \
                         mk_1910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2386[k] = f_13 * lk_1584[k]
                    + pb_y[k] * mk_1908[k];

        t_2387[k] = pa_y[k] * ll_1982[k];

        t_2388[k] = f_14 * lk_1585[k]
                    + pa_y[k] * ll_1983[k];

        t_2389[k] = f_13 * lk_1586[k]
                    + pb_y[k] * mk_1910[k];

        t_2390[k] = pa_y[k] * ll_1985[k];
    }

#pragma omp simd aligned(t_2391, t_2392, t_2393, t_2394, pa_y, pb_y, pb_z, lk_1551, lk_1587, \
                         lk_1589, ll_1986, ll_1989, mk_1911, mk_1913 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2391[k] = f_15 * lk_1587[k]
                    + pa_y[k] * ll_1986[k];

        t_2392[k] = f_19 * lk_1551[k]
                    + pb_z[k] * mk_1911[k];

        t_2393[k] = f_13 * lk_1589[k]
                    + pb_y[k] * mk_1913[k];

        t_2394[k] = pa_y[k] * ll_1989[k];
    }

#pragma omp simd aligned(t_2395, t_2396, t_2397, t_2398, pa_y, pb_y, pb_z, lk_1554, lk_1590, \
                         lk_1592, lk_1593, ll_1990, ll_1992, mk_1914, \
                         mk_1917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2395[k] = f_16 * lk_1590[k]
                    + pa_y[k] * ll_1990[k];

        t_2396[k] = f_19 * lk_1554[k]
                    + pb_z[k] * mk_1914[k];

        t_2397[k] = f_14 * lk_1592[k]
                    + pa_y[k] * ll_1992[k];

        t_2398[k] = f_13 * lk_1593[k]
                    + pb_y[k] * mk_1917[k];
    }

#pragma omp simd aligned(t_2399, t_2400, t_2401, t_2402, t_2403, pa_y, pb_z, lk_1558, lk_1594, \
                         lk_1596, lk_1597, ll_1994, ll_1995, ll_1997, ll_1998, \
                         mk_1918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2399[k] = pa_y[k] * ll_1994[k];

        t_2400[k] = f_17 * lk_1594[k]
                    + pa_y[k] * ll_1995[k];

        t_2401[k] = f_19 * lk_1558[k]
                    + pb_z[k] * mk_1918[k];

        t_2402[k] = f_15 * lk_1596[k]
                    + pa_y[k] * ll_1997[k];

        t_2403[k] = f_14 * lk_1597[k]
                    + pa_y[k] * ll_1998[k];
    }

#pragma omp simd aligned(t_2404, t_2405, t_2406, t_2407, pa_y, pb_y, pb_z, lk_1563, lk_1598, \
                         lk_1599, ll_2000, ll_2001, mk_1922, mk_1923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2404[k] = f_13 * lk_1598[k]
                    + pb_y[k] * mk_1922[k];

        t_2405[k] = pa_y[k] * ll_2000[k];

        t_2406[k] = f_18 * lk_1599[k]
                    + pa_y[k] * ll_2001[k];

        t_2407[k] = f_19 * lk_1563[k]
                    + pb_z[k] * mk_1923[k];
    }

#pragma omp simd aligned(t_2408, t_2409, t_2410, t_2411, t_2412, pa_y, pb_y, lk_1601, lk_1602, \
                         lk_1603, lk_1604, ll_2003, ll_2004, ll_2005, ll_2007, \
                         mk_1928 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2408[k] = f_16 * lk_1601[k]
                    + pa_y[k] * ll_2003[k];

        t_2409[k] = f_15 * lk_1602[k]
                    + pa_y[k] * ll_2004[k];

        t_2410[k] = f_14 * lk_1603[k]
                    + pa_y[k] * ll_2005[k];

        t_2411[k] = f_13 * lk_1604[k]
                    + pb_y[k] * mk_1928[k];

        t_2412[k] = pa_y[k] * ll_2007[k];
    }

#pragma omp simd aligned(t_2413, t_2414, t_2415, t_2416, t_2417, t_2418, t_2419, pb_x, \
                         mk_1936, mk_1937, mk_1938, mk_1939, mk_1940, mk_1941, \
                         mk_1942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2413[k] = pb_x[k] * mk_1936[k];

        t_2414[k] = pb_x[k] * mk_1937[k];

        t_2415[k] = pb_x[k] * mk_1938[k];

        t_2416[k] = pb_x[k] * mk_1939[k];

        t_2417[k] = pb_x[k] * mk_1940[k];

        t_2418[k] = pb_x[k] * mk_1941[k];

        t_2419[k] = pb_x[k] * mk_1942[k];
    }

#pragma omp simd aligned(t_2420, t_2421, t_2422, t_2423, pa_y, pb_x, pb_z, lk_1576, lk_1612, \
                         lk_1614, ll_2016, ll_2018, mk_1936, mk_1943 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2420[k] = pb_x[k] * mk_1943[k];

        t_2421[k] = f_19 * lk_1612[k]
                    + pa_y[k] * ll_2016[k];

        t_2422[k] = f_19 * lk_1576[k]
                    + pb_z[k] * mk_1936[k];

        t_2423[k] = f_18 * lk_1614[k]
                    + pa_y[k] * ll_2018[k];
    }

#pragma omp simd aligned(t_2424, t_2425, t_2426, t_2427, pa_y, lk_1615, lk_1616, lk_1617, \
                         lk_1618, ll_2019, ll_2020, ll_2021, ll_2022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2424[k] = f_17 * lk_1615[k]
                    + pa_y[k] * ll_2019[k];

        t_2425[k] = f_16 * lk_1616[k]
                    + pa_y[k] * ll_2020[k];

        t_2426[k] = f_15 * lk_1617[k]
                    + pa_y[k] * ll_2021[k];

        t_2427[k] = f_14 * lk_1618[k]
                    + pa_y[k] * ll_2022[k];
    }

#pragma omp simd aligned(t_2428, t_2429, t_2430, t_2431, t_2432, pa_y, pb_x, pb_y, pb_z, \
                         lk_1584, lk_1619, ll_2024, mi0_1512, mi1_1512, mk_1943, \
                         mk_1944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2428[k] = f_13 * lk_1619[k]
                    + pb_y[k] * mk_1943[k];

        t_2429[k] = pa_y[k] * ll_2024[k];

        t_2430[k] = f_1 * mi0_1512[k]
                    - f_2 * mi1_1512[k]
                    + pb_x[k] * mk_1944[k];

        t_2431[k] = pb_y[k] * mk_1944[k];

        t_2432[k] = f_0 * lk_1584[k]
                    + pb_z[k] * mk_1944[k];
    }

#pragma omp simd aligned(t_2433, t_2434, t_2435, t_2436, pb_x, pb_y, mi0_1515, mi0_1517, \
                         mi0_1518, mi1_1515, mi1_1517, mi1_1518, mk_1946, mk_1947, mk_1949, \
                         mk_1950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2433[k] = f_11 * mi0_1515[k]
                    - f_12 * mi1_1515[k]
                    + pb_x[k] * mk_1947[k];

        t_2434[k] = pb_y[k] * mk_1946[k];

        t_2435[k] = f_11 * mi0_1517[k]
                    - f_12 * mi1_1517[k]
                    + pb_x[k] * mk_1949[k];

        t_2436[k] = f_9 * mi0_1518[k]
                    - f_10 * mi1_1518[k]
                    + pb_x[k] * mk_1950[k];
    }

#pragma omp simd aligned(t_2437, t_2438, t_2439, t_2440, pb_x, pb_y, pb_z, lk_1587, mi0_1521, \
                         mi0_1522, mi1_1521, mi1_1522, mk_1947, mk_1949, mk_1953, \
                         mk_1954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2437[k] = f_0 * lk_1587[k]
                    + pb_z[k] * mk_1947[k];

        t_2438[k] = pb_y[k] * mk_1949[k];

        t_2439[k] = f_9 * mi0_1521[k]
                    - f_10 * mi1_1521[k]
                    + pb_x[k] * mk_1953[k];

        t_2440[k] = f_7 * mi0_1522[k]
                    - f_8 * mi1_1522[k]
                    + pb_x[k] * mk_1954[k];
    }

#pragma omp simd aligned(t_2441, t_2442, t_2443, t_2444, pb_x, pb_y, pb_z, lk_1590, mi0_1524, \
                         mi0_1526, mi1_1524, mi1_1526, mk_1950, mk_1953, mk_1956, \
                         mk_1958 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2441[k] = f_0 * lk_1590[k]
                    + pb_z[k] * mk_1950[k];

        t_2442[k] = f_7 * mi0_1524[k]
                    - f_8 * mi1_1524[k]
                    + pb_x[k] * mk_1956[k];

        t_2443[k] = pb_y[k] * mk_1953[k];

        t_2444[k] = f_7 * mi0_1526[k]
                    - f_8 * mi1_1526[k]
                    + pb_x[k] * mk_1958[k];
    }

#pragma omp simd aligned(t_2445, t_2446, t_2447, pb_x, pb_z, lk_1594, mi0_1527, mi0_1529, \
                         mi1_1527, mi1_1529, mk_1954, mk_1959, \
                         mk_1961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2445[k] = f_5 * mi0_1527[k]
                    - f_6 * mi1_1527[k]
                    + pb_x[k] * mk_1959[k];

        t_2446[k] = f_0 * lk_1594[k]
                    + pb_z[k] * mk_1954[k];

        t_2447[k] = f_5 * mi0_1529[k]
                    - f_6 * mi1_1529[k]
                    + pb_x[k] * mk_1961[k];
    }

#pragma omp simd aligned(t_2448, t_2449, t_2450, t_2451, pb_x, pb_y, mi0_1530, mi0_1532, \
                         mi0_1533, mi1_1530, mi1_1532, mi1_1533, mk_1958, mk_1962, mk_1964, \
                         mk_1965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2448[k] = f_5 * mi0_1530[k]
                    - f_6 * mi1_1530[k]
                    + pb_x[k] * mk_1962[k];

        t_2449[k] = pb_y[k] * mk_1958[k];

        t_2450[k] = f_5 * mi0_1532[k]
                    - f_6 * mi1_1532[k]
                    + pb_x[k] * mk_1964[k];

        t_2451[k] = f_3 * mi0_1533[k]
                    - f_4 * mi1_1533[k]
                    + pb_x[k] * mk_1965[k];
    }
}

static auto
compute_prim_ml_electron_repulsion_0_piece19(CSimdMatrix &buffer, const size_t target,
                                             const size_t pb, const size_t lk, const size_t mi0,
                                             const size_t mi1, const size_t mk,
                                             const size_t ncols, const double alpha,
                                             const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
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

    auto *t_2452 = buffer.data(target + 2452);
    auto *t_2453 = buffer.data(target + 2453);
    auto *t_2454 = buffer.data(target + 2454);
    auto *t_2455 = buffer.data(target + 2455);
    auto *t_2456 = buffer.data(target + 2456);
    auto *t_2457 = buffer.data(target + 2457);
    auto *t_2458 = buffer.data(target + 2458);
    auto *t_2459 = buffer.data(target + 2459);
    auto *t_2460 = buffer.data(target + 2460);
    auto *t_2461 = buffer.data(target + 2461);
    auto *t_2462 = buffer.data(target + 2462);
    auto *t_2463 = buffer.data(target + 2463);
    auto *t_2464 = buffer.data(target + 2464);
    auto *t_2465 = buffer.data(target + 2465);
    auto *t_2466 = buffer.data(target + 2466);
    auto *t_2467 = buffer.data(target + 2467);
    auto *t_2468 = buffer.data(target + 2468);
    auto *t_2469 = buffer.data(target + 2469);
    auto *t_2470 = buffer.data(target + 2470);
    auto *t_2471 = buffer.data(target + 2471);
    auto *t_2472 = buffer.data(target + 2472);
    auto *t_2473 = buffer.data(target + 2473);
    auto *t_2474 = buffer.data(target + 2474);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *lk_1599 = buffer.data(lk + 1599);
    const auto *lk_1612 = buffer.data(lk + 1612);
    const auto *lk_1619 = buffer.data(lk + 1619);

    const auto *mi0_1533 = buffer.data(mi0 + 1533);
    const auto *mi0_1535 = buffer.data(mi0 + 1535);
    const auto *mi0_1536 = buffer.data(mi0 + 1536);
    const auto *mi0_1537 = buffer.data(mi0 + 1537);
    const auto *mi0_1538 = buffer.data(mi0 + 1538);
    const auto *mi0_1539 = buffer.data(mi0 + 1539);

    const auto *mi1_1533 = buffer.data(mi1 + 1533);
    const auto *mi1_1535 = buffer.data(mi1 + 1535);
    const auto *mi1_1536 = buffer.data(mi1 + 1536);
    const auto *mi1_1537 = buffer.data(mi1 + 1537);
    const auto *mi1_1538 = buffer.data(mi1 + 1538);
    const auto *mi1_1539 = buffer.data(mi1 + 1539);

    const auto *mk_1959 = buffer.data(mk + 1959);
    const auto *mk_1964 = buffer.data(mk + 1964);
    const auto *mk_1967 = buffer.data(mk + 1967);
    const auto *mk_1968 = buffer.data(mk + 1968);
    const auto *mk_1969 = buffer.data(mk + 1969);
    const auto *mk_1971 = buffer.data(mk + 1971);
    const auto *mk_1972 = buffer.data(mk + 1972);
    const auto *mk_1973 = buffer.data(mk + 1973);
    const auto *mk_1974 = buffer.data(mk + 1974);
    const auto *mk_1975 = buffer.data(mk + 1975);
    const auto *mk_1976 = buffer.data(mk + 1976);
    const auto *mk_1977 = buffer.data(mk + 1977);
    const auto *mk_1978 = buffer.data(mk + 1978);
    const auto *mk_1979 = buffer.data(mk + 1979);

#pragma omp simd aligned(t_2452, t_2453, t_2454, pb_x, pb_z, lk_1599, mi0_1535, mi0_1536, \
                         mi1_1535, mi1_1536, mk_1959, mk_1967, \
                         mk_1968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2452[k] = f_0 * lk_1599[k]
                    + pb_z[k] * mk_1959[k];

        t_2453[k] = f_3 * mi0_1535[k]
                    - f_4 * mi1_1535[k]
                    + pb_x[k] * mk_1967[k];

        t_2454[k] = f_3 * mi0_1536[k]
                    - f_4 * mi1_1536[k]
                    + pb_x[k] * mk_1968[k];
    }

#pragma omp simd aligned(t_2455, t_2456, t_2457, t_2458, t_2459, pb_x, pb_y, mi0_1537, \
                         mi0_1539, mi1_1537, mi1_1539, mk_1964, mk_1969, mk_1971, mk_1972, \
                         mk_1973 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2455[k] = f_3 * mi0_1537[k]
                    - f_4 * mi1_1537[k]
                    + pb_x[k] * mk_1969[k];

        t_2456[k] = pb_y[k] * mk_1964[k];

        t_2457[k] = f_3 * mi0_1539[k]
                    - f_4 * mi1_1539[k]
                    + pb_x[k] * mk_1971[k];

        t_2458[k] = pb_x[k] * mk_1972[k];

        t_2459[k] = pb_x[k] * mk_1973[k];
    }

#pragma omp simd aligned(t_2460, t_2461, t_2462, t_2463, t_2464, t_2465, pb_x, mk_1974, \
                         mk_1975, mk_1976, mk_1977, mk_1978, mk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2460[k] = pb_x[k] * mk_1974[k];

        t_2461[k] = pb_x[k] * mk_1975[k];

        t_2462[k] = pb_x[k] * mk_1976[k];

        t_2463[k] = pb_x[k] * mk_1977[k];

        t_2464[k] = pb_x[k] * mk_1978[k];

        t_2465[k] = pb_x[k] * mk_1979[k];
    }

#pragma omp simd aligned(t_2466, t_2467, t_2468, t_2469, pb_y, pb_z, lk_1612, mi0_1533, \
                         mi0_1535, mi0_1536, mi1_1533, mi1_1535, mi1_1536, mk_1972, mk_1974, \
                         mk_1975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2466[k] = f_1 * mi0_1533[k]
                    - f_2 * mi1_1533[k]
                    + pb_y[k] * mk_1972[k];

        t_2467[k] = f_0 * lk_1612[k]
                    + pb_z[k] * mk_1972[k];

        t_2468[k] = f_11 * mi0_1535[k]
                    - f_12 * mi1_1535[k]
                    + pb_y[k] * mk_1974[k];

        t_2469[k] = f_9 * mi0_1536[k]
                    - f_10 * mi1_1536[k]
                    + pb_y[k] * mk_1975[k];
    }

#pragma omp simd aligned(t_2470, t_2471, t_2472, t_2473, pb_y, mi0_1537, mi0_1538, mi0_1539, \
                         mi1_1537, mi1_1538, mi1_1539, mk_1976, mk_1977, mk_1978, \
                         mk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2470[k] = f_7 * mi0_1537[k]
                    - f_8 * mi1_1537[k]
                    + pb_y[k] * mk_1976[k];

        t_2471[k] = f_5 * mi0_1538[k]
                    - f_6 * mi1_1538[k]
                    + pb_y[k] * mk_1977[k];

        t_2472[k] = f_3 * mi0_1539[k]
                    - f_4 * mi1_1539[k]
                    + pb_y[k] * mk_1978[k];

        t_2473[k] = pb_y[k] * mk_1979[k];
    }

#pragma omp simd aligned(t_2474, pb_z, lk_1619, mi0_1539, mi1_1539, \
                         mk_1979 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2474[k] = f_0 * lk_1619[k]
                    + f_1 * mi0_1539[k]
                    - f_2 * mi1_1539[k]
                    + pb_z[k] * mk_1979[k];
    }
}

auto
compute_prim_ml_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t kl0, const size_t kl1,
                                     const size_t lk, const size_t ll, const size_t mi0,
                                     const size_t mi1, const size_t mk, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_ml_electron_repulsion_0_piece0(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece1(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece2(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece3(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece4(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece5(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece6(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece7(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece8(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece9(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece10(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                 mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece11(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                 mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece12(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                 mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece13(buffer, target, pa, pb, lk, ll, mk, ncols, p);

    compute_prim_ml_electron_repulsion_0_piece14(buffer, target, pa, pb, lk, ll, mk, ncols, p);

    compute_prim_ml_electron_repulsion_0_piece15(buffer, target, pa, pb, lk, ll, mi0, mi1, mk,
                                                 ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece16(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                 mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece17(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                 mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece18(buffer, target, pa, pb, kl0, kl1, lk, ll, mi0,
                                                 mi1, mk, ncols, alpha, beta, p);

    compute_prim_ml_electron_repulsion_0_piece19(buffer, target, pb, lk, mi0, mi1, mk, ncols,
                                                 alpha, beta, p);
}

}  // namespace simdt2ceri
