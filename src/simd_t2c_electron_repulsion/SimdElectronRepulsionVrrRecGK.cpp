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


#include "SimdElectronRepulsionVrrRecGK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_gk_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t dk0,
                                            const size_t dk1, const size_t fi, const size_t fk,
                                            const size_t gh0, const size_t gh1, const size_t gi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_136 = buffer.data(dk0 + 136);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_136 = buffer.data(dk1 + 136);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_28 = buffer.data(fk + 28);
    const auto *fk_30 = buffer.data(fk + 30);
    const auto *fk_31 = buffer.data(fk + 31);
    const auto *fk_32 = buffer.data(fk + 32);
    const auto *fk_33 = buffer.data(fk + 33);
    const auto *fk_35 = buffer.data(fk + 35);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_136 = buffer.data(fk + 136);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_1 = buffer.data(gh0 + 1);
    const auto *gh0_2 = buffer.data(gh0 + 2);
    const auto *gh0_3 = buffer.data(gh0 + 3);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);
    const auto *gh0_63 = buffer.data(gh0 + 63);
    const auto *gh0_65 = buffer.data(gh0 + 65);
    const auto *gh0_66 = buffer.data(gh0 + 66);
    const auto *gh0_68 = buffer.data(gh0 + 68);
    const auto *gh0_69 = buffer.data(gh0 + 69);
    const auto *gh0_70 = buffer.data(gh0 + 70);
    const auto *gh0_72 = buffer.data(gh0 + 72);
    const auto *gh0_73 = buffer.data(gh0 + 73);
    const auto *gh0_78 = buffer.data(gh0 + 78);
    const auto *gh0_79 = buffer.data(gh0 + 79);
    const auto *gh0_80 = buffer.data(gh0 + 80);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_1 = buffer.data(gh1 + 1);
    const auto *gh1_2 = buffer.data(gh1 + 2);
    const auto *gh1_3 = buffer.data(gh1 + 3);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);
    const auto *gh1_63 = buffer.data(gh1 + 63);
    const auto *gh1_65 = buffer.data(gh1 + 65);
    const auto *gh1_66 = buffer.data(gh1 + 66);
    const auto *gh1_68 = buffer.data(gh1 + 68);
    const auto *gh1_69 = buffer.data(gh1 + 69);
    const auto *gh1_70 = buffer.data(gh1 + 70);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_73 = buffer.data(gh1 + 73);
    const auto *gh1_78 = buffer.data(gh1 + 78);
    const auto *gh1_79 = buffer.data(gh1 + 79);
    const auto *gh1_80 = buffer.data(gh1 + 80);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fi_0, gh0_0, gh1_0, \
                         gi_0, gi_1, gi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fi_0[k]
                 + f_1 * gh0_0[k]
                 - f_2 * gh1_0[k]
                 + pb_x[k] * gi_0[k];

        t_1[k] = pb_y[k] * gi_0[k];

        t_2[k] = pb_z[k] * gi_0[k];

        t_3[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_y[k] * gi_1[k];

        t_4[k] = pb_y[k] * gi_2[k];

        t_5[k] = f_3 * gh0_0[k]
                 - f_4 * gh1_0[k]
                 + pb_z[k] * gi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, gh0_1, gh0_2, gh0_3, gh1_1, \
                         gh1_2, gh1_3, gi_3, gi_5, gi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * gh0_1[k]
                 - f_6 * gh1_1[k]
                 + pb_y[k] * gi_3[k];

        t_7[k] = pb_z[k] * gi_3[k];

        t_8[k] = pb_y[k] * gi_5[k];

        t_9[k] = f_5 * gh0_2[k]
                 - f_6 * gh1_2[k]
                 + pb_z[k] * gi_5[k];

        t_10[k] = f_7 * gh0_3[k]
                  - f_8 * gh1_3[k]
                  + pb_y[k] * gi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, gh0_5, gh0_6, gh1_5, \
                         gh1_6, gi_6, gi_8, gi_9, gi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * gi_6[k];

        t_12[k] = f_3 * gh0_5[k]
                  - f_4 * gh1_5[k]
                  + pb_y[k] * gi_8[k];

        t_13[k] = pb_y[k] * gi_9[k];

        t_14[k] = f_7 * gh0_5[k]
                  - f_8 * gh1_5[k]
                  + pb_z[k] * gi_9[k];

        t_15[k] = f_9 * gh0_6[k]
                  - f_10 * gh1_6[k]
                  + pb_y[k] * gi_10[k];

        t_16[k] = pb_z[k] * gi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, gh0_8, gh0_9, gh1_8, gh1_9, \
                         gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * gh0_8[k]
                  - f_6 * gh1_8[k]
                  + pb_y[k] * gi_12[k];

        t_18[k] = f_3 * gh0_9[k]
                  - f_4 * gh1_9[k]
                  + pb_y[k] * gi_13[k];

        t_19[k] = pb_y[k] * gi_14[k];

        t_20[k] = f_9 * gh0_9[k]
                  - f_10 * gh1_9[k]
                  + pb_z[k] * gi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, fi_21, fi_23, fi_24, fi_25, \
                         gi_15, gi_21, gi_23, gi_24, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * fi_21[k]
                  + pb_x[k] * gi_21[k];

        t_22[k] = pb_z[k] * gi_15[k];

        t_23[k] = f_0 * fi_23[k]
                  + pb_x[k] * gi_23[k];

        t_24[k] = f_0 * fi_24[k]
                  + pb_x[k] * gi_24[k];

        t_25[k] = f_0 * fi_25[k]
                  + pb_x[k] * gi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, fi_27, gh0_15, gh1_15, \
                         gi_20, gi_21, gi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * gi_20[k];

        t_27[k] = f_0 * fi_27[k]
                  + pb_x[k] * gi_27[k];

        t_28[k] = f_1 * gh0_15[k]
                  - f_2 * gh1_15[k]
                  + pb_y[k] * gi_21[k];

        t_29[k] = pb_z[k] * gi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, gh0_17, gh0_18, gh0_19, gh1_17, gh1_18, \
                         gh1_19, gi_23, gi_24, gi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * gh0_17[k]
                  - f_10 * gh1_17[k]
                  + pb_y[k] * gi_23[k];

        t_31[k] = f_7 * gh0_18[k]
                  - f_8 * gh1_18[k]
                  + pb_y[k] * gi_24[k];

        t_32[k] = f_5 * gh0_19[k]
                  - f_6 * gh1_19[k]
                  + pb_y[k] * gi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, fi_0, fk_0, \
                         gh0_20, gh1_20, gi_26, gi_27, gi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gh0_20[k]
                  - f_4 * gh1_20[k]
                  + pb_y[k] * gi_26[k];

        t_34[k] = pb_y[k] * gi_27[k];

        t_35[k] = f_1 * gh0_20[k]
                  - f_2 * gh1_20[k]
                  + pb_z[k] * gi_27[k];

        t_36[k] = pa_y[k] * fk_0[k];

        t_37[k] = f_11 * fi_0[k]
                  + pb_y[k] * gi_28[k];

        t_38[k] = pb_z[k] * gi_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, fi_1, fi_3, fk_3, fk_5, \
                         fk_6, gi_29, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * fi_1[k]
                  + pa_y[k] * fk_3[k];

        t_40[k] = pb_z[k] * gi_29[k];

        t_41[k] = pa_y[k] * fk_5[k];

        t_42[k] = f_13 * fi_3[k]
                  + pa_y[k] * fk_6[k];

        t_43[k] = pb_z[k] * gi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, fi_5, fi_6, fi_8, \
                         fk_9, fk_10, fk_12, gi_33, gi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * fi_5[k]
                  + pb_y[k] * gi_33[k];

        t_45[k] = pa_y[k] * fk_9[k];

        t_46[k] = f_0 * fi_6[k]
                  + pa_y[k] * fk_10[k];

        t_47[k] = pb_z[k] * gi_34[k];

        t_48[k] = f_12 * fi_8[k]
                  + pa_y[k] * fk_12[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, fi_9, fi_10, fi_12, \
                         fk_14, fk_15, fk_17, gi_37, gi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * fi_9[k]
                  + pb_y[k] * gi_37[k];

        t_50[k] = pa_y[k] * fk_14[k];

        t_51[k] = f_14 * fi_10[k]
                  + pa_y[k] * fk_15[k];

        t_52[k] = pb_z[k] * gi_38[k];

        t_53[k] = f_13 * fi_12[k]
                  + pa_y[k] * fk_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, fi_13, fi_14, fi_49, fk_18, \
                         fk_20, gi_42, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * fi_13[k]
                  + pa_y[k] * fk_18[k];

        t_55[k] = f_11 * fi_14[k]
                  + pb_y[k] * gi_42[k];

        t_56[k] = pa_y[k] * fk_20[k];

        t_57[k] = f_13 * fi_49[k]
                  + pb_x[k] * gi_49[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, fi_51, fi_52, fi_53, fi_54, \
                         gi_43, gi_51, gi_52, gi_53, gi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * gi_43[k];

        t_59[k] = f_13 * fi_51[k]
                  + pb_x[k] * gi_51[k];

        t_60[k] = f_13 * fi_52[k]
                  + pb_x[k] * gi_52[k];

        t_61[k] = f_13 * fi_53[k]
                  + pb_x[k] * gi_53[k];

        t_62[k] = f_13 * fi_54[k]
                  + pb_x[k] * gi_54[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, fi_21, fi_23, fi_24, fk_27, \
                         fk_28, fk_30, fk_31, gi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * fk_27[k];

        t_64[k] = f_15 * fi_21[k]
                  + pa_y[k] * fk_28[k];

        t_65[k] = pb_z[k] * gi_49[k];

        t_66[k] = f_14 * fi_23[k]
                  + pa_y[k] * fk_30[k];

        t_67[k] = f_0 * fi_24[k]
                  + pa_y[k] * fk_31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, fi_25, fi_26, fi_27, \
                         fk_0, fk_32, fk_33, fk_35, gi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_13 * fi_25[k]
                  + pa_y[k] * fk_32[k];

        t_69[k] = f_12 * fi_26[k]
                  + pa_y[k] * fk_33[k];

        t_70[k] = f_11 * fi_27[k]
                  + pb_y[k] * gi_55[k];

        t_71[k] = pa_y[k] * fk_35[k];

        t_72[k] = pa_z[k] * fk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, fi_0, fi_2, \
                         fk_3, fk_5, fk_6, gi_56, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * gi_56[k];

        t_74[k] = f_11 * fi_0[k]
                  + pb_z[k] * gi_56[k];

        t_75[k] = pa_z[k] * fk_3[k];

        t_76[k] = pb_y[k] * gi_58[k];

        t_77[k] = f_12 * fi_2[k]
                  + pa_z[k] * fk_5[k];

        t_78[k] = pa_z[k] * fk_6[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, fi_3, fi_5, fi_6, \
                         fk_9, fk_10, gi_59, gi_61, gi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * fi_3[k]
                  + pb_z[k] * gi_59[k];

        t_80[k] = pb_y[k] * gi_61[k];

        t_81[k] = f_13 * fi_5[k]
                  + pa_z[k] * fk_9[k];

        t_82[k] = pa_z[k] * fk_10[k];

        t_83[k] = f_11 * fi_6[k]
                  + pb_z[k] * gi_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, fi_7, fi_9, fi_10, \
                         fk_12, fk_14, fk_15, gi_65, gi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * fi_7[k]
                  + pa_z[k] * fk_12[k];

        t_85[k] = pb_y[k] * gi_65[k];

        t_86[k] = f_0 * fi_9[k]
                  + pa_z[k] * fk_14[k];

        t_87[k] = pa_z[k] * fk_15[k];

        t_88[k] = f_11 * fi_10[k]
                  + pb_z[k] * gi_66[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, fi_11, fi_12, fi_14, fk_17, \
                         fk_18, fk_20, fk_21, gi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * fi_11[k]
                  + pa_z[k] * fk_17[k];

        t_90[k] = f_13 * fi_12[k]
                  + pa_z[k] * fk_18[k];

        t_91[k] = pb_y[k] * gi_70[k];

        t_92[k] = f_14 * fi_14[k]
                  + pa_z[k] * fk_20[k];

        t_93[k] = pa_z[k] * fk_21[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, fi_78, fi_79, fi_80, fi_81, \
                         gi_76, gi_78, gi_79, gi_80, gi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_13 * fi_78[k]
                  + pb_x[k] * gi_78[k];

        t_95[k] = f_13 * fi_79[k]
                  + pb_x[k] * gi_79[k];

        t_96[k] = f_13 * fi_80[k]
                  + pb_x[k] * gi_80[k];

        t_97[k] = f_13 * fi_81[k]
                  + pb_x[k] * gi_81[k];

        t_98[k] = pb_y[k] * gi_76[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, fi_21, fi_22, fi_83, \
                         fk_28, fk_30, gi_77, gi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * fi_83[k]
                  + pb_x[k] * gi_83[k];

        t_100[k] = pa_z[k] * fk_28[k];

        t_101[k] = f_11 * fi_21[k]
                   + pb_z[k] * gi_77[k];

        t_102[k] = f_12 * fi_22[k]
                   + pa_z[k] * fk_30[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, fi_23, fi_24, fi_25, \
                         fi_27, fk_31, fk_32, fk_33, fk_35, gi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_13 * fi_23[k]
                   + pa_z[k] * fk_31[k];

        t_104[k] = f_0 * fi_24[k]
                   + pa_z[k] * fk_32[k];

        t_105[k] = f_14 * fi_25[k]
                   + pa_z[k] * fk_33[k];

        t_106[k] = pb_y[k] * gi_83[k];

        t_107[k] = f_15 * fi_27[k]
                   + pa_z[k] * fk_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_y, pb_z, dk0_0, dk1_0, fi_28, fk_36, \
                         gi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_16 * dk0_0[k]
                   - f_17 * dk1_0[k]
                   + pa_y[k] * fk_36[k];

        t_109[k] = f_12 * fi_28[k]
                   + pb_y[k] * gi_84[k];

        t_110[k] = pb_z[k] * gi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pb_z, fi_87, gh0_63, gh0_66, gh1_63, \
                         gh1_66, gi_85, gi_86, gi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * fi_87[k]
                   + f_9 * gh0_66[k]
                   - f_10 * gh1_66[k]
                   + pb_x[k] * gi_87[k];

        t_112[k] = pb_z[k] * gi_85[k];

        t_113[k] = f_3 * gh0_63[k]
                   - f_4 * gh1_63[k]
                   + pb_z[k] * gi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pb_y, pb_z, fi_33, fi_90, gh0_65, \
                         gh0_69, gh1_65, gh1_69, gi_87, gi_89, gi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_12 * fi_90[k]
                   + f_7 * gh0_69[k]
                   - f_8 * gh1_69[k]
                   + pb_x[k] * gi_90[k];

        t_115[k] = pb_z[k] * gi_87[k];

        t_116[k] = f_12 * fi_33[k]
                   + pb_y[k] * gi_89[k];

        t_117[k] = f_5 * gh0_65[k]
                   - f_6 * gh1_65[k]
                   + pb_z[k] * gi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_z, fi_94, gh0_66, gh0_73, gh1_66, \
                         gh1_73, gi_90, gi_91, gi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * fi_94[k]
                   + f_5 * gh0_73[k]
                   - f_6 * gh1_73[k]
                   + pb_x[k] * gi_94[k];

        t_119[k] = pb_z[k] * gi_90[k];

        t_120[k] = f_3 * gh0_66[k]
                   - f_4 * gh1_66[k]
                   + pb_z[k] * gi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, fi_37, fi_99, gh0_68, \
                         gh0_78, gh1_68, gh1_78, gi_93, gi_94, gi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * fi_37[k]
                   + pb_y[k] * gi_93[k];

        t_122[k] = f_7 * gh0_68[k]
                   - f_8 * gh1_68[k]
                   + pb_z[k] * gi_93[k];

        t_123[k] = f_12 * fi_99[k]
                   + f_3 * gh0_78[k]
                   - f_4 * gh1_78[k]
                   + pb_x[k] * gi_99[k];

        t_124[k] = pb_z[k] * gi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, fi_42, gh0_69, gh0_70, \
                         gh0_72, gh1_69, gh1_70, gh1_72, gi_95, gi_96, \
                         gi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * gh0_69[k]
                   - f_4 * gh1_69[k]
                   + pb_z[k] * gi_95[k];

        t_126[k] = f_5 * gh0_70[k]
                   - f_6 * gh1_70[k]
                   + pb_z[k] * gi_96[k];

        t_127[k] = f_12 * fi_42[k]
                   + pb_y[k] * gi_98[k];

        t_128[k] = f_9 * gh0_72[k]
                   - f_10 * gh1_72[k]
                   + pb_z[k] * gi_98[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, fi_105, fi_107, \
                         fi_108, fi_109, gi_99, gi_105, gi_107, gi_108, \
                         gi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * fi_105[k]
                   + pb_x[k] * gi_105[k];

        t_130[k] = pb_z[k] * gi_99[k];

        t_131[k] = f_12 * fi_107[k]
                   + pb_x[k] * gi_107[k];

        t_132[k] = f_12 * fi_108[k]
                   + pb_x[k] * gi_108[k];

        t_133[k] = f_12 * fi_109[k]
                   + pb_x[k] * gi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_z, dk0_136, dk1_136, \
                         fi_110, fi_111, fk_136, gi_105, gi_110, \
                         gi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_12 * fi_110[k]
                   + pb_x[k] * gi_110[k];

        t_135[k] = f_12 * fi_111[k]
                   + pb_x[k] * gi_111[k];

        t_136[k] = f_16 * dk0_136[k]
                   - f_17 * dk1_136[k]
                   + pa_x[k] * fk_136[k];

        t_137[k] = pb_z[k] * gi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_z, gh0_78, gh0_79, gh0_80, gh1_78, gh1_79, \
                         gh1_80, gi_106, gi_107, gi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * gh0_78[k]
                   - f_4 * gh1_78[k]
                   + pb_z[k] * gi_106[k];

        t_139[k] = f_5 * gh0_79[k]
                   - f_6 * gh1_79[k]
                   + pb_z[k] * gi_107[k];

        t_140[k] = f_7 * gh0_80[k]
                   - f_8 * gh1_80[k]
                   + pb_z[k] * gi_108[k];
    }
}

static auto
compute_prim_gk_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t dk0,
                                            const size_t dk1, const size_t fi, const size_t fk,
                                            const size_t gh0, const size_t gh1, const size_t gi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_0 = buffer.data(dk0 + 0);
    const auto *dk0_215 = buffer.data(dk0 + 215);

    const auto *dk1_0 = buffer.data(dk1 + 0);
    const auto *dk1_215 = buffer.data(dk1 + 215);

    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);

    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_64 = buffer.data(fk + 64);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_102 = buffer.data(fk + 102);
    const auto *fk_103 = buffer.data(fk + 103);
    const auto *fk_104 = buffer.data(fk + 104);
    const auto *fk_105 = buffer.data(fk + 105);
    const auto *fk_107 = buffer.data(fk + 107);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_215 = buffer.data(fk + 215);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);

    const auto *gh0_81 = buffer.data(gh0 + 81);
    const auto *gh0_83 = buffer.data(gh0 + 83);
    const auto *gh0_105 = buffer.data(gh0 + 105);
    const auto *gh0_106 = buffer.data(gh0 + 106);
    const auto *gh0_108 = buffer.data(gh0 + 108);
    const auto *gh0_110 = buffer.data(gh0 + 110);
    const auto *gh0_111 = buffer.data(gh0 + 111);
    const auto *gh0_113 = buffer.data(gh0 + 113);
    const auto *gh0_114 = buffer.data(gh0 + 114);
    const auto *gh0_119 = buffer.data(gh0 + 119);
    const auto *gh0_120 = buffer.data(gh0 + 120);
    const auto *gh0_122 = buffer.data(gh0 + 122);
    const auto *gh0_123 = buffer.data(gh0 + 123);
    const auto *gh0_124 = buffer.data(gh0 + 124);
    const auto *gh0_125 = buffer.data(gh0 + 125);

    const auto *gh1_81 = buffer.data(gh1 + 81);
    const auto *gh1_83 = buffer.data(gh1 + 83);
    const auto *gh1_105 = buffer.data(gh1 + 105);
    const auto *gh1_106 = buffer.data(gh1 + 106);
    const auto *gh1_108 = buffer.data(gh1 + 108);
    const auto *gh1_110 = buffer.data(gh1 + 110);
    const auto *gh1_111 = buffer.data(gh1 + 111);
    const auto *gh1_113 = buffer.data(gh1 + 113);
    const auto *gh1_114 = buffer.data(gh1 + 114);
    const auto *gh1_119 = buffer.data(gh1 + 119);
    const auto *gh1_120 = buffer.data(gh1 + 120);
    const auto *gh1_122 = buffer.data(gh1 + 122);
    const auto *gh1_123 = buffer.data(gh1 + 123);
    const auto *gh1_124 = buffer.data(gh1 + 124);
    const auto *gh1_125 = buffer.data(gh1 + 125);

    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, fi_55, fk_72, gh0_81, \
                         gh0_83, gh1_81, gh1_83, gi_109, gi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * gh0_81[k]
                   - f_10 * gh1_81[k]
                   + pb_z[k] * gi_109[k];

        t_142[k] = f_12 * fi_55[k]
                   + pb_y[k] * gi_111[k];

        t_143[k] = f_1 * gh0_83[k]
                   - f_2 * gh1_83[k]
                   + pb_z[k] * gi_111[k];

        t_144[k] = pa_y[k] * fk_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, fi_58, \
                         fk_37, fk_39, fk_42, fk_74, fk_77, gi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * fk_37[k];

        t_146[k] = pa_y[k] * fk_74[k];

        t_147[k] = pa_z[k] * fk_39[k];

        t_148[k] = f_11 * fi_58[k]
                   + pb_y[k] * gi_114[k];

        t_149[k] = pa_y[k] * fk_77[k];

        t_150[k] = pa_z[k] * fk_42[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, fi_31, fi_61, \
                         fk_46, fk_81, gi_115, gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * fi_31[k]
                   + pb_z[k] * gi_115[k];

        t_152[k] = f_11 * fi_61[k]
                   + pb_y[k] * gi_117[k];

        t_153[k] = pa_y[k] * fk_81[k];

        t_154[k] = pa_z[k] * fk_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, fi_34, fi_64, fi_65, \
                         fk_84, fk_86, gi_118, gi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * fi_34[k]
                   + pb_z[k] * gi_118[k];

        t_156[k] = f_12 * fi_64[k]
                   + pa_y[k] * fk_84[k];

        t_157[k] = f_11 * fi_65[k]
                   + pb_y[k] * gi_121[k];

        t_158[k] = pa_y[k] * fk_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, fi_38, fi_68, fi_69, \
                         fk_51, fk_89, fk_90, gi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * fk_51[k];

        t_160[k] = f_11 * fi_38[k]
                   + pb_z[k] * gi_122[k];

        t_161[k] = f_13 * fi_68[k]
                   + pa_y[k] * fk_89[k];

        t_162[k] = f_12 * fi_69[k]
                   + pa_y[k] * fk_90[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, fi_70, fi_134, \
                         fk_57, fk_92, gi_126, gi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * fi_70[k]
                   + pb_y[k] * gi_126[k];

        t_164[k] = pa_y[k] * fk_92[k];

        t_165[k] = pa_z[k] * fk_57[k];

        t_166[k] = f_12 * fi_134[k]
                   + pb_x[k] * gi_134[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, fi_135, fi_136, \
                         fi_137, fi_138, fk_99, gi_135, gi_136, gi_137, \
                         gi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_12 * fi_135[k]
                   + pb_x[k] * gi_135[k];

        t_168[k] = f_12 * fi_136[k]
                   + pb_x[k] * gi_136[k];

        t_169[k] = f_12 * fi_137[k]
                   + pb_x[k] * gi_137[k];

        t_170[k] = f_12 * fi_138[k]
                   + pb_x[k] * gi_138[k];

        t_171[k] = pa_y[k] * fk_99[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pa_z, pb_z, fi_49, fi_79, fi_80, \
                         fk_64, fk_102, fk_103, gi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_z[k] * fk_64[k];

        t_173[k] = f_11 * fi_49[k]
                   + pb_z[k] * gi_133[k];

        t_174[k] = f_14 * fi_79[k]
                   + pa_y[k] * fk_102[k];

        t_175[k] = f_0 * fi_80[k]
                   + pa_y[k] * fk_103[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, fi_81, fi_82, fi_83, fk_104, \
                         fk_105, fk_107, gi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * fi_81[k]
                   + pa_y[k] * fk_104[k];

        t_177[k] = f_12 * fi_82[k]
                   + pa_y[k] * fk_105[k];

        t_178[k] = f_11 * fi_83[k]
                   + pb_y[k] * gi_139[k];

        t_179[k] = pa_y[k] * fk_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_y, pb_z, dk0_0, dk1_0, fi_56, \
                         fk_72, gh0_105, gh1_105, gi_140, gi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_16 * dk0_0[k]
                   - f_17 * dk1_0[k]
                   + pa_z[k] * fk_72[k];

        t_181[k] = pb_y[k] * gi_140[k];

        t_182[k] = f_12 * fi_56[k]
                   + pb_z[k] * gi_140[k];

        t_183[k] = f_3 * gh0_105[k]
                   - f_4 * gh1_105[k]
                   + pb_y[k] * gi_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pb_y, pb_z, fi_59, fi_145, gh0_106, \
                         gh0_110, gh1_106, gh1_110, gi_142, gi_143, \
                         gi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * gi_142[k];

        t_185[k] = f_12 * fi_145[k]
                   + f_9 * gh0_110[k]
                   - f_10 * gh1_110[k]
                   + pb_x[k] * gi_145[k];

        t_186[k] = f_5 * gh0_106[k]
                   - f_6 * gh1_106[k]
                   + pb_y[k] * gi_143[k];

        t_187[k] = f_12 * fi_59[k]
                   + pb_z[k] * gi_143[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pb_x, pb_y, pb_z, fi_62, fi_149, gh0_108, \
                         gh0_114, gh1_108, gh1_114, gi_145, gi_146, \
                         gi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_y[k] * gi_145[k];

        t_189[k] = f_12 * fi_149[k]
                   + f_7 * gh0_114[k]
                   - f_8 * gh1_114[k]
                   + pb_x[k] * gi_149[k];

        t_190[k] = f_7 * gh0_108[k]
                   - f_8 * gh1_108[k]
                   + pb_y[k] * gi_146[k];

        t_191[k] = f_12 * fi_62[k]
                   + pb_z[k] * gi_146[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, fi_154, gh0_110, gh0_119, gh1_110, \
                         gh1_119, gi_148, gi_149, gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * gh0_110[k]
                   - f_4 * gh1_110[k]
                   + pb_y[k] * gi_148[k];

        t_193[k] = pb_y[k] * gi_149[k];

        t_194[k] = f_12 * fi_154[k]
                   + f_5 * gh0_119[k]
                   - f_6 * gh1_119[k]
                   + pb_x[k] * gi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, fi_66, gh0_111, gh0_113, \
                         gh0_114, gh1_111, gh1_113, gh1_114, gi_150, gi_152, \
                         gi_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_9 * gh0_111[k]
                   - f_10 * gh1_111[k]
                   + pb_y[k] * gi_150[k];

        t_196[k] = f_12 * fi_66[k]
                   + pb_z[k] * gi_150[k];

        t_197[k] = f_5 * gh0_113[k]
                   - f_6 * gh1_113[k]
                   + pb_y[k] * gi_152[k];

        t_198[k] = f_3 * gh0_114[k]
                   - f_4 * gh1_114[k]
                   + pb_y[k] * gi_153[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_x, pb_y, fi_160, fi_161, fi_162, \
                         gh0_125, gh1_125, gi_154, gi_160, gi_161, \
                         gi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * gi_154[k];

        t_200[k] = f_12 * fi_160[k]
                   + f_3 * gh0_125[k]
                   - f_4 * gh1_125[k]
                   + pb_x[k] * gi_160[k];

        t_201[k] = f_12 * fi_161[k]
                   + pb_x[k] * gi_161[k];

        t_202[k] = f_12 * fi_162[k]
                   + pb_x[k] * gi_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, pb_y, fi_163, fi_164, \
                         fi_165, fi_167, gi_160, gi_163, gi_164, gi_165, \
                         gi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_12 * fi_163[k]
                   + pb_x[k] * gi_163[k];

        t_204[k] = f_12 * fi_164[k]
                   + pb_x[k] * gi_164[k];

        t_205[k] = f_12 * fi_165[k]
                   + pb_x[k] * gi_165[k];

        t_206[k] = pb_y[k] * gi_160[k];

        t_207[k] = f_12 * fi_167[k]
                   + pb_x[k] * gi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pb_z, fi_77, gh0_120, gh0_122, \
                         gh0_123, gh1_120, gh1_122, gh1_123, gi_161, gi_163, \
                         gi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * gh0_120[k]
                   - f_2 * gh1_120[k]
                   + pb_y[k] * gi_161[k];

        t_209[k] = f_12 * fi_77[k]
                   + pb_z[k] * gi_161[k];

        t_210[k] = f_9 * gh0_122[k]
                   - f_10 * gh1_122[k]
                   + pb_y[k] * gi_163[k];

        t_211[k] = f_7 * gh0_123[k]
                   - f_8 * gh1_123[k]
                   + pb_y[k] * gi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pb_y, dk0_215, dk1_215, fk_215, \
                         gh0_124, gh0_125, gh1_124, gh1_125, gi_165, gi_166, \
                         gi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * gh0_124[k]
                   - f_6 * gh1_124[k]
                   + pb_y[k] * gi_165[k];

        t_213[k] = f_3 * gh0_125[k]
                   - f_4 * gh1_125[k]
                   + pb_y[k] * gi_166[k];

        t_214[k] = pb_y[k] * gi_167[k];

        t_215[k] = f_16 * dk0_215[k]
                   - f_17 * dk1_215[k]
                   + pa_x[k] * fk_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pa_x, pb_y, pb_z, fi_84, fi_168, \
                         fi_171, fk_216, fk_219, gi_168, gi_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_15 * fi_168[k]
                   + pa_x[k] * fk_216[k];

        t_217[k] = f_13 * fi_84[k]
                   + pb_y[k] * gi_168[k];

        t_218[k] = pb_z[k] * gi_168[k];

        t_219[k] = f_14 * fi_171[k]
                   + pa_x[k] * fk_219[k];

        t_220[k] = pb_z[k] * gi_169[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_x, pb_y, pb_z, fi_89, fi_173, fi_174, \
                         fk_221, fk_222, gi_171, gi_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_14 * fi_173[k]
                   + pa_x[k] * fk_221[k];

        t_222[k] = f_0 * fi_174[k]
                   + pa_x[k] * fk_222[k];

        t_223[k] = pb_z[k] * gi_171[k];

        t_224[k] = f_13 * fi_89[k]
                   + pb_y[k] * gi_173[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pb_z, fi_177, fi_178, fi_180, \
                         fk_225, fk_226, fk_228, gi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * fi_177[k]
                   + pa_x[k] * fk_225[k];

        t_226[k] = f_13 * fi_178[k]
                   + pa_x[k] * fk_226[k];

        t_227[k] = pb_z[k] * gi_174[k];

        t_228[k] = f_13 * fi_180[k]
                   + pa_x[k] * fk_228[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_x, pb_y, pb_z, fi_93, fi_182, fi_183, \
                         fk_230, fk_231, gi_177, gi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * fi_93[k]
                   + pb_y[k] * gi_177[k];

        t_230[k] = f_13 * fi_182[k]
                   + pa_x[k] * fk_230[k];

        t_231[k] = f_12 * fi_183[k]
                   + pa_x[k] * fk_231[k];

        t_232[k] = pb_z[k] * gi_178[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pb_y, fi_98, fi_185, fi_186, \
                         fi_188, fk_233, fk_234, fk_236, gi_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_12 * fi_185[k]
                   + pa_x[k] * fk_233[k];

        t_234[k] = f_12 * fi_186[k]
                   + pa_x[k] * fk_234[k];

        t_235[k] = f_13 * fi_98[k]
                   + pb_y[k] * gi_182[k];

        t_236[k] = f_12 * fi_188[k]
                   + pa_x[k] * fk_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_x, pb_z, fi_189, fi_191, \
                         fi_192, fi_193, gi_183, gi_189, gi_191, gi_192, \
                         gi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_11 * fi_189[k]
                   + pb_x[k] * gi_189[k];

        t_238[k] = pb_z[k] * gi_183[k];

        t_239[k] = f_11 * fi_191[k]
                   + pb_x[k] * gi_191[k];

        t_240[k] = f_11 * fi_192[k]
                   + pb_x[k] * gi_192[k];

        t_241[k] = f_11 * fi_193[k]
                   + pb_x[k] * gi_193[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pa_x, pb_x, pb_z, fi_194, fi_195, \
                         fk_244, fk_246, gi_189, gi_194, gi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_11 * fi_194[k]
                   + pb_x[k] * gi_194[k];

        t_243[k] = f_11 * fi_195[k]
                   + pb_x[k] * gi_195[k];

        t_244[k] = pa_x[k] * fk_244[k];

        t_245[k] = pb_z[k] * gi_189[k];

        t_246[k] = pa_x[k] * fk_246[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, t_252, t_253, pa_x, pa_z, fk_108, \
                         fk_109, fk_247, fk_248, fk_249, fk_250, \
                         fk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = pa_x[k] * fk_247[k];

        t_248[k] = pa_x[k] * fk_248[k];

        t_249[k] = pa_x[k] * fk_249[k];

        t_250[k] = pa_x[k] * fk_250[k];

        t_251[k] = pa_x[k] * fk_251[k];

        t_252[k] = pa_z[k] * fk_108[k];

        t_253[k] = pa_z[k] * fk_109[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pa_x, pa_z, pb_y, pb_z, fi_84, fi_114, \
                         fi_201, fk_111, fk_257, gi_196, gi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_11 * fi_84[k]
                   + pb_z[k] * gi_196[k];

        t_255[k] = pa_z[k] * fk_111[k];

        t_256[k] = f_12 * fi_114[k]
                   + pb_y[k] * gi_198[k];

        t_257[k] = f_14 * fi_201[k]
                   + pa_x[k] * fk_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pa_x, pa_z, pb_y, pb_z, fi_87, fi_117, \
                         fi_205, fk_114, fk_261, gi_199, gi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_z[k] * fk_114[k];

        t_259[k] = f_11 * fi_87[k]
                   + pb_z[k] * gi_199[k];

        t_260[k] = f_12 * fi_117[k]
                   + pb_y[k] * gi_201[k];

        t_261[k] = f_0 * fi_205[k]
                   + pa_x[k] * fk_261[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, pa_x, pa_z, pb_y, pb_z, fi_90, fi_121, \
                         fi_208, fk_118, fk_264, gi_202, gi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = pa_z[k] * fk_118[k];

        t_263[k] = f_11 * fi_90[k]
                   + pb_z[k] * gi_202[k];

        t_264[k] = f_13 * fi_208[k]
                   + pa_x[k] * fk_264[k];

        t_265[k] = f_12 * fi_121[k]
                   + pb_y[k] * gi_205[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_x, pa_z, pb_z, fi_94, fi_210, fi_213, \
                         fk_123, fk_266, fk_269, gi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_13 * fi_210[k]
                   + pa_x[k] * fk_266[k];

        t_267[k] = pa_z[k] * fk_123[k];

        t_268[k] = f_11 * fi_94[k]
                   + pb_z[k] * gi_206[k];

        t_269[k] = f_12 * fi_213[k]
                   + pa_x[k] * fk_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pa_x, pa_z, pb_y, fi_126, fi_214, fi_216, \
                         fk_129, fk_270, fk_272, gi_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_12 * fi_214[k]
                   + pa_x[k] * fk_270[k];

        t_271[k] = f_12 * fi_126[k]
                   + pb_y[k] * gi_210[k];

        t_272[k] = f_12 * fi_216[k]
                   + pa_x[k] * fk_272[k];

        t_273[k] = pa_z[k] * fk_129[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_x, fi_218, fi_219, fi_220, \
                         fi_221, fi_222, gi_218, gi_219, gi_220, gi_221, \
                         gi_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_11 * fi_218[k]
                   + pb_x[k] * gi_218[k];

        t_275[k] = f_11 * fi_219[k]
                   + pb_x[k] * gi_219[k];

        t_276[k] = f_11 * fi_220[k]
                   + pb_x[k] * gi_220[k];

        t_277[k] = f_11 * fi_221[k]
                   + pb_x[k] * gi_221[k];

        t_278[k] = f_11 * fi_222[k]
                   + pb_x[k] * gi_222[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, t_283, t_284, pa_x, pb_x, fi_223, fk_280, \
                         fk_281, fk_282, fk_283, fk_284, gi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_11 * fi_223[k]
                   + pb_x[k] * gi_223[k];

        t_280[k] = pa_x[k] * fk_280[k];

        t_281[k] = pa_x[k] * fk_281[k];

        t_282[k] = pa_x[k] * fk_282[k];

        t_283[k] = pa_x[k] * fk_283[k];

        t_284[k] = pa_x[k] * fk_284[k];
    }
}

static auto
compute_prim_gk_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t fi,
                                            const size_t fk, const size_t gh0, const size_t gh1,
                                            const size_t gi, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_279 = buffer.data(fi + 279);

    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_251 = buffer.data(fk + 251);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_287 = buffer.data(fk + 287);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_359 = buffer.data(fk + 359);

    const auto *gh0_210 = buffer.data(gh0 + 210);
    const auto *gh0_213 = buffer.data(gh0 + 213);
    const auto *gh0_215 = buffer.data(gh0 + 215);
    const auto *gh0_216 = buffer.data(gh0 + 216);
    const auto *gh0_219 = buffer.data(gh0 + 219);
    const auto *gh0_220 = buffer.data(gh0 + 220);
    const auto *gh0_222 = buffer.data(gh0 + 222);
    const auto *gh0_224 = buffer.data(gh0 + 224);
    const auto *gh0_225 = buffer.data(gh0 + 225);
    const auto *gh0_226 = buffer.data(gh0 + 226);
    const auto *gh0_227 = buffer.data(gh0 + 227);
    const auto *gh0_228 = buffer.data(gh0 + 228);
    const auto *gh0_230 = buffer.data(gh0 + 230);
    const auto *gh0_252 = buffer.data(gh0 + 252);

    const auto *gh1_210 = buffer.data(gh1 + 210);
    const auto *gh1_213 = buffer.data(gh1 + 213);
    const auto *gh1_215 = buffer.data(gh1 + 215);
    const auto *gh1_216 = buffer.data(gh1 + 216);
    const auto *gh1_219 = buffer.data(gh1 + 219);
    const auto *gh1_220 = buffer.data(gh1 + 220);
    const auto *gh1_222 = buffer.data(gh1 + 222);
    const auto *gh1_224 = buffer.data(gh1 + 224);
    const auto *gh1_225 = buffer.data(gh1 + 225);
    const auto *gh1_226 = buffer.data(gh1 + 226);
    const auto *gh1_227 = buffer.data(gh1 + 227);
    const auto *gh1_228 = buffer.data(gh1 + 228);
    const auto *gh1_230 = buffer.data(gh1 + 230);
    const auto *gh1_252 = buffer.data(gh1 + 252);

    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, t_290, pa_x, pa_y, pb_y, fi_140, \
                         fk_180, fk_182, fk_285, fk_286, fk_287, \
                         gi_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = pa_x[k] * fk_285[k];

        t_286[k] = pa_x[k] * fk_286[k];

        t_287[k] = pa_x[k] * fk_287[k];

        t_288[k] = pa_y[k] * fk_180[k];

        t_289[k] = f_11 * fi_140[k]
                   + pb_y[k] * gi_224[k];

        t_290[k] = pa_y[k] * fk_182[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pa_x, pa_y, pb_y, fi_142, fi_227, fi_230, \
                         fk_185, fk_291, fk_294, gi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_14 * fi_227[k]
                   + pa_x[k] * fk_291[k];

        t_292[k] = f_11 * fi_142[k]
                   + pb_y[k] * gi_226[k];

        t_293[k] = pa_y[k] * fk_185[k];

        t_294[k] = f_0 * fi_230[k]
                   + pa_x[k] * fk_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pa_x, pa_y, pb_y, pb_z, fi_115, fi_145, \
                         fi_234, fk_189, fk_298, gi_227, gi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_12 * fi_115[k]
                   + pb_z[k] * gi_227[k];

        t_296[k] = f_11 * fi_145[k]
                   + pb_y[k] * gi_229[k];

        t_297[k] = pa_y[k] * fk_189[k];

        t_298[k] = f_13 * fi_234[k]
                   + pa_x[k] * fk_298[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pa_x, pa_y, pb_y, pb_z, fi_118, fi_149, \
                         fi_236, fk_194, fk_300, gi_230, gi_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_12 * fi_118[k]
                   + pb_z[k] * gi_230[k];

        t_300[k] = f_13 * fi_236[k]
                   + pa_x[k] * fk_300[k];

        t_301[k] = f_11 * fi_149[k]
                   + pb_y[k] * gi_233[k];

        t_302[k] = pa_y[k] * fk_194[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pa_x, pb_z, fi_122, fi_239, fi_241, \
                         fi_242, fk_303, fk_305, fk_306, gi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_12 * fi_239[k]
                   + pa_x[k] * fk_303[k];

        t_304[k] = f_12 * fi_122[k]
                   + pb_z[k] * gi_234[k];

        t_305[k] = f_12 * fi_241[k]
                   + pa_x[k] * fk_305[k];

        t_306[k] = f_12 * fi_242[k]
                   + pa_x[k] * fk_306[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_y, pb_x, pb_y, fi_154, fi_245, fi_246, \
                         fk_200, gi_238, gi_245, gi_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_11 * fi_154[k]
                   + pb_y[k] * gi_238[k];

        t_308[k] = pa_y[k] * fk_200[k];

        t_309[k] = f_11 * fi_245[k]
                   + pb_x[k] * gi_245[k];

        t_310[k] = f_11 * fi_246[k]
                   + pb_x[k] * gi_246[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, t_315, pa_y, pb_x, fi_247, fi_248, \
                         fi_249, fi_250, fk_207, gi_247, gi_248, gi_249, \
                         gi_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_11 * fi_247[k]
                   + pb_x[k] * gi_247[k];

        t_312[k] = f_11 * fi_248[k]
                   + pb_x[k] * gi_248[k];

        t_313[k] = f_11 * fi_249[k]
                   + pb_x[k] * gi_249[k];

        t_314[k] = f_11 * fi_250[k]
                   + pb_x[k] * gi_250[k];

        t_315[k] = pa_y[k] * fk_207[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, t_320, t_321, t_322, pa_x, fk_316, \
                         fk_317, fk_318, fk_319, fk_320, fk_321, \
                         fk_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_x[k] * fk_316[k];

        t_317[k] = pa_x[k] * fk_317[k];

        t_318[k] = pa_x[k] * fk_318[k];

        t_319[k] = pa_x[k] * fk_319[k];

        t_320[k] = pa_x[k] * fk_320[k];

        t_321[k] = pa_x[k] * fk_321[k];

        t_322[k] = pa_x[k] * fk_322[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, t_327, pa_x, pb_y, pb_z, fi_140, fi_252, \
                         fi_255, fk_323, fk_324, fk_327, gi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = pa_x[k] * fk_323[k];

        t_324[k] = f_15 * fi_252[k]
                   + pa_x[k] * fk_324[k];

        t_325[k] = pb_y[k] * gi_252[k];

        t_326[k] = f_13 * fi_140[k]
                   + pb_z[k] * gi_252[k];

        t_327[k] = f_14 * fi_255[k]
                   + pa_x[k] * fk_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, pa_x, pb_y, pb_z, fi_143, fi_257, \
                         fi_258, fk_329, fk_330, gi_254, gi_255, \
                         gi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pb_y[k] * gi_254[k];

        t_329[k] = f_14 * fi_257[k]
                   + pa_x[k] * fk_329[k];

        t_330[k] = f_0 * fi_258[k]
                   + pa_x[k] * fk_330[k];

        t_331[k] = f_13 * fi_143[k]
                   + pb_z[k] * gi_255[k];

        t_332[k] = pb_y[k] * gi_257[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pa_x, pb_z, fi_146, fi_261, fi_262, \
                         fi_264, fk_333, fk_334, fk_336, gi_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_0 * fi_261[k]
                   + pa_x[k] * fk_333[k];

        t_334[k] = f_13 * fi_262[k]
                   + pa_x[k] * fk_334[k];

        t_335[k] = f_13 * fi_146[k]
                   + pb_z[k] * gi_258[k];

        t_336[k] = f_13 * fi_264[k]
                   + pa_x[k] * fk_336[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_x, pb_y, pb_z, fi_150, fi_266, fi_267, \
                         fk_338, fk_339, gi_261, gi_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = pb_y[k] * gi_261[k];

        t_338[k] = f_13 * fi_266[k]
                   + pa_x[k] * fk_338[k];

        t_339[k] = f_12 * fi_267[k]
                   + pa_x[k] * fk_339[k];

        t_340[k] = f_13 * fi_150[k]
                   + pb_z[k] * gi_262[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_x, pb_y, fi_269, fi_270, fi_272, \
                         fk_341, fk_342, fk_344, gi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_12 * fi_269[k]
                   + pa_x[k] * fk_341[k];

        t_342[k] = f_12 * fi_270[k]
                   + pa_x[k] * fk_342[k];

        t_343[k] = pb_y[k] * gi_266[k];

        t_344[k] = f_12 * fi_272[k]
                   + pa_x[k] * fk_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pb_x, fi_273, fi_274, fi_275, \
                         fi_276, fi_277, gi_273, gi_274, gi_275, gi_276, \
                         gi_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_11 * fi_273[k]
                   + pb_x[k] * gi_273[k];

        t_346[k] = f_11 * fi_274[k]
                   + pb_x[k] * gi_274[k];

        t_347[k] = f_11 * fi_275[k]
                   + pb_x[k] * gi_275[k];

        t_348[k] = f_11 * fi_276[k]
                   + pb_x[k] * gi_276[k];

        t_349[k] = f_11 * fi_277[k]
                   + pb_x[k] * gi_277[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, t_355, pa_x, pb_x, pb_y, fi_279, \
                         fk_352, fk_353, fk_354, fk_355, gi_272, \
                         gi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_y[k] * gi_272[k];

        t_351[k] = f_11 * fi_279[k]
                   + pb_x[k] * gi_279[k];

        t_352[k] = pa_x[k] * fk_352[k];

        t_353[k] = pa_x[k] * fk_353[k];

        t_354[k] = pa_x[k] * fk_354[k];

        t_355[k] = pa_x[k] * fk_355[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, t_360, pa_x, pb_x, pb_y, fk_356, fk_357, \
                         fk_359, gh0_210, gh1_210, gi_279, gi_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * fk_356[k];

        t_357[k] = pa_x[k] * fk_357[k];

        t_358[k] = pb_y[k] * gi_279[k];

        t_359[k] = pa_x[k] * fk_359[k];

        t_360[k] = f_1 * gh0_210[k]
                   - f_2 * gh1_210[k]
                   + pb_x[k] * gi_280[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, t_364, pb_x, pb_y, pb_z, fi_168, gh0_213, \
                         gh1_213, gi_280, gi_281, gi_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = f_0 * fi_168[k]
                   + pb_y[k] * gi_280[k];

        t_362[k] = pb_z[k] * gi_280[k];

        t_363[k] = f_9 * gh0_213[k]
                   - f_10 * gh1_213[k]
                   + pb_x[k] * gi_283[k];

        t_364[k] = pb_z[k] * gi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pb_x, pb_y, pb_z, fi_173, gh0_215, \
                         gh0_216, gh1_215, gh1_216, gi_283, gi_285, \
                         gi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_9 * gh0_215[k]
                   - f_10 * gh1_215[k]
                   + pb_x[k] * gi_285[k];

        t_366[k] = f_7 * gh0_216[k]
                   - f_8 * gh1_216[k]
                   + pb_x[k] * gi_286[k];

        t_367[k] = pb_z[k] * gi_283[k];

        t_368[k] = f_0 * fi_173[k]
                   + pb_y[k] * gi_285[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pb_x, pb_z, gh0_219, gh0_220, gh0_222, \
                         gh1_219, gh1_220, gh1_222, gi_286, gi_289, gi_290, \
                         gi_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_7 * gh0_219[k]
                   - f_8 * gh1_219[k]
                   + pb_x[k] * gi_289[k];

        t_370[k] = f_5 * gh0_220[k]
                   - f_6 * gh1_220[k]
                   + pb_x[k] * gi_290[k];

        t_371[k] = pb_z[k] * gi_286[k];

        t_372[k] = f_5 * gh0_222[k]
                   - f_6 * gh1_222[k]
                   + pb_x[k] * gi_292[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_x, pb_y, pb_z, fi_177, gh0_224, \
                         gh0_225, gh1_224, gh1_225, gi_289, gi_290, gi_294, \
                         gi_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_0 * fi_177[k]
                   + pb_y[k] * gi_289[k];

        t_374[k] = f_5 * gh0_224[k]
                   - f_6 * gh1_224[k]
                   + pb_x[k] * gi_294[k];

        t_375[k] = f_3 * gh0_225[k]
                   - f_4 * gh1_225[k]
                   + pb_x[k] * gi_295[k];

        t_376[k] = pb_z[k] * gi_290[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pb_x, pb_y, fi_182, gh0_227, gh0_228, gh1_227, \
                         gh1_228, gi_294, gi_297, gi_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_3 * gh0_227[k]
                   - f_4 * gh1_227[k]
                   + pb_x[k] * gi_297[k];

        t_378[k] = f_3 * gh0_228[k]
                   - f_4 * gh1_228[k]
                   + pb_x[k] * gi_298[k];

        t_379[k] = f_0 * fi_182[k]
                   + pb_y[k] * gi_294[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, t_385, pb_x, gh0_230, gh1_230, \
                         gi_300, gi_301, gi_302, gi_303, gi_304, \
                         gi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_3 * gh0_230[k]
                   - f_4 * gh1_230[k]
                   + pb_x[k] * gi_300[k];

        t_381[k] = pb_x[k] * gi_301[k];

        t_382[k] = pb_x[k] * gi_302[k];

        t_383[k] = pb_x[k] * gi_303[k];

        t_384[k] = pb_x[k] * gi_304[k];

        t_385[k] = pb_x[k] * gi_305[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, pb_x, pb_y, pb_z, fi_189, gh0_225, \
                         gh1_225, gi_301, gi_302, gi_306, gi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = pb_x[k] * gi_306[k];

        t_387[k] = pb_x[k] * gi_307[k];

        t_388[k] = f_0 * fi_189[k]
                   + f_1 * gh0_225[k]
                   - f_2 * gh1_225[k]
                   + pb_y[k] * gi_301[k];

        t_389[k] = pb_z[k] * gi_301[k];

        t_390[k] = f_3 * gh0_225[k]
                   - f_4 * gh1_225[k]
                   + pb_z[k] * gi_302[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pb_z, gh0_226, gh0_227, gh0_228, gh1_226, \
                         gh1_227, gh1_228, gi_303, gi_304, gi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_5 * gh0_226[k]
                   - f_6 * gh1_226[k]
                   + pb_z[k] * gi_303[k];

        t_392[k] = f_7 * gh0_227[k]
                   - f_8 * gh1_227[k]
                   + pb_z[k] * gi_304[k];

        t_393[k] = f_9 * gh0_228[k]
                   - f_10 * gh1_228[k]
                   + pb_z[k] * gi_305[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pa_z, pb_y, pb_z, fi_168, fi_195, \
                         fk_216, fk_217, gh0_230, gh1_230, gi_307, \
                         gi_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_0 * fi_195[k]
                   + pb_y[k] * gi_307[k];

        t_395[k] = f_1 * gh0_230[k]
                   - f_2 * gh1_230[k]
                   + pb_z[k] * gi_307[k];

        t_396[k] = pa_z[k] * fk_216[k];

        t_397[k] = pa_z[k] * fk_217[k];

        t_398[k] = f_11 * fi_168[k]
                   + pb_z[k] * gi_308[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, pa_z, pb_y, pb_z, fi_170, fi_171, \
                         fi_198, fk_219, fk_221, fk_222, gi_310, \
                         gi_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pa_z[k] * fk_219[k];

        t_400[k] = f_13 * fi_198[k]
                   + pb_y[k] * gi_310[k];

        t_401[k] = f_12 * fi_170[k]
                   + pa_z[k] * fk_221[k];

        t_402[k] = pa_z[k] * fk_222[k];

        t_403[k] = f_11 * fi_171[k]
                   + pb_z[k] * gi_311[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_z, pb_y, pb_z, fi_173, fi_174, fi_201, \
                         fk_225, fk_226, gi_313, gi_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_13 * fi_201[k]
                   + pb_y[k] * gi_313[k];

        t_405[k] = f_13 * fi_173[k]
                   + pa_z[k] * fk_225[k];

        t_406[k] = pa_z[k] * fk_226[k];

        t_407[k] = f_11 * fi_174[k]
                   + pb_z[k] * gi_314[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pa_z, pb_y, fi_175, fi_177, fi_205, \
                         fk_228, fk_230, fk_231, gi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_12 * fi_175[k]
                   + pa_z[k] * fk_228[k];

        t_409[k] = f_13 * fi_205[k]
                   + pb_y[k] * gi_317[k];

        t_410[k] = f_0 * fi_177[k]
                   + pa_z[k] * fk_230[k];

        t_411[k] = pa_z[k] * fk_231[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_z, pb_y, pb_z, fi_178, fi_179, fi_180, \
                         fi_210, fk_233, fk_234, gi_318, gi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_11 * fi_178[k]
                   + pb_z[k] * gi_318[k];

        t_413[k] = f_12 * fi_179[k]
                   + pa_z[k] * fk_233[k];

        t_414[k] = f_13 * fi_180[k]
                   + pa_z[k] * fk_234[k];

        t_415[k] = f_13 * fi_210[k]
                   + pb_y[k] * gi_322[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, t_420, t_421, pa_z, pb_x, fi_182, fk_236, \
                         gi_329, gi_330, gi_331, gi_332, gi_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_14 * fi_182[k]
                   + pa_z[k] * fk_236[k];

        t_417[k] = pb_x[k] * gi_329[k];

        t_418[k] = pb_x[k] * gi_330[k];

        t_419[k] = pb_x[k] * gi_331[k];

        t_420[k] = pb_x[k] * gi_332[k];

        t_421[k] = pb_x[k] * gi_333[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, pa_z, pb_x, pb_z, fi_189, fi_190, \
                         fk_244, fk_246, gi_329, gi_334, gi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = pb_x[k] * gi_334[k];

        t_423[k] = pb_x[k] * gi_335[k];

        t_424[k] = pa_z[k] * fk_244[k];

        t_425[k] = f_11 * fi_189[k]
                   + pb_z[k] * gi_329[k];

        t_426[k] = f_12 * fi_190[k]
                   + pa_z[k] * fk_246[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pa_z, pb_y, fi_191, fi_192, fi_193, \
                         fi_223, fk_247, fk_248, fk_249, gi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_13 * fi_191[k]
                   + pa_z[k] * fk_247[k];

        t_428[k] = f_0 * fi_192[k]
                   + pa_z[k] * fk_248[k];

        t_429[k] = f_14 * fi_193[k]
                   + pa_z[k] * fk_249[k];

        t_430[k] = f_13 * fi_223[k]
                   + pb_y[k] * gi_335[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pa_z, pb_x, pb_y, pb_z, fi_195, fi_196, \
                         fi_224, fk_251, gh0_252, gh1_252, gi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_15 * fi_195[k]
                   + pa_z[k] * fk_251[k];

        t_432[k] = f_1 * gh0_252[k]
                   - f_2 * gh1_252[k]
                   + pb_x[k] * gi_336[k];

        t_433[k] = f_12 * fi_224[k]
                   + pb_y[k] * gi_336[k];

        t_434[k] = f_12 * fi_196[k]
                   + pb_z[k] * gi_336[k];
    }
}

static auto
compute_prim_gk_electron_repulsion_0_piece3(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t dk0,
                                            const size_t dk1, const size_t fi, const size_t fk,
                                            const size_t gh0, const size_t gh1, const size_t gi,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.5 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dk0_136 = buffer.data(dk0 + 136);
    const auto *dk0_215 = buffer.data(dk0 + 215);

    const auto *dk1_136 = buffer.data(dk1 + 136);
    const auto *dk1_215 = buffer.data(dk1 + 215);

    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_323 = buffer.data(fk + 323);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_359 = buffer.data(fk + 359);

    const auto *gh0_255 = buffer.data(gh0 + 255);
    const auto *gh0_257 = buffer.data(gh0 + 257);
    const auto *gh0_258 = buffer.data(gh0 + 258);
    const auto *gh0_261 = buffer.data(gh0 + 261);
    const auto *gh0_262 = buffer.data(gh0 + 262);
    const auto *gh0_264 = buffer.data(gh0 + 264);
    const auto *gh0_266 = buffer.data(gh0 + 266);
    const auto *gh0_267 = buffer.data(gh0 + 267);
    const auto *gh0_269 = buffer.data(gh0 + 269);
    const auto *gh0_270 = buffer.data(gh0 + 270);
    const auto *gh0_271 = buffer.data(gh0 + 271);
    const auto *gh0_272 = buffer.data(gh0 + 272);
    const auto *gh0_294 = buffer.data(gh0 + 294);
    const auto *gh0_297 = buffer.data(gh0 + 297);
    const auto *gh0_299 = buffer.data(gh0 + 299);
    const auto *gh0_300 = buffer.data(gh0 + 300);
    const auto *gh0_303 = buffer.data(gh0 + 303);
    const auto *gh0_304 = buffer.data(gh0 + 304);
    const auto *gh0_306 = buffer.data(gh0 + 306);
    const auto *gh0_308 = buffer.data(gh0 + 308);
    const auto *gh0_309 = buffer.data(gh0 + 309);
    const auto *gh0_311 = buffer.data(gh0 + 311);
    const auto *gh0_312 = buffer.data(gh0 + 312);
    const auto *gh0_313 = buffer.data(gh0 + 313);
    const auto *gh0_314 = buffer.data(gh0 + 314);

    const auto *gh1_255 = buffer.data(gh1 + 255);
    const auto *gh1_257 = buffer.data(gh1 + 257);
    const auto *gh1_258 = buffer.data(gh1 + 258);
    const auto *gh1_261 = buffer.data(gh1 + 261);
    const auto *gh1_262 = buffer.data(gh1 + 262);
    const auto *gh1_264 = buffer.data(gh1 + 264);
    const auto *gh1_266 = buffer.data(gh1 + 266);
    const auto *gh1_267 = buffer.data(gh1 + 267);
    const auto *gh1_269 = buffer.data(gh1 + 269);
    const auto *gh1_270 = buffer.data(gh1 + 270);
    const auto *gh1_271 = buffer.data(gh1 + 271);
    const auto *gh1_272 = buffer.data(gh1 + 272);
    const auto *gh1_294 = buffer.data(gh1 + 294);
    const auto *gh1_297 = buffer.data(gh1 + 297);
    const auto *gh1_299 = buffer.data(gh1 + 299);
    const auto *gh1_300 = buffer.data(gh1 + 300);
    const auto *gh1_303 = buffer.data(gh1 + 303);
    const auto *gh1_304 = buffer.data(gh1 + 304);
    const auto *gh1_306 = buffer.data(gh1 + 306);
    const auto *gh1_308 = buffer.data(gh1 + 308);
    const auto *gh1_309 = buffer.data(gh1 + 309);
    const auto *gh1_311 = buffer.data(gh1 + 311);
    const auto *gh1_312 = buffer.data(gh1 + 312);
    const auto *gh1_313 = buffer.data(gh1 + 313);
    const auto *gh1_314 = buffer.data(gh1 + 314);

    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

#pragma omp simd aligned(t_435, t_436, t_437, pb_x, pb_y, fi_226, gh0_255, gh0_257, gh1_255, \
                         gh1_257, gi_338, gi_339, gi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_9 * gh0_255[k]
                   - f_10 * gh1_255[k]
                   + pb_x[k] * gi_339[k];

        t_436[k] = f_12 * fi_226[k]
                   + pb_y[k] * gi_338[k];

        t_437[k] = f_9 * gh0_257[k]
                   - f_10 * gh1_257[k]
                   + pb_x[k] * gi_341[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pb_x, pb_y, pb_z, fi_199, fi_229, gh0_258, \
                         gh1_258, gi_339, gi_341, gi_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_7 * gh0_258[k]
                   - f_8 * gh1_258[k]
                   + pb_x[k] * gi_342[k];

        t_439[k] = f_12 * fi_199[k]
                   + pb_z[k] * gi_339[k];

        t_440[k] = f_12 * fi_229[k]
                   + pb_y[k] * gi_341[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pb_x, pb_z, fi_202, gh0_261, gh0_262, gh1_261, \
                         gh1_262, gi_342, gi_345, gi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_7 * gh0_261[k]
                   - f_8 * gh1_261[k]
                   + pb_x[k] * gi_345[k];

        t_442[k] = f_5 * gh0_262[k]
                   - f_6 * gh1_262[k]
                   + pb_x[k] * gi_346[k];

        t_443[k] = f_12 * fi_202[k]
                   + pb_z[k] * gi_342[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_x, pb_y, fi_233, gh0_264, gh0_266, gh1_264, \
                         gh1_266, gi_345, gi_348, gi_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_5 * gh0_264[k]
                   - f_6 * gh1_264[k]
                   + pb_x[k] * gi_348[k];

        t_445[k] = f_12 * fi_233[k]
                   + pb_y[k] * gi_345[k];

        t_446[k] = f_5 * gh0_266[k]
                   - f_6 * gh1_266[k]
                   + pb_x[k] * gi_350[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pb_x, pb_z, fi_206, gh0_267, gh0_269, gh1_267, \
                         gh1_269, gi_346, gi_351, gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_3 * gh0_267[k]
                   - f_4 * gh1_267[k]
                   + pb_x[k] * gi_351[k];

        t_448[k] = f_12 * fi_206[k]
                   + pb_z[k] * gi_346[k];

        t_449[k] = f_3 * gh0_269[k]
                   - f_4 * gh1_269[k]
                   + pb_x[k] * gi_353[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pb_x, pb_y, fi_238, gh0_270, gh0_272, \
                         gh1_270, gh1_272, gi_350, gi_354, gi_356, \
                         gi_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_3 * gh0_270[k]
                   - f_4 * gh1_270[k]
                   + pb_x[k] * gi_354[k];

        t_451[k] = f_12 * fi_238[k]
                   + pb_y[k] * gi_350[k];

        t_452[k] = f_3 * gh0_272[k]
                   - f_4 * gh1_272[k]
                   + pb_x[k] * gi_356[k];

        t_453[k] = pb_x[k] * gi_357[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, t_459, pb_x, gi_358, gi_359, \
                         gi_360, gi_361, gi_362, gi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pb_x[k] * gi_358[k];

        t_455[k] = pb_x[k] * gi_359[k];

        t_456[k] = pb_x[k] * gi_360[k];

        t_457[k] = pb_x[k] * gi_361[k];

        t_458[k] = pb_x[k] * gi_362[k];

        t_459[k] = pb_x[k] * gi_363[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pa_z, pb_y, pb_z, dk0_136, dk1_136, fi_217, \
                         fi_247, fk_280, gh0_269, gh1_269, gi_357, \
                         gi_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_16 * dk0_136[k]
                   - f_17 * dk1_136[k]
                   + pa_z[k] * fk_280[k];

        t_461[k] = f_12 * fi_217[k]
                   + pb_z[k] * gi_357[k];

        t_462[k] = f_12 * fi_247[k]
                   + f_9 * gh0_269[k]
                   - f_10 * gh1_269[k]
                   + pb_y[k] * gi_359[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pb_y, fi_248, fi_249, fi_250, gh0_270, gh0_271, \
                         gh0_272, gh1_270, gh1_271, gh1_272, gi_360, gi_361, \
                         gi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_12 * fi_248[k]
                   + f_7 * gh0_270[k]
                   - f_8 * gh1_270[k]
                   + pb_y[k] * gi_360[k];

        t_464[k] = f_12 * fi_249[k]
                   + f_5 * gh0_271[k]
                   - f_6 * gh1_271[k]
                   + pb_y[k] * gi_361[k];

        t_465[k] = f_12 * fi_250[k]
                   + f_3 * gh0_272[k]
                   - f_4 * gh1_272[k]
                   + pb_y[k] * gi_362[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, pa_y, pb_y, dk0_215, dk1_215, \
                         fi_251, fi_252, fk_323, fk_324, fk_326, gi_363, \
                         gi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_12 * fi_251[k]
                   + pb_y[k] * gi_363[k];

        t_467[k] = f_16 * dk0_215[k]
                   - f_17 * dk1_215[k]
                   + pa_y[k] * fk_323[k];

        t_468[k] = pa_y[k] * fk_324[k];

        t_469[k] = f_11 * fi_252[k]
                   + pb_y[k] * gi_364[k];

        t_470[k] = pa_y[k] * fk_326[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pa_y, pb_y, fi_253, fi_254, fi_255, \
                         fk_327, fk_329, fk_330, gi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_12 * fi_253[k]
                   + pa_y[k] * fk_327[k];

        t_472[k] = f_11 * fi_254[k]
                   + pb_y[k] * gi_366[k];

        t_473[k] = pa_y[k] * fk_329[k];

        t_474[k] = f_13 * fi_255[k]
                   + pa_y[k] * fk_330[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_y, pb_y, pb_z, fi_227, fi_257, fi_258, \
                         fk_333, fk_334, gi_367, gi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_13 * fi_227[k]
                   + pb_z[k] * gi_367[k];

        t_476[k] = f_11 * fi_257[k]
                   + pb_y[k] * gi_369[k];

        t_477[k] = pa_y[k] * fk_333[k];

        t_478[k] = f_0 * fi_258[k]
                   + pa_y[k] * fk_334[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_y, pb_y, pb_z, fi_230, fi_260, fi_261, \
                         fk_336, fk_338, gi_370, gi_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * fi_230[k]
                   + pb_z[k] * gi_370[k];

        t_480[k] = f_12 * fi_260[k]
                   + pa_y[k] * fk_336[k];

        t_481[k] = f_11 * fi_261[k]
                   + pb_y[k] * gi_373[k];

        t_482[k] = pa_y[k] * fk_338[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_y, pb_z, fi_234, fi_262, fi_264, \
                         fi_265, fk_339, fk_341, fk_342, gi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_14 * fi_262[k]
                   + pa_y[k] * fk_339[k];

        t_484[k] = f_13 * fi_234[k]
                   + pb_z[k] * gi_374[k];

        t_485[k] = f_13 * fi_264[k]
                   + pa_y[k] * fk_341[k];

        t_486[k] = f_12 * fi_265[k]
                   + pa_y[k] * fk_342[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, t_492, pa_y, pb_x, pb_y, fi_266, \
                         fk_344, gi_378, gi_385, gi_386, gi_387, \
                         gi_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_11 * fi_266[k]
                   + pb_y[k] * gi_378[k];

        t_488[k] = pa_y[k] * fk_344[k];

        t_489[k] = pb_x[k] * gi_385[k];

        t_490[k] = pb_x[k] * gi_386[k];

        t_491[k] = pb_x[k] * gi_387[k];

        t_492[k] = pb_x[k] * gi_388[k];
    }

#pragma omp simd aligned(t_493, t_494, t_495, t_496, t_497, pa_y, pb_x, pb_z, fi_245, fi_273, \
                         fk_352, gi_385, gi_389, gi_390, gi_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = pb_x[k] * gi_389[k];

        t_494[k] = pb_x[k] * gi_390[k];

        t_495[k] = pb_x[k] * gi_391[k];

        t_496[k] = f_15 * fi_273[k]
                   + pa_y[k] * fk_352[k];

        t_497[k] = f_13 * fi_245[k]
                   + pb_z[k] * gi_385[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pa_y, fi_275, fi_276, fi_277, fi_278, \
                         fk_354, fk_355, fk_356, fk_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_14 * fi_275[k]
                   + pa_y[k] * fk_354[k];

        t_499[k] = f_0 * fi_276[k]
                   + pa_y[k] * fk_355[k];

        t_500[k] = f_13 * fi_277[k]
                   + pa_y[k] * fk_356[k];

        t_501[k] = f_12 * fi_278[k]
                   + pa_y[k] * fk_357[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, pa_y, pb_x, pb_y, pb_z, fi_252, \
                         fi_279, fk_359, gh0_294, gh1_294, gi_391, \
                         gi_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_11 * fi_279[k]
                   + pb_y[k] * gi_391[k];

        t_503[k] = pa_y[k] * fk_359[k];

        t_504[k] = f_1 * gh0_294[k]
                   - f_2 * gh1_294[k]
                   + pb_x[k] * gi_392[k];

        t_505[k] = pb_y[k] * gi_392[k];

        t_506[k] = f_0 * fi_252[k]
                   + pb_z[k] * gi_392[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pb_x, pb_y, gh0_297, gh0_299, gh0_300, \
                         gh1_297, gh1_299, gh1_300, gi_394, gi_395, gi_397, \
                         gi_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_9 * gh0_297[k]
                   - f_10 * gh1_297[k]
                   + pb_x[k] * gi_395[k];

        t_508[k] = pb_y[k] * gi_394[k];

        t_509[k] = f_9 * gh0_299[k]
                   - f_10 * gh1_299[k]
                   + pb_x[k] * gi_397[k];

        t_510[k] = f_7 * gh0_300[k]
                   - f_8 * gh1_300[k]
                   + pb_x[k] * gi_398[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pb_x, pb_y, pb_z, fi_255, gh0_303, \
                         gh0_304, gh1_303, gh1_304, gi_395, gi_397, gi_401, \
                         gi_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_0 * fi_255[k]
                   + pb_z[k] * gi_395[k];

        t_512[k] = pb_y[k] * gi_397[k];

        t_513[k] = f_7 * gh0_303[k]
                   - f_8 * gh1_303[k]
                   + pb_x[k] * gi_401[k];

        t_514[k] = f_5 * gh0_304[k]
                   - f_6 * gh1_304[k]
                   + pb_x[k] * gi_402[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pb_x, pb_y, pb_z, fi_258, gh0_306, \
                         gh0_308, gh1_306, gh1_308, gi_398, gi_401, gi_404, \
                         gi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_0 * fi_258[k]
                   + pb_z[k] * gi_398[k];

        t_516[k] = f_5 * gh0_306[k]
                   - f_6 * gh1_306[k]
                   + pb_x[k] * gi_404[k];

        t_517[k] = pb_y[k] * gi_401[k];

        t_518[k] = f_5 * gh0_308[k]
                   - f_6 * gh1_308[k]
                   + pb_x[k] * gi_406[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pb_x, pb_z, fi_262, gh0_309, gh0_311, gh1_309, \
                         gh1_311, gi_402, gi_407, gi_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_3 * gh0_309[k]
                   - f_4 * gh1_309[k]
                   + pb_x[k] * gi_407[k];

        t_520[k] = f_0 * fi_262[k]
                   + pb_z[k] * gi_402[k];

        t_521[k] = f_3 * gh0_311[k]
                   - f_4 * gh1_311[k]
                   + pb_x[k] * gi_409[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pb_x, pb_y, gh0_312, gh0_314, \
                         gh1_312, gh1_314, gi_406, gi_410, gi_412, gi_413, \
                         gi_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_3 * gh0_312[k]
                   - f_4 * gh1_312[k]
                   + pb_x[k] * gi_410[k];

        t_523[k] = pb_y[k] * gi_406[k];

        t_524[k] = f_3 * gh0_314[k]
                   - f_4 * gh1_314[k]
                   + pb_x[k] * gi_412[k];

        t_525[k] = pb_x[k] * gi_413[k];

        t_526[k] = pb_x[k] * gi_414[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, t_531, t_532, pb_x, pb_y, gh0_309, \
                         gh1_309, gi_413, gi_415, gi_416, gi_417, gi_418, \
                         gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = pb_x[k] * gi_415[k];

        t_528[k] = pb_x[k] * gi_416[k];

        t_529[k] = pb_x[k] * gi_417[k];

        t_530[k] = pb_x[k] * gi_418[k];

        t_531[k] = pb_x[k] * gi_419[k];

        t_532[k] = f_1 * gh0_309[k]
                   - f_2 * gh1_309[k]
                   + pb_y[k] * gi_413[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pb_y, pb_z, fi_273, gh0_311, gh0_312, gh1_311, \
                         gh1_312, gi_413, gi_415, gi_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_0 * fi_273[k]
                   + pb_z[k] * gi_413[k];

        t_534[k] = f_9 * gh0_311[k]
                   - f_10 * gh1_311[k]
                   + pb_y[k] * gi_415[k];

        t_535[k] = f_7 * gh0_312[k]
                   - f_8 * gh1_312[k]
                   + pb_y[k] * gi_416[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pb_y, pb_z, fi_279, gh0_313, gh0_314, \
                         gh1_313, gh1_314, gi_417, gi_418, gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_5 * gh0_313[k]
                   - f_6 * gh1_313[k]
                   + pb_y[k] * gi_417[k];

        t_537[k] = f_3 * gh0_314[k]
                   - f_4 * gh1_314[k]
                   + pb_y[k] * gi_418[k];

        t_538[k] = pb_y[k] * gi_419[k];

        t_539[k] = f_0 * fi_279[k]
                   + f_1 * gh0_314[k]
                   - f_2 * gh1_314[k]
                   + pb_z[k] * gi_419[k];
    }
}

auto
compute_prim_gk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t dk0, const size_t dk1,
                                     const size_t fi, const size_t fk, const size_t gh0,
                                     const size_t gh1, const size_t gi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_gk_electron_repulsion_0_piece0(buffer, target, pa, pb, dk0, dk1, fi, fk, gh0,
                                                gh1, gi, ncols, alpha, beta, p);

    compute_prim_gk_electron_repulsion_0_piece1(buffer, target, pa, pb, dk0, dk1, fi, fk, gh0,
                                                gh1, gi, ncols, alpha, beta, p);

    compute_prim_gk_electron_repulsion_0_piece2(buffer, target, pa, pb, fi, fk, gh0, gh1, gi,
                                                ncols, alpha, beta, p);

    compute_prim_gk_electron_repulsion_0_piece3(buffer, target, pa, pb, dk0, dk1, fi, fk, gh0,
                                                gh1, gi, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
