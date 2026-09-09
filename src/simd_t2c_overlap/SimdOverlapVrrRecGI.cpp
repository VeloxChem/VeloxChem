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


#include "SimdOverlapVrrRecGI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

static auto
compute_prim_gi_overlap_0_piece0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t di, const size_t fh,
                                 const size_t fi, const size_t gg, const size_t gh,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_114 = buffer.data(fh + 114);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_84 = buffer.data(gg + 84);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fh_0, gg_0, gh_0, \
                         gh_1, gh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fh_0[k]
                 + f_1 * gg_0[k]
                 + pb_x[k] * gh_0[k];

        t_1[k] = pb_y[k] * gh_0[k];

        t_2[k] = pb_z[k] * gh_0[k];

        t_3[k] = f_2 * gg_0[k]
                 + pb_y[k] * gh_1[k];

        t_4[k] = pb_y[k] * gh_2[k];

        t_5[k] = f_2 * gg_0[k]
                 + pb_z[k] * gh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, gg_1, gg_2, gg_3, gh_3, \
                         gh_5, gh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gg_1[k]
                 + pb_y[k] * gh_3[k];

        t_7[k] = pb_z[k] * gh_3[k];

        t_8[k] = pb_y[k] * gh_5[k];

        t_9[k] = f_3 * gg_2[k]
                 + pb_z[k] * gh_5[k];

        t_10[k] = f_4 * gg_3[k]
                  + pb_y[k] * gh_6[k];

        t_11[k] = pb_z[k] * gh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, fh_15, gg_5, gh_8, \
                         gh_9, gh_10, gh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * gg_5[k]
                  + pb_y[k] * gh_8[k];

        t_13[k] = pb_y[k] * gh_9[k];

        t_14[k] = f_4 * gg_5[k]
                  + pb_z[k] * gh_9[k];

        t_15[k] = f_0 * fh_15[k]
                  + pb_x[k] * gh_15[k];

        t_16[k] = pb_z[k] * gh_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, fh_17, fh_18, fh_20, gg_10, \
                         gh_14, gh_15, gh_17, gh_18, gh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fh_17[k]
                  + pb_x[k] * gh_17[k];

        t_18[k] = f_0 * fh_18[k]
                  + pb_x[k] * gh_18[k];

        t_19[k] = pb_y[k] * gh_14[k];

        t_20[k] = f_0 * fh_20[k]
                  + pb_x[k] * gh_20[k];

        t_21[k] = f_1 * gg_10[k]
                  + pb_y[k] * gh_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, gg_12, gg_13, gg_14, \
                         gh_15, gh_17, gh_18, gh_19, gh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * gh_15[k];

        t_23[k] = f_4 * gg_12[k]
                  + pb_y[k] * gh_17[k];

        t_24[k] = f_3 * gg_13[k]
                  + pb_y[k] * gh_18[k];

        t_25[k] = f_2 * gg_14[k]
                  + pb_y[k] * gh_19[k];

        t_26[k] = pb_y[k] * gh_20[k];

        t_27[k] = f_1 * gg_14[k]
                  + pb_z[k] * gh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pa_y, pb_y, pb_z, fh_0, fh_1, \
                         fi_0, fi_3, fi_5, gh_21, gh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * fi_0[k];

        t_29[k] = f_2 * fh_0[k]
                  + pb_y[k] * gh_21[k];

        t_30[k] = pb_z[k] * gh_21[k];

        t_31[k] = f_3 * fh_1[k]
                  + pa_y[k] * fi_3[k];

        t_32[k] = pb_z[k] * gh_22[k];

        t_33[k] = pa_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, fh_3, fh_5, fh_6, \
                         fi_6, fi_9, fi_10, gh_24, gh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_4 * fh_3[k]
                  + pa_y[k] * fi_6[k];

        t_35[k] = pb_z[k] * gh_24[k];

        t_36[k] = f_2 * fh_5[k]
                  + pb_y[k] * gh_26[k];

        t_37[k] = pa_y[k] * fi_9[k];

        t_38[k] = f_0 * fh_6[k]
                  + pa_y[k] * fi_10[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_y, pb_z, fh_9, fi_14, gg_18, gh_27, \
                         gh_28, gh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_z[k] * gh_27[k];

        t_40[k] = f_2 * gg_18[k]
                  + pb_z[k] * gh_28[k];

        t_41[k] = f_2 * fh_9[k]
                  + pb_y[k] * gh_30[k];

        t_42[k] = pa_y[k] * fi_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_z, fh_36, fh_38, fh_39, fh_40, \
                         gh_31, gh_36, gh_38, gh_39, gh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_4 * fh_36[k]
                  + pb_x[k] * gh_36[k];

        t_44[k] = pb_z[k] * gh_31[k];

        t_45[k] = f_4 * fh_38[k]
                  + pb_x[k] * gh_38[k];

        t_46[k] = f_4 * fh_39[k]
                  + pb_x[k] * gh_39[k];

        t_47[k] = f_4 * fh_40[k]
                  + pb_x[k] * gh_40[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pa_y, pb_z, di_49, fi_20, fi_49, \
                         gg_25, gg_26, gh_36, gh_37, gh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * fi_20[k];

        t_49[k] = f_3 * di_49[k]
                  + pa_x[k] * fi_49[k];

        t_50[k] = pb_z[k] * gh_36[k];

        t_51[k] = f_2 * gg_25[k]
                  + pb_z[k] * gh_37[k];

        t_52[k] = f_3 * gg_26[k]
                  + pb_z[k] * gh_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pa_z, pb_y, pb_z, fh_20, fi_0, \
                         fi_27, gg_27, gh_39, gh_41, gh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * gg_27[k]
                  + pb_z[k] * gh_39[k];

        t_54[k] = f_2 * fh_20[k]
                  + pb_y[k] * gh_41[k];

        t_55[k] = pa_y[k] * fi_27[k];

        t_56[k] = pa_z[k] * fi_0[k];

        t_57[k] = pb_y[k] * gh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, fh_0, fh_2, fi_3, \
                         fi_5, fi_6, gh_42, gh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * fh_0[k]
                  + pb_z[k] * gh_42[k];

        t_59[k] = pa_z[k] * fi_3[k];

        t_60[k] = pb_y[k] * gh_44[k];

        t_61[k] = f_3 * fh_2[k]
                  + pa_z[k] * fi_5[k];

        t_62[k] = pa_z[k] * fi_6[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_z, pb_y, fh_5, fi_9, fi_10, gg_32, \
                         gg_34, gh_46, gh_47, gh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * gg_32[k]
                  + pb_y[k] * gh_46[k];

        t_64[k] = pb_y[k] * gh_47[k];

        t_65[k] = f_4 * fh_5[k]
                  + pa_z[k] * fi_9[k];

        t_66[k] = pa_z[k] * fi_10[k];

        t_67[k] = f_3 * gg_34[k]
                  + pb_y[k] * gh_49[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_x, pb_y, fh_9, fh_58, fi_14, \
                         fi_15, gg_35, gh_50, gh_51, gh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * gg_35[k]
                  + pb_y[k] * gh_50[k];

        t_69[k] = pb_y[k] * gh_51[k];

        t_70[k] = f_0 * fh_9[k]
                  + pa_z[k] * fi_14[k];

        t_71[k] = pa_z[k] * fi_15[k];

        t_72[k] = f_4 * fh_58[k]
                  + pb_x[k] * gh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, fh_59, fh_60, fh_62, \
                         fi_21, gh_56, gh_59, gh_60, gh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * fh_59[k]
                  + pb_x[k] * gh_59[k];

        t_74[k] = f_4 * fh_60[k]
                  + pb_x[k] * gh_60[k];

        t_75[k] = pb_y[k] * gh_56[k];

        t_76[k] = f_4 * fh_62[k]
                  + pb_x[k] * gh_62[k];

        t_77[k] = pa_z[k] * fi_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_y, gg_41, gg_42, gg_43, gg_44, \
                         gh_58, gh_59, gh_60, gh_61, gh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * gg_41[k]
                  + pb_y[k] * gh_58[k];

        t_79[k] = f_4 * gg_42[k]
                  + pb_y[k] * gh_59[k];

        t_80[k] = f_3 * gg_43[k]
                  + pb_y[k] * gh_60[k];

        t_81[k] = f_2 * gg_44[k]
                  + pb_y[k] * gh_61[k];

        t_82[k] = pb_y[k] * gh_62[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pa_y, pb_y, pb_z, di_0, di_83, fh_21, \
                         fi_28, fi_83, gh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * di_83[k]
                  + pa_x[k] * fi_83[k];

        t_84[k] = f_2 * di_0[k]
                  + pa_y[k] * fi_28[k];

        t_85[k] = f_3 * fh_21[k]
                  + pb_y[k] * gh_63[k];

        t_86[k] = pb_z[k] * gh_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pb_x, pb_z, fh_66, fh_69, gg_45, gg_48, \
                         gg_51, gh_64, gh_65, gh_66, gh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * fh_66[k]
                  + f_4 * gg_48[k]
                  + pb_x[k] * gh_66[k];

        t_88[k] = pb_z[k] * gh_64[k];

        t_89[k] = f_2 * gg_45[k]
                  + pb_z[k] * gh_65[k];

        t_90[k] = f_3 * fh_69[k]
                  + f_3 * gg_51[k]
                  + pb_x[k] * gh_69[k];

        t_91[k] = pb_z[k] * gh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, fh_26, fh_73, gg_47, gg_55, \
                         gh_68, gh_69, gh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_3 * fh_26[k]
                  + pb_y[k] * gh_68[k];

        t_93[k] = f_3 * gg_47[k]
                  + pb_z[k] * gh_68[k];

        t_94[k] = f_3 * fh_73[k]
                  + f_2 * gg_55[k]
                  + pb_x[k] * gh_73[k];

        t_95[k] = pb_z[k] * gh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, fh_30, fh_78, gg_48, \
                         gg_50, gh_70, gh_72, gh_73, gh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_2 * gg_48[k]
                  + pb_z[k] * gh_70[k];

        t_97[k] = f_3 * fh_30[k]
                  + pb_y[k] * gh_72[k];

        t_98[k] = f_4 * gg_50[k]
                  + pb_z[k] * gh_72[k];

        t_99[k] = f_3 * fh_78[k]
                  + pb_x[k] * gh_78[k];

        t_100[k] = pb_z[k] * gh_73[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, fh_80, fh_81, fh_82, fh_83, gh_80, \
                         gh_81, gh_82, gh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_3 * fh_80[k]
                   + pb_x[k] * gh_80[k];

        t_102[k] = f_3 * fh_81[k]
                   + pb_x[k] * gh_81[k];

        t_103[k] = f_3 * fh_82[k]
                   + pb_x[k] * gh_82[k];

        t_104[k] = f_3 * fh_83[k]
                   + pb_x[k] * gh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_x, pb_z, di_105, fi_105, gg_55, \
                         gg_56, gg_57, gh_78, gh_79, gh_80, gh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_2 * di_105[k]
                   + pa_x[k] * fi_105[k];

        t_106[k] = pb_z[k] * gh_78[k];

        t_107[k] = f_2 * gg_55[k]
                   + pb_z[k] * gh_79[k];

        t_108[k] = f_3 * gg_56[k]
                   + pb_z[k] * gh_80[k];

        t_109[k] = f_4 * gg_57[k]
                   + pb_z[k] * gh_81[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_y, pa_z, pb_y, pb_z, fh_41, \
                         fi_29, fi_56, fi_58, gg_59, gh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * fh_41[k]
                   + pb_y[k] * gh_83[k];

        t_111[k] = f_1 * gg_59[k]
                   + pb_z[k] * gh_83[k];

        t_112[k] = pa_y[k] * fi_56[k];

        t_113[k] = pa_z[k] * fi_29[k];

        t_114[k] = pa_y[k] * fi_58[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pa_z, pb_y, pb_z, fh_24, \
                         fh_44, fi_31, fi_34, fi_61, gh_86, gh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * fi_31[k];

        t_116[k] = f_2 * fh_44[k]
                   + pb_y[k] * gh_86[k];

        t_117[k] = pa_y[k] * fi_61[k];

        t_118[k] = pa_z[k] * fi_34[k];

        t_119[k] = f_2 * fh_24[k]
                   + pb_z[k] * gh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pa_z, pb_y, pb_z, fh_27, fh_47, \
                         fi_38, fi_65, gh_89, gh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_2 * fh_47[k]
                   + pb_y[k] * gh_89[k];

        t_121[k] = pa_y[k] * fi_65[k];

        t_122[k] = pa_z[k] * fi_38[k];

        t_123[k] = f_2 * fh_27[k]
                   + pb_z[k] * gh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pa_z, pb_y, fh_50, fh_51, fi_43, \
                         fi_68, fi_70, gh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_3 * fh_50[k]
                   + pa_y[k] * fi_68[k];

        t_125[k] = f_2 * fh_51[k]
                   + pb_y[k] * gh_93[k];

        t_126[k] = pa_y[k] * fi_70[k];

        t_127[k] = pa_z[k] * fi_43[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, fh_100, fh_101, \
                         fh_102, fh_103, fi_76, gh_100, gh_101, gh_102, \
                         gh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * fh_100[k]
                   + pb_x[k] * gh_100[k];

        t_129[k] = f_3 * fh_101[k]
                   + pb_x[k] * gh_101[k];

        t_130[k] = f_3 * fh_102[k]
                   + pb_x[k] * gh_102[k];

        t_131[k] = f_3 * fh_103[k]
                   + pb_x[k] * gh_103[k];

        t_132[k] = pa_y[k] * fi_76[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_x, pa_z, pb_z, di_135, di_136, fh_36, \
                         fi_49, fi_135, fi_136, gh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * fi_49[k];

        t_134[k] = f_2 * fh_36[k]
                   + pb_z[k] * gh_99[k];

        t_135[k] = f_2 * di_135[k]
                   + pa_x[k] * fi_135[k];

        t_136[k] = f_2 * di_136[k]
                   + pa_x[k] * fi_136[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_x, pa_y, pa_z, pb_y, di_0, di_137, \
                         fh_62, fi_56, fi_83, fi_137, gh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * di_137[k]
                   + pa_x[k] * fi_137[k];

        t_138[k] = f_2 * fh_62[k]
                   + pb_y[k] * gh_104[k];

        t_139[k] = pa_y[k] * fi_83[k];

        t_140[k] = f_2 * di_0[k]
                   + pa_z[k] * fi_56[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, fh_42, fh_110, \
                         gg_75, gg_80, gh_105, gh_106, gh_107, gh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_y[k] * gh_105[k];

        t_142[k] = f_3 * fh_42[k]
                   + pb_z[k] * gh_105[k];

        t_143[k] = f_2 * gg_75[k]
                   + pb_y[k] * gh_106[k];

        t_144[k] = pb_y[k] * gh_107[k];

        t_145[k] = f_3 * fh_110[k]
                   + f_4 * gg_80[k]
                   + pb_x[k] * gh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_x, pb_y, fh_114, gg_76, gg_77, gg_84, \
                         gh_108, gh_109, gh_110, gh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_3 * gg_76[k]
                   + pb_y[k] * gh_108[k];

        t_147[k] = f_2 * gg_77[k]
                   + pb_y[k] * gh_109[k];

        t_148[k] = pb_y[k] * gh_110[k];

        t_149[k] = f_3 * fh_114[k]
                   + f_3 * gg_84[k]
                   + pb_x[k] * gh_114[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, gg_78, gg_79, gg_80, gh_111, \
                         gh_112, gh_113, gh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_4 * gg_78[k]
                   + pb_y[k] * gh_111[k];

        t_151[k] = f_3 * gg_79[k]
                   + pb_y[k] * gh_112[k];

        t_152[k] = f_2 * gg_80[k]
                   + pb_y[k] * gh_113[k];

        t_153[k] = pb_y[k] * gh_114[k];
    }
}

static auto
compute_prim_gi_overlap_0_piece1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t di, const size_t fh,
                                 const size_t fi, const size_t gg, const size_t gh,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_167 = buffer.data(di + 167);

    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_178 = buffer.data(fi + 178);
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
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_279 = buffer.data(fi + 279);

    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_136 = buffer.data(gg + 136);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_138 = buffer.data(gg + 138);
    const auto *gg_139 = buffer.data(gg + 139);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_169 = buffer.data(gg + 169);
    const auto *gg_170 = buffer.data(gg + 170);

    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pb_x, fh_119, fh_120, fh_121, fh_122, \
                         gg_89, gh_119, gh_120, gh_121, gh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_3 * fh_119[k]
                   + f_2 * gg_89[k]
                   + pb_x[k] * gh_119[k];

        t_155[k] = f_3 * fh_120[k]
                   + pb_x[k] * gh_120[k];

        t_156[k] = f_3 * fh_121[k]
                   + pb_x[k] * gh_121[k];

        t_157[k] = f_3 * fh_122[k]
                   + pb_x[k] * gh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pb_x, pb_y, fh_123, fh_125, gg_85, \
                         gg_86, gh_119, gh_120, gh_121, gh_123, \
                         gh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * fh_123[k]
                   + pb_x[k] * gh_123[k];

        t_159[k] = pb_y[k] * gh_119[k];

        t_160[k] = f_3 * fh_125[k]
                   + pb_x[k] * gh_125[k];

        t_161[k] = f_1 * gg_85[k]
                   + pb_y[k] * gh_120[k];

        t_162[k] = f_0 * gg_86[k]
                   + pb_y[k] * gh_121[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_x, pb_y, di_167, fi_167, gg_87, \
                         gg_88, gg_89, gh_122, gh_123, gh_124, gh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_4 * gg_87[k]
                   + pb_y[k] * gh_122[k];

        t_164[k] = f_3 * gg_88[k]
                   + pb_y[k] * gh_123[k];

        t_165[k] = f_2 * gg_89[k]
                   + pb_y[k] * gh_124[k];

        t_166[k] = pb_y[k] * gh_125[k];

        t_167[k] = f_2 * di_167[k]
                   + pa_x[k] * fi_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pa_x, pb_y, pb_z, fh_63, fh_126, \
                         fh_129, fi_168, fi_171, gh_126, gh_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_5 * fh_126[k]
                   + pa_x[k] * fi_168[k];

        t_169[k] = f_4 * fh_63[k]
                   + pb_y[k] * gh_126[k];

        t_170[k] = pb_z[k] * gh_126[k];

        t_171[k] = f_0 * fh_129[k]
                   + pa_x[k] * fi_171[k];

        t_172[k] = pb_z[k] * gh_127[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pa_x, pb_y, pb_z, fh_68, fh_132, \
                         fi_174, gg_90, gg_92, gh_128, gh_129, gh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_2 * gg_90[k]
                   + pb_z[k] * gh_128[k];

        t_174[k] = f_4 * fh_132[k]
                   + pa_x[k] * fi_174[k];

        t_175[k] = pb_z[k] * gh_129[k];

        t_176[k] = f_4 * fh_68[k]
                   + pb_y[k] * gh_131[k];

        t_177[k] = f_3 * gg_92[k]
                   + pb_z[k] * gh_131[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_x, pb_y, pb_z, fh_72, fh_136, \
                         fi_178, gg_93, gg_95, gh_132, gh_133, gh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * fh_136[k]
                   + pa_x[k] * fi_178[k];

        t_179[k] = pb_z[k] * gh_132[k];

        t_180[k] = f_2 * gg_93[k]
                   + pb_z[k] * gh_133[k];

        t_181[k] = f_4 * fh_72[k]
                   + pb_y[k] * gh_135[k];

        t_182[k] = f_4 * gg_95[k]
                   + pb_z[k] * gh_135[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pb_x, pb_z, fh_141, fh_143, \
                         fh_144, fh_145, gh_136, gh_141, gh_143, gh_144, \
                         gh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_2 * fh_141[k]
                   + pb_x[k] * gh_141[k];

        t_184[k] = pb_z[k] * gh_136[k];

        t_185[k] = f_2 * fh_143[k]
                   + pb_x[k] * gh_143[k];

        t_186[k] = f_2 * fh_144[k]
                   + pb_x[k] * gh_144[k];

        t_187[k] = f_2 * fh_145[k]
                   + pb_x[k] * gh_145[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, t_193, pa_x, pb_x, pb_z, fh_146, \
                         fi_189, fi_191, fi_192, fi_193, gh_141, \
                         gh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_2 * fh_146[k]
                   + pb_x[k] * gh_146[k];

        t_189[k] = pa_x[k] * fi_189[k];

        t_190[k] = pb_z[k] * gh_141[k];

        t_191[k] = pa_x[k] * fi_191[k];

        t_192[k] = pa_x[k] * fi_192[k];

        t_193[k] = pa_x[k] * fi_193[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, t_199, pa_x, pa_z, pb_z, fh_63, \
                         fi_84, fi_85, fi_87, fi_194, fi_195, gh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * fi_194[k];

        t_195[k] = pa_x[k] * fi_195[k];

        t_196[k] = pa_z[k] * fi_84[k];

        t_197[k] = pa_z[k] * fi_85[k];

        t_198[k] = f_2 * fh_63[k]
                   + pb_z[k] * gh_147[k];

        t_199[k] = pa_z[k] * fi_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_x, pa_z, pb_y, pb_z, fh_66, fh_86, \
                         fh_152, fi_90, fi_201, gh_149, gh_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_3 * fh_86[k]
                   + pb_y[k] * gh_149[k];

        t_201[k] = f_0 * fh_152[k]
                   + pa_x[k] * fi_201[k];

        t_202[k] = pa_z[k] * fi_90[k];

        t_203[k] = f_2 * fh_66[k]
                   + pb_z[k] * gh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, pa_z, pb_y, pb_z, fh_69, fh_89, \
                         fh_156, fi_94, fi_205, gh_152, gh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_3 * fh_89[k]
                   + pb_y[k] * gh_152[k];

        t_205[k] = f_4 * fh_156[k]
                   + pa_x[k] * fi_205[k];

        t_206[k] = pa_z[k] * fi_94[k];

        t_207[k] = f_2 * fh_69[k]
                   + pb_z[k] * gh_153[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pa_z, pb_y, fh_93, fh_159, fh_161, \
                         fi_99, fi_208, fi_210, gh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * fh_159[k]
                   + pa_x[k] * fi_208[k];

        t_209[k] = f_3 * fh_93[k]
                   + pb_y[k] * gh_156[k];

        t_210[k] = f_3 * fh_161[k]
                   + pa_x[k] * fi_210[k];

        t_211[k] = pa_z[k] * fi_99[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pb_x, fh_163, fh_164, fh_165, \
                         fh_166, fh_167, gh_163, gh_164, gh_165, gh_166, \
                         gh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_2 * fh_163[k]
                   + pb_x[k] * gh_163[k];

        t_213[k] = f_2 * fh_164[k]
                   + pb_x[k] * gh_164[k];

        t_214[k] = f_2 * fh_165[k]
                   + pb_x[k] * gh_165[k];

        t_215[k] = f_2 * fh_166[k]
                   + pb_x[k] * gh_166[k];

        t_216[k] = f_2 * fh_167[k]
                   + pb_x[k] * gh_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, t_222, t_223, pa_x, fi_217, \
                         fi_218, fi_219, fi_220, fi_221, fi_222, \
                         fi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_x[k] * fi_217[k];

        t_218[k] = pa_x[k] * fi_218[k];

        t_219[k] = pa_x[k] * fi_219[k];

        t_220[k] = pa_x[k] * fi_220[k];

        t_221[k] = pa_x[k] * fi_221[k];

        t_222[k] = pa_x[k] * fi_222[k];

        t_223[k] = pa_x[k] * fi_223[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_x, pa_y, pb_y, fh_105, fh_107, \
                         fh_171, fi_140, fi_142, fi_227, gh_168, \
                         gh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * fi_140[k];

        t_225[k] = f_2 * fh_105[k]
                   + pb_y[k] * gh_168[k];

        t_226[k] = pa_y[k] * fi_142[k];

        t_227[k] = f_0 * fh_171[k]
                   + pa_x[k] * fi_227[k];

        t_228[k] = f_2 * fh_107[k]
                   + pb_y[k] * gh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_x, pa_y, pb_y, pb_z, fh_87, fh_110, \
                         fh_174, fi_145, fi_230, gh_171, gh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * fi_145[k];

        t_230[k] = f_4 * fh_174[k]
                   + pa_x[k] * fi_230[k];

        t_231[k] = f_3 * fh_87[k]
                   + pb_z[k] * gh_171[k];

        t_232[k] = f_2 * fh_110[k]
                   + pb_y[k] * gh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pa_y, pb_z, fh_90, fh_178, fh_180, \
                         fi_149, fi_234, fi_236, gh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * fi_149[k];

        t_234[k] = f_3 * fh_178[k]
                   + pa_x[k] * fi_234[k];

        t_235[k] = f_3 * fh_90[k]
                   + pb_z[k] * gh_174[k];

        t_236[k] = f_3 * fh_180[k]
                   + pa_x[k] * fi_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_y, pb_x, pb_y, fh_114, fh_183, fh_184, \
                         fi_154, gh_177, gh_183, gh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_2 * fh_114[k]
                   + pb_y[k] * gh_177[k];

        t_238[k] = pa_y[k] * fi_154[k];

        t_239[k] = f_2 * fh_183[k]
                   + pb_x[k] * gh_183[k];

        t_240[k] = f_2 * fh_184[k]
                   + pb_x[k] * gh_184[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pa_x, pa_y, pb_x, fh_185, fh_186, \
                         fh_187, fi_160, fi_245, gh_185, gh_186, \
                         gh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_2 * fh_185[k]
                   + pb_x[k] * gh_185[k];

        t_242[k] = f_2 * fh_186[k]
                   + pb_x[k] * gh_186[k];

        t_243[k] = f_2 * fh_187[k]
                   + pb_x[k] * gh_187[k];

        t_244[k] = pa_y[k] * fi_160[k];

        t_245[k] = pa_x[k] * fi_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, pa_x, fh_189, \
                         fi_246, fi_247, fi_248, fi_249, fi_250, fi_251, \
                         fi_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_x[k] * fi_246[k];

        t_247[k] = pa_x[k] * fi_247[k];

        t_248[k] = pa_x[k] * fi_248[k];

        t_249[k] = pa_x[k] * fi_249[k];

        t_250[k] = pa_x[k] * fi_250[k];

        t_251[k] = pa_x[k] * fi_251[k];

        t_252[k] = f_5 * fh_189[k]
                   + pa_x[k] * fi_252[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, pa_x, pb_y, pb_z, fh_105, fh_194, \
                         fi_257, gg_135, gh_189, gh_190, gh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pb_y[k] * gh_189[k];

        t_254[k] = f_4 * fh_105[k]
                   + pb_z[k] * gh_189[k];

        t_255[k] = f_2 * gg_135[k]
                   + pb_y[k] * gh_190[k];

        t_256[k] = pb_y[k] * gh_191[k];

        t_257[k] = f_0 * fh_194[k]
                   + pa_x[k] * fi_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pa_x, pb_y, fh_198, fi_261, \
                         gg_136, gg_137, gg_138, gh_192, gh_193, gh_194, \
                         gh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_3 * gg_136[k]
                   + pb_y[k] * gh_192[k];

        t_259[k] = f_2 * gg_137[k]
                   + pb_y[k] * gh_193[k];

        t_260[k] = pb_y[k] * gh_194[k];

        t_261[k] = f_4 * fh_198[k]
                   + pa_x[k] * fi_261[k];

        t_262[k] = f_4 * gg_138[k]
                   + pb_y[k] * gh_195[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pa_x, pb_y, fh_203, fi_266, gg_139, \
                         gg_140, gh_196, gh_197, gh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_3 * gg_139[k]
                   + pb_y[k] * gh_196[k];

        t_264[k] = f_2 * gg_140[k]
                   + pb_y[k] * gh_197[k];

        t_265[k] = pb_y[k] * gh_198[k];

        t_266[k] = f_3 * fh_203[k]
                   + pa_x[k] * fi_266[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pb_x, pb_y, fh_204, fh_205, \
                         fh_206, fh_207, gh_203, gh_204, gh_205, gh_206, \
                         gh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_2 * fh_204[k]
                   + pb_x[k] * gh_204[k];

        t_268[k] = f_2 * fh_205[k]
                   + pb_x[k] * gh_205[k];

        t_269[k] = f_2 * fh_206[k]
                   + pb_x[k] * gh_206[k];

        t_270[k] = f_2 * fh_207[k]
                   + pb_x[k] * gh_207[k];

        t_271[k] = pb_y[k] * gh_203[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, t_277, pa_x, pb_x, fh_209, fi_273, \
                         fi_274, fi_275, fi_276, fi_277, gh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_2 * fh_209[k]
                   + pb_x[k] * gh_209[k];

        t_273[k] = pa_x[k] * fi_273[k];

        t_274[k] = pa_x[k] * fi_274[k];

        t_275[k] = pa_x[k] * fi_275[k];

        t_276[k] = pa_x[k] * fi_276[k];

        t_277[k] = pa_x[k] * fi_277[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, pa_x, pb_x, pb_y, pb_z, fi_279, \
                         gg_150, gg_151, gh_209, gh_210, gh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pb_y[k] * gh_209[k];

        t_279[k] = pa_x[k] * fi_279[k];

        t_280[k] = f_1 * gg_150[k]
                   + pb_x[k] * gh_210[k];

        t_281[k] = f_0 * gg_151[k]
                   + pb_x[k] * gh_211[k];

        t_282[k] = pb_z[k] * gh_210[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pb_x, pb_z, gg_153, gg_155, \
                         gg_156, gh_211, gh_213, gh_215, gh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_4 * gg_153[k]
                   + pb_x[k] * gh_213[k];

        t_284[k] = pb_z[k] * gh_211[k];

        t_285[k] = f_4 * gg_155[k]
                   + pb_x[k] * gh_215[k];

        t_286[k] = f_3 * gg_156[k]
                   + pb_x[k] * gh_216[k];

        t_287[k] = pb_z[k] * gh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pb_x, pb_z, gg_158, gg_159, \
                         gg_160, gg_162, gh_216, gh_218, gh_219, gh_220, \
                         gh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_3 * gg_158[k]
                   + pb_x[k] * gh_218[k];

        t_289[k] = f_3 * gg_159[k]
                   + pb_x[k] * gh_219[k];

        t_290[k] = f_2 * gg_160[k]
                   + pb_x[k] * gh_220[k];

        t_291[k] = pb_z[k] * gh_216[k];

        t_292[k] = f_2 * gg_162[k]
                   + pb_x[k] * gh_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, t_298, pb_x, gg_163, gg_164, \
                         gh_223, gh_224, gh_225, gh_226, gh_227, \
                         gh_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_2 * gg_163[k]
                   + pb_x[k] * gh_223[k];

        t_294[k] = f_2 * gg_164[k]
                   + pb_x[k] * gh_224[k];

        t_295[k] = pb_x[k] * gh_225[k];

        t_296[k] = pb_x[k] * gh_226[k];

        t_297[k] = pb_x[k] * gh_227[k];

        t_298[k] = pb_x[k] * gh_228[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pb_x, pb_y, pb_z, fh_141, gg_160, \
                         gh_225, gh_226, gh_229, gh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pb_x[k] * gh_229[k];

        t_300[k] = pb_x[k] * gh_230[k];

        t_301[k] = f_0 * fh_141[k]
                   + f_1 * gg_160[k]
                   + pb_y[k] * gh_225[k];

        t_302[k] = pb_z[k] * gh_225[k];

        t_303[k] = f_2 * gg_160[k]
                   + pb_z[k] * gh_226[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pa_z, pb_y, pb_z, fh_146, fi_168, \
                         gg_161, gg_162, gg_164, gh_227, gh_228, \
                         gh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * gg_161[k]
                   + pb_z[k] * gh_227[k];

        t_305[k] = f_4 * gg_162[k]
                   + pb_z[k] * gh_228[k];

        t_306[k] = f_0 * fh_146[k]
                   + pb_y[k] * gh_230[k];

        t_307[k] = f_1 * gg_164[k]
                   + pb_z[k] * gh_230[k];

        t_308[k] = pa_z[k] * fi_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pa_z, pb_x, fi_169, fi_171, \
                         gg_167, gg_169, gg_170, gh_233, gh_235, \
                         gh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_z[k] * fi_169[k];

        t_310[k] = f_0 * gg_167[k]
                   + pb_x[k] * gh_233[k];

        t_311[k] = pa_z[k] * fi_171[k];

        t_312[k] = f_4 * gg_169[k]
                   + pb_x[k] * gh_235[k];

        t_313[k] = f_4 * gg_170[k]
                   + pb_x[k] * gh_236[k];
    }
}

static auto
compute_prim_gi_overlap_0_piece2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t di, const size_t fh,
                                 const size_t fi, const size_t gg, const size_t gh,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_105 = buffer.data(di + 105);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_167 = buffer.data(di + 167);

    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_279 = buffer.data(fi + 279);

    const auto *gg_172 = buffer.data(gg + 172);
    const auto *gg_173 = buffer.data(gg + 173);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_196 = buffer.data(gg + 196);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_199 = buffer.data(gg + 199);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_202 = buffer.data(gg + 202);
    const auto *gg_203 = buffer.data(gg + 203);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, pa_z, pb_x, fi_174, fi_178, \
                         gg_172, gg_173, gg_174, gh_238, gh_239, \
                         gh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_z[k] * fi_174[k];

        t_315[k] = f_3 * gg_172[k]
                   + pb_x[k] * gh_238[k];

        t_316[k] = f_3 * gg_173[k]
                   + pb_x[k] * gh_239[k];

        t_317[k] = f_3 * gg_174[k]
                   + pb_x[k] * gh_240[k];

        t_318[k] = pa_z[k] * fi_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, pb_x, gg_176, gg_177, gg_178, \
                         gg_179, gh_242, gh_243, gh_244, gh_245, \
                         gh_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_2 * gg_176[k]
                   + pb_x[k] * gh_242[k];

        t_320[k] = f_2 * gg_177[k]
                   + pb_x[k] * gh_243[k];

        t_321[k] = f_2 * gg_178[k]
                   + pb_x[k] * gh_244[k];

        t_322[k] = f_2 * gg_179[k]
                   + pb_x[k] * gh_245[k];

        t_323[k] = pb_x[k] * gh_246[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pa_z, pb_x, fi_189, gh_247, \
                         gh_248, gh_249, gh_250, gh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pb_x[k] * gh_247[k];

        t_325[k] = pb_x[k] * gh_248[k];

        t_326[k] = pb_x[k] * gh_249[k];

        t_327[k] = pb_x[k] * gh_250[k];

        t_328[k] = pb_x[k] * gh_251[k];

        t_329[k] = pa_z[k] * fi_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pa_z, pb_z, fh_141, fh_142, fh_143, \
                         fh_144, fi_191, fi_192, fi_193, gh_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_2 * fh_141[k]
                   + pb_z[k] * gh_246[k];

        t_331[k] = f_3 * fh_142[k]
                   + pa_z[k] * fi_191[k];

        t_332[k] = f_4 * fh_143[k]
                   + pa_z[k] * fi_192[k];

        t_333[k] = f_0 * fh_144[k]
                   + pa_z[k] * fi_193[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_y, pb_x, pb_y, di_139, fh_167, fi_223, \
                         gg_180, gg_181, gh_251, gh_252, gh_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_4 * fh_167[k]
                   + pb_y[k] * gh_251[k];

        t_335[k] = f_3 * di_139[k]
                   + pa_y[k] * fi_223[k];

        t_336[k] = f_1 * gg_180[k]
                   + pb_x[k] * gh_252[k];

        t_337[k] = f_0 * gg_181[k]
                   + pb_x[k] * gh_253[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, pb_x, gg_182, gg_183, gg_184, \
                         gg_185, gg_186, gh_254, gh_255, gh_256, gh_257, \
                         gh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_0 * gg_182[k]
                   + pb_x[k] * gh_254[k];

        t_339[k] = f_4 * gg_183[k]
                   + pb_x[k] * gh_255[k];

        t_340[k] = f_4 * gg_184[k]
                   + pb_x[k] * gh_256[k];

        t_341[k] = f_4 * gg_185[k]
                   + pb_x[k] * gh_257[k];

        t_342[k] = f_3 * gg_186[k]
                   + pb_x[k] * gh_258[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, pb_x, gg_187, gg_188, gg_189, \
                         gg_190, gg_191, gh_259, gh_260, gh_261, gh_262, \
                         gh_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * gg_187[k]
                   + pb_x[k] * gh_259[k];

        t_344[k] = f_3 * gg_188[k]
                   + pb_x[k] * gh_260[k];

        t_345[k] = f_3 * gg_189[k]
                   + pb_x[k] * gh_261[k];

        t_346[k] = f_2 * gg_190[k]
                   + pb_x[k] * gh_262[k];

        t_347[k] = f_2 * gg_191[k]
                   + pb_x[k] * gh_263[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, t_353, pb_x, gg_192, gg_193, \
                         gg_194, gh_264, gh_265, gh_266, gh_267, gh_268, \
                         gh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_2 * gg_192[k]
                   + pb_x[k] * gh_264[k];

        t_349[k] = f_2 * gg_193[k]
                   + pb_x[k] * gh_265[k];

        t_350[k] = f_2 * gg_194[k]
                   + pb_x[k] * gh_266[k];

        t_351[k] = pb_x[k] * gh_267[k];

        t_352[k] = pb_x[k] * gh_268[k];

        t_353[k] = pb_x[k] * gh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pa_z, pb_x, pb_z, di_105, fh_162, \
                         fi_217, gh_267, gh_270, gh_271, gh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * gh_270[k];

        t_355[k] = pb_x[k] * gh_271[k];

        t_356[k] = pb_x[k] * gh_272[k];

        t_357[k] = f_2 * di_105[k]
                   + pa_z[k] * fi_217[k];

        t_358[k] = f_3 * fh_162[k]
                   + pb_z[k] * gh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, fh_185, fh_186, fh_187, fh_188, \
                         gg_192, gg_193, gg_194, gh_269, gh_270, gh_271, \
                         gh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_3 * fh_185[k]
                   + f_4 * gg_192[k]
                   + pb_y[k] * gh_269[k];

        t_360[k] = f_3 * fh_186[k]
                   + f_3 * gg_193[k]
                   + pb_y[k] * gh_270[k];

        t_361[k] = f_3 * fh_187[k]
                   + f_2 * gg_194[k]
                   + pb_y[k] * gh_271[k];

        t_362[k] = f_3 * fh_188[k]
                   + pb_y[k] * gh_272[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, pa_y, pb_x, di_167, fi_251, \
                         fi_252, fi_254, gg_196, gg_198, gh_274, \
                         gh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_2 * di_167[k]
                   + pa_y[k] * fi_251[k];

        t_364[k] = pa_y[k] * fi_252[k];

        t_365[k] = f_0 * gg_196[k]
                   + pb_x[k] * gh_274[k];

        t_366[k] = pa_y[k] * fi_254[k];

        t_367[k] = f_4 * gg_198[k]
                   + pb_x[k] * gh_276[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, pa_y, pb_x, fi_257, gg_199, \
                         gg_201, gg_202, gg_203, gh_277, gh_279, gh_280, \
                         gh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_4 * gg_199[k]
                   + pb_x[k] * gh_277[k];

        t_369[k] = pa_y[k] * fi_257[k];

        t_370[k] = f_3 * gg_201[k]
                   + pb_x[k] * gh_279[k];

        t_371[k] = f_3 * gg_202[k]
                   + pb_x[k] * gh_280[k];

        t_372[k] = f_3 * gg_203[k]
                   + pb_x[k] * gh_281[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, pa_y, pb_x, fi_261, gg_205, \
                         gg_206, gg_207, gg_208, gh_283, gh_284, gh_285, \
                         gh_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = pa_y[k] * fi_261[k];

        t_374[k] = f_2 * gg_205[k]
                   + pb_x[k] * gh_283[k];

        t_375[k] = f_2 * gg_206[k]
                   + pb_x[k] * gh_284[k];

        t_376[k] = f_2 * gg_207[k]
                   + pb_x[k] * gh_285[k];

        t_377[k] = f_2 * gg_208[k]
                   + pb_x[k] * gh_286[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, t_383, t_384, pa_y, pb_x, fi_266, \
                         gh_288, gh_289, gh_290, gh_291, gh_292, \
                         gh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * fi_266[k];

        t_379[k] = pb_x[k] * gh_288[k];

        t_380[k] = pb_x[k] * gh_289[k];

        t_381[k] = pb_x[k] * gh_290[k];

        t_382[k] = pb_x[k] * gh_291[k];

        t_383[k] = pb_x[k] * gh_292[k];

        t_384[k] = pb_x[k] * gh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pb_z, fh_183, fh_204, fh_206, \
                         fh_207, fi_273, fi_275, fi_276, gh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_5 * fh_204[k]
                   + pa_y[k] * fi_273[k];

        t_386[k] = f_4 * fh_183[k]
                   + pb_z[k] * gh_288[k];

        t_387[k] = f_0 * fh_206[k]
                   + pa_y[k] * fi_275[k];

        t_388[k] = f_4 * fh_207[k]
                   + pa_y[k] * fi_276[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pa_y, pb_x, pb_y, fh_208, fh_209, \
                         fi_277, fi_279, gg_210, gh_293, gh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_3 * fh_208[k]
                   + pa_y[k] * fi_277[k];

        t_390[k] = f_2 * fh_209[k]
                   + pb_y[k] * gh_293[k];

        t_391[k] = pa_y[k] * fi_279[k];

        t_392[k] = f_1 * gg_210[k]
                   + pb_x[k] * gh_294[k];

        t_393[k] = pb_y[k] * gh_294[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pb_x, pb_y, gg_212, gg_213, \
                         gg_215, gg_216, gh_296, gh_297, gh_299, \
                         gh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_0 * gg_212[k]
                   + pb_x[k] * gh_296[k];

        t_395[k] = f_4 * gg_213[k]
                   + pb_x[k] * gh_297[k];

        t_396[k] = pb_y[k] * gh_296[k];

        t_397[k] = f_4 * gg_215[k]
                   + pb_x[k] * gh_299[k];

        t_398[k] = f_3 * gg_216[k]
                   + pb_x[k] * gh_300[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, pb_x, pb_y, gg_217, gg_219, \
                         gg_220, gg_221, gh_299, gh_301, gh_303, gh_304, \
                         gh_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_3 * gg_217[k]
                   + pb_x[k] * gh_301[k];

        t_400[k] = pb_y[k] * gh_299[k];

        t_401[k] = f_3 * gg_219[k]
                   + pb_x[k] * gh_303[k];

        t_402[k] = f_2 * gg_220[k]
                   + pb_x[k] * gh_304[k];

        t_403[k] = f_2 * gg_221[k]
                   + pb_x[k] * gh_305[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, t_409, pb_x, pb_y, gg_222, gg_224, \
                         gh_303, gh_306, gh_308, gh_309, gh_310, \
                         gh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_2 * gg_222[k]
                   + pb_x[k] * gh_306[k];

        t_405[k] = pb_y[k] * gh_303[k];

        t_406[k] = f_2 * gg_224[k]
                   + pb_x[k] * gh_308[k];

        t_407[k] = pb_x[k] * gh_309[k];

        t_408[k] = pb_x[k] * gh_310[k];

        t_409[k] = pb_x[k] * gh_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pb_x, pb_y, gg_220, gg_221, \
                         gh_309, gh_310, gh_312, gh_313, gh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * gh_312[k];

        t_411[k] = pb_x[k] * gh_313[k];

        t_412[k] = pb_x[k] * gh_314[k];

        t_413[k] = f_1 * gg_220[k]
                   + pb_y[k] * gh_309[k];

        t_414[k] = f_0 * gg_221[k]
                   + pb_y[k] * gh_310[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, pb_y, pb_z, fh_209, gg_222, \
                         gg_223, gg_224, gh_311, gh_312, gh_313, \
                         gh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_4 * gg_222[k]
                   + pb_y[k] * gh_311[k];

        t_416[k] = f_3 * gg_223[k]
                   + pb_y[k] * gh_312[k];

        t_417[k] = f_2 * gg_224[k]
                   + pb_y[k] * gh_313[k];

        t_418[k] = pb_y[k] * gh_314[k];

        t_419[k] = f_0 * fh_209[k]
                   + f_1 * gg_224[k]
                   + pb_z[k] * gh_314[k];
    }
}

auto
compute_prim_gi_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t di, const size_t fh, const size_t fi,
                          const size_t gg, const size_t gh, const size_t ncols,
                          const double p) -> void
{
    compute_prim_gi_overlap_0_piece0(buffer, target, pa, pb, di, fh, fi, gg, gh, ncols, p);

    compute_prim_gi_overlap_0_piece1(buffer, target, pa, pb, di, fh, fi, gg, gh, ncols, p);

    compute_prim_gi_overlap_0_piece2(buffer, target, pa, pb, di, fh, fi, gg, gh, ncols, p);
}

}  // namespace simdovl
