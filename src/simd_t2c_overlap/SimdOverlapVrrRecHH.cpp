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


#include "SimdOverlapVrrRecHH.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

static auto
compute_prim_hh_overlap_0_piece0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fh, const size_t gg,
                                 const size_t gh, const size_t hf, const size_t hg,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_141 = buffer.data(fh + 141);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_141 = buffer.data(gh + 141);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_69 = buffer.data(hf + 69);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_107 = buffer.data(hg + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gg_0, hf_0, hg_0, \
                         hg_1, hg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_0[k]
                 + f_1 * hf_0[k]
                 + pb_x[k] * hg_0[k];

        t_1[k] = pb_y[k] * hg_0[k];

        t_2[k] = pb_z[k] * hg_0[k];

        t_3[k] = f_2 * hf_0[k]
                 + pb_y[k] * hg_1[k];

        t_4[k] = pb_y[k] * hg_2[k];

        t_5[k] = f_2 * hf_0[k]
                 + pb_z[k] * hg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, gg_10, hf_1, hf_2, \
                         hg_3, hg_5, hg_6, hg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hf_1[k]
                 + pb_y[k] * hg_3[k];

        t_7[k] = pb_z[k] * hg_3[k];

        t_8[k] = pb_y[k] * hg_5[k];

        t_9[k] = f_3 * hf_2[k]
                 + pb_z[k] * hg_5[k];

        t_10[k] = f_0 * gg_10[k]
                  + pb_x[k] * hg_10[k];

        t_11[k] = pb_z[k] * hg_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, gg_12, gg_14, hf_6, \
                         hg_9, hg_10, hg_12, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * gg_12[k]
                  + pb_x[k] * hg_12[k];

        t_13[k] = pb_y[k] * hg_9[k];

        t_14[k] = f_0 * gg_14[k]
                  + pb_x[k] * hg_14[k];

        t_15[k] = f_1 * hf_6[k]
                  + pb_y[k] * hg_10[k];

        t_16[k] = pb_z[k] * hg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_y, pb_z, gh_0, hf_8, hf_9, \
                         hg_12, hg_13, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * hf_8[k]
                  + pb_y[k] * hg_12[k];

        t_18[k] = f_2 * hf_9[k]
                  + pb_y[k] * hg_13[k];

        t_19[k] = pb_y[k] * hg_14[k];

        t_20[k] = f_1 * hf_9[k]
                  + pb_z[k] * hg_14[k];

        t_21[k] = pa_y[k] * gh_0[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, pb_z, gg_0, gg_1, gh_3, \
                         gh_5, hg_15, hg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * gg_0[k]
                  + pb_y[k] * hg_15[k];

        t_23[k] = pb_z[k] * hg_15[k];

        t_24[k] = f_3 * gg_1[k]
                  + pa_y[k] * gh_3[k];

        t_25[k] = pb_z[k] * hg_16[k];

        t_26[k] = pa_y[k] * gh_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pb_y, pb_z, gg_3, gg_5, gh_6, gh_9, \
                         hg_18, hg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * gg_3[k]
                  + pa_y[k] * gh_6[k];

        t_28[k] = pb_z[k] * hg_18[k];

        t_29[k] = f_2 * gg_5[k]
                  + pb_y[k] * hg_20[k];

        t_30[k] = pa_y[k] * gh_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_y, pb_x, pb_z, gg_25, gg_27, gg_28, \
                         gh_14, hg_21, hg_25, hg_27, hg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * gg_25[k]
                  + pb_x[k] * hg_25[k];

        t_32[k] = pb_z[k] * hg_21[k];

        t_33[k] = f_1 * gg_27[k]
                  + pb_x[k] * hg_27[k];

        t_34[k] = f_1 * gg_28[k]
                  + pb_x[k] * hg_28[k];

        t_35[k] = pa_y[k] * gh_14[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_z, fh_36, gh_36, hf_16, hf_17, \
                         hg_25, hg_26, hg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * fh_36[k]
                  + pa_x[k] * gh_36[k];

        t_37[k] = pb_z[k] * hg_25[k];

        t_38[k] = f_2 * hf_16[k]
                  + pb_z[k] * hg_26[k];

        t_39[k] = f_3 * hf_17[k]
                  + pb_z[k] * hg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, pb_z, gg_0, gg_14, \
                         gh_0, gh_20, hg_29, hg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * gg_14[k]
                  + pb_y[k] * hg_29[k];

        t_41[k] = pa_y[k] * gh_20[k];

        t_42[k] = pa_z[k] * gh_0[k];

        t_43[k] = pb_y[k] * hg_30[k];

        t_44[k] = f_2 * gg_0[k]
                  + pb_z[k] * hg_30[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pa_z, pb_y, gg_2, gh_3, gh_5, \
                         gh_6, hf_22, hg_32, hg_34, hg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * gh_3[k];

        t_46[k] = pb_y[k] * hg_32[k];

        t_47[k] = f_3 * gg_2[k]
                  + pa_z[k] * gh_5[k];

        t_48[k] = pa_z[k] * gh_6[k];

        t_49[k] = f_2 * hf_22[k]
                  + pb_y[k] * hg_34[k];

        t_50[k] = pb_y[k] * hg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_z, pb_x, pb_y, gg_5, gg_41, gg_42, \
                         gh_9, gh_10, hg_39, hg_41, hg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_4 * gg_5[k]
                  + pa_z[k] * gh_9[k];

        t_52[k] = pa_z[k] * gh_10[k];

        t_53[k] = f_1 * gg_41[k]
                  + pb_x[k] * hg_41[k];

        t_54[k] = f_1 * gg_42[k]
                  + pb_x[k] * hg_42[k];

        t_55[k] = pb_y[k] * hg_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, pb_y, gg_44, gh_15, hf_27, hf_28, \
                         hg_41, hg_42, hg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * gg_44[k]
                  + pb_x[k] * hg_44[k];

        t_57[k] = pa_z[k] * gh_15[k];

        t_58[k] = f_4 * hf_27[k]
                  + pb_y[k] * hg_41[k];

        t_59[k] = f_3 * hf_28[k]
                  + pb_y[k] * hg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pa_y, pb_y, fh_0, fh_62, gh_21, gh_62, \
                         hf_29, hg_43, hg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * hf_29[k]
                  + pb_y[k] * hg_43[k];

        t_61[k] = pb_y[k] * hg_44[k];

        t_62[k] = f_4 * fh_62[k]
                  + pa_x[k] * gh_62[k];

        t_63[k] = f_2 * fh_0[k]
                  + pa_y[k] * gh_21[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pb_x, pb_y, pb_z, gg_15, gg_48, hf_30, \
                         hf_33, hg_45, hg_46, hg_47, hg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * gg_15[k]
                  + pb_y[k] * hg_45[k];

        t_65[k] = pb_z[k] * hg_45[k];

        t_66[k] = f_4 * gg_48[k]
                  + f_3 * hf_33[k]
                  + pb_x[k] * hg_48[k];

        t_67[k] = pb_z[k] * hg_46[k];

        t_68[k] = f_2 * hf_30[k]
                  + pb_z[k] * hg_47[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, gg_20, gg_51, hf_32, hf_36, \
                         hg_48, hg_50, hg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_4 * gg_51[k]
                  + f_2 * hf_36[k]
                  + pb_x[k] * hg_51[k];

        t_70[k] = pb_z[k] * hg_48[k];

        t_71[k] = f_3 * gg_20[k]
                  + pb_y[k] * hg_50[k];

        t_72[k] = f_3 * hf_32[k]
                  + pb_z[k] * hg_50[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, gg_55, gg_57, gg_58, gg_59, \
                         hg_51, hg_55, hg_57, hg_58, hg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * gg_55[k]
                  + pb_x[k] * hg_55[k];

        t_74[k] = pb_z[k] * hg_51[k];

        t_75[k] = f_4 * gg_57[k]
                  + pb_x[k] * hg_57[k];

        t_76[k] = f_4 * gg_58[k]
                  + pb_x[k] * hg_58[k];

        t_77[k] = f_4 * gg_59[k]
                  + pb_x[k] * hg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, fh_78, gh_78, hf_36, hf_37, \
                         hg_55, hg_56, hg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_3 * fh_78[k]
                  + pa_x[k] * gh_78[k];

        t_79[k] = pb_z[k] * hg_55[k];

        t_80[k] = f_2 * hf_36[k]
                  + pb_z[k] * hg_56[k];

        t_81[k] = f_3 * hf_37[k]
                  + pb_z[k] * hg_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, gg_29, gh_22, \
                         gh_42, gh_44, hf_39, hg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * gg_29[k]
                  + pb_y[k] * hg_59[k];

        t_83[k] = f_1 * hf_39[k]
                  + pb_z[k] * hg_59[k];

        t_84[k] = pa_y[k] * gh_42[k];

        t_85[k] = pa_z[k] * gh_22[k];

        t_86[k] = pa_y[k] * gh_44[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, gg_18, gg_32, \
                         gh_24, gh_27, gh_47, hg_62, hg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * gh_24[k];

        t_88[k] = f_2 * gg_32[k]
                  + pb_y[k] * hg_62[k];

        t_89[k] = pa_y[k] * gh_47[k];

        t_90[k] = pa_z[k] * gh_27[k];

        t_91[k] = f_2 * gg_18[k]
                  + pb_z[k] * hg_63[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, gg_35, gg_71, gh_31, \
                         gh_51, hg_65, hg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * gg_35[k]
                  + pb_y[k] * hg_65[k];

        t_93[k] = pa_y[k] * gh_51[k];

        t_94[k] = pa_z[k] * gh_31[k];

        t_95[k] = f_4 * gg_71[k]
                  + pb_x[k] * hg_71[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pa_z, pb_x, gg_72, gg_73, gh_36, gh_56, \
                         hg_72, hg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * gg_72[k]
                  + pb_x[k] * hg_72[k];

        t_97[k] = f_4 * gg_73[k]
                  + pb_x[k] * hg_73[k];

        t_98[k] = pa_y[k] * gh_56[k];

        t_99[k] = pa_z[k] * gh_36[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, fh_101, fh_102, gg_25, \
                         gg_44, gh_101, gh_102, hg_70, hg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_2 * gg_25[k]
                   + pb_z[k] * hg_70[k];

        t_101[k] = f_3 * fh_101[k]
                   + pa_x[k] * gh_101[k];

        t_102[k] = f_3 * fh_102[k]
                   + pa_x[k] * gh_102[k];

        t_103[k] = f_2 * gg_44[k]
                   + pb_y[k] * hg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_y, pa_z, pb_y, pb_z, fh_0, \
                         gg_30, gh_42, gh_62, hf_50, hg_75, hg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * gh_62[k];

        t_105[k] = f_2 * fh_0[k]
                   + pa_z[k] * gh_42[k];

        t_106[k] = pb_y[k] * hg_75[k];

        t_107[k] = f_3 * gg_30[k]
                   + pb_z[k] * hg_75[k];

        t_108[k] = f_2 * hf_50[k]
                   + pb_y[k] * hg_76[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, gg_80, hf_51, hf_52, \
                         hf_55, hg_77, hg_78, hg_79, hg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_y[k] * hg_77[k];

        t_110[k] = f_4 * gg_80[k]
                   + f_3 * hf_55[k]
                   + pb_x[k] * hg_80[k];

        t_111[k] = f_3 * hf_51[k]
                   + pb_y[k] * hg_78[k];

        t_112[k] = f_2 * hf_52[k]
                   + pb_y[k] * hg_79[k];

        t_113[k] = pb_y[k] * hg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pb_x, pb_y, gg_84, gg_85, gg_86, \
                         gg_87, hf_59, hg_84, hg_85, hg_86, hg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_4 * gg_84[k]
                   + f_2 * hf_59[k]
                   + pb_x[k] * hg_84[k];

        t_115[k] = f_4 * gg_85[k]
                   + pb_x[k] * hg_85[k];

        t_116[k] = f_4 * gg_86[k]
                   + pb_x[k] * hg_86[k];

        t_117[k] = f_4 * gg_87[k]
                   + pb_x[k] * hg_87[k];

        t_118[k] = pb_y[k] * hg_84[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, gg_89, hf_56, hf_57, hf_58, \
                         hg_85, hg_86, hg_87, hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_4 * gg_89[k]
                   + pb_x[k] * hg_89[k];

        t_120[k] = f_1 * hf_56[k]
                   + pb_y[k] * hg_85[k];

        t_121[k] = f_4 * hf_57[k]
                   + pb_y[k] * hg_86[k];

        t_122[k] = f_3 * hf_58[k]
                   + pb_y[k] * hg_87[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pa_y, pb_y, fh_21, fh_125, gh_63, \
                         gh_125, hf_59, hg_88, hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_2 * hf_59[k]
                   + pb_y[k] * hg_88[k];

        t_124[k] = pb_y[k] * hg_89[k];

        t_125[k] = f_3 * fh_125[k]
                   + pa_x[k] * gh_125[k];

        t_126[k] = f_3 * fh_21[k]
                   + pa_y[k] * gh_63[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pb_x, pb_y, pb_z, gg_45, gg_93, \
                         hf_60, hf_63, hg_90, hg_91, hg_92, hg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_4 * gg_45[k]
                   + pb_y[k] * hg_90[k];

        t_128[k] = pb_z[k] * hg_90[k];

        t_129[k] = f_3 * gg_93[k]
                   + f_3 * hf_63[k]
                   + pb_x[k] * hg_93[k];

        t_130[k] = pb_z[k] * hg_91[k];

        t_131[k] = f_2 * hf_60[k]
                   + pb_z[k] * hg_92[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_y, pb_z, gg_50, gg_96, hf_62, \
                         hf_66, hg_93, hg_95, hg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_3 * gg_96[k]
                   + f_2 * hf_66[k]
                   + pb_x[k] * hg_96[k];

        t_133[k] = pb_z[k] * hg_93[k];

        t_134[k] = f_4 * gg_50[k]
                   + pb_y[k] * hg_95[k];

        t_135[k] = f_3 * hf_62[k]
                   + pb_z[k] * hg_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, pb_z, gg_100, gg_102, \
                         gg_103, gg_104, hg_96, hg_100, hg_102, hg_103, \
                         hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_3 * gg_100[k]
                   + pb_x[k] * hg_100[k];

        t_137[k] = pb_z[k] * hg_96[k];

        t_138[k] = f_3 * gg_102[k]
                   + pb_x[k] * hg_102[k];

        t_139[k] = f_3 * gg_103[k]
                   + pb_x[k] * hg_103[k];

        t_140[k] = f_3 * gg_104[k]
                   + pb_x[k] * hg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, pb_z, fh_141, gh_141, hf_66, hf_67, \
                         hg_100, hg_101, hg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_2 * fh_141[k]
                   + pa_x[k] * gh_141[k];

        t_142[k] = pb_z[k] * hg_100[k];

        t_143[k] = f_2 * hf_66[k]
                   + pb_z[k] * hg_101[k];

        t_144[k] = f_3 * hf_67[k]
                   + pb_z[k] * hg_102[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, gg_45, gg_59, \
                         gh_63, gh_64, hf_69, hg_104, hg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_4 * gg_59[k]
                   + pb_y[k] * hg_104[k];

        t_146[k] = f_1 * hf_69[k]
                   + pb_z[k] * hg_104[k];

        t_147[k] = pa_z[k] * gh_63[k];

        t_148[k] = pa_z[k] * gh_64[k];

        t_149[k] = f_2 * gg_45[k]
                   + pb_z[k] * hg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pa_z, pb_y, fh_47, gg_62, gh_66, \
                         gh_69, gh_89, hg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * gh_66[k];

        t_151[k] = f_3 * gg_62[k]
                   + pb_y[k] * hg_107[k];

        t_152[k] = f_2 * fh_47[k]
                   + pa_y[k] * gh_89[k];

        t_153[k] = pa_z[k] * gh_69[k];
    }
}

static auto
compute_prim_hh_overlap_0_piece1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fh, const size_t gg,
                                 const size_t gh, const size_t hf, const size_t hg,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_219 = buffer.data(gg + 219);

    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_303 = buffer.data(gh + 303);

    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_141 = buffer.data(hf + 141);
    const auto *hf_142 = buffer.data(hf + 142);

    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_211 = buffer.data(hg + 211);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_213 = buffer.data(hg + 213);
    const auto *hg_214 = buffer.data(hg + 214);
    const auto *hg_215 = buffer.data(hg + 215);

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_y, pa_z, pb_y, pb_z, fh_51, gg_48, \
                         gg_65, gh_73, gh_93, hg_108, hg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_2 * gg_48[k]
                   + pb_z[k] * hg_108[k];

        t_155[k] = f_3 * gg_65[k]
                   + pb_y[k] * hg_110[k];

        t_156[k] = f_2 * fh_51[k]
                   + pa_y[k] * gh_93[k];

        t_157[k] = pa_z[k] * gh_73[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_x, gg_116, gg_117, \
                         gg_118, gg_119, gh_78, hg_116, hg_117, hg_118, \
                         hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * gg_116[k]
                   + pb_x[k] * hg_116[k];

        t_159[k] = f_3 * gg_117[k]
                   + pb_x[k] * hg_117[k];

        t_160[k] = f_3 * gg_118[k]
                   + pb_x[k] * hg_118[k];

        t_161[k] = f_3 * gg_119[k]
                   + pb_x[k] * hg_119[k];

        t_162[k] = pa_z[k] * gh_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_x, pb_y, pb_z, fh_164, fh_165, gg_55, \
                         gg_74, gh_164, gh_165, hg_115, hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_2 * gg_55[k]
                   + pb_z[k] * hg_115[k];

        t_164[k] = f_2 * fh_164[k]
                   + pa_x[k] * gh_164[k];

        t_165[k] = f_2 * fh_165[k]
                   + pa_x[k] * gh_165[k];

        t_166[k] = f_3 * gg_74[k]
                   + pb_y[k] * hg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_x, pa_y, pb_y, fh_167, gg_75, \
                         gg_76, gh_105, gh_107, gh_108, gh_167, \
                         hg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_2 * fh_167[k]
                   + pa_x[k] * gh_167[k];

        t_168[k] = pa_y[k] * gh_105[k];

        t_169[k] = f_2 * gg_75[k]
                   + pb_y[k] * hg_120[k];

        t_170[k] = pa_y[k] * gh_107[k];

        t_171[k] = f_3 * gg_76[k]
                   + pa_y[k] * gh_108[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_y, pb_z, gg_63, gg_77, gg_78, \
                         gh_110, gh_111, hg_122, hg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_2 * gg_77[k]
                   + pb_y[k] * hg_122[k];

        t_173[k] = pa_y[k] * gh_110[k];

        t_174[k] = f_4 * gg_78[k]
                   + pa_y[k] * gh_111[k];

        t_175[k] = f_3 * gg_63[k]
                   + pb_z[k] * hg_123[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, pb_y, gg_80, gg_130, gg_131, \
                         gh_114, hg_125, hg_130, hg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_2 * gg_80[k]
                   + pb_y[k] * hg_125[k];

        t_177[k] = pa_y[k] * gh_114[k];

        t_178[k] = f_3 * gg_130[k]
                   + pb_x[k] * hg_130[k];

        t_179[k] = f_3 * gg_131[k]
                   + pb_x[k] * hg_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pa_y, pb_x, fh_183, gg_132, gg_133, \
                         gh_119, gh_183, hg_132, hg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_3 * gg_132[k]
                   + pb_x[k] * hg_132[k];

        t_181[k] = f_3 * gg_133[k]
                   + pb_x[k] * hg_133[k];

        t_182[k] = pa_y[k] * gh_119[k];

        t_183[k] = f_2 * fh_183[k]
                   + pa_x[k] * gh_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pb_y, pb_z, fh_185, fh_186, gg_70, \
                         gg_89, gh_185, gh_186, hg_130, hg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * gg_70[k]
                   + pb_z[k] * hg_130[k];

        t_185[k] = f_2 * fh_185[k]
                   + pa_x[k] * gh_185[k];

        t_186[k] = f_2 * fh_186[k]
                   + pa_x[k] * gh_186[k];

        t_187[k] = f_2 * gg_89[k]
                   + pb_y[k] * hg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pa_y, pa_z, pb_y, pb_z, fh_42, \
                         gg_75, gh_105, gh_125, hf_90, hg_135, hg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * gh_125[k];

        t_189[k] = f_3 * fh_42[k]
                   + pa_z[k] * gh_105[k];

        t_190[k] = pb_y[k] * hg_135[k];

        t_191[k] = f_4 * gg_75[k]
                   + pb_z[k] * hg_135[k];

        t_192[k] = f_2 * hf_90[k]
                   + pb_y[k] * hg_136[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pb_x, pb_y, gg_140, hf_91, hf_92, \
                         hf_95, hg_137, hg_138, hg_139, hg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_y[k] * hg_137[k];

        t_194[k] = f_3 * gg_140[k]
                   + f_3 * hf_95[k]
                   + pb_x[k] * hg_140[k];

        t_195[k] = f_3 * hf_91[k]
                   + pb_y[k] * hg_138[k];

        t_196[k] = f_2 * hf_92[k]
                   + pb_y[k] * hg_139[k];

        t_197[k] = pb_y[k] * hg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pb_x, pb_y, gg_144, gg_145, \
                         gg_146, gg_147, hf_99, hg_144, hg_145, hg_146, \
                         hg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * gg_144[k]
                   + f_2 * hf_99[k]
                   + pb_x[k] * hg_144[k];

        t_199[k] = f_3 * gg_145[k]
                   + pb_x[k] * hg_145[k];

        t_200[k] = f_3 * gg_146[k]
                   + pb_x[k] * hg_146[k];

        t_201[k] = f_3 * gg_147[k]
                   + pb_x[k] * hg_147[k];

        t_202[k] = pb_y[k] * hg_144[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pb_x, pb_y, gg_149, hf_96, hf_97, hf_98, \
                         hg_145, hg_146, hg_147, hg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * gg_149[k]
                   + pb_x[k] * hg_149[k];

        t_204[k] = f_1 * hf_96[k]
                   + pb_y[k] * hg_145[k];

        t_205[k] = f_4 * hf_97[k]
                   + pb_y[k] * hg_146[k];

        t_206[k] = f_3 * hf_98[k]
                   + pb_y[k] * hg_147[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_x, pb_y, fh_209, gg_90, gg_150, \
                         gh_209, gh_210, hf_99, hg_148, hg_149, \
                         hg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_2 * hf_99[k]
                   + pb_y[k] * hg_148[k];

        t_208[k] = pb_y[k] * hg_149[k];

        t_209[k] = f_2 * fh_209[k]
                   + pa_x[k] * gh_209[k];

        t_210[k] = f_0 * gg_150[k]
                   + pa_x[k] * gh_210[k];

        t_211[k] = f_1 * gg_90[k]
                   + pb_y[k] * hg_150[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pb_z, gg_153, gg_156, \
                         gh_213, gh_216, hf_100, hg_150, hg_151, \
                         hg_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_z[k] * hg_150[k];

        t_213[k] = f_4 * gg_153[k]
                   + pa_x[k] * gh_213[k];

        t_214[k] = pb_z[k] * hg_151[k];

        t_215[k] = f_2 * hf_100[k]
                   + pb_z[k] * hg_152[k];

        t_216[k] = f_3 * gg_156[k]
                   + pa_x[k] * gh_216[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pb_x, pb_y, pb_z, gg_95, gg_160, \
                         hf_102, hg_153, hg_155, hg_156, hg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_z[k] * hg_153[k];

        t_218[k] = f_1 * gg_95[k]
                   + pb_y[k] * hg_155[k];

        t_219[k] = f_3 * hf_102[k]
                   + pb_z[k] * hg_155[k];

        t_220[k] = f_2 * gg_160[k]
                   + pb_x[k] * hg_160[k];

        t_221[k] = pb_z[k] * hg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pa_x, pb_x, pb_z, gg_162, gg_163, \
                         gg_164, gh_225, hg_160, hg_162, hg_163, \
                         hg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_2 * gg_162[k]
                   + pb_x[k] * hg_162[k];

        t_223[k] = f_2 * gg_163[k]
                   + pb_x[k] * hg_163[k];

        t_224[k] = f_2 * gg_164[k]
                   + pb_x[k] * hg_164[k];

        t_225[k] = pa_x[k] * gh_225[k];

        t_226[k] = pb_z[k] * hg_160[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, t_232, pa_x, pa_z, gh_126, gh_127, \
                         gh_227, gh_228, gh_229, gh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pa_x[k] * gh_227[k];

        t_228[k] = pa_x[k] * gh_228[k];

        t_229[k] = pa_x[k] * gh_229[k];

        t_230[k] = pa_x[k] * gh_230[k];

        t_231[k] = pa_z[k] * gh_126[k];

        t_232[k] = pa_z[k] * gh_127[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pa_z, pb_y, pb_z, gg_90, gg_107, \
                         gg_170, gh_129, gh_236, hg_165, hg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_2 * gg_90[k]
                   + pb_z[k] * hg_165[k];

        t_234[k] = pa_z[k] * gh_129[k];

        t_235[k] = f_4 * gg_107[k]
                   + pb_y[k] * hg_167[k];

        t_236[k] = f_4 * gg_170[k]
                   + pa_x[k] * gh_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_x, pa_z, pb_y, pb_z, gg_93, gg_110, \
                         gg_174, gh_132, gh_240, hg_168, hg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * gh_132[k];

        t_238[k] = f_2 * gg_93[k]
                   + pb_z[k] * hg_168[k];

        t_239[k] = f_4 * gg_110[k]
                   + pb_y[k] * hg_170[k];

        t_240[k] = f_3 * gg_174[k]
                   + pa_x[k] * gh_240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pa_z, pb_x, gg_176, gg_177, \
                         gg_178, gg_179, gh_136, hg_176, hg_177, hg_178, \
                         hg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = pa_z[k] * gh_136[k];

        t_242[k] = f_2 * gg_176[k]
                   + pb_x[k] * hg_176[k];

        t_243[k] = f_2 * gg_177[k]
                   + pb_x[k] * hg_177[k];

        t_244[k] = f_2 * gg_178[k]
                   + pb_x[k] * hg_178[k];

        t_245[k] = f_2 * gg_179[k]
                   + pb_x[k] * hg_179[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, pa_x, gg_180, \
                         gh_246, gh_247, gh_248, gh_249, gh_250, gh_251, \
                         gh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_x[k] * gh_246[k];

        t_247[k] = pa_x[k] * gh_247[k];

        t_248[k] = pa_x[k] * gh_248[k];

        t_249[k] = pa_x[k] * gh_249[k];

        t_250[k] = pa_x[k] * gh_250[k];

        t_251[k] = pa_x[k] * gh_251[k];

        t_252[k] = f_0 * gg_180[k]
                   + pa_x[k] * gh_252[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_x, pb_y, pb_z, gg_105, gg_120, gg_122, \
                         gg_183, gh_255, hg_180, hg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_3 * gg_120[k]
                   + pb_y[k] * hg_180[k];

        t_254[k] = f_3 * gg_105[k]
                   + pb_z[k] * hg_180[k];

        t_255[k] = f_4 * gg_183[k]
                   + pa_x[k] * gh_255[k];

        t_256[k] = f_3 * gg_122[k]
                   + pb_y[k] * hg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_x, pb_y, pb_z, gg_108, gg_125, gg_185, \
                         gg_186, gh_257, gh_258, hg_183, hg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_4 * gg_185[k]
                   + pa_x[k] * gh_257[k];

        t_258[k] = f_3 * gg_186[k]
                   + pa_x[k] * gh_258[k];

        t_259[k] = f_3 * gg_108[k]
                   + pb_z[k] * hg_183[k];

        t_260[k] = f_3 * gg_125[k]
                   + pb_y[k] * hg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_x, pb_x, gg_189, gg_190, gg_191, \
                         gg_192, gh_261, hg_190, hg_191, hg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * gg_189[k]
                   + pa_x[k] * gh_261[k];

        t_262[k] = f_2 * gg_190[k]
                   + pb_x[k] * hg_190[k];

        t_263[k] = f_2 * gg_191[k]
                   + pb_x[k] * hg_191[k];

        t_264[k] = f_2 * gg_192[k]
                   + pb_x[k] * hg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, pa_x, pb_x, gg_193, gg_194, \
                         gh_267, gh_268, gh_269, gh_270, hg_193, \
                         hg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_2 * gg_193[k]
                   + pb_x[k] * hg_193[k];

        t_266[k] = f_2 * gg_194[k]
                   + pb_x[k] * hg_194[k];

        t_267[k] = pa_x[k] * gh_267[k];

        t_268[k] = pa_x[k] * gh_268[k];

        t_269[k] = pa_x[k] * gh_269[k];

        t_270[k] = pa_x[k] * gh_270[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, pa_x, pa_y, pb_y, gg_135, gh_189, \
                         gh_191, gh_271, gh_272, hg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = pa_x[k] * gh_271[k];

        t_272[k] = pa_x[k] * gh_272[k];

        t_273[k] = pa_y[k] * gh_189[k];

        t_274[k] = f_2 * gg_135[k]
                   + pb_y[k] * hg_195[k];

        t_275[k] = pa_y[k] * gh_191[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pa_y, pb_y, gg_137, gg_198, gg_201, \
                         gh_194, gh_276, gh_279, hg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_4 * gg_198[k]
                   + pa_x[k] * gh_276[k];

        t_277[k] = f_2 * gg_137[k]
                   + pb_y[k] * hg_197[k];

        t_278[k] = pa_y[k] * gh_194[k];

        t_279[k] = f_3 * gg_201[k]
                   + pa_x[k] * gh_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_y, pb_z, gg_123, gg_140, \
                         gg_205, gh_198, hg_198, hg_200, hg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_4 * gg_123[k]
                   + pb_z[k] * hg_198[k];

        t_281[k] = f_2 * gg_140[k]
                   + pb_y[k] * hg_200[k];

        t_282[k] = pa_y[k] * gh_198[k];

        t_283[k] = f_2 * gg_205[k]
                   + pb_x[k] * hg_205[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_x, pa_y, pb_x, gg_206, gg_207, \
                         gg_208, gh_203, gh_288, hg_206, hg_207, \
                         hg_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_2 * gg_206[k]
                   + pb_x[k] * hg_206[k];

        t_285[k] = f_2 * gg_207[k]
                   + pb_x[k] * hg_207[k];

        t_286[k] = f_2 * gg_208[k]
                   + pb_x[k] * hg_208[k];

        t_287[k] = pa_y[k] * gh_203[k];

        t_288[k] = pa_x[k] * gh_288[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, t_294, pa_x, gg_210, gh_289, \
                         gh_290, gh_291, gh_292, gh_293, gh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = pa_x[k] * gh_289[k];

        t_290[k] = pa_x[k] * gh_290[k];

        t_291[k] = pa_x[k] * gh_291[k];

        t_292[k] = pa_x[k] * gh_292[k];

        t_293[k] = pa_x[k] * gh_293[k];

        t_294[k] = f_0 * gg_210[k]
                   + pa_x[k] * gh_294[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, pa_x, pb_y, pb_z, gg_135, gg_215, \
                         gh_299, hf_140, hg_210, hg_211, hg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = pb_y[k] * hg_210[k];

        t_296[k] = f_1 * gg_135[k]
                   + pb_z[k] * hg_210[k];

        t_297[k] = f_2 * hf_140[k]
                   + pb_y[k] * hg_211[k];

        t_298[k] = pb_y[k] * hg_212[k];

        t_299[k] = f_4 * gg_215[k]
                   + pa_x[k] * gh_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pb_y, gg_219, gh_303, hf_141, \
                         hf_142, hg_213, hg_214, hg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_3 * hf_141[k]
                   + pb_y[k] * hg_213[k];

        t_301[k] = f_2 * hf_142[k]
                   + pb_y[k] * hg_214[k];

        t_302[k] = pb_y[k] * hg_215[k];

        t_303[k] = f_3 * gg_219[k]
                   + pa_x[k] * gh_303[k];
    }
}

static auto
compute_prim_hh_overlap_0_piece2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fh, const size_t gg,
                                 const size_t gh, const size_t hf, const size_t hg,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_164 = buffer.data(hf + 164);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_191 = buffer.data(hf + 191);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_194 = buffer.data(hf + 194);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_220 = buffer.data(hg + 220);
    const auto *hg_221 = buffer.data(hg + 221);
    const auto *hg_222 = buffer.data(hg + 222);
    const auto *hg_224 = buffer.data(hg + 224);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_231 = buffer.data(hg + 231);
    const auto *hg_233 = buffer.data(hg + 233);
    const auto *hg_234 = buffer.data(hg + 234);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_236 = buffer.data(hg + 236);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_242 = buffer.data(hg + 242);
    const auto *hg_244 = buffer.data(hg + 244);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_247 = buffer.data(hg + 247);
    const auto *hg_248 = buffer.data(hg + 248);
    const auto *hg_249 = buffer.data(hg + 249);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_256 = buffer.data(hg + 256);
    const auto *hg_257 = buffer.data(hg + 257);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_259 = buffer.data(hg + 259);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_261 = buffer.data(hg + 261);
    const auto *hg_262 = buffer.data(hg + 262);
    const auto *hg_263 = buffer.data(hg + 263);
    const auto *hg_264 = buffer.data(hg + 264);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_271 = buffer.data(hg + 271);
    const auto *hg_272 = buffer.data(hg + 272);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_274 = buffer.data(hg + 274);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_276 = buffer.data(hg + 276);
    const auto *hg_277 = buffer.data(hg + 277);
    const auto *hg_278 = buffer.data(hg + 278);
    const auto *hg_279 = buffer.data(hg + 279);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_286 = buffer.data(hg + 286);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_289 = buffer.data(hg + 289);
    const auto *hg_291 = buffer.data(hg + 291);
    const auto *hg_292 = buffer.data(hg + 292);
    const auto *hg_293 = buffer.data(hg + 293);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_306 = buffer.data(hg + 306);
    const auto *hg_307 = buffer.data(hg + 307);
    const auto *hg_309 = buffer.data(hg + 309);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_313 = buffer.data(hg + 313);
    const auto *hg_314 = buffer.data(hg + 314);

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, gg_220, gg_221, \
                         gg_222, gg_224, hg_219, hg_220, hg_221, hg_222, \
                         hg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_2 * gg_220[k]
                   + pb_x[k] * hg_220[k];

        t_305[k] = f_2 * gg_221[k]
                   + pb_x[k] * hg_221[k];

        t_306[k] = f_2 * gg_222[k]
                   + pb_x[k] * hg_222[k];

        t_307[k] = pb_y[k] * hg_219[k];

        t_308[k] = f_2 * gg_224[k]
                   + pb_x[k] * hg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, t_314, pa_x, pb_y, gh_309, gh_310, \
                         gh_311, gh_312, gh_314, hg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_x[k] * gh_309[k];

        t_310[k] = pa_x[k] * gh_310[k];

        t_311[k] = pa_x[k] * gh_311[k];

        t_312[k] = pa_x[k] * gh_312[k];

        t_313[k] = pb_y[k] * hg_224[k];

        t_314[k] = pa_x[k] * gh_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, pb_x, pb_z, hf_150, hf_151, \
                         hf_153, hf_155, hg_225, hg_226, hg_228, \
                         hg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_1 * hf_150[k]
                   + pb_x[k] * hg_225[k];

        t_316[k] = f_4 * hf_151[k]
                   + pb_x[k] * hg_226[k];

        t_317[k] = pb_z[k] * hg_225[k];

        t_318[k] = f_3 * hf_153[k]
                   + pb_x[k] * hg_228[k];

        t_319[k] = pb_z[k] * hg_226[k];

        t_320[k] = f_3 * hf_155[k]
                   + pb_x[k] * hg_230[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, t_325, pb_x, pb_z, hf_156, hf_158, \
                         hf_159, hg_228, hg_231, hg_233, hg_234, \
                         hg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_2 * hf_156[k]
                   + pb_x[k] * hg_231[k];

        t_322[k] = pb_z[k] * hg_228[k];

        t_323[k] = f_2 * hf_158[k]
                   + pb_x[k] * hg_233[k];

        t_324[k] = f_2 * hf_159[k]
                   + pb_x[k] * hg_234[k];

        t_325[k] = pb_x[k] * hg_235[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, pb_x, pb_y, pb_z, gg_160, \
                         hf_156, hg_235, hg_236, hg_237, hg_238, \
                         hg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = pb_x[k] * hg_236[k];

        t_327[k] = pb_x[k] * hg_237[k];

        t_328[k] = pb_x[k] * hg_238[k];

        t_329[k] = pb_x[k] * hg_239[k];

        t_330[k] = f_0 * gg_160[k]
                   + f_1 * hf_156[k]
                   + pb_y[k] * hg_235[k];

        t_331[k] = pb_z[k] * hg_235[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pa_z, pb_y, pb_z, gg_164, gh_210, \
                         hf_156, hf_157, hf_159, hg_236, hg_237, \
                         hg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_2 * hf_156[k]
                   + pb_z[k] * hg_236[k];

        t_333[k] = f_3 * hf_157[k]
                   + pb_z[k] * hg_237[k];

        t_334[k] = f_0 * gg_164[k]
                   + pb_y[k] * hg_239[k];

        t_335[k] = f_1 * hf_159[k]
                   + pb_z[k] * hg_239[k];

        t_336[k] = pa_z[k] * gh_210[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, pa_z, pb_x, gh_211, gh_213, \
                         hf_162, hf_164, hf_165, hg_242, hg_244, \
                         hg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = pa_z[k] * gh_211[k];

        t_338[k] = f_4 * hf_162[k]
                   + pb_x[k] * hg_242[k];

        t_339[k] = pa_z[k] * gh_213[k];

        t_340[k] = f_3 * hf_164[k]
                   + pb_x[k] * hg_244[k];

        t_341[k] = f_3 * hf_165[k]
                   + pb_x[k] * hg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, pa_z, pb_x, gh_216, hf_167, \
                         hf_168, hf_169, hg_247, hg_248, hg_249, \
                         hg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * gh_216[k];

        t_343[k] = f_2 * hf_167[k]
                   + pb_x[k] * hg_247[k];

        t_344[k] = f_2 * hf_168[k]
                   + pb_x[k] * hg_248[k];

        t_345[k] = f_2 * hf_169[k]
                   + pb_x[k] * hg_249[k];

        t_346[k] = pb_x[k] * hg_250[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, t_352, pa_z, pb_x, pb_z, gg_160, \
                         gh_225, hg_250, hg_251, hg_252, hg_253, \
                         hg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = pb_x[k] * hg_251[k];

        t_348[k] = pb_x[k] * hg_252[k];

        t_349[k] = pb_x[k] * hg_253[k];

        t_350[k] = pb_x[k] * hg_254[k];

        t_351[k] = pa_z[k] * gh_225[k];

        t_352[k] = f_2 * gg_160[k]
                   + pb_z[k] * hg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pa_y, pa_z, pb_y, fh_167, gg_161, gg_162, \
                         gg_179, gh_227, gh_228, gh_251, hg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_3 * gg_161[k]
                   + pa_z[k] * gh_227[k];

        t_354[k] = f_4 * gg_162[k]
                   + pa_z[k] * gh_228[k];

        t_355[k] = f_1 * gg_179[k]
                   + pb_y[k] * hg_254[k];

        t_356[k] = f_4 * fh_167[k]
                   + pa_y[k] * gh_251[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, pb_x, hf_170, hf_171, hf_172, \
                         hf_173, hf_174, hg_255, hg_256, hg_257, hg_258, \
                         hg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_1 * hf_170[k]
                   + pb_x[k] * hg_255[k];

        t_358[k] = f_4 * hf_171[k]
                   + pb_x[k] * hg_256[k];

        t_359[k] = f_4 * hf_172[k]
                   + pb_x[k] * hg_257[k];

        t_360[k] = f_3 * hf_173[k]
                   + pb_x[k] * hg_258[k];

        t_361[k] = f_3 * hf_174[k]
                   + pb_x[k] * hg_259[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, pb_x, hf_175, hf_176, hf_177, \
                         hf_178, hf_179, hg_260, hg_261, hg_262, hg_263, \
                         hg_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_3 * hf_175[k]
                   + pb_x[k] * hg_260[k];

        t_363[k] = f_2 * hf_176[k]
                   + pb_x[k] * hg_261[k];

        t_364[k] = f_2 * hf_177[k]
                   + pb_x[k] * hg_262[k];

        t_365[k] = f_2 * hf_178[k]
                   + pb_x[k] * hg_263[k];

        t_366[k] = f_2 * hf_179[k]
                   + pb_x[k] * hg_264[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, t_372, pa_z, pb_x, fh_141, gh_246, \
                         hg_265, hg_266, hg_267, hg_268, hg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pb_x[k] * hg_265[k];

        t_368[k] = pb_x[k] * hg_266[k];

        t_369[k] = pb_x[k] * hg_267[k];

        t_370[k] = pb_x[k] * hg_268[k];

        t_371[k] = pb_x[k] * hg_269[k];

        t_372[k] = f_2 * fh_141[k]
                   + pa_z[k] * gh_246[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_y, pb_z, gg_175, gg_192, gg_193, \
                         gg_194, hf_178, hf_179, hg_265, hg_267, hg_268, \
                         hg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_3 * gg_175[k]
                   + pb_z[k] * hg_265[k];

        t_374[k] = f_4 * gg_192[k]
                   + f_3 * hf_178[k]
                   + pb_y[k] * hg_267[k];

        t_375[k] = f_4 * gg_193[k]
                   + f_2 * hf_179[k]
                   + pb_y[k] * hg_268[k];

        t_376[k] = f_4 * gg_194[k]
                   + pb_y[k] * hg_269[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pb_x, fh_188, gh_272, hf_180, \
                         hf_181, hf_182, hg_270, hg_271, hg_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_3 * fh_188[k]
                   + pa_y[k] * gh_272[k];

        t_378[k] = f_1 * hf_180[k]
                   + pb_x[k] * hg_270[k];

        t_379[k] = f_4 * hf_181[k]
                   + pb_x[k] * hg_271[k];

        t_380[k] = f_4 * hf_182[k]
                   + pb_x[k] * hg_272[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, pb_x, hf_183, hf_184, hf_185, \
                         hf_186, hf_187, hg_273, hg_274, hg_275, hg_276, \
                         hg_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_3 * hf_183[k]
                   + pb_x[k] * hg_273[k];

        t_382[k] = f_3 * hf_184[k]
                   + pb_x[k] * hg_274[k];

        t_383[k] = f_3 * hf_185[k]
                   + pb_x[k] * hg_275[k];

        t_384[k] = f_2 * hf_186[k]
                   + pb_x[k] * hg_276[k];

        t_385[k] = f_2 * hf_187[k]
                   + pb_x[k] * hg_277[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, t_391, pb_x, hf_188, hf_189, \
                         hg_278, hg_279, hg_280, hg_281, hg_282, \
                         hg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_2 * hf_188[k]
                   + pb_x[k] * hg_278[k];

        t_387[k] = f_2 * hf_189[k]
                   + pb_x[k] * hg_279[k];

        t_388[k] = pb_x[k] * hg_280[k];

        t_389[k] = pb_x[k] * hg_281[k];

        t_390[k] = pb_x[k] * hg_282[k];

        t_391[k] = pb_x[k] * hg_283[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pa_z, pb_x, pb_y, pb_z, fh_162, gg_190, \
                         gg_207, gh_267, hf_188, hg_280, hg_282, \
                         hg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = pb_x[k] * hg_284[k];

        t_393[k] = f_3 * fh_162[k]
                   + pa_z[k] * gh_267[k];

        t_394[k] = f_4 * gg_190[k]
                   + pb_z[k] * hg_280[k];

        t_395[k] = f_3 * gg_207[k]
                   + f_3 * hf_188[k]
                   + pb_y[k] * hg_282[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pa_y, pb_y, fh_209, gg_208, gg_209, \
                         gh_293, gh_294, hf_189, hg_283, hg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_3 * gg_208[k]
                   + f_2 * hf_189[k]
                   + pb_y[k] * hg_283[k];

        t_397[k] = f_3 * gg_209[k]
                   + pb_y[k] * hg_284[k];

        t_398[k] = f_2 * fh_209[k]
                   + pa_y[k] * gh_293[k];

        t_399[k] = pa_y[k] * gh_294[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, pa_y, pb_x, gh_296, gh_299, \
                         hf_191, hf_193, hf_194, hg_286, hg_288, \
                         hg_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_4 * hf_191[k]
                   + pb_x[k] * hg_286[k];

        t_401[k] = pa_y[k] * gh_296[k];

        t_402[k] = f_3 * hf_193[k]
                   + pb_x[k] * hg_288[k];

        t_403[k] = f_3 * hf_194[k]
                   + pb_x[k] * hg_289[k];

        t_404[k] = pa_y[k] * gh_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pa_y, pb_x, gh_303, hf_196, \
                         hf_197, hf_198, hg_291, hg_292, hg_293, \
                         hg_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_2 * hf_196[k]
                   + pb_x[k] * hg_291[k];

        t_406[k] = f_2 * hf_197[k]
                   + pb_x[k] * hg_292[k];

        t_407[k] = f_2 * hf_198[k]
                   + pb_x[k] * hg_293[k];

        t_408[k] = pa_y[k] * gh_303[k];

        t_409[k] = pb_x[k] * hg_295[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pa_y, pb_x, gg_220, gh_309, \
                         hg_296, hg_297, hg_298, hg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * hg_296[k];

        t_411[k] = pb_x[k] * hg_297[k];

        t_412[k] = pb_x[k] * hg_298[k];

        t_413[k] = pb_x[k] * hg_299[k];

        t_414[k] = f_0 * gg_220[k]
                   + pa_y[k] * gh_309[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_y, pb_y, pb_z, gg_205, gg_222, gg_223, \
                         gg_224, gh_311, gh_312, hg_295, hg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_1 * gg_205[k]
                   + pb_z[k] * hg_295[k];

        t_416[k] = f_4 * gg_222[k]
                   + pa_y[k] * gh_311[k];

        t_417[k] = f_3 * gg_223[k]
                   + pa_y[k] * gh_312[k];

        t_418[k] = f_2 * gg_224[k]
                   + pb_y[k] * hg_299[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, t_424, pa_y, pb_x, pb_y, gh_314, \
                         hf_200, hf_202, hf_203, hg_300, hg_302, \
                         hg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pa_y[k] * gh_314[k];

        t_420[k] = f_1 * hf_200[k]
                   + pb_x[k] * hg_300[k];

        t_421[k] = pb_y[k] * hg_300[k];

        t_422[k] = f_4 * hf_202[k]
                   + pb_x[k] * hg_302[k];

        t_423[k] = f_3 * hf_203[k]
                   + pb_x[k] * hg_303[k];

        t_424[k] = pb_y[k] * hg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, pb_x, pb_y, hf_205, hf_206, \
                         hf_207, hf_209, hg_305, hg_306, hg_307, \
                         hg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_3 * hf_205[k]
                   + pb_x[k] * hg_305[k];

        t_426[k] = f_2 * hf_206[k]
                   + pb_x[k] * hg_306[k];

        t_427[k] = f_2 * hf_207[k]
                   + pb_x[k] * hg_307[k];

        t_428[k] = pb_y[k] * hg_305[k];

        t_429[k] = f_2 * hf_209[k]
                   + pb_x[k] * hg_309[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, t_435, t_436, pb_x, pb_y, hf_206, \
                         hf_207, hg_310, hg_311, hg_312, hg_313, \
                         hg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pb_x[k] * hg_310[k];

        t_431[k] = pb_x[k] * hg_311[k];

        t_432[k] = pb_x[k] * hg_312[k];

        t_433[k] = pb_x[k] * hg_313[k];

        t_434[k] = pb_x[k] * hg_314[k];

        t_435[k] = f_1 * hf_206[k]
                   + pb_y[k] * hg_310[k];

        t_436[k] = f_4 * hf_207[k]
                   + pb_y[k] * hg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pb_y, pb_z, gg_224, hf_208, hf_209, \
                         hg_312, hg_313, hg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_3 * hf_208[k]
                   + pb_y[k] * hg_312[k];

        t_438[k] = f_2 * hf_209[k]
                   + pb_y[k] * hg_313[k];

        t_439[k] = pb_y[k] * hg_314[k];

        t_440[k] = f_0 * gg_224[k]
                   + f_1 * hf_209[k]
                   + pb_z[k] * hg_314[k];
    }
}

auto
compute_prim_hh_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fh, const size_t gg, const size_t gh,
                          const size_t hf, const size_t hg, const size_t ncols,
                          const double p) -> void
{
    compute_prim_hh_overlap_0_piece0(buffer, target, pa, pb, fh, gg, gh, hf, hg, ncols, p);

    compute_prim_hh_overlap_0_piece1(buffer, target, pa, pb, fh, gg, gh, hf, hg, ncols, p);

    compute_prim_hh_overlap_0_piece2(buffer, target, pa, pb, fh, gg, gh, hf, hg, ncols, p);
}

}  // namespace simdovl
