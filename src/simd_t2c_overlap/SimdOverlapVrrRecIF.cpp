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


#include "SimdOverlapVrrRecIF.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

static auto
compute_prim_if_overlap_0_piece0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gf, const size_t hd,
                                 const size_t hf, const size_t ip, const size_t id,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_89 = buffer.data(hd + 89);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_149 = buffer.data(hf + 149);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hd_0, hd_3, ip_0, id_0, \
                         id_2, id_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 + f_1 * ip_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = pb_y[k] * id_0[k];

        t_2[k] = pb_z[k] * id_0[k];

        t_3[k] = f_0 * hd_3[k]
                 + pb_x[k] * id_3[k];

        t_4[k] = pb_y[k] * id_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pb_x, pb_y, pb_z, hd_5, hf_0, \
                         ip_1, ip_2, id_3, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * hd_5[k]
                 + pb_x[k] * id_5[k];

        t_6[k] = f_1 * ip_1[k]
                 + pb_y[k] * id_3[k];

        t_7[k] = pb_z[k] * id_3[k];

        t_8[k] = pb_y[k] * id_5[k];

        t_9[k] = f_1 * ip_2[k]
                 + pb_z[k] * id_5[k];

        t_10[k] = pa_y[k] * hf_0[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, pb_z, hd_0, hd_9, \
                         hf_5, id_6, id_7, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * hd_0[k]
                  + pb_y[k] * id_6[k];

        t_12[k] = pb_z[k] * id_6[k];

        t_13[k] = f_3 * hd_9[k]
                  + pb_x[k] * id_9[k];

        t_14[k] = pb_z[k] * id_7[k];

        t_15[k] = pa_y[k] * hf_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, gf_16, hd_5, hf_9, \
                         hf_16, id_9, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * gf_16[k]
                  + pa_x[k] * hf_16[k];

        t_17[k] = pb_z[k] * id_9[k];

        t_18[k] = f_2 * hd_5[k]
                  + pb_y[k] * id_11[k];

        t_19[k] = pa_y[k] * hf_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, hd_0, hf_0, hf_3, \
                         id_12, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * hf_0[k];

        t_21[k] = pb_y[k] * id_12[k];

        t_22[k] = f_2 * hd_0[k]
                  + pb_z[k] * id_12[k];

        t_23[k] = pa_z[k] * hf_3[k];

        t_24[k] = pb_y[k] * id_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pa_z, pb_x, pb_y, gf_29, hd_17, \
                         hf_6, hf_29, ip_8, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * hd_17[k]
                  + pb_x[k] * id_17[k];

        t_26[k] = pa_z[k] * hf_6[k];

        t_27[k] = f_2 * ip_8[k]
                  + pb_y[k] * id_16[k];

        t_28[k] = pb_y[k] * id_17[k];

        t_29[k] = f_4 * gf_29[k]
                  + pa_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_y, pb_x, pb_y, pb_z, gf_0, hd_6, \
                         hd_21, hf_10, id_18, id_19, id_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * gf_0[k]
                  + pa_y[k] * hf_10[k];

        t_31[k] = f_1 * hd_6[k]
                  + pb_y[k] * id_18[k];

        t_32[k] = pb_z[k] * id_18[k];

        t_33[k] = f_4 * hd_21[k]
                  + pb_x[k] * id_21[k];

        t_34[k] = pb_z[k] * id_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pb_x, pb_y, pb_z, gf_36, hd_11, \
                         hd_23, hf_36, ip_11, id_21, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_4 * hd_23[k]
                  + pb_x[k] * id_23[k];

        t_36[k] = f_5 * gf_36[k]
                  + pa_x[k] * hf_36[k];

        t_37[k] = pb_z[k] * id_21[k];

        t_38[k] = f_1 * hd_11[k]
                  + pb_y[k] * id_23[k];

        t_39[k] = f_1 * ip_11[k]
                  + pb_z[k] * id_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, hd_28, hf_11, \
                         hf_13, hf_20, hf_22, hf_25, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * hf_20[k];

        t_41[k] = pa_z[k] * hf_11[k];

        t_42[k] = pa_y[k] * hf_22[k];

        t_43[k] = pa_z[k] * hf_13[k];

        t_44[k] = f_4 * hd_28[k]
                  + pb_x[k] * id_28[k];

        t_45[k] = pa_y[k] * hf_25[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_y, pb_z, hd_9, hd_17, hf_16, \
                         hf_29, id_27, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * hf_16[k];

        t_47[k] = f_2 * hd_9[k]
                  + pb_z[k] * id_27[k];

        t_48[k] = f_2 * hd_17[k]
                  + pb_y[k] * id_29[k];

        t_49[k] = pa_y[k] * hf_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_z, pb_x, pb_y, pb_z, gf_0, hd_12, \
                         hd_33, hf_20, id_30, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * gf_0[k]
                  + pa_z[k] * hf_20[k];

        t_51[k] = pb_y[k] * id_30[k];

        t_52[k] = f_1 * hd_12[k]
                  + pb_z[k] * id_30[k];

        t_53[k] = f_4 * hd_33[k]
                  + pb_x[k] * id_33[k];

        t_54[k] = pb_y[k] * id_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_x, pb_y, gf_59, hd_35, hf_59, \
                         ip_16, ip_17, id_33, id_34, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * hd_35[k]
                  + pb_x[k] * id_35[k];

        t_56[k] = f_1 * ip_16[k]
                  + pb_y[k] * id_33[k];

        t_57[k] = f_2 * ip_17[k]
                  + pb_y[k] * id_34[k];

        t_58[k] = pb_y[k] * id_35[k];

        t_59[k] = f_5 * gf_59[k]
                  + pa_x[k] * hf_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_y, pb_x, pb_y, pb_z, gf_10, hd_18, \
                         hd_39, hf_30, id_36, id_37, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * gf_10[k]
                  + pa_y[k] * hf_30[k];

        t_61[k] = f_5 * hd_18[k]
                  + pb_y[k] * id_36[k];

        t_62[k] = pb_z[k] * id_36[k];

        t_63[k] = f_5 * hd_39[k]
                  + pb_x[k] * id_39[k];

        t_64[k] = pb_z[k] * id_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pb_x, pb_y, pb_z, gf_66, hd_23, \
                         hd_41, hf_66, ip_20, id_39, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_5 * hd_41[k]
                  + pb_x[k] * id_41[k];

        t_66[k] = f_1 * gf_66[k]
                  + pa_x[k] * hf_66[k];

        t_67[k] = pb_z[k] * id_39[k];

        t_68[k] = f_5 * hd_23[k]
                  + pb_y[k] * id_41[k];

        t_69[k] = f_1 * ip_20[k]
                  + pb_z[k] * id_41[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_x, pb_z, hd_18, hd_46, hf_30, \
                         hf_31, hf_33, id_42, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * hf_30[k];

        t_71[k] = pa_z[k] * hf_31[k];

        t_72[k] = f_2 * hd_18[k]
                  + pb_z[k] * id_42[k];

        t_73[k] = pa_z[k] * hf_33[k];

        t_74[k] = f_5 * hd_46[k]
                  + pb_x[k] * id_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_x, pb_y, pb_z, hd_21, hd_29, hd_47, \
                         hf_36, id_45, id_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * hd_47[k]
                  + pb_x[k] * id_47[k];

        t_76[k] = pa_z[k] * hf_36[k];

        t_77[k] = f_2 * hd_21[k]
                  + pb_z[k] * id_45[k];

        t_78[k] = f_1 * hd_29[k]
                  + pb_y[k] * id_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_x, pa_y, pb_y, gf_79, hd_30, hf_50, hf_52, \
                         hf_79, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * gf_79[k]
                  + pa_x[k] * hf_79[k];

        t_80[k] = pa_y[k] * hf_50[k];

        t_81[k] = f_2 * hd_30[k]
                  + pb_y[k] * id_48[k];

        t_82[k] = pa_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pa_y, pb_x, gf_86, hd_51, hd_52, hf_55, \
                         hf_86, id_51, id_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_5 * hd_51[k]
                  + pb_x[k] * id_51[k];

        t_84[k] = f_5 * hd_52[k]
                  + pb_x[k] * id_52[k];

        t_85[k] = pa_y[k] * hf_55[k];

        t_86[k] = f_1 * gf_86[k]
                  + pa_x[k] * hf_86[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_y, pb_z, gf_20, hd_27, hd_35, \
                         hf_50, hf_59, id_51, id_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_1 * hd_27[k]
                  + pb_z[k] * id_51[k];

        t_88[k] = f_2 * hd_35[k]
                  + pb_y[k] * id_53[k];

        t_89[k] = pa_y[k] * hf_59[k];

        t_90[k] = f_1 * gf_20[k]
                  + pa_z[k] * hf_50[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, hd_30, hd_57, hd_59, \
                         id_54, id_56, id_57, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_y[k] * id_54[k];

        t_92[k] = f_5 * hd_30[k]
                  + pb_z[k] * id_54[k];

        t_93[k] = f_5 * hd_57[k]
                  + pb_x[k] * id_57[k];

        t_94[k] = pb_y[k] * id_56[k];

        t_95[k] = f_5 * hd_59[k]
                  + pb_x[k] * id_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, gf_99, hf_99, ip_28, ip_29, \
                         id_57, id_58, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * ip_28[k]
                  + pb_y[k] * id_57[k];

        t_97[k] = f_2 * ip_29[k]
                  + pb_y[k] * id_58[k];

        t_98[k] = pb_y[k] * id_59[k];

        t_99[k] = f_1 * gf_99[k]
                  + pa_x[k] * hf_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pa_y, pb_x, pb_y, pb_z, gf_30, \
                         hd_36, hd_63, hf_60, id_60, id_61, id_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * gf_30[k]
                   + pa_y[k] * hf_60[k];

        t_101[k] = f_4 * hd_36[k]
                   + pb_y[k] * id_60[k];

        t_102[k] = pb_z[k] * id_60[k];

        t_103[k] = f_1 * hd_63[k]
                   + pb_x[k] * id_63[k];

        t_104[k] = pb_z[k] * id_61[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_x, pb_x, pb_y, pb_z, gf_106, \
                         hd_41, hd_65, hf_106, ip_32, id_63, id_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_1 * hd_65[k]
                   + pb_x[k] * id_65[k];

        t_106[k] = f_2 * gf_106[k]
                   + pa_x[k] * hf_106[k];

        t_107[k] = pb_z[k] * id_63[k];

        t_108[k] = f_4 * hd_41[k]
                   + pb_y[k] * id_65[k];

        t_109[k] = f_1 * ip_32[k]
                   + pb_z[k] * id_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_z, pb_x, pb_z, hd_36, hd_70, \
                         hf_60, hf_61, hf_63, id_66, id_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * hf_60[k];

        t_111[k] = pa_z[k] * hf_61[k];

        t_112[k] = f_2 * hd_36[k]
                   + pb_z[k] * id_66[k];

        t_113[k] = pa_z[k] * hf_63[k];

        t_114[k] = f_1 * hd_70[k]
                   + pb_x[k] * id_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_z, pb_x, pb_y, pb_z, hd_39, hd_47, \
                         hd_71, hf_66, id_69, id_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_1 * hd_71[k]
                   + pb_x[k] * id_71[k];

        t_116[k] = pa_z[k] * hf_66[k];

        t_117[k] = f_2 * hd_39[k]
                   + pb_z[k] * id_69[k];

        t_118[k] = f_5 * hd_47[k]
                   + pb_y[k] * id_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_x, pa_y, pb_y, pb_z, gf_50, gf_119, \
                         hd_42, hd_48, hf_80, hf_119, id_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * gf_119[k]
                   + pa_x[k] * hf_119[k];

        t_120[k] = f_2 * gf_50[k]
                   + pa_y[k] * hf_80[k];

        t_121[k] = f_1 * hd_48[k]
                   + pb_y[k] * id_72[k];

        t_122[k] = f_1 * hd_42[k]
                   + pb_z[k] * id_72[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pb_x, gf_126, hd_75, hd_76, hd_77, \
                         hf_126, id_75, id_76, id_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_1 * hd_75[k]
                   + pb_x[k] * id_75[k];

        t_124[k] = f_1 * hd_76[k]
                   + pb_x[k] * id_76[k];

        t_125[k] = f_1 * hd_77[k]
                   + pb_x[k] * id_77[k];

        t_126[k] = f_2 * gf_126[k]
                   + pa_x[k] * hf_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pa_y, pb_y, pb_z, gf_129, hd_45, \
                         hd_53, hf_90, hf_129, id_75, id_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_1 * hd_45[k]
                   + pb_z[k] * id_75[k];

        t_128[k] = f_1 * hd_53[k]
                   + pb_y[k] * id_77[k];

        t_129[k] = f_2 * gf_129[k]
                   + pa_x[k] * hf_129[k];

        t_130[k] = pa_y[k] * hf_90[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_y, pb_x, pb_y, hd_54, hd_81, \
                         hd_82, hf_92, hf_95, id_78, id_81, id_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_2 * hd_54[k]
                   + pb_y[k] * id_78[k];

        t_132[k] = pa_y[k] * hf_92[k];

        t_133[k] = f_1 * hd_81[k]
                   + pb_x[k] * id_81[k];

        t_134[k] = f_1 * hd_82[k]
                   + pb_x[k] * id_82[k];

        t_135[k] = pa_y[k] * hf_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, pa_y, pb_y, pb_z, gf_136, hd_51, \
                         hd_59, hf_99, hf_136, id_81, id_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * gf_136[k]
                   + pa_x[k] * hf_136[k];

        t_137[k] = f_5 * hd_51[k]
                   + pb_z[k] * id_81[k];

        t_138[k] = f_2 * hd_59[k]
                   + pb_y[k] * id_83[k];

        t_139[k] = pa_y[k] * hf_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pa_z, pb_x, pb_y, pb_z, gf_50, \
                         hd_54, hd_87, hf_90, id_84, id_86, id_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_5 * gf_50[k]
                   + pa_z[k] * hf_90[k];

        t_141[k] = pb_y[k] * id_84[k];

        t_142[k] = f_4 * hd_54[k]
                   + pb_z[k] * id_84[k];

        t_143[k] = f_1 * hd_87[k]
                   + pb_x[k] * id_87[k];

        t_144[k] = pb_y[k] * id_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_x, pb_x, pb_y, gf_149, hd_89, \
                         hf_149, ip_43, ip_44, id_87, id_88, id_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_1 * hd_89[k]
                   + pb_x[k] * id_89[k];

        t_146[k] = f_1 * ip_43[k]
                   + pb_y[k] * id_87[k];

        t_147[k] = f_2 * ip_44[k]
                   + pb_y[k] * id_88[k];

        t_148[k] = pb_y[k] * id_89[k];

        t_149[k] = f_2 * gf_149[k]
                   + pa_x[k] * hf_149[k];
    }
}

static auto
compute_prim_if_overlap_0_piece1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gf, const size_t hd,
                                 const size_t hf, const size_t ip, const size_t id,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 2.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_x, pb_x, pb_y, pb_z, hd_60, \
                         hd_90, hd_93, hf_150, id_90, id_91, id_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * hd_90[k]
                   + pa_x[k] * hf_150[k];

        t_151[k] = f_3 * hd_60[k]
                   + pb_y[k] * id_90[k];

        t_152[k] = pb_z[k] * id_90[k];

        t_153[k] = f_2 * hd_93[k]
                   + pb_x[k] * id_93[k];

        t_154[k] = pb_z[k] * id_91[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_x, pb_x, pb_z, hd_95, hf_156, \
                         hf_158, hf_159, id_93, id_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_2 * hd_95[k]
                   + pb_x[k] * id_95[k];

        t_156[k] = pa_x[k] * hf_156[k];

        t_157[k] = pb_z[k] * id_93[k];

        t_158[k] = pa_x[k] * hf_158[k];

        t_159[k] = pa_x[k] * hf_159[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pa_z, pb_x, pb_z, hd_60, hd_100, \
                         hf_100, hf_101, hf_103, id_96, id_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_z[k] * hf_100[k];

        t_161[k] = pa_z[k] * hf_101[k];

        t_162[k] = f_2 * hd_60[k]
                   + pb_z[k] * id_96[k];

        t_163[k] = pa_z[k] * hf_103[k];

        t_164[k] = f_2 * hd_100[k]
                   + pb_x[k] * id_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, t_170, pa_x, pb_x, hd_101, hd_102, \
                         hf_166, hf_167, hf_168, hf_169, hf_170, \
                         id_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_2 * hd_101[k]
                   + pb_x[k] * id_101[k];

        t_166[k] = pa_x[k] * hf_166[k];

        t_167[k] = pa_x[k] * hf_167[k];

        t_168[k] = pa_x[k] * hf_168[k];

        t_169[k] = pa_x[k] * hf_169[k];

        t_170[k] = f_5 * hd_102[k]
                   + pa_x[k] * hf_170[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pb_x, pb_y, pb_z, hd_66, hd_72, hd_105, \
                         hd_106, id_102, id_105, id_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_5 * hd_72[k]
                   + pb_y[k] * id_102[k];

        t_172[k] = f_1 * hd_66[k]
                   + pb_z[k] * id_102[k];

        t_173[k] = f_2 * hd_105[k]
                   + pb_x[k] * id_105[k];

        t_174[k] = f_2 * hd_106[k]
                   + pb_x[k] * id_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, pa_x, pb_x, hd_107, hd_108, \
                         hf_176, hf_177, hf_178, hf_179, hf_180, \
                         id_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_2 * hd_107[k]
                   + pb_x[k] * id_107[k];

        t_176[k] = pa_x[k] * hf_176[k];

        t_177[k] = pa_x[k] * hf_177[k];

        t_178[k] = pa_x[k] * hf_178[k];

        t_179[k] = pa_x[k] * hf_179[k];

        t_180[k] = f_5 * hd_108[k]
                   + pa_x[k] * hf_180[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_x, pb_y, pb_z, hd_72, hd_78, hd_111, \
                         hd_112, id_108, id_111, id_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_1 * hd_78[k]
                   + pb_y[k] * id_108[k];

        t_182[k] = f_5 * hd_72[k]
                   + pb_z[k] * id_108[k];

        t_183[k] = f_2 * hd_111[k]
                   + pb_x[k] * id_111[k];

        t_184[k] = f_2 * hd_112[k]
                   + pb_x[k] * id_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, pa_x, pa_y, pb_x, hd_113, \
                         hf_140, hf_186, hf_187, hf_188, hf_189, \
                         id_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_2 * hd_113[k]
                   + pb_x[k] * id_113[k];

        t_186[k] = pa_x[k] * hf_186[k];

        t_187[k] = pa_x[k] * hf_187[k];

        t_188[k] = pa_x[k] * hf_188[k];

        t_189[k] = pa_x[k] * hf_189[k];

        t_190[k] = pa_y[k] * hf_140[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pa_y, pb_x, pb_y, hd_84, hd_117, \
                         hd_118, hf_142, hf_145, id_114, id_117, \
                         id_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_2 * hd_84[k]
                   + pb_y[k] * id_114[k];

        t_192[k] = pa_y[k] * hf_142[k];

        t_193[k] = f_2 * hd_117[k]
                   + pb_x[k] * id_117[k];

        t_194[k] = f_2 * hd_118[k]
                   + pb_x[k] * id_118[k];

        t_195[k] = pa_y[k] * hf_145[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, t_201, pa_x, pb_y, hd_120, hf_196, \
                         hf_197, hf_198, hf_199, hf_200, id_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * hf_196[k];

        t_197[k] = pa_x[k] * hf_197[k];

        t_198[k] = pa_x[k] * hf_198[k];

        t_199[k] = pa_x[k] * hf_199[k];

        t_200[k] = f_5 * hd_120[k]
                   + pa_x[k] * hf_200[k];

        t_201[k] = pb_y[k] * id_120[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, hd_84, hd_123, hd_125, \
                         id_120, id_122, id_123, id_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * hd_84[k]
                   + pb_z[k] * id_120[k];

        t_203[k] = f_2 * hd_123[k]
                   + pb_x[k] * id_123[k];

        t_204[k] = pb_y[k] * id_122[k];

        t_205[k] = f_2 * hd_125[k]
                   + pb_x[k] * id_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, pa_x, pb_x, pb_y, hf_206, hf_207, \
                         hf_209, ip_63, id_125, id_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_x[k] * hf_206[k];

        t_207[k] = pa_x[k] * hf_207[k];

        t_208[k] = pb_y[k] * id_125[k];

        t_209[k] = pa_x[k] * hf_209[k];

        t_210[k] = f_1 * ip_63[k]
                   + pb_x[k] * id_126[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, t_216, pb_x, pb_y, pb_z, hd_93, \
                         ip_64, id_126, id_127, id_129, id_130, \
                         id_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_2 * ip_64[k]
                   + pb_x[k] * id_127[k];

        t_212[k] = pb_z[k] * id_126[k];

        t_213[k] = pb_x[k] * id_129[k];

        t_214[k] = pb_x[k] * id_130[k];

        t_215[k] = pb_x[k] * id_131[k];

        t_216[k] = f_0 * hd_93[k]
                   + f_1 * ip_64[k]
                   + pb_y[k] * id_129[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_y, pb_z, hd_95, hf_150, \
                         hf_151, ip_65, id_129, id_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_z[k] * id_129[k];

        t_218[k] = f_0 * hd_95[k]
                   + pb_y[k] * id_131[k];

        t_219[k] = f_1 * ip_65[k]
                   + pb_z[k] * id_131[k];

        t_220[k] = pa_z[k] * hf_150[k];

        t_221[k] = pa_z[k] * hf_151[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, hd_93, \
                         hf_156, ip_68, id_134, id_135, id_136, \
                         id_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_2 * ip_68[k]
                   + pb_x[k] * id_134[k];

        t_223[k] = pb_x[k] * id_135[k];

        t_224[k] = pb_x[k] * id_136[k];

        t_225[k] = pb_x[k] * id_137[k];

        t_226[k] = pa_z[k] * hf_156[k];

        t_227[k] = f_2 * hd_93[k]
                   + pb_z[k] * id_135[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, pb_x, pb_y, gf_119, hd_101, hf_169, \
                         ip_69, ip_70, id_137, id_138, id_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_3 * hd_101[k]
                   + pb_y[k] * id_137[k];

        t_229[k] = f_4 * gf_119[k]
                   + pa_y[k] * hf_169[k];

        t_230[k] = f_1 * ip_69[k]
                   + pb_x[k] * id_138[k];

        t_231[k] = f_2 * ip_70[k]
                   + pb_x[k] * id_139[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pa_z, pb_x, gf_106, hf_166, ip_71, \
                         id_140, id_141, id_142, id_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_2 * ip_71[k]
                   + pb_x[k] * id_140[k];

        t_233[k] = pb_x[k] * id_141[k];

        t_234[k] = pb_x[k] * id_142[k];

        t_235[k] = pb_x[k] * id_143[k];

        t_236[k] = f_2 * gf_106[k]
                   + pa_z[k] * hf_166[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_y, pb_x, pb_y, pb_z, gf_129, hd_99, \
                         hd_107, hf_179, ip_72, id_141, id_143, \
                         id_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_1 * hd_99[k]
                   + pb_z[k] * id_141[k];

        t_238[k] = f_4 * hd_107[k]
                   + pb_y[k] * id_143[k];

        t_239[k] = f_5 * gf_129[k]
                   + pa_y[k] * hf_179[k];

        t_240[k] = f_1 * ip_72[k]
                   + pb_x[k] * id_144[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pb_x, ip_73, ip_74, id_145, \
                         id_146, id_147, id_148, id_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_2 * ip_73[k]
                   + pb_x[k] * id_145[k];

        t_242[k] = f_2 * ip_74[k]
                   + pb_x[k] * id_146[k];

        t_243[k] = pb_x[k] * id_147[k];

        t_244[k] = pb_x[k] * id_148[k];

        t_245[k] = pb_x[k] * id_149[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pa_y, pa_z, pb_y, pb_z, gf_116, gf_139, \
                         hd_105, hd_113, hf_176, hf_189, id_147, \
                         id_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_1 * gf_116[k]
                   + pa_z[k] * hf_176[k];

        t_247[k] = f_5 * hd_105[k]
                   + pb_z[k] * id_147[k];

        t_248[k] = f_5 * hd_113[k]
                   + pb_y[k] * id_149[k];

        t_249[k] = f_1 * gf_139[k]
                   + pa_y[k] * hf_189[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, t_255, pb_x, ip_75, ip_76, ip_77, \
                         id_150, id_151, id_152, id_153, id_154, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_1 * ip_75[k]
                   + pb_x[k] * id_150[k];

        t_251[k] = f_2 * ip_76[k]
                   + pb_x[k] * id_151[k];

        t_252[k] = f_2 * ip_77[k]
                   + pb_x[k] * id_152[k];

        t_253[k] = pb_x[k] * id_153[k];

        t_254[k] = pb_x[k] * id_154[k];

        t_255[k] = pb_x[k] * id_155[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_y, pa_z, pb_y, pb_z, gf_126, gf_149, \
                         hd_111, hd_119, hf_186, hf_199, id_153, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_5 * gf_126[k]
                   + pa_z[k] * hf_186[k];

        t_257[k] = f_4 * hd_111[k]
                   + pb_z[k] * id_153[k];

        t_258[k] = f_1 * hd_119[k]
                   + pb_y[k] * id_155[k];

        t_259[k] = f_2 * gf_149[k]
                   + pa_y[k] * hf_199[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, t_265, pa_y, pb_x, hf_200, hf_202, \
                         ip_79, id_157, id_159, id_160, id_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * hf_200[k];

        t_261[k] = f_2 * ip_79[k]
                   + pb_x[k] * id_157[k];

        t_262[k] = pa_y[k] * hf_202[k];

        t_263[k] = pb_x[k] * id_159[k];

        t_264[k] = pb_x[k] * id_160[k];

        t_265[k] = pb_x[k] * id_161[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pa_y, pb_y, pb_z, hd_117, hd_123, hd_125, \
                         hf_206, hf_209, id_159, id_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_5 * hd_123[k]
                   + pa_y[k] * hf_206[k];

        t_267[k] = f_3 * hd_117[k]
                   + pb_z[k] * id_159[k];

        t_268[k] = f_2 * hd_125[k]
                   + pb_y[k] * id_161[k];

        t_269[k] = pa_y[k] * hf_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, pb_x, pb_y, ip_81, ip_83, \
                         id_162, id_164, id_165, id_166, id_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * ip_81[k]
                   + pb_x[k] * id_162[k];

        t_271[k] = pb_y[k] * id_162[k];

        t_272[k] = f_2 * ip_83[k]
                   + pb_x[k] * id_164[k];

        t_273[k] = pb_x[k] * id_165[k];

        t_274[k] = pb_x[k] * id_166[k];

        t_275[k] = pb_x[k] * id_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_y, pb_z, hd_125, ip_82, ip_83, id_165, \
                         id_166, id_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * ip_82[k]
                   + pb_y[k] * id_165[k];

        t_277[k] = f_2 * ip_83[k]
                   + pb_y[k] * id_166[k];

        t_278[k] = pb_y[k] * id_167[k];

        t_279[k] = f_0 * hd_125[k]
                   + f_1 * ip_83[k]
                   + pb_z[k] * id_167[k];
    }
}

auto
compute_prim_if_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t gf, const size_t hd, const size_t hf,
                          const size_t ip, const size_t id, const size_t ncols,
                          const double p) -> void
{
    compute_prim_if_overlap_0_piece0(buffer, target, pa, pb, gf, hd, hf, ip, id, ncols, p);

    compute_prim_if_overlap_0_piece1(buffer, target, pa, pb, gf, hd, hf, ip, id, ncols, p);
}

}  // namespace simdovl
