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


#include "SimdOverlapVrrRecFI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

static auto
compute_prim_fi_overlap_0_piece0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pi, const size_t dh,
                                 const size_t di, const size_t fg, const size_t fh,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_83 = buffer.data(pi + 83);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_96 = buffer.data(dh + 96);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_3 = buffer.data(di + 3);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_14 = buffer.data(di + 14);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_27 = buffer.data(di + 27);
    const auto *di_29 = buffer.data(di + 29);
    const auto *di_31 = buffer.data(di + 31);
    const auto *di_34 = buffer.data(di + 34);
    const auto *di_38 = buffer.data(di + 38);
    const auto *di_43 = buffer.data(di + 43);
    const auto *di_49 = buffer.data(di + 49);
    const auto *di_56 = buffer.data(di + 56);
    const auto *di_58 = buffer.data(di + 58);
    const auto *di_61 = buffer.data(di + 61);
    const auto *di_65 = buffer.data(di + 65);
    const auto *di_70 = buffer.data(di + 70);
    const auto *di_76 = buffer.data(di + 76);
    const auto *di_83 = buffer.data(di + 83);
    const auto *di_84 = buffer.data(di + 84);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_110 = buffer.data(di + 110);
    const auto *di_111 = buffer.data(di + 111);
    const auto *di_124 = buffer.data(di + 124);
    const auto *di_133 = buffer.data(di + 133);
    const auto *di_134 = buffer.data(di + 134);
    const auto *di_135 = buffer.data(di + 135);
    const auto *di_136 = buffer.data(di + 136);
    const auto *di_137 = buffer.data(di + 137);
    const auto *di_138 = buffer.data(di + 138);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_154 = buffer.data(di + 154);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_34 = buffer.data(fg + 34);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_47 = buffer.data(fg + 47);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_76 = buffer.data(fg + 76);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_78 = buffer.data(fg + 78);
    const auto *fg_79 = buffer.data(fg + 79);
    const auto *fg_80 = buffer.data(fg + 80);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_108 = buffer.data(fh + 108);
    const auto *fh_109 = buffer.data(fh + 109);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_111 = buffer.data(fh + 111);
    const auto *fh_112 = buffer.data(fh + 112);
    const auto *fh_113 = buffer.data(fh + 113);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, dh_0, fg_0, fh_0, \
                         fh_1, fh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dh_0[k]
                 + f_1 * fg_0[k]
                 + pb_x[k] * fh_0[k];

        t_1[k] = pb_y[k] * fh_0[k];

        t_2[k] = pb_z[k] * fh_0[k];

        t_3[k] = f_2 * fg_0[k]
                 + pb_y[k] * fh_1[k];

        t_4[k] = pb_y[k] * fh_2[k];

        t_5[k] = f_2 * fg_0[k]
                 + pb_z[k] * fh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, fg_1, fg_2, fg_3, fh_3, \
                         fh_5, fh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * fg_1[k]
                 + pb_y[k] * fh_3[k];

        t_7[k] = pb_z[k] * fh_3[k];

        t_8[k] = pb_y[k] * fh_5[k];

        t_9[k] = f_3 * fg_2[k]
                 + pb_z[k] * fh_5[k];

        t_10[k] = f_0 * fg_3[k]
                  + pb_y[k] * fh_6[k];

        t_11[k] = pb_z[k] * fh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, dh_15, fg_5, fh_8, \
                         fh_9, fh_10, fh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * fg_5[k]
                  + pb_y[k] * fh_8[k];

        t_13[k] = pb_y[k] * fh_9[k];

        t_14[k] = f_0 * fg_5[k]
                  + pb_z[k] * fh_9[k];

        t_15[k] = f_0 * dh_15[k]
                  + pb_x[k] * fh_15[k];

        t_16[k] = pb_z[k] * fh_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, dh_17, dh_18, dh_20, fg_10, \
                         fh_14, fh_15, fh_17, fh_18, fh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * dh_17[k]
                  + pb_x[k] * fh_17[k];

        t_18[k] = f_0 * dh_18[k]
                  + pb_x[k] * fh_18[k];

        t_19[k] = pb_y[k] * fh_14[k];

        t_20[k] = f_0 * dh_20[k]
                  + pb_x[k] * fh_20[k];

        t_21[k] = f_1 * fg_10[k]
                  + pb_y[k] * fh_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, fg_12, fg_13, fg_14, \
                         fh_15, fh_17, fh_18, fh_19, fh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * fh_15[k];

        t_23[k] = f_0 * fg_12[k]
                  + pb_y[k] * fh_17[k];

        t_24[k] = f_3 * fg_13[k]
                  + pb_y[k] * fh_18[k];

        t_25[k] = f_2 * fg_14[k]
                  + pb_y[k] * fh_19[k];

        t_26[k] = pb_y[k] * fh_20[k];

        t_27[k] = f_1 * fg_14[k]
                  + pb_z[k] * fh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pa_y, pb_y, pb_z, dh_0, dh_1, \
                         di_0, di_3, di_5, fh_21, fh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * di_0[k];

        t_29[k] = f_2 * dh_0[k]
                  + pb_y[k] * fh_21[k];

        t_30[k] = pb_z[k] * fh_21[k];

        t_31[k] = f_3 * dh_1[k]
                  + pa_y[k] * di_3[k];

        t_32[k] = pb_z[k] * fh_22[k];

        t_33[k] = pa_y[k] * di_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, dh_3, dh_5, dh_6, \
                         di_6, di_9, di_10, fh_24, fh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_0 * dh_3[k]
                  + pa_y[k] * di_6[k];

        t_35[k] = pb_z[k] * fh_24[k];

        t_36[k] = f_2 * dh_5[k]
                  + pb_y[k] * fh_26[k];

        t_37[k] = pa_y[k] * di_9[k];

        t_38[k] = f_4 * dh_6[k]
                  + pa_y[k] * di_10[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_y, pb_z, dh_9, di_14, fg_18, fh_27, \
                         fh_28, fh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_z[k] * fh_27[k];

        t_40[k] = f_2 * fg_18[k]
                  + pb_z[k] * fh_28[k];

        t_41[k] = f_2 * dh_9[k]
                  + pb_y[k] * fh_30[k];

        t_42[k] = pa_y[k] * di_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_z, dh_36, dh_38, dh_39, dh_40, \
                         fh_31, fh_36, fh_38, fh_39, fh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * dh_36[k]
                  + pb_x[k] * fh_36[k];

        t_44[k] = pb_z[k] * fh_31[k];

        t_45[k] = f_3 * dh_38[k]
                  + pb_x[k] * fh_38[k];

        t_46[k] = f_3 * dh_39[k]
                  + pb_x[k] * fh_39[k];

        t_47[k] = f_3 * dh_40[k]
                  + pb_x[k] * fh_40[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pa_y, pb_z, pi_49, di_20, di_49, \
                         fg_25, fg_26, fh_36, fh_37, fh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * di_20[k];

        t_49[k] = f_2 * pi_49[k]
                  + pa_x[k] * di_49[k];

        t_50[k] = pb_z[k] * fh_36[k];

        t_51[k] = f_2 * fg_25[k]
                  + pb_z[k] * fh_37[k];

        t_52[k] = f_3 * fg_26[k]
                  + pb_z[k] * fh_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pa_z, pb_y, pb_z, dh_20, di_0, \
                         di_27, fg_27, fh_39, fh_41, fh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * fg_27[k]
                  + pb_z[k] * fh_39[k];

        t_54[k] = f_2 * dh_20[k]
                  + pb_y[k] * fh_41[k];

        t_55[k] = pa_y[k] * di_27[k];

        t_56[k] = pa_z[k] * di_0[k];

        t_57[k] = pb_y[k] * fh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, dh_0, dh_2, di_3, \
                         di_5, di_6, fh_42, fh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * dh_0[k]
                  + pb_z[k] * fh_42[k];

        t_59[k] = pa_z[k] * di_3[k];

        t_60[k] = pb_y[k] * fh_44[k];

        t_61[k] = f_3 * dh_2[k]
                  + pa_z[k] * di_5[k];

        t_62[k] = pa_z[k] * di_6[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_z, pb_y, dh_5, di_9, di_10, fg_32, \
                         fg_34, fh_46, fh_47, fh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * fg_32[k]
                  + pb_y[k] * fh_46[k];

        t_64[k] = pb_y[k] * fh_47[k];

        t_65[k] = f_0 * dh_5[k]
                  + pa_z[k] * di_9[k];

        t_66[k] = pa_z[k] * di_10[k];

        t_67[k] = f_3 * fg_34[k]
                  + pb_y[k] * fh_49[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_x, pb_y, dh_9, dh_58, di_14, \
                         di_15, fg_35, fh_50, fh_51, fh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * fg_35[k]
                  + pb_y[k] * fh_50[k];

        t_69[k] = pb_y[k] * fh_51[k];

        t_70[k] = f_4 * dh_9[k]
                  + pa_z[k] * di_14[k];

        t_71[k] = pa_z[k] * di_15[k];

        t_72[k] = f_3 * dh_58[k]
                  + pb_x[k] * fh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, dh_59, dh_60, dh_62, \
                         di_21, fh_56, fh_59, fh_60, fh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * dh_59[k]
                  + pb_x[k] * fh_59[k];

        t_74[k] = f_3 * dh_60[k]
                  + pb_x[k] * fh_60[k];

        t_75[k] = pb_y[k] * fh_56[k];

        t_76[k] = f_3 * dh_62[k]
                  + pb_x[k] * fh_62[k];

        t_77[k] = pa_z[k] * di_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_y, fg_41, fg_42, fg_43, fg_44, \
                         fh_58, fh_59, fh_60, fh_61, fh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * fg_41[k]
                  + pb_y[k] * fh_58[k];

        t_79[k] = f_0 * fg_42[k]
                  + pb_y[k] * fh_59[k];

        t_80[k] = f_3 * fg_43[k]
                  + pb_y[k] * fh_60[k];

        t_81[k] = f_2 * fg_44[k]
                  + pb_y[k] * fh_61[k];

        t_82[k] = pb_y[k] * fh_62[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_x, pb_y, pb_z, pi_83, dh_21, dh_63, \
                         dh_66, di_83, di_84, di_87, fh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * pi_83[k]
                  + pa_x[k] * di_83[k];

        t_84[k] = f_5 * dh_63[k]
                  + pa_x[k] * di_84[k];

        t_85[k] = f_3 * dh_21[k]
                  + pb_y[k] * fh_63[k];

        t_86[k] = pb_z[k] * fh_63[k];

        t_87[k] = f_4 * dh_66[k]
                  + pa_x[k] * di_87[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_x, pb_y, pb_z, dh_26, dh_69, di_90, \
                         fg_45, fh_64, fh_65, fh_66, fh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_z[k] * fh_64[k];

        t_89[k] = f_2 * fg_45[k]
                  + pb_z[k] * fh_65[k];

        t_90[k] = f_0 * dh_69[k]
                  + pa_x[k] * di_90[k];

        t_91[k] = pb_z[k] * fh_66[k];

        t_92[k] = f_3 * dh_26[k]
                  + pb_y[k] * fh_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_x, pb_z, dh_73, di_94, fg_47, fg_48, \
                         fh_68, fh_69, fh_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_3 * fg_47[k]
                  + pb_z[k] * fh_68[k];

        t_94[k] = f_3 * dh_73[k]
                  + pa_x[k] * di_94[k];

        t_95[k] = pb_z[k] * fh_69[k];

        t_96[k] = f_2 * fg_48[k]
                  + pb_z[k] * fh_70[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pb_x, pb_y, pb_z, dh_30, dh_78, \
                         dh_80, fg_50, fh_72, fh_73, fh_78, fh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * dh_30[k]
                  + pb_y[k] * fh_72[k];

        t_98[k] = f_0 * fg_50[k]
                  + pb_z[k] * fh_72[k];

        t_99[k] = f_2 * dh_78[k]
                  + pb_x[k] * fh_78[k];

        t_100[k] = pb_z[k] * fh_73[k];

        t_101[k] = f_2 * dh_80[k]
                   + pb_x[k] * fh_80[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pa_x, pb_x, pb_z, dh_81, dh_82, \
                         dh_83, di_105, fh_78, fh_81, fh_82, fh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_2 * dh_81[k]
                   + pb_x[k] * fh_81[k];

        t_103[k] = f_2 * dh_82[k]
                   + pb_x[k] * fh_82[k];

        t_104[k] = f_2 * dh_83[k]
                   + pb_x[k] * fh_83[k];

        t_105[k] = pa_x[k] * di_105[k];

        t_106[k] = pb_z[k] * fh_78[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, pa_x, pa_y, di_56, di_107, \
                         di_108, di_109, di_110, di_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * di_107[k];

        t_108[k] = pa_x[k] * di_108[k];

        t_109[k] = pa_x[k] * di_109[k];

        t_110[k] = pa_x[k] * di_110[k];

        t_111[k] = pa_x[k] * di_111[k];

        t_112[k] = pa_y[k] * di_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_y, pa_z, pb_y, dh_44, \
                         di_29, di_31, di_34, di_58, di_61, fh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * di_29[k];

        t_114[k] = pa_y[k] * di_58[k];

        t_115[k] = pa_z[k] * di_31[k];

        t_116[k] = f_2 * dh_44[k]
                   + pb_y[k] * fh_86[k];

        t_117[k] = pa_y[k] * di_61[k];

        t_118[k] = pa_z[k] * di_34[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_y, pa_z, pb_y, pb_z, dh_24, dh_47, \
                         di_38, di_65, fh_87, fh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * dh_24[k]
                   + pb_z[k] * fh_87[k];

        t_120[k] = f_2 * dh_47[k]
                   + pb_y[k] * fh_89[k];

        t_121[k] = pa_y[k] * di_65[k];

        t_122[k] = pa_z[k] * di_38[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pa_y, pb_y, pb_z, dh_27, dh_51, \
                         dh_96, di_70, di_124, fh_90, fh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_2 * dh_27[k]
                   + pb_z[k] * fh_90[k];

        t_124[k] = f_3 * dh_96[k]
                   + pa_x[k] * di_124[k];

        t_125[k] = f_2 * dh_51[k]
                   + pb_y[k] * fh_93[k];

        t_126[k] = pa_y[k] * di_70[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_z, pb_x, dh_100, dh_101, \
                         dh_102, dh_103, di_43, fh_100, fh_101, fh_102, \
                         fh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_z[k] * di_43[k];

        t_128[k] = f_2 * dh_100[k]
                   + pb_x[k] * fh_100[k];

        t_129[k] = f_2 * dh_101[k]
                   + pb_x[k] * fh_101[k];

        t_130[k] = f_2 * dh_102[k]
                   + pb_x[k] * fh_102[k];

        t_131[k] = f_2 * dh_103[k]
                   + pb_x[k] * fh_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, t_137, t_138, pa_x, pa_y, di_76, \
                         di_133, di_134, di_135, di_136, di_137, \
                         di_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * di_76[k];

        t_133[k] = pa_x[k] * di_133[k];

        t_134[k] = pa_x[k] * di_134[k];

        t_135[k] = pa_x[k] * di_135[k];

        t_136[k] = pa_x[k] * di_136[k];

        t_137[k] = pa_x[k] * di_137[k];

        t_138[k] = pa_x[k] * di_138[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pa_x, pb_y, pb_z, dh_42, dh_105, \
                         di_139, di_140, fg_75, fh_105, fh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_x[k] * di_139[k];

        t_140[k] = f_5 * dh_105[k]
                   + pa_x[k] * di_140[k];

        t_141[k] = pb_y[k] * fh_105[k];

        t_142[k] = f_3 * dh_42[k]
                   + pb_z[k] * fh_105[k];

        t_143[k] = f_2 * fg_75[k]
                   + pb_y[k] * fh_106[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pa_x, pb_y, dh_110, di_145, fg_76, \
                         fg_77, fh_107, fh_108, fh_109, fh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_y[k] * fh_107[k];

        t_145[k] = f_4 * dh_110[k]
                   + pa_x[k] * di_145[k];

        t_146[k] = f_3 * fg_76[k]
                   + pb_y[k] * fh_108[k];

        t_147[k] = f_2 * fg_77[k]
                   + pb_y[k] * fh_109[k];

        t_148[k] = pb_y[k] * fh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, pa_x, pb_y, dh_114, di_149, fg_78, \
                         fg_79, fg_80, fh_111, fh_112, fh_113, fh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * dh_114[k]
                   + pa_x[k] * di_149[k];

        t_150[k] = f_0 * fg_78[k]
                   + pb_y[k] * fh_111[k];

        t_151[k] = f_3 * fg_79[k]
                   + pb_y[k] * fh_112[k];

        t_152[k] = f_2 * fg_80[k]
                   + pb_y[k] * fh_113[k];

        t_153[k] = pb_y[k] * fh_114[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, dh_119, dh_120, dh_121, \
                         dh_122, di_154, fh_120, fh_121, fh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_3 * dh_119[k]
                   + pa_x[k] * di_154[k];

        t_155[k] = f_2 * dh_120[k]
                   + pb_x[k] * fh_120[k];

        t_156[k] = f_2 * dh_121[k]
                   + pb_x[k] * fh_121[k];

        t_157[k] = f_2 * dh_122[k]
                   + pb_x[k] * fh_122[k];
    }
}

static auto
compute_prim_fi_overlap_0_piece1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pi, const size_t dh,
                                 const size_t di, const size_t fg, const size_t fh,
                                 const size_t ncols, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pi_83 = buffer.data(pi + 83);

    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *di_84 = buffer.data(di + 84);
    const auto *di_85 = buffer.data(di + 85);
    const auto *di_87 = buffer.data(di + 87);
    const auto *di_90 = buffer.data(di + 90);
    const auto *di_94 = buffer.data(di + 94);
    const auto *di_105 = buffer.data(di + 105);
    const auto *di_107 = buffer.data(di + 107);
    const auto *di_108 = buffer.data(di + 108);
    const auto *di_109 = buffer.data(di + 109);
    const auto *di_139 = buffer.data(di + 139);
    const auto *di_140 = buffer.data(di + 140);
    const auto *di_142 = buffer.data(di + 142);
    const auto *di_145 = buffer.data(di + 145);
    const auto *di_149 = buffer.data(di + 149);
    const auto *di_154 = buffer.data(di + 154);
    const auto *di_161 = buffer.data(di + 161);
    const auto *di_162 = buffer.data(di + 162);
    const auto *di_163 = buffer.data(di + 163);
    const auto *di_164 = buffer.data(di + 164);
    const auto *di_165 = buffer.data(di + 165);
    const auto *di_167 = buffer.data(di + 167);

    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_98 = buffer.data(fg + 98);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_107 = buffer.data(fg + 107);
    const auto *fg_109 = buffer.data(fg + 109);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_112 = buffer.data(fg + 112);
    const auto *fg_113 = buffer.data(fg + 113);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_121 = buffer.data(fg + 121);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_124 = buffer.data(fg + 124);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_127 = buffer.data(fg + 127);
    const auto *fg_128 = buffer.data(fg + 128);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_158 = buffer.data(fh + 158);
    const auto *fh_159 = buffer.data(fh + 159);
    const auto *fh_160 = buffer.data(fh + 160);
    const auto *fh_161 = buffer.data(fh + 161);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_x, pb_x, pb_y, dh_123, dh_125, \
                         di_161, di_162, fh_119, fh_123, fh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_2 * dh_123[k]
                   + pb_x[k] * fh_123[k];

        t_159[k] = pb_y[k] * fh_119[k];

        t_160[k] = f_2 * dh_125[k]
                   + pb_x[k] * fh_125[k];

        t_161[k] = pa_x[k] * di_161[k];

        t_162[k] = pa_x[k] * di_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pb_x, pb_y, di_163, \
                         di_164, di_165, di_167, fg_90, fh_125, \
                         fh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * di_163[k];

        t_164[k] = pa_x[k] * di_164[k];

        t_165[k] = pa_x[k] * di_165[k];

        t_166[k] = pb_y[k] * fh_125[k];

        t_167[k] = pa_x[k] * di_167[k];

        t_168[k] = f_1 * fg_90[k]
                   + pb_x[k] * fh_126[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pb_x, pb_z, fg_91, fg_93, fg_95, \
                         fh_126, fh_127, fh_129, fh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_4 * fg_91[k]
                   + pb_x[k] * fh_127[k];

        t_170[k] = pb_z[k] * fh_126[k];

        t_171[k] = f_0 * fg_93[k]
                   + pb_x[k] * fh_129[k];

        t_172[k] = pb_z[k] * fh_127[k];

        t_173[k] = f_0 * fg_95[k]
                   + pb_x[k] * fh_131[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pb_x, pb_z, fg_96, fg_98, fg_99, \
                         fg_100, fh_129, fh_132, fh_134, fh_135, \
                         fh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_3 * fg_96[k]
                   + pb_x[k] * fh_132[k];

        t_175[k] = pb_z[k] * fh_129[k];

        t_176[k] = f_3 * fg_98[k]
                   + pb_x[k] * fh_134[k];

        t_177[k] = f_3 * fg_99[k]
                   + pb_x[k] * fh_135[k];

        t_178[k] = f_2 * fg_100[k]
                   + pb_x[k] * fh_136[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, t_183, pb_x, pb_z, fg_102, fg_103, \
                         fg_104, fh_132, fh_138, fh_139, fh_140, \
                         fh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = pb_z[k] * fh_132[k];

        t_180[k] = f_2 * fg_102[k]
                   + pb_x[k] * fh_138[k];

        t_181[k] = f_2 * fg_103[k]
                   + pb_x[k] * fh_139[k];

        t_182[k] = f_2 * fg_104[k]
                   + pb_x[k] * fh_140[k];

        t_183[k] = pb_x[k] * fh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, t_189, pb_x, pb_y, dh_78, fg_100, \
                         fh_141, fh_142, fh_143, fh_144, fh_145, \
                         fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_x[k] * fh_142[k];

        t_185[k] = pb_x[k] * fh_143[k];

        t_186[k] = pb_x[k] * fh_144[k];

        t_187[k] = pb_x[k] * fh_145[k];

        t_188[k] = pb_x[k] * fh_146[k];

        t_189[k] = f_0 * dh_78[k]
                   + f_1 * fg_100[k]
                   + pb_y[k] * fh_141[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pb_y, pb_z, dh_83, fg_100, fg_101, \
                         fg_102, fh_141, fh_142, fh_143, fh_144, \
                         fh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pb_z[k] * fh_141[k];

        t_191[k] = f_2 * fg_100[k]
                   + pb_z[k] * fh_142[k];

        t_192[k] = f_3 * fg_101[k]
                   + pb_z[k] * fh_143[k];

        t_193[k] = f_0 * fg_102[k]
                   + pb_z[k] * fh_144[k];

        t_194[k] = f_0 * dh_83[k]
                   + pb_y[k] * fh_146[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pa_z, pb_x, pb_z, di_84, di_85, \
                         di_87, fg_104, fg_107, fh_146, fh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * fg_104[k]
                   + pb_z[k] * fh_146[k];

        t_196[k] = pa_z[k] * di_84[k];

        t_197[k] = pa_z[k] * di_85[k];

        t_198[k] = f_4 * fg_107[k]
                   + pb_x[k] * fh_149[k];

        t_199[k] = pa_z[k] * di_87[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pa_z, pb_x, di_90, fg_109, fg_110, \
                         fg_112, fg_113, fh_151, fh_152, fh_154, \
                         fh_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_0 * fg_109[k]
                   + pb_x[k] * fh_151[k];

        t_201[k] = f_0 * fg_110[k]
                   + pb_x[k] * fh_152[k];

        t_202[k] = pa_z[k] * di_90[k];

        t_203[k] = f_3 * fg_112[k]
                   + pb_x[k] * fh_154[k];

        t_204[k] = f_3 * fg_113[k]
                   + pb_x[k] * fh_155[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, pa_z, pb_x, di_94, fg_114, fg_116, \
                         fg_117, fg_118, fh_156, fh_158, fh_159, \
                         fh_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_3 * fg_114[k]
                   + pb_x[k] * fh_156[k];

        t_206[k] = pa_z[k] * di_94[k];

        t_207[k] = f_2 * fg_116[k]
                   + pb_x[k] * fh_158[k];

        t_208[k] = f_2 * fg_117[k]
                   + pb_x[k] * fh_159[k];

        t_209[k] = f_2 * fg_118[k]
                   + pb_x[k] * fh_160[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, t_215, t_216, pb_x, fg_119, \
                         fh_161, fh_162, fh_163, fh_164, fh_165, fh_166, \
                         fh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_2 * fg_119[k]
                   + pb_x[k] * fh_161[k];

        t_211[k] = pb_x[k] * fh_162[k];

        t_212[k] = pb_x[k] * fh_163[k];

        t_213[k] = pb_x[k] * fh_164[k];

        t_214[k] = pb_x[k] * fh_165[k];

        t_215[k] = pb_x[k] * fh_166[k];

        t_216[k] = pb_x[k] * fh_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_z, dh_78, dh_79, dh_80, \
                         dh_81, di_105, di_107, di_108, di_109, \
                         fh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_z[k] * di_105[k];

        t_218[k] = f_2 * dh_78[k]
                   + pb_z[k] * fh_162[k];

        t_219[k] = f_3 * dh_79[k]
                   + pa_z[k] * di_107[k];

        t_220[k] = f_0 * dh_80[k]
                   + pa_z[k] * di_108[k];

        t_221[k] = f_4 * dh_81[k]
                   + pa_z[k] * di_109[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pa_y, pb_x, pb_y, pi_83, dh_104, \
                         di_139, di_140, di_142, fg_121, fh_167, \
                         fh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_3 * dh_104[k]
                   + pb_y[k] * fh_167[k];

        t_223[k] = f_2 * pi_83[k]
                   + pa_y[k] * di_139[k];

        t_224[k] = pa_y[k] * di_140[k];

        t_225[k] = f_4 * fg_121[k]
                   + pb_x[k] * fh_169[k];

        t_226[k] = pa_y[k] * di_142[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, pa_y, pb_x, di_145, fg_123, \
                         fg_124, fg_126, fg_127, fh_171, fh_172, fh_174, \
                         fh_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_0 * fg_123[k]
                   + pb_x[k] * fh_171[k];

        t_228[k] = f_0 * fg_124[k]
                   + pb_x[k] * fh_172[k];

        t_229[k] = pa_y[k] * di_145[k];

        t_230[k] = f_3 * fg_126[k]
                   + pb_x[k] * fh_174[k];

        t_231[k] = f_3 * fg_127[k]
                   + pb_x[k] * fh_175[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pa_y, pb_x, di_149, fg_128, \
                         fg_130, fg_131, fg_132, fh_176, fh_178, fh_179, \
                         fh_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * fg_128[k]
                   + pb_x[k] * fh_176[k];

        t_233[k] = pa_y[k] * di_149[k];

        t_234[k] = f_2 * fg_130[k]
                   + pb_x[k] * fh_178[k];

        t_235[k] = f_2 * fg_131[k]
                   + pb_x[k] * fh_179[k];

        t_236[k] = f_2 * fg_132[k]
                   + pb_x[k] * fh_180[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, t_242, pa_y, pb_x, di_154, fg_133, \
                         fh_181, fh_183, fh_184, fh_185, fh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_2 * fg_133[k]
                   + pb_x[k] * fh_181[k];

        t_238[k] = pa_y[k] * di_154[k];

        t_239[k] = pb_x[k] * fh_183[k];

        t_240[k] = pb_x[k] * fh_184[k];

        t_241[k] = pb_x[k] * fh_185[k];

        t_242[k] = pb_x[k] * fh_186[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, pa_y, pb_x, pb_z, dh_99, dh_120, \
                         dh_122, di_161, di_163, fh_183, fh_187, \
                         fh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pb_x[k] * fh_187[k];

        t_244[k] = pb_x[k] * fh_188[k];

        t_245[k] = f_5 * dh_120[k]
                   + pa_y[k] * di_161[k];

        t_246[k] = f_3 * dh_99[k]
                   + pb_z[k] * fh_183[k];

        t_247[k] = f_4 * dh_122[k]
                   + pa_y[k] * di_163[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_y, pb_y, dh_123, dh_124, dh_125, \
                         di_164, di_165, di_167, fh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_0 * dh_123[k]
                   + pa_y[k] * di_164[k];

        t_249[k] = f_3 * dh_124[k]
                   + pa_y[k] * di_165[k];

        t_250[k] = f_2 * dh_125[k]
                   + pb_y[k] * fh_188[k];

        t_251[k] = pa_y[k] * di_167[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, t_257, pb_x, pb_y, fg_135, fg_137, \
                         fg_138, fg_140, fh_189, fh_191, fh_192, \
                         fh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * fg_135[k]
                   + pb_x[k] * fh_189[k];

        t_253[k] = pb_y[k] * fh_189[k];

        t_254[k] = f_4 * fg_137[k]
                   + pb_x[k] * fh_191[k];

        t_255[k] = f_0 * fg_138[k]
                   + pb_x[k] * fh_192[k];

        t_256[k] = pb_y[k] * fh_191[k];

        t_257[k] = f_0 * fg_140[k]
                   + pb_x[k] * fh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pb_x, pb_y, fg_141, fg_142, \
                         fg_144, fg_145, fh_194, fh_195, fh_196, fh_198, \
                         fh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_3 * fg_141[k]
                   + pb_x[k] * fh_195[k];

        t_259[k] = f_3 * fg_142[k]
                   + pb_x[k] * fh_196[k];

        t_260[k] = pb_y[k] * fh_194[k];

        t_261[k] = f_3 * fg_144[k]
                   + pb_x[k] * fh_198[k];

        t_262[k] = f_2 * fg_145[k]
                   + pb_x[k] * fh_199[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, pb_x, pb_y, fg_146, fg_147, \
                         fg_149, fh_198, fh_200, fh_201, fh_203, \
                         fh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_2 * fg_146[k]
                   + pb_x[k] * fh_200[k];

        t_264[k] = f_2 * fg_147[k]
                   + pb_x[k] * fh_201[k];

        t_265[k] = pb_y[k] * fh_198[k];

        t_266[k] = f_2 * fg_149[k]
                   + pb_x[k] * fh_203[k];

        t_267[k] = pb_x[k] * fh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, t_273, pb_x, pb_y, fg_145, fh_204, \
                         fh_205, fh_206, fh_207, fh_208, fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = pb_x[k] * fh_205[k];

        t_269[k] = pb_x[k] * fh_206[k];

        t_270[k] = pb_x[k] * fh_207[k];

        t_271[k] = pb_x[k] * fh_208[k];

        t_272[k] = pb_x[k] * fh_209[k];

        t_273[k] = f_1 * fg_145[k]
                   + pb_y[k] * fh_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_y, fg_146, fg_147, fg_148, \
                         fg_149, fh_205, fh_206, fh_207, fh_208, \
                         fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_4 * fg_146[k]
                   + pb_y[k] * fh_205[k];

        t_275[k] = f_0 * fg_147[k]
                   + pb_y[k] * fh_206[k];

        t_276[k] = f_3 * fg_148[k]
                   + pb_y[k] * fh_207[k];

        t_277[k] = f_2 * fg_149[k]
                   + pb_y[k] * fh_208[k];

        t_278[k] = pb_y[k] * fh_209[k];
    }

#pragma omp simd aligned(t_279, pb_z, dh_125, fg_149, fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_0 * dh_125[k]
                   + f_1 * fg_149[k]
                   + pb_z[k] * fh_209[k];
    }
}

auto
compute_prim_fi_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t pi, const size_t dh, const size_t di,
                          const size_t fg, const size_t fh, const size_t ncols,
                          const double p) -> void
{
    compute_prim_fi_overlap_0_piece0(buffer, target, pa, pb, pi, dh, di, fg, fh, ncols, p);

    compute_prim_fi_overlap_0_piece1(buffer, target, pa, pb, pi, dh, di, fg, fh, ncols, p);
}

}  // namespace simdovl
