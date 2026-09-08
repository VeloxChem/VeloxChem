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


#include "SimdTransferHD.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t hp, const size_t ip, const size_t nmax) -> void
{
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_48 = buffer.data(hp + 48);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_57 = buffer.data(hp + 57);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_61 = buffer.data(ip + 61);
    const auto *ip_62 = buffer.data(ip + 62);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_67 = buffer.data(ip + 67);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_83 = buffer.data(ip + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, hp_0, hp_1, hp_2, ip_0, ip_1, \
                         ip_2, ip_4, ip_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * hp_0[k]
                 + ip_0[k];

        t_1[k] = ab_x[k] * hp_1[k]
                 + ip_1[k];

        t_2[k] = ab_x[k] * hp_2[k]
                 + ip_2[k];

        t_3[k] = ab_y[k] * hp_1[k]
                 + ip_4[k];

        t_4[k] = ab_y[k] * hp_2[k]
                 + ip_5[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, hp_2, hp_3, hp_4, hp_5, ip_3, ip_4, \
                         ip_5, ip_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_z[k] * hp_2[k]
                 + ip_8[k];

        t_6[k] = ab_x[k] * hp_3[k]
                 + ip_3[k];

        t_7[k] = ab_x[k] * hp_4[k]
                 + ip_4[k];

        t_8[k] = ab_x[k] * hp_5[k]
                 + ip_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, hp_4, hp_5, hp_6, ip_6, \
                         ip_10, ip_11, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = ab_y[k] * hp_4[k]
                 + ip_10[k];

        t_10[k] = ab_y[k] * hp_5[k]
                  + ip_11[k];

        t_11[k] = ab_z[k] * hp_5[k]
                  + ip_14[k];

        t_12[k] = ab_x[k] * hp_6[k]
                  + ip_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, hp_7, hp_8, ip_7, \
                         ip_8, ip_13, ip_14, ip_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = ab_x[k] * hp_7[k]
                  + ip_7[k];

        t_14[k] = ab_x[k] * hp_8[k]
                  + ip_8[k];

        t_15[k] = ab_y[k] * hp_7[k]
                  + ip_13[k];

        t_16[k] = ab_y[k] * hp_8[k]
                  + ip_14[k];

        t_17[k] = ab_z[k] * hp_8[k]
                  + ip_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, hp_9, hp_10, hp_11, ip_9, \
                         ip_10, ip_11, ip_19, ip_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_18[k] = ab_x[k] * hp_9[k]
                  + ip_9[k];

        t_19[k] = ab_x[k] * hp_10[k]
                  + ip_10[k];

        t_20[k] = ab_x[k] * hp_11[k]
                  + ip_11[k];

        t_21[k] = ab_y[k] * hp_10[k]
                  + ip_19[k];

        t_22[k] = ab_y[k] * hp_11[k]
                  + ip_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, hp_11, hp_12, hp_13, hp_14, \
                         ip_12, ip_13, ip_14, ip_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_23[k] = ab_z[k] * hp_11[k]
                  + ip_23[k];

        t_24[k] = ab_x[k] * hp_12[k]
                  + ip_12[k];

        t_25[k] = ab_x[k] * hp_13[k]
                  + ip_13[k];

        t_26[k] = ab_x[k] * hp_14[k]
                  + ip_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, hp_13, hp_14, hp_15, ip_15, \
                         ip_22, ip_23, ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_y[k] * hp_13[k]
                  + ip_22[k];

        t_28[k] = ab_y[k] * hp_14[k]
                  + ip_23[k];

        t_29[k] = ab_z[k] * hp_14[k]
                  + ip_26[k];

        t_30[k] = ab_x[k] * hp_15[k]
                  + ip_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, hp_16, hp_17, ip_16, \
                         ip_17, ip_25, ip_26, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_31[k] = ab_x[k] * hp_16[k]
                  + ip_16[k];

        t_32[k] = ab_x[k] * hp_17[k]
                  + ip_17[k];

        t_33[k] = ab_y[k] * hp_16[k]
                  + ip_25[k];

        t_34[k] = ab_y[k] * hp_17[k]
                  + ip_26[k];

        t_35[k] = ab_z[k] * hp_17[k]
                  + ip_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, hp_18, hp_19, hp_20, ip_18, \
                         ip_19, ip_20, ip_31, ip_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_36[k] = ab_x[k] * hp_18[k]
                  + ip_18[k];

        t_37[k] = ab_x[k] * hp_19[k]
                  + ip_19[k];

        t_38[k] = ab_x[k] * hp_20[k]
                  + ip_20[k];

        t_39[k] = ab_y[k] * hp_19[k]
                  + ip_31[k];

        t_40[k] = ab_y[k] * hp_20[k]
                  + ip_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, hp_20, hp_21, hp_22, hp_23, \
                         ip_21, ip_22, ip_23, ip_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_41[k] = ab_z[k] * hp_20[k]
                  + ip_35[k];

        t_42[k] = ab_x[k] * hp_21[k]
                  + ip_21[k];

        t_43[k] = ab_x[k] * hp_22[k]
                  + ip_22[k];

        t_44[k] = ab_x[k] * hp_23[k]
                  + ip_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, hp_22, hp_23, hp_24, ip_24, \
                         ip_34, ip_35, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_y[k] * hp_22[k]
                  + ip_34[k];

        t_46[k] = ab_y[k] * hp_23[k]
                  + ip_35[k];

        t_47[k] = ab_z[k] * hp_23[k]
                  + ip_38[k];

        t_48[k] = ab_x[k] * hp_24[k]
                  + ip_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, hp_25, hp_26, ip_25, \
                         ip_26, ip_37, ip_38, ip_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_x[k] * hp_25[k]
                  + ip_25[k];

        t_50[k] = ab_x[k] * hp_26[k]
                  + ip_26[k];

        t_51[k] = ab_y[k] * hp_25[k]
                  + ip_37[k];

        t_52[k] = ab_y[k] * hp_26[k]
                  + ip_38[k];

        t_53[k] = ab_z[k] * hp_26[k]
                  + ip_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, hp_27, hp_28, hp_29, ip_27, \
                         ip_28, ip_29, ip_40, ip_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * hp_27[k]
                  + ip_27[k];

        t_55[k] = ab_x[k] * hp_28[k]
                  + ip_28[k];

        t_56[k] = ab_x[k] * hp_29[k]
                  + ip_29[k];

        t_57[k] = ab_y[k] * hp_28[k]
                  + ip_40[k];

        t_58[k] = ab_y[k] * hp_29[k]
                  + ip_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, hp_29, hp_30, hp_31, hp_32, \
                         ip_30, ip_31, ip_32, ip_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_59[k] = ab_z[k] * hp_29[k]
                  + ip_44[k];

        t_60[k] = ab_x[k] * hp_30[k]
                  + ip_30[k];

        t_61[k] = ab_x[k] * hp_31[k]
                  + ip_31[k];

        t_62[k] = ab_x[k] * hp_32[k]
                  + ip_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, hp_31, hp_32, hp_33, ip_33, \
                         ip_46, ip_47, ip_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_63[k] = ab_y[k] * hp_31[k]
                  + ip_46[k];

        t_64[k] = ab_y[k] * hp_32[k]
                  + ip_47[k];

        t_65[k] = ab_z[k] * hp_32[k]
                  + ip_50[k];

        t_66[k] = ab_x[k] * hp_33[k]
                  + ip_33[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, hp_34, hp_35, ip_34, \
                         ip_35, ip_49, ip_50, ip_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_67[k] = ab_x[k] * hp_34[k]
                  + ip_34[k];

        t_68[k] = ab_x[k] * hp_35[k]
                  + ip_35[k];

        t_69[k] = ab_y[k] * hp_34[k]
                  + ip_49[k];

        t_70[k] = ab_y[k] * hp_35[k]
                  + ip_50[k];

        t_71[k] = ab_z[k] * hp_35[k]
                  + ip_53[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, hp_36, hp_37, hp_38, ip_36, \
                         ip_37, ip_38, ip_52, ip_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_72[k] = ab_x[k] * hp_36[k]
                  + ip_36[k];

        t_73[k] = ab_x[k] * hp_37[k]
                  + ip_37[k];

        t_74[k] = ab_x[k] * hp_38[k]
                  + ip_38[k];

        t_75[k] = ab_y[k] * hp_37[k]
                  + ip_52[k];

        t_76[k] = ab_y[k] * hp_38[k]
                  + ip_53[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, hp_38, hp_39, hp_40, hp_41, \
                         ip_39, ip_40, ip_41, ip_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * hp_38[k]
                  + ip_56[k];

        t_78[k] = ab_x[k] * hp_39[k]
                  + ip_39[k];

        t_79[k] = ab_x[k] * hp_40[k]
                  + ip_40[k];

        t_80[k] = ab_x[k] * hp_41[k]
                  + ip_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, hp_40, hp_41, hp_42, ip_42, \
                         ip_55, ip_56, ip_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_81[k] = ab_y[k] * hp_40[k]
                  + ip_55[k];

        t_82[k] = ab_y[k] * hp_41[k]
                  + ip_56[k];

        t_83[k] = ab_z[k] * hp_41[k]
                  + ip_59[k];

        t_84[k] = ab_x[k] * hp_42[k]
                  + ip_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, hp_43, hp_44, ip_43, \
                         ip_44, ip_58, ip_59, ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * hp_43[k]
                  + ip_43[k];

        t_86[k] = ab_x[k] * hp_44[k]
                  + ip_44[k];

        t_87[k] = ab_y[k] * hp_43[k]
                  + ip_58[k];

        t_88[k] = ab_y[k] * hp_44[k]
                  + ip_59[k];

        t_89[k] = ab_z[k] * hp_44[k]
                  + ip_62[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, ab_y, hp_45, hp_46, hp_47, ip_45, \
                         ip_46, ip_47, ip_64, ip_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * hp_45[k]
                  + ip_45[k];

        t_91[k] = ab_x[k] * hp_46[k]
                  + ip_46[k];

        t_92[k] = ab_x[k] * hp_47[k]
                  + ip_47[k];

        t_93[k] = ab_y[k] * hp_46[k]
                  + ip_64[k];

        t_94[k] = ab_y[k] * hp_47[k]
                  + ip_65[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, ab_x, ab_z, hp_47, hp_48, hp_49, hp_50, \
                         ip_48, ip_49, ip_50, ip_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_z[k] * hp_47[k]
                  + ip_68[k];

        t_96[k] = ab_x[k] * hp_48[k]
                  + ip_48[k];

        t_97[k] = ab_x[k] * hp_49[k]
                  + ip_49[k];

        t_98[k] = ab_x[k] * hp_50[k]
                  + ip_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, ab_x, ab_y, ab_z, hp_49, hp_50, hp_51, \
                         ip_51, ip_67, ip_68, ip_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = ab_y[k] * hp_49[k]
                  + ip_67[k];

        t_100[k] = ab_y[k] * hp_50[k]
                   + ip_68[k];

        t_101[k] = ab_z[k] * hp_50[k]
                   + ip_71[k];

        t_102[k] = ab_x[k] * hp_51[k]
                   + ip_51[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ab_x, ab_y, ab_z, hp_52, hp_53, \
                         ip_52, ip_53, ip_70, ip_71, ip_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_103[k] = ab_x[k] * hp_52[k]
                   + ip_52[k];

        t_104[k] = ab_x[k] * hp_53[k]
                   + ip_53[k];

        t_105[k] = ab_y[k] * hp_52[k]
                   + ip_70[k];

        t_106[k] = ab_y[k] * hp_53[k]
                   + ip_71[k];

        t_107[k] = ab_z[k] * hp_53[k]
                   + ip_74[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_x, ab_y, hp_54, hp_55, hp_56, \
                         ip_54, ip_55, ip_56, ip_73, ip_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_108[k] = ab_x[k] * hp_54[k]
                   + ip_54[k];

        t_109[k] = ab_x[k] * hp_55[k]
                   + ip_55[k];

        t_110[k] = ab_x[k] * hp_56[k]
                   + ip_56[k];

        t_111[k] = ab_y[k] * hp_55[k]
                   + ip_73[k];

        t_112[k] = ab_y[k] * hp_56[k]
                   + ip_74[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, ab_x, ab_z, hp_56, hp_57, hp_58, hp_59, \
                         ip_57, ip_58, ip_59, ip_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_113[k] = ab_z[k] * hp_56[k]
                   + ip_77[k];

        t_114[k] = ab_x[k] * hp_57[k]
                   + ip_57[k];

        t_115[k] = ab_x[k] * hp_58[k]
                   + ip_58[k];

        t_116[k] = ab_x[k] * hp_59[k]
                   + ip_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, ab_x, ab_y, ab_z, hp_58, hp_59, hp_60, \
                         ip_60, ip_76, ip_77, ip_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_117[k] = ab_y[k] * hp_58[k]
                   + ip_76[k];

        t_118[k] = ab_y[k] * hp_59[k]
                   + ip_77[k];

        t_119[k] = ab_z[k] * hp_59[k]
                   + ip_80[k];

        t_120[k] = ab_x[k] * hp_60[k]
                   + ip_60[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, ab_x, ab_y, ab_z, hp_61, hp_62, \
                         ip_61, ip_62, ip_79, ip_80, ip_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_121[k] = ab_x[k] * hp_61[k]
                   + ip_61[k];

        t_122[k] = ab_x[k] * hp_62[k]
                   + ip_62[k];

        t_123[k] = ab_y[k] * hp_61[k]
                   + ip_79[k];

        t_124[k] = ab_y[k] * hp_62[k]
                   + ip_80[k];

        t_125[k] = ab_z[k] * hp_62[k]
                   + ip_83[k];
    }
}

}  // namespace simdtrf
