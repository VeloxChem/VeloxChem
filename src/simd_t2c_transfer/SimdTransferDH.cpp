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


#include "SimdTransferDH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dh(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t ph, const size_t pi, const size_t nmax) -> void
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

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_16 = buffer.data(ph + 16);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_23 = buffer.data(ph + 23);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_26 = buffer.data(ph + 26);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_28 = buffer.data(ph + 28);
    const auto *ph_29 = buffer.data(ph + 29);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_32 = buffer.data(ph + 32);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_34 = buffer.data(ph + 34);
    const auto *ph_35 = buffer.data(ph + 35);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_37 = buffer.data(ph + 37);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_41 = buffer.data(ph + 41);
    const auto *ph_42 = buffer.data(ph + 42);
    const auto *ph_43 = buffer.data(ph + 43);
    const auto *ph_44 = buffer.data(ph + 44);
    const auto *ph_45 = buffer.data(ph + 45);
    const auto *ph_46 = buffer.data(ph + 46);
    const auto *ph_47 = buffer.data(ph + 47);
    const auto *ph_48 = buffer.data(ph + 48);
    const auto *ph_49 = buffer.data(ph + 49);
    const auto *ph_50 = buffer.data(ph + 50);
    const auto *ph_51 = buffer.data(ph + 51);
    const auto *ph_52 = buffer.data(ph + 52);
    const auto *ph_53 = buffer.data(ph + 53);
    const auto *ph_54 = buffer.data(ph + 54);
    const auto *ph_55 = buffer.data(ph + 55);
    const auto *ph_56 = buffer.data(ph + 56);
    const auto *ph_57 = buffer.data(ph + 57);
    const auto *ph_58 = buffer.data(ph + 58);
    const auto *ph_59 = buffer.data(ph + 59);
    const auto *ph_60 = buffer.data(ph + 60);
    const auto *ph_61 = buffer.data(ph + 61);
    const auto *ph_62 = buffer.data(ph + 62);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_1 = buffer.data(pi + 1);
    const auto *pi_2 = buffer.data(pi + 2);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_4 = buffer.data(pi + 4);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_7 = buffer.data(pi + 7);
    const auto *pi_8 = buffer.data(pi + 8);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_11 = buffer.data(pi + 11);
    const auto *pi_12 = buffer.data(pi + 12);
    const auto *pi_13 = buffer.data(pi + 13);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_16 = buffer.data(pi + 16);
    const auto *pi_17 = buffer.data(pi + 17);
    const auto *pi_18 = buffer.data(pi + 18);
    const auto *pi_19 = buffer.data(pi + 19);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_28 = buffer.data(pi + 28);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_30 = buffer.data(pi + 30);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_32 = buffer.data(pi + 32);
    const auto *pi_33 = buffer.data(pi + 33);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_35 = buffer.data(pi + 35);
    const auto *pi_36 = buffer.data(pi + 36);
    const auto *pi_37 = buffer.data(pi + 37);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_39 = buffer.data(pi + 39);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_41 = buffer.data(pi + 41);
    const auto *pi_42 = buffer.data(pi + 42);
    const auto *pi_43 = buffer.data(pi + 43);
    const auto *pi_44 = buffer.data(pi + 44);
    const auto *pi_45 = buffer.data(pi + 45);
    const auto *pi_46 = buffer.data(pi + 46);
    const auto *pi_47 = buffer.data(pi + 47);
    const auto *pi_48 = buffer.data(pi + 48);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_50 = buffer.data(pi + 50);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_57 = buffer.data(pi + 57);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_59 = buffer.data(pi + 59);
    const auto *pi_60 = buffer.data(pi + 60);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_62 = buffer.data(pi + 62);
    const auto *pi_63 = buffer.data(pi + 63);
    const auto *pi_64 = buffer.data(pi + 64);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_66 = buffer.data(pi + 66);
    const auto *pi_67 = buffer.data(pi + 67);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_69 = buffer.data(pi + 69);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_71 = buffer.data(pi + 71);
    const auto *pi_72 = buffer.data(pi + 72);
    const auto *pi_73 = buffer.data(pi + 73);
    const auto *pi_74 = buffer.data(pi + 74);
    const auto *pi_75 = buffer.data(pi + 75);
    const auto *pi_76 = buffer.data(pi + 76);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_82 = buffer.data(pi + 82);
    const auto *pi_83 = buffer.data(pi + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ph_0, ph_1, ph_2, ph_3, ph_4, pi_0, \
                         pi_1, pi_2, pi_3, pi_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * ph_0[k]
                 + pi_0[k];

        t_1[k] = -ab_x[k] * ph_1[k]
                 + pi_1[k];

        t_2[k] = -ab_x[k] * ph_2[k]
                 + pi_2[k];

        t_3[k] = -ab_x[k] * ph_3[k]
                 + pi_3[k];

        t_4[k] = -ab_x[k] * ph_4[k]
                 + pi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ph_5, ph_6, ph_7, ph_8, ph_9, pi_5, \
                         pi_6, pi_7, pi_8, pi_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * ph_5[k]
                 + pi_5[k];

        t_6[k] = -ab_x[k] * ph_6[k]
                 + pi_6[k];

        t_7[k] = -ab_x[k] * ph_7[k]
                 + pi_7[k];

        t_8[k] = -ab_x[k] * ph_8[k]
                 + pi_8[k];

        t_9[k] = -ab_x[k] * ph_9[k]
                 + pi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, ph_10, ph_11, ph_12, ph_13, \
                         ph_14, pi_10, pi_11, pi_12, pi_13, pi_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * ph_10[k]
                  + pi_10[k];

        t_11[k] = -ab_x[k] * ph_11[k]
                  + pi_11[k];

        t_12[k] = -ab_x[k] * ph_12[k]
                  + pi_12[k];

        t_13[k] = -ab_x[k] * ph_13[k]
                  + pi_13[k];

        t_14[k] = -ab_x[k] * ph_14[k]
                  + pi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ph_15, ph_16, ph_17, ph_18, \
                         ph_19, pi_15, pi_16, pi_17, pi_18, pi_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * ph_15[k]
                  + pi_15[k];

        t_16[k] = -ab_x[k] * ph_16[k]
                  + pi_16[k];

        t_17[k] = -ab_x[k] * ph_17[k]
                  + pi_17[k];

        t_18[k] = -ab_x[k] * ph_18[k]
                  + pi_18[k];

        t_19[k] = -ab_x[k] * ph_19[k]
                  + pi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, ph_20, ph_21, ph_22, ph_23, \
                         ph_24, pi_20, pi_28, pi_29, pi_30, pi_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * ph_20[k]
                  + pi_20[k];

        t_21[k] = -ab_x[k] * ph_21[k]
                  + pi_28[k];

        t_22[k] = -ab_x[k] * ph_22[k]
                  + pi_29[k];

        t_23[k] = -ab_x[k] * ph_23[k]
                  + pi_30[k];

        t_24[k] = -ab_x[k] * ph_24[k]
                  + pi_31[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ph_25, ph_26, ph_27, ph_28, \
                         ph_29, pi_32, pi_33, pi_34, pi_35, pi_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * ph_25[k]
                  + pi_32[k];

        t_26[k] = -ab_x[k] * ph_26[k]
                  + pi_33[k];

        t_27[k] = -ab_x[k] * ph_27[k]
                  + pi_34[k];

        t_28[k] = -ab_x[k] * ph_28[k]
                  + pi_35[k];

        t_29[k] = -ab_x[k] * ph_29[k]
                  + pi_36[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, ph_30, ph_31, ph_32, ph_33, \
                         ph_34, pi_37, pi_38, pi_39, pi_40, pi_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * ph_30[k]
                  + pi_37[k];

        t_31[k] = -ab_x[k] * ph_31[k]
                  + pi_38[k];

        t_32[k] = -ab_x[k] * ph_32[k]
                  + pi_39[k];

        t_33[k] = -ab_x[k] * ph_33[k]
                  + pi_40[k];

        t_34[k] = -ab_x[k] * ph_34[k]
                  + pi_41[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ph_35, ph_36, ph_37, ph_38, \
                         ph_39, pi_42, pi_43, pi_44, pi_45, pi_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * ph_35[k]
                  + pi_42[k];

        t_36[k] = -ab_x[k] * ph_36[k]
                  + pi_43[k];

        t_37[k] = -ab_x[k] * ph_37[k]
                  + pi_44[k];

        t_38[k] = -ab_x[k] * ph_38[k]
                  + pi_45[k];

        t_39[k] = -ab_x[k] * ph_39[k]
                  + pi_46[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, ph_40, ph_41, ph_42, ph_43, \
                         ph_44, pi_47, pi_48, pi_56, pi_57, pi_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * ph_40[k]
                  + pi_47[k];

        t_41[k] = -ab_x[k] * ph_41[k]
                  + pi_48[k];

        t_42[k] = -ab_x[k] * ph_42[k]
                  + pi_56[k];

        t_43[k] = -ab_x[k] * ph_43[k]
                  + pi_57[k];

        t_44[k] = -ab_x[k] * ph_44[k]
                  + pi_58[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ph_45, ph_46, ph_47, ph_48, \
                         ph_49, pi_59, pi_60, pi_61, pi_62, pi_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * ph_45[k]
                  + pi_59[k];

        t_46[k] = -ab_x[k] * ph_46[k]
                  + pi_60[k];

        t_47[k] = -ab_x[k] * ph_47[k]
                  + pi_61[k];

        t_48[k] = -ab_x[k] * ph_48[k]
                  + pi_62[k];

        t_49[k] = -ab_x[k] * ph_49[k]
                  + pi_63[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, ph_50, ph_51, ph_52, ph_53, \
                         ph_54, pi_64, pi_65, pi_66, pi_67, pi_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * ph_50[k]
                  + pi_64[k];

        t_51[k] = -ab_x[k] * ph_51[k]
                  + pi_65[k];

        t_52[k] = -ab_x[k] * ph_52[k]
                  + pi_66[k];

        t_53[k] = -ab_x[k] * ph_53[k]
                  + pi_67[k];

        t_54[k] = -ab_x[k] * ph_54[k]
                  + pi_68[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ph_55, ph_56, ph_57, ph_58, \
                         ph_59, pi_69, pi_70, pi_71, pi_72, pi_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * ph_55[k]
                  + pi_69[k];

        t_56[k] = -ab_x[k] * ph_56[k]
                  + pi_70[k];

        t_57[k] = -ab_x[k] * ph_57[k]
                  + pi_71[k];

        t_58[k] = -ab_x[k] * ph_58[k]
                  + pi_72[k];

        t_59[k] = -ab_x[k] * ph_59[k]
                  + pi_73[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, ab_x, ab_y, ph_21, ph_60, ph_61, ph_62, \
                         pi_29, pi_74, pi_75, pi_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * ph_60[k]
                  + pi_74[k];

        t_61[k] = -ab_x[k] * ph_61[k]
                  + pi_75[k];

        t_62[k] = -ab_x[k] * ph_62[k]
                  + pi_76[k];

        t_63[k] = -ab_y[k] * ph_21[k]
                  + pi_29[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, ab_y, ph_22, ph_23, ph_24, ph_25, \
                         ph_26, pi_31, pi_32, pi_34, pi_35, pi_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_64[k] = -ab_y[k] * ph_22[k]
                  + pi_31[k];

        t_65[k] = -ab_y[k] * ph_23[k]
                  + pi_32[k];

        t_66[k] = -ab_y[k] * ph_24[k]
                  + pi_34[k];

        t_67[k] = -ab_y[k] * ph_25[k]
                  + pi_35[k];

        t_68[k] = -ab_y[k] * ph_26[k]
                  + pi_36[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, ab_y, ph_27, ph_28, ph_29, ph_30, \
                         ph_31, pi_38, pi_39, pi_40, pi_41, pi_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_69[k] = -ab_y[k] * ph_27[k]
                  + pi_38[k];

        t_70[k] = -ab_y[k] * ph_28[k]
                  + pi_39[k];

        t_71[k] = -ab_y[k] * ph_29[k]
                  + pi_40[k];

        t_72[k] = -ab_y[k] * ph_30[k]
                  + pi_41[k];

        t_73[k] = -ab_y[k] * ph_31[k]
                  + pi_43[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, ab_y, ph_32, ph_33, ph_34, ph_35, \
                         ph_36, pi_44, pi_45, pi_46, pi_47, pi_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_74[k] = -ab_y[k] * ph_32[k]
                  + pi_44[k];

        t_75[k] = -ab_y[k] * ph_33[k]
                  + pi_45[k];

        t_76[k] = -ab_y[k] * ph_34[k]
                  + pi_46[k];

        t_77[k] = -ab_y[k] * ph_35[k]
                  + pi_47[k];

        t_78[k] = -ab_y[k] * ph_36[k]
                  + pi_49[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, ab_y, ph_37, ph_38, ph_39, ph_40, \
                         ph_41, pi_50, pi_51, pi_52, pi_53, pi_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_79[k] = -ab_y[k] * ph_37[k]
                  + pi_50[k];

        t_80[k] = -ab_y[k] * ph_38[k]
                  + pi_51[k];

        t_81[k] = -ab_y[k] * ph_39[k]
                  + pi_52[k];

        t_82[k] = -ab_y[k] * ph_40[k]
                  + pi_53[k];

        t_83[k] = -ab_y[k] * ph_41[k]
                  + pi_54[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ab_y, ph_42, ph_43, ph_44, ph_45, \
                         ph_46, pi_57, pi_59, pi_60, pi_62, pi_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_84[k] = -ab_y[k] * ph_42[k]
                  + pi_57[k];

        t_85[k] = -ab_y[k] * ph_43[k]
                  + pi_59[k];

        t_86[k] = -ab_y[k] * ph_44[k]
                  + pi_60[k];

        t_87[k] = -ab_y[k] * ph_45[k]
                  + pi_62[k];

        t_88[k] = -ab_y[k] * ph_46[k]
                  + pi_63[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ab_y, ph_47, ph_48, ph_49, ph_50, \
                         ph_51, pi_64, pi_66, pi_67, pi_68, pi_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_89[k] = -ab_y[k] * ph_47[k]
                  + pi_64[k];

        t_90[k] = -ab_y[k] * ph_48[k]
                  + pi_66[k];

        t_91[k] = -ab_y[k] * ph_49[k]
                  + pi_67[k];

        t_92[k] = -ab_y[k] * ph_50[k]
                  + pi_68[k];

        t_93[k] = -ab_y[k] * ph_51[k]
                  + pi_69[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ab_y, ph_52, ph_53, ph_54, ph_55, \
                         ph_56, pi_71, pi_72, pi_73, pi_74, pi_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_94[k] = -ab_y[k] * ph_52[k]
                  + pi_71[k];

        t_95[k] = -ab_y[k] * ph_53[k]
                  + pi_72[k];

        t_96[k] = -ab_y[k] * ph_54[k]
                  + pi_73[k];

        t_97[k] = -ab_y[k] * ph_55[k]
                  + pi_74[k];

        t_98[k] = -ab_y[k] * ph_56[k]
                  + pi_75[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ab_y, ph_57, ph_58, ph_59, ph_60, \
                         ph_61, pi_77, pi_78, pi_79, pi_80, pi_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_99[k] = -ab_y[k] * ph_57[k]
                  + pi_77[k];

        t_100[k] = -ab_y[k] * ph_58[k]
                   + pi_78[k];

        t_101[k] = -ab_y[k] * ph_59[k]
                   + pi_79[k];

        t_102[k] = -ab_y[k] * ph_60[k]
                   + pi_80[k];

        t_103[k] = -ab_y[k] * ph_61[k]
                   + pi_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, ab_y, ab_z, ph_42, ph_43, ph_44, ph_62, \
                         pi_58, pi_60, pi_61, pi_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_104[k] = -ab_y[k] * ph_62[k]
                   + pi_82[k];

        t_105[k] = -ab_z[k] * ph_42[k]
                   + pi_58[k];

        t_106[k] = -ab_z[k] * ph_43[k]
                   + pi_60[k];

        t_107[k] = -ab_z[k] * ph_44[k]
                   + pi_61[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, ab_z, ph_45, ph_46, ph_47, ph_48, \
                         ph_49, pi_63, pi_64, pi_65, pi_67, pi_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_108[k] = -ab_z[k] * ph_45[k]
                   + pi_63[k];

        t_109[k] = -ab_z[k] * ph_46[k]
                   + pi_64[k];

        t_110[k] = -ab_z[k] * ph_47[k]
                   + pi_65[k];

        t_111[k] = -ab_z[k] * ph_48[k]
                   + pi_67[k];

        t_112[k] = -ab_z[k] * ph_49[k]
                   + pi_68[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, ab_z, ph_50, ph_51, ph_52, ph_53, \
                         ph_54, pi_69, pi_70, pi_72, pi_73, pi_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_113[k] = -ab_z[k] * ph_50[k]
                   + pi_69[k];

        t_114[k] = -ab_z[k] * ph_51[k]
                   + pi_70[k];

        t_115[k] = -ab_z[k] * ph_52[k]
                   + pi_72[k];

        t_116[k] = -ab_z[k] * ph_53[k]
                   + pi_73[k];

        t_117[k] = -ab_z[k] * ph_54[k]
                   + pi_74[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, ab_z, ph_55, ph_56, ph_57, ph_58, \
                         ph_59, pi_75, pi_76, pi_78, pi_79, pi_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_118[k] = -ab_z[k] * ph_55[k]
                   + pi_75[k];

        t_119[k] = -ab_z[k] * ph_56[k]
                   + pi_76[k];

        t_120[k] = -ab_z[k] * ph_57[k]
                   + pi_78[k];

        t_121[k] = -ab_z[k] * ph_58[k]
                   + pi_79[k];

        t_122[k] = -ab_z[k] * ph_59[k]
                   + pi_80[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, ab_z, ph_60, ph_61, ph_62, pi_81, pi_82, \
                         pi_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_123[k] = -ab_z[k] * ph_60[k]
                   + pi_81[k];

        t_124[k] = -ab_z[k] * ph_61[k]
                   + pi_82[k];

        t_125[k] = -ab_z[k] * ph_62[k]
                   + pi_83[k];
    }
}

}  // namespace simdtrf
