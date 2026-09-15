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


#include "SimdElectronRepulsionGeom10VrrRecDH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_dh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ph, const size_t fh,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
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
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ph_0, ph_1, ph_2, ph_3, ph_4, fh_0, fh_1, \
                         fh_2, fh_3, fh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -2.0 * ph_0[k]
                 + f_0 * fh_0[k];

        t_1[k] = -2.0 * ph_1[k]
                 + f_0 * fh_1[k];

        t_2[k] = -2.0 * ph_2[k]
                 + f_0 * fh_2[k];

        t_3[k] = -2.0 * ph_3[k]
                 + f_0 * fh_3[k];

        t_4[k] = -2.0 * ph_4[k]
                 + f_0 * fh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ph_5, ph_6, ph_7, ph_8, ph_9, fh_5, fh_6, \
                         fh_7, fh_8, fh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -2.0 * ph_5[k]
                 + f_0 * fh_5[k];

        t_6[k] = -2.0 * ph_6[k]
                 + f_0 * fh_6[k];

        t_7[k] = -2.0 * ph_7[k]
                 + f_0 * fh_7[k];

        t_8[k] = -2.0 * ph_8[k]
                 + f_0 * fh_8[k];

        t_9[k] = -2.0 * ph_9[k]
                 + f_0 * fh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ph_10, ph_11, ph_12, ph_13, ph_14, \
                         fh_10, fh_11, fh_12, fh_13, fh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -2.0 * ph_10[k]
                  + f_0 * fh_10[k];

        t_11[k] = -2.0 * ph_11[k]
                  + f_0 * fh_11[k];

        t_12[k] = -2.0 * ph_12[k]
                  + f_0 * fh_12[k];

        t_13[k] = -2.0 * ph_13[k]
                  + f_0 * fh_13[k];

        t_14[k] = -2.0 * ph_14[k]
                  + f_0 * fh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ph_15, ph_16, ph_17, ph_18, ph_19, \
                         fh_15, fh_16, fh_17, fh_18, fh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -2.0 * ph_15[k]
                  + f_0 * fh_15[k];

        t_16[k] = -2.0 * ph_16[k]
                  + f_0 * fh_16[k];

        t_17[k] = -2.0 * ph_17[k]
                  + f_0 * fh_17[k];

        t_18[k] = -2.0 * ph_18[k]
                  + f_0 * fh_18[k];

        t_19[k] = -2.0 * ph_19[k]
                  + f_0 * fh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ph_20, ph_21, ph_22, ph_23, ph_24, \
                         fh_20, fh_21, fh_22, fh_23, fh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * ph_20[k]
                  + f_0 * fh_20[k];

        t_21[k] = -ph_21[k]
                  + f_0 * fh_21[k];

        t_22[k] = -ph_22[k]
                  + f_0 * fh_22[k];

        t_23[k] = -ph_23[k]
                  + f_0 * fh_23[k];

        t_24[k] = -ph_24[k]
                  + f_0 * fh_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ph_25, ph_26, ph_27, ph_28, ph_29, \
                         fh_25, fh_26, fh_27, fh_28, fh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -ph_25[k]
                  + f_0 * fh_25[k];

        t_26[k] = -ph_26[k]
                  + f_0 * fh_26[k];

        t_27[k] = -ph_27[k]
                  + f_0 * fh_27[k];

        t_28[k] = -ph_28[k]
                  + f_0 * fh_28[k];

        t_29[k] = -ph_29[k]
                  + f_0 * fh_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ph_30, ph_31, ph_32, ph_33, ph_34, \
                         fh_30, fh_31, fh_32, fh_33, fh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -ph_30[k]
                  + f_0 * fh_30[k];

        t_31[k] = -ph_31[k]
                  + f_0 * fh_31[k];

        t_32[k] = -ph_32[k]
                  + f_0 * fh_32[k];

        t_33[k] = -ph_33[k]
                  + f_0 * fh_33[k];

        t_34[k] = -ph_34[k]
                  + f_0 * fh_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ph_35, ph_36, ph_37, ph_38, ph_39, \
                         fh_35, fh_36, fh_37, fh_38, fh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -ph_35[k]
                  + f_0 * fh_35[k];

        t_36[k] = -ph_36[k]
                  + f_0 * fh_36[k];

        t_37[k] = -ph_37[k]
                  + f_0 * fh_37[k];

        t_38[k] = -ph_38[k]
                  + f_0 * fh_38[k];

        t_39[k] = -ph_39[k]
                  + f_0 * fh_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ph_40, ph_41, ph_42, ph_43, ph_44, \
                         fh_40, fh_41, fh_42, fh_43, fh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -ph_40[k]
                  + f_0 * fh_40[k];

        t_41[k] = -ph_41[k]
                  + f_0 * fh_41[k];

        t_42[k] = -ph_42[k]
                  + f_0 * fh_42[k];

        t_43[k] = -ph_43[k]
                  + f_0 * fh_43[k];

        t_44[k] = -ph_44[k]
                  + f_0 * fh_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ph_45, ph_46, ph_47, ph_48, ph_49, \
                         fh_45, fh_46, fh_47, fh_48, fh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -ph_45[k]
                  + f_0 * fh_45[k];

        t_46[k] = -ph_46[k]
                  + f_0 * fh_46[k];

        t_47[k] = -ph_47[k]
                  + f_0 * fh_47[k];

        t_48[k] = -ph_48[k]
                  + f_0 * fh_48[k];

        t_49[k] = -ph_49[k]
                  + f_0 * fh_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ph_50, ph_51, ph_52, ph_53, ph_54, \
                         fh_50, fh_51, fh_52, fh_53, fh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -ph_50[k]
                  + f_0 * fh_50[k];

        t_51[k] = -ph_51[k]
                  + f_0 * fh_51[k];

        t_52[k] = -ph_52[k]
                  + f_0 * fh_52[k];

        t_53[k] = -ph_53[k]
                  + f_0 * fh_53[k];

        t_54[k] = -ph_54[k]
                  + f_0 * fh_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ph_55, ph_56, ph_57, ph_58, ph_59, \
                         fh_55, fh_56, fh_57, fh_58, fh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -ph_55[k]
                  + f_0 * fh_55[k];

        t_56[k] = -ph_56[k]
                  + f_0 * fh_56[k];

        t_57[k] = -ph_57[k]
                  + f_0 * fh_57[k];

        t_58[k] = -ph_58[k]
                  + f_0 * fh_58[k];

        t_59[k] = -ph_59[k]
                  + f_0 * fh_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, ph_60, ph_61, ph_62, fh_60, \
                         fh_61, fh_62, fh_63, fh_64, fh_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -ph_60[k]
                  + f_0 * fh_60[k];

        t_61[k] = -ph_61[k]
                  + f_0 * fh_61[k];

        t_62[k] = -ph_62[k]
                  + f_0 * fh_62[k];

        t_63[k] = f_0 * fh_63[k];

        t_64[k] = f_0 * fh_64[k];

        t_65[k] = f_0 * fh_65[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, t_71, t_72, t_73, fh_66, fh_67, fh_68, \
                         fh_69, fh_70, fh_71, fh_72, fh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * fh_66[k];

        t_67[k] = f_0 * fh_67[k];

        t_68[k] = f_0 * fh_68[k];

        t_69[k] = f_0 * fh_69[k];

        t_70[k] = f_0 * fh_70[k];

        t_71[k] = f_0 * fh_71[k];

        t_72[k] = f_0 * fh_72[k];

        t_73[k] = f_0 * fh_73[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, t_79, t_80, t_81, fh_74, fh_75, fh_76, \
                         fh_77, fh_78, fh_79, fh_80, fh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * fh_74[k];

        t_75[k] = f_0 * fh_75[k];

        t_76[k] = f_0 * fh_76[k];

        t_77[k] = f_0 * fh_77[k];

        t_78[k] = f_0 * fh_78[k];

        t_79[k] = f_0 * fh_79[k];

        t_80[k] = f_0 * fh_80[k];

        t_81[k] = f_0 * fh_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, t_88, t_89, fh_82, fh_83, fh_84, \
                         fh_85, fh_86, fh_87, fh_88, fh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * fh_82[k];

        t_83[k] = f_0 * fh_83[k];

        t_84[k] = f_0 * fh_84[k];

        t_85[k] = f_0 * fh_85[k];

        t_86[k] = f_0 * fh_86[k];

        t_87[k] = f_0 * fh_87[k];

        t_88[k] = f_0 * fh_88[k];

        t_89[k] = f_0 * fh_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, fh_90, fh_91, fh_92, \
                         fh_93, fh_94, fh_95, fh_96, fh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * fh_90[k];

        t_91[k] = f_0 * fh_91[k];

        t_92[k] = f_0 * fh_92[k];

        t_93[k] = f_0 * fh_93[k];

        t_94[k] = f_0 * fh_94[k];

        t_95[k] = f_0 * fh_95[k];

        t_96[k] = f_0 * fh_96[k];

        t_97[k] = f_0 * fh_97[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, fh_98, fh_99, \
                         fh_100, fh_101, fh_102, fh_103, fh_104, \
                         fh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * fh_98[k];

        t_99[k] = f_0 * fh_99[k];

        t_100[k] = f_0 * fh_100[k];

        t_101[k] = f_0 * fh_101[k];

        t_102[k] = f_0 * fh_102[k];

        t_103[k] = f_0 * fh_103[k];

        t_104[k] = f_0 * fh_104[k];

        t_105[k] = f_0 * fh_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, t_112, t_113, fh_106, \
                         fh_107, fh_108, fh_109, fh_110, fh_111, fh_112, \
                         fh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * fh_106[k];

        t_107[k] = f_0 * fh_107[k];

        t_108[k] = f_0 * fh_108[k];

        t_109[k] = f_0 * fh_109[k];

        t_110[k] = f_0 * fh_110[k];

        t_111[k] = f_0 * fh_111[k];

        t_112[k] = f_0 * fh_112[k];

        t_113[k] = f_0 * fh_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, fh_114, \
                         fh_115, fh_116, fh_117, fh_118, fh_119, fh_120, \
                         fh_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * fh_114[k];

        t_115[k] = f_0 * fh_115[k];

        t_116[k] = f_0 * fh_116[k];

        t_117[k] = f_0 * fh_117[k];

        t_118[k] = f_0 * fh_118[k];

        t_119[k] = f_0 * fh_119[k];

        t_120[k] = f_0 * fh_120[k];

        t_121[k] = f_0 * fh_121[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, fh_122, fh_123, fh_124, \
                         fh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * fh_122[k];

        t_123[k] = f_0 * fh_123[k];

        t_124[k] = f_0 * fh_124[k];

        t_125[k] = f_0 * fh_125[k];
    }
}

auto
compute_prim_geom_10_dh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ph, const size_t fh,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_65 = buffer.data(fh + 65);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_68 = buffer.data(fh + 68);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_70 = buffer.data(fh + 70);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_74 = buffer.data(fh + 74);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_128 = buffer.data(fh + 128);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_136 = buffer.data(fh + 136);
    const auto *fh_137 = buffer.data(fh + 137);
    const auto *fh_138 = buffer.data(fh + 138);
    const auto *fh_139 = buffer.data(fh + 139);
    const auto *fh_140 = buffer.data(fh + 140);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_142 = buffer.data(fh + 142);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
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
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fh_21, fh_22, fh_23, fh_24, \
                         fh_25, fh_26, fh_27, fh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fh_21[k];

        t_1[k] = f_0 * fh_22[k];

        t_2[k] = f_0 * fh_23[k];

        t_3[k] = f_0 * fh_24[k];

        t_4[k] = f_0 * fh_25[k];

        t_5[k] = f_0 * fh_26[k];

        t_6[k] = f_0 * fh_27[k];

        t_7[k] = f_0 * fh_28[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fh_29, fh_30, fh_31, \
                         fh_32, fh_33, fh_34, fh_35, fh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fh_29[k];

        t_9[k] = f_0 * fh_30[k];

        t_10[k] = f_0 * fh_31[k];

        t_11[k] = f_0 * fh_32[k];

        t_12[k] = f_0 * fh_33[k];

        t_13[k] = f_0 * fh_34[k];

        t_14[k] = f_0 * fh_35[k];

        t_15[k] = f_0 * fh_36[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, ph_0, ph_1, fh_37, fh_38, \
                         fh_39, fh_40, fh_41, fh_63, fh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fh_37[k];

        t_17[k] = f_0 * fh_38[k];

        t_18[k] = f_0 * fh_39[k];

        t_19[k] = f_0 * fh_40[k];

        t_20[k] = f_0 * fh_41[k];

        t_21[k] = -ph_0[k]
                  + f_0 * fh_63[k];

        t_22[k] = -ph_1[k]
                  + f_0 * fh_64[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, ph_2, ph_3, ph_4, ph_5, ph_6, fh_65, \
                         fh_66, fh_67, fh_68, fh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -ph_2[k]
                  + f_0 * fh_65[k];

        t_24[k] = -ph_3[k]
                  + f_0 * fh_66[k];

        t_25[k] = -ph_4[k]
                  + f_0 * fh_67[k];

        t_26[k] = -ph_5[k]
                  + f_0 * fh_68[k];

        t_27[k] = -ph_6[k]
                  + f_0 * fh_69[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, ph_7, ph_8, ph_9, ph_10, ph_11, fh_70, \
                         fh_71, fh_72, fh_73, fh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -ph_7[k]
                  + f_0 * fh_70[k];

        t_29[k] = -ph_8[k]
                  + f_0 * fh_71[k];

        t_30[k] = -ph_9[k]
                  + f_0 * fh_72[k];

        t_31[k] = -ph_10[k]
                  + f_0 * fh_73[k];

        t_32[k] = -ph_11[k]
                  + f_0 * fh_74[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, ph_12, ph_13, ph_14, ph_15, ph_16, \
                         fh_75, fh_76, fh_77, fh_78, fh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -ph_12[k]
                  + f_0 * fh_75[k];

        t_34[k] = -ph_13[k]
                  + f_0 * fh_76[k];

        t_35[k] = -ph_14[k]
                  + f_0 * fh_77[k];

        t_36[k] = -ph_15[k]
                  + f_0 * fh_78[k];

        t_37[k] = -ph_16[k]
                  + f_0 * fh_79[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ph_17, ph_18, ph_19, ph_20, \
                         fh_80, fh_81, fh_82, fh_83, fh_84, fh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = -ph_17[k]
                  + f_0 * fh_80[k];

        t_39[k] = -ph_18[k]
                  + f_0 * fh_81[k];

        t_40[k] = -ph_19[k]
                  + f_0 * fh_82[k];

        t_41[k] = -ph_20[k]
                  + f_0 * fh_83[k];

        t_42[k] = f_0 * fh_84[k];

        t_43[k] = f_0 * fh_85[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, t_49, t_50, t_51, fh_86, fh_87, fh_88, \
                         fh_89, fh_90, fh_91, fh_92, fh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * fh_86[k];

        t_45[k] = f_0 * fh_87[k];

        t_46[k] = f_0 * fh_88[k];

        t_47[k] = f_0 * fh_89[k];

        t_48[k] = f_0 * fh_90[k];

        t_49[k] = f_0 * fh_91[k];

        t_50[k] = f_0 * fh_92[k];

        t_51[k] = f_0 * fh_93[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, t_57, t_58, t_59, fh_94, fh_95, fh_96, \
                         fh_97, fh_98, fh_99, fh_100, fh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * fh_94[k];

        t_53[k] = f_0 * fh_95[k];

        t_54[k] = f_0 * fh_96[k];

        t_55[k] = f_0 * fh_97[k];

        t_56[k] = f_0 * fh_98[k];

        t_57[k] = f_0 * fh_99[k];

        t_58[k] = f_0 * fh_100[k];

        t_59[k] = f_0 * fh_101[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, ph_21, ph_22, ph_23, fh_102, \
                         fh_103, fh_104, fh_126, fh_127, fh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * fh_102[k];

        t_61[k] = f_0 * fh_103[k];

        t_62[k] = f_0 * fh_104[k];

        t_63[k] = -2.0 * ph_21[k]
                  + f_0 * fh_126[k];

        t_64[k] = -2.0 * ph_22[k]
                  + f_0 * fh_127[k];

        t_65[k] = -2.0 * ph_23[k]
                  + f_0 * fh_128[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ph_24, ph_25, ph_26, ph_27, ph_28, \
                         fh_129, fh_130, fh_131, fh_132, fh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -2.0 * ph_24[k]
                  + f_0 * fh_129[k];

        t_67[k] = -2.0 * ph_25[k]
                  + f_0 * fh_130[k];

        t_68[k] = -2.0 * ph_26[k]
                  + f_0 * fh_131[k];

        t_69[k] = -2.0 * ph_27[k]
                  + f_0 * fh_132[k];

        t_70[k] = -2.0 * ph_28[k]
                  + f_0 * fh_133[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ph_29, ph_30, ph_31, ph_32, ph_33, \
                         fh_134, fh_135, fh_136, fh_137, fh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -2.0 * ph_29[k]
                  + f_0 * fh_134[k];

        t_72[k] = -2.0 * ph_30[k]
                  + f_0 * fh_135[k];

        t_73[k] = -2.0 * ph_31[k]
                  + f_0 * fh_136[k];

        t_74[k] = -2.0 * ph_32[k]
                  + f_0 * fh_137[k];

        t_75[k] = -2.0 * ph_33[k]
                  + f_0 * fh_138[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, ph_34, ph_35, ph_36, ph_37, ph_38, \
                         fh_139, fh_140, fh_141, fh_142, fh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * ph_34[k]
                  + f_0 * fh_139[k];

        t_77[k] = -2.0 * ph_35[k]
                  + f_0 * fh_140[k];

        t_78[k] = -2.0 * ph_36[k]
                  + f_0 * fh_141[k];

        t_79[k] = -2.0 * ph_37[k]
                  + f_0 * fh_142[k];

        t_80[k] = -2.0 * ph_38[k]
                  + f_0 * fh_143[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, ph_39, ph_40, ph_41, ph_42, ph_43, \
                         fh_144, fh_145, fh_146, fh_147, fh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -2.0 * ph_39[k]
                  + f_0 * fh_144[k];

        t_82[k] = -2.0 * ph_40[k]
                  + f_0 * fh_145[k];

        t_83[k] = -2.0 * ph_41[k]
                  + f_0 * fh_146[k];

        t_84[k] = -ph_42[k]
                  + f_0 * fh_147[k];

        t_85[k] = -ph_43[k]
                  + f_0 * fh_148[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, ph_44, ph_45, ph_46, ph_47, ph_48, \
                         fh_149, fh_150, fh_151, fh_152, fh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -ph_44[k]
                  + f_0 * fh_149[k];

        t_87[k] = -ph_45[k]
                  + f_0 * fh_150[k];

        t_88[k] = -ph_46[k]
                  + f_0 * fh_151[k];

        t_89[k] = -ph_47[k]
                  + f_0 * fh_152[k];

        t_90[k] = -ph_48[k]
                  + f_0 * fh_153[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, ph_49, ph_50, ph_51, ph_52, ph_53, \
                         fh_154, fh_155, fh_156, fh_157, fh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = -ph_49[k]
                  + f_0 * fh_154[k];

        t_92[k] = -ph_50[k]
                  + f_0 * fh_155[k];

        t_93[k] = -ph_51[k]
                  + f_0 * fh_156[k];

        t_94[k] = -ph_52[k]
                  + f_0 * fh_157[k];

        t_95[k] = -ph_53[k]
                  + f_0 * fh_158[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, ph_54, ph_55, ph_56, ph_57, ph_58, \
                         fh_159, fh_160, fh_161, fh_162, fh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -ph_54[k]
                  + f_0 * fh_159[k];

        t_97[k] = -ph_55[k]
                  + f_0 * fh_160[k];

        t_98[k] = -ph_56[k]
                  + f_0 * fh_161[k];

        t_99[k] = -ph_57[k]
                  + f_0 * fh_162[k];

        t_100[k] = -ph_58[k]
                   + f_0 * fh_163[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, ph_59, ph_60, ph_61, ph_62, \
                         fh_164, fh_165, fh_166, fh_167, fh_168, \
                         fh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -ph_59[k]
                   + f_0 * fh_164[k];

        t_102[k] = -ph_60[k]
                   + f_0 * fh_165[k];

        t_103[k] = -ph_61[k]
                   + f_0 * fh_166[k];

        t_104[k] = -ph_62[k]
                   + f_0 * fh_167[k];

        t_105[k] = f_0 * fh_168[k];

        t_106[k] = f_0 * fh_169[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, t_113, t_114, fh_170, \
                         fh_171, fh_172, fh_173, fh_174, fh_175, fh_176, \
                         fh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * fh_170[k];

        t_108[k] = f_0 * fh_171[k];

        t_109[k] = f_0 * fh_172[k];

        t_110[k] = f_0 * fh_173[k];

        t_111[k] = f_0 * fh_174[k];

        t_112[k] = f_0 * fh_175[k];

        t_113[k] = f_0 * fh_176[k];

        t_114[k] = f_0 * fh_177[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, t_121, t_122, fh_178, \
                         fh_179, fh_180, fh_181, fh_182, fh_183, fh_184, \
                         fh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * fh_178[k];

        t_116[k] = f_0 * fh_179[k];

        t_117[k] = f_0 * fh_180[k];

        t_118[k] = f_0 * fh_181[k];

        t_119[k] = f_0 * fh_182[k];

        t_120[k] = f_0 * fh_183[k];

        t_121[k] = f_0 * fh_184[k];

        t_122[k] = f_0 * fh_185[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, fh_186, fh_187, fh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_0 * fh_186[k];

        t_124[k] = f_0 * fh_187[k];

        t_125[k] = f_0 * fh_188[k];
    }
}

auto
compute_prim_geom_10_dh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ph, const size_t fh,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);
    const auto *fh_103 = buffer.data(fh + 103);
    const auto *fh_104 = buffer.data(fh + 104);
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
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);
    const auto *fh_121 = buffer.data(fh + 121);
    const auto *fh_122 = buffer.data(fh + 122);
    const auto *fh_123 = buffer.data(fh + 123);
    const auto *fh_124 = buffer.data(fh + 124);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_147 = buffer.data(fh + 147);
    const auto *fh_148 = buffer.data(fh + 148);
    const auto *fh_149 = buffer.data(fh + 149);
    const auto *fh_150 = buffer.data(fh + 150);
    const auto *fh_151 = buffer.data(fh + 151);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_153 = buffer.data(fh + 153);
    const auto *fh_154 = buffer.data(fh + 154);
    const auto *fh_155 = buffer.data(fh + 155);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_157 = buffer.data(fh + 157);
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
    const auto *fh_168 = buffer.data(fh + 168);
    const auto *fh_169 = buffer.data(fh + 169);
    const auto *fh_170 = buffer.data(fh + 170);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_172 = buffer.data(fh + 172);
    const auto *fh_173 = buffer.data(fh + 173);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_175 = buffer.data(fh + 175);
    const auto *fh_176 = buffer.data(fh + 176);
    const auto *fh_177 = buffer.data(fh + 177);
    const auto *fh_178 = buffer.data(fh + 178);
    const auto *fh_179 = buffer.data(fh + 179);
    const auto *fh_180 = buffer.data(fh + 180);
    const auto *fh_181 = buffer.data(fh + 181);
    const auto *fh_182 = buffer.data(fh + 182);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_199 = buffer.data(fh + 199);
    const auto *fh_200 = buffer.data(fh + 200);
    const auto *fh_201 = buffer.data(fh + 201);
    const auto *fh_202 = buffer.data(fh + 202);
    const auto *fh_203 = buffer.data(fh + 203);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_208 = buffer.data(fh + 208);
    const auto *fh_209 = buffer.data(fh + 209);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, fh_42, fh_43, fh_44, fh_45, \
                         fh_46, fh_47, fh_48, fh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fh_42[k];

        t_1[k] = f_0 * fh_43[k];

        t_2[k] = f_0 * fh_44[k];

        t_3[k] = f_0 * fh_45[k];

        t_4[k] = f_0 * fh_46[k];

        t_5[k] = f_0 * fh_47[k];

        t_6[k] = f_0 * fh_48[k];

        t_7[k] = f_0 * fh_49[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, t_14, t_15, fh_50, fh_51, fh_52, \
                         fh_53, fh_54, fh_55, fh_56, fh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * fh_50[k];

        t_9[k] = f_0 * fh_51[k];

        t_10[k] = f_0 * fh_52[k];

        t_11[k] = f_0 * fh_53[k];

        t_12[k] = f_0 * fh_54[k];

        t_13[k] = f_0 * fh_55[k];

        t_14[k] = f_0 * fh_56[k];

        t_15[k] = f_0 * fh_57[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, fh_58, fh_59, fh_60, \
                         fh_61, fh_62, fh_84, fh_85, fh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * fh_58[k];

        t_17[k] = f_0 * fh_59[k];

        t_18[k] = f_0 * fh_60[k];

        t_19[k] = f_0 * fh_61[k];

        t_20[k] = f_0 * fh_62[k];

        t_21[k] = f_0 * fh_84[k];

        t_22[k] = f_0 * fh_85[k];

        t_23[k] = f_0 * fh_86[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, t_29, t_30, t_31, fh_87, fh_88, fh_89, \
                         fh_90, fh_91, fh_92, fh_93, fh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * fh_87[k];

        t_25[k] = f_0 * fh_88[k];

        t_26[k] = f_0 * fh_89[k];

        t_27[k] = f_0 * fh_90[k];

        t_28[k] = f_0 * fh_91[k];

        t_29[k] = f_0 * fh_92[k];

        t_30[k] = f_0 * fh_93[k];

        t_31[k] = f_0 * fh_94[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, t_38, t_39, fh_95, fh_96, fh_97, \
                         fh_98, fh_99, fh_100, fh_101, fh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_0 * fh_95[k];

        t_33[k] = f_0 * fh_96[k];

        t_34[k] = f_0 * fh_97[k];

        t_35[k] = f_0 * fh_98[k];

        t_36[k] = f_0 * fh_99[k];

        t_37[k] = f_0 * fh_100[k];

        t_38[k] = f_0 * fh_101[k];

        t_39[k] = f_0 * fh_102[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, t_45, ph_0, ph_1, ph_2, ph_3, fh_103, \
                         fh_104, fh_105, fh_106, fh_107, fh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * fh_103[k];

        t_41[k] = f_0 * fh_104[k];

        t_42[k] = -ph_0[k]
                  + f_0 * fh_105[k];

        t_43[k] = -ph_1[k]
                  + f_0 * fh_106[k];

        t_44[k] = -ph_2[k]
                  + f_0 * fh_107[k];

        t_45[k] = -ph_3[k]
                  + f_0 * fh_108[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, ph_4, ph_5, ph_6, ph_7, ph_8, fh_109, \
                         fh_110, fh_111, fh_112, fh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -ph_4[k]
                  + f_0 * fh_109[k];

        t_47[k] = -ph_5[k]
                  + f_0 * fh_110[k];

        t_48[k] = -ph_6[k]
                  + f_0 * fh_111[k];

        t_49[k] = -ph_7[k]
                  + f_0 * fh_112[k];

        t_50[k] = -ph_8[k]
                  + f_0 * fh_113[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, ph_9, ph_10, ph_11, ph_12, ph_13, \
                         fh_114, fh_115, fh_116, fh_117, fh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -ph_9[k]
                  + f_0 * fh_114[k];

        t_52[k] = -ph_10[k]
                  + f_0 * fh_115[k];

        t_53[k] = -ph_11[k]
                  + f_0 * fh_116[k];

        t_54[k] = -ph_12[k]
                  + f_0 * fh_117[k];

        t_55[k] = -ph_13[k]
                  + f_0 * fh_118[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, ph_14, ph_15, ph_16, ph_17, ph_18, \
                         fh_119, fh_120, fh_121, fh_122, fh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -ph_14[k]
                  + f_0 * fh_119[k];

        t_57[k] = -ph_15[k]
                  + f_0 * fh_120[k];

        t_58[k] = -ph_16[k]
                  + f_0 * fh_121[k];

        t_59[k] = -ph_17[k]
                  + f_0 * fh_122[k];

        t_60[k] = -ph_18[k]
                  + f_0 * fh_123[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, ph_19, ph_20, fh_124, \
                         fh_125, fh_147, fh_148, fh_149, fh_150, \
                         fh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -ph_19[k]
                  + f_0 * fh_124[k];

        t_62[k] = -ph_20[k]
                  + f_0 * fh_125[k];

        t_63[k] = f_0 * fh_147[k];

        t_64[k] = f_0 * fh_148[k];

        t_65[k] = f_0 * fh_149[k];

        t_66[k] = f_0 * fh_150[k];

        t_67[k] = f_0 * fh_151[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, t_73, t_74, t_75, fh_152, fh_153, \
                         fh_154, fh_155, fh_156, fh_157, fh_158, \
                         fh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * fh_152[k];

        t_69[k] = f_0 * fh_153[k];

        t_70[k] = f_0 * fh_154[k];

        t_71[k] = f_0 * fh_155[k];

        t_72[k] = f_0 * fh_156[k];

        t_73[k] = f_0 * fh_157[k];

        t_74[k] = f_0 * fh_158[k];

        t_75[k] = f_0 * fh_159[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, t_83, fh_160, fh_161, \
                         fh_162, fh_163, fh_164, fh_165, fh_166, \
                         fh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_0 * fh_160[k];

        t_77[k] = f_0 * fh_161[k];

        t_78[k] = f_0 * fh_162[k];

        t_79[k] = f_0 * fh_163[k];

        t_80[k] = f_0 * fh_164[k];

        t_81[k] = f_0 * fh_165[k];

        t_82[k] = f_0 * fh_166[k];

        t_83[k] = f_0 * fh_167[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, ph_21, ph_22, ph_23, ph_24, ph_25, \
                         fh_168, fh_169, fh_170, fh_171, fh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -ph_21[k]
                  + f_0 * fh_168[k];

        t_85[k] = -ph_22[k]
                  + f_0 * fh_169[k];

        t_86[k] = -ph_23[k]
                  + f_0 * fh_170[k];

        t_87[k] = -ph_24[k]
                  + f_0 * fh_171[k];

        t_88[k] = -ph_25[k]
                  + f_0 * fh_172[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, ph_26, ph_27, ph_28, ph_29, ph_30, \
                         fh_173, fh_174, fh_175, fh_176, fh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -ph_26[k]
                  + f_0 * fh_173[k];

        t_90[k] = -ph_27[k]
                  + f_0 * fh_174[k];

        t_91[k] = -ph_28[k]
                  + f_0 * fh_175[k];

        t_92[k] = -ph_29[k]
                  + f_0 * fh_176[k];

        t_93[k] = -ph_30[k]
                  + f_0 * fh_177[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, ph_31, ph_32, ph_33, ph_34, ph_35, \
                         fh_178, fh_179, fh_180, fh_181, fh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = -ph_31[k]
                  + f_0 * fh_178[k];

        t_95[k] = -ph_32[k]
                  + f_0 * fh_179[k];

        t_96[k] = -ph_33[k]
                  + f_0 * fh_180[k];

        t_97[k] = -ph_34[k]
                  + f_0 * fh_181[k];

        t_98[k] = -ph_35[k]
                  + f_0 * fh_182[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, ph_36, ph_37, ph_38, ph_39, ph_40, \
                         fh_183, fh_184, fh_185, fh_186, fh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = -ph_36[k]
                  + f_0 * fh_183[k];

        t_100[k] = -ph_37[k]
                   + f_0 * fh_184[k];

        t_101[k] = -ph_38[k]
                   + f_0 * fh_185[k];

        t_102[k] = -ph_39[k]
                   + f_0 * fh_186[k];

        t_103[k] = -ph_40[k]
                   + f_0 * fh_187[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, ph_41, ph_42, ph_43, ph_44, ph_45, \
                         fh_188, fh_189, fh_190, fh_191, fh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = -ph_41[k]
                   + f_0 * fh_188[k];

        t_105[k] = -2.0 * ph_42[k]
                   + f_0 * fh_189[k];

        t_106[k] = -2.0 * ph_43[k]
                   + f_0 * fh_190[k];

        t_107[k] = -2.0 * ph_44[k]
                   + f_0 * fh_191[k];

        t_108[k] = -2.0 * ph_45[k]
                   + f_0 * fh_192[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, ph_46, ph_47, ph_48, ph_49, ph_50, \
                         fh_193, fh_194, fh_195, fh_196, fh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = -2.0 * ph_46[k]
                   + f_0 * fh_193[k];

        t_110[k] = -2.0 * ph_47[k]
                   + f_0 * fh_194[k];

        t_111[k] = -2.0 * ph_48[k]
                   + f_0 * fh_195[k];

        t_112[k] = -2.0 * ph_49[k]
                   + f_0 * fh_196[k];

        t_113[k] = -2.0 * ph_50[k]
                   + f_0 * fh_197[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, ph_51, ph_52, ph_53, ph_54, ph_55, \
                         fh_198, fh_199, fh_200, fh_201, fh_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = -2.0 * ph_51[k]
                   + f_0 * fh_198[k];

        t_115[k] = -2.0 * ph_52[k]
                   + f_0 * fh_199[k];

        t_116[k] = -2.0 * ph_53[k]
                   + f_0 * fh_200[k];

        t_117[k] = -2.0 * ph_54[k]
                   + f_0 * fh_201[k];

        t_118[k] = -2.0 * ph_55[k]
                   + f_0 * fh_202[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, ph_56, ph_57, ph_58, ph_59, ph_60, \
                         fh_203, fh_204, fh_205, fh_206, fh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -2.0 * ph_56[k]
                   + f_0 * fh_203[k];

        t_120[k] = -2.0 * ph_57[k]
                   + f_0 * fh_204[k];

        t_121[k] = -2.0 * ph_58[k]
                   + f_0 * fh_205[k];

        t_122[k] = -2.0 * ph_59[k]
                   + f_0 * fh_206[k];

        t_123[k] = -2.0 * ph_60[k]
                   + f_0 * fh_207[k];
    }

#pragma omp simd aligned(t_124, t_125, ph_61, ph_62, fh_208, fh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = -2.0 * ph_61[k]
                   + f_0 * fh_208[k];

        t_125[k] = -2.0 * ph_62[k]
                   + f_0 * fh_209[k];
    }
}

}  // namespace simdt2ceri
