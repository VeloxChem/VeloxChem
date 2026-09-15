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


#include "SimdElectronRepulsionGeom10VrrRecHD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_hd_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gd, const size_t id,
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

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_62 = buffer.data(gd + 62);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_85 = buffer.data(gd + 85);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gd_0, gd_1, gd_2, gd_3, gd_4, id_0, id_1, \
                         id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gd_0[k]
                 + f_0 * id_0[k];

        t_1[k] = -5.0 * gd_1[k]
                 + f_0 * id_1[k];

        t_2[k] = -5.0 * gd_2[k]
                 + f_0 * id_2[k];

        t_3[k] = -5.0 * gd_3[k]
                 + f_0 * id_3[k];

        t_4[k] = -5.0 * gd_4[k]
                 + f_0 * id_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gd_5, gd_6, gd_7, gd_8, gd_9, id_5, id_6, \
                         id_7, id_8, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -5.0 * gd_5[k]
                 + f_0 * id_5[k];

        t_6[k] = -4.0 * gd_6[k]
                 + f_0 * id_6[k];

        t_7[k] = -4.0 * gd_7[k]
                 + f_0 * id_7[k];

        t_8[k] = -4.0 * gd_8[k]
                 + f_0 * id_8[k];

        t_9[k] = -4.0 * gd_9[k]
                 + f_0 * id_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gd_10, gd_11, gd_12, gd_13, gd_14, \
                         id_10, id_11, id_12, id_13, id_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -4.0 * gd_10[k]
                  + f_0 * id_10[k];

        t_11[k] = -4.0 * gd_11[k]
                  + f_0 * id_11[k];

        t_12[k] = -4.0 * gd_12[k]
                  + f_0 * id_12[k];

        t_13[k] = -4.0 * gd_13[k]
                  + f_0 * id_13[k];

        t_14[k] = -4.0 * gd_14[k]
                  + f_0 * id_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gd_15, gd_16, gd_17, gd_18, gd_19, \
                         id_15, id_16, id_17, id_18, id_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -4.0 * gd_15[k]
                  + f_0 * id_15[k];

        t_16[k] = -4.0 * gd_16[k]
                  + f_0 * id_16[k];

        t_17[k] = -4.0 * gd_17[k]
                  + f_0 * id_17[k];

        t_18[k] = -3.0 * gd_18[k]
                  + f_0 * id_18[k];

        t_19[k] = -3.0 * gd_19[k]
                  + f_0 * id_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gd_20, gd_21, gd_22, gd_23, gd_24, \
                         id_20, id_21, id_22, id_23, id_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -3.0 * gd_20[k]
                  + f_0 * id_20[k];

        t_21[k] = -3.0 * gd_21[k]
                  + f_0 * id_21[k];

        t_22[k] = -3.0 * gd_22[k]
                  + f_0 * id_22[k];

        t_23[k] = -3.0 * gd_23[k]
                  + f_0 * id_23[k];

        t_24[k] = -3.0 * gd_24[k]
                  + f_0 * id_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gd_25, gd_26, gd_27, gd_28, gd_29, \
                         id_25, id_26, id_27, id_28, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -3.0 * gd_25[k]
                  + f_0 * id_25[k];

        t_26[k] = -3.0 * gd_26[k]
                  + f_0 * id_26[k];

        t_27[k] = -3.0 * gd_27[k]
                  + f_0 * id_27[k];

        t_28[k] = -3.0 * gd_28[k]
                  + f_0 * id_28[k];

        t_29[k] = -3.0 * gd_29[k]
                  + f_0 * id_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gd_30, gd_31, gd_32, gd_33, gd_34, \
                         id_30, id_31, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * gd_30[k]
                  + f_0 * id_30[k];

        t_31[k] = -3.0 * gd_31[k]
                  + f_0 * id_31[k];

        t_32[k] = -3.0 * gd_32[k]
                  + f_0 * id_32[k];

        t_33[k] = -3.0 * gd_33[k]
                  + f_0 * id_33[k];

        t_34[k] = -3.0 * gd_34[k]
                  + f_0 * id_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gd_35, gd_36, gd_37, gd_38, gd_39, \
                         id_35, id_36, id_37, id_38, id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * gd_35[k]
                  + f_0 * id_35[k];

        t_36[k] = -2.0 * gd_36[k]
                  + f_0 * id_36[k];

        t_37[k] = -2.0 * gd_37[k]
                  + f_0 * id_37[k];

        t_38[k] = -2.0 * gd_38[k]
                  + f_0 * id_38[k];

        t_39[k] = -2.0 * gd_39[k]
                  + f_0 * id_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gd_40, gd_41, gd_42, gd_43, gd_44, \
                         id_40, id_41, id_42, id_43, id_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -2.0 * gd_40[k]
                  + f_0 * id_40[k];

        t_41[k] = -2.0 * gd_41[k]
                  + f_0 * id_41[k];

        t_42[k] = -2.0 * gd_42[k]
                  + f_0 * id_42[k];

        t_43[k] = -2.0 * gd_43[k]
                  + f_0 * id_43[k];

        t_44[k] = -2.0 * gd_44[k]
                  + f_0 * id_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gd_45, gd_46, gd_47, gd_48, gd_49, \
                         id_45, id_46, id_47, id_48, id_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * gd_45[k]
                  + f_0 * id_45[k];

        t_46[k] = -2.0 * gd_46[k]
                  + f_0 * id_46[k];

        t_47[k] = -2.0 * gd_47[k]
                  + f_0 * id_47[k];

        t_48[k] = -2.0 * gd_48[k]
                  + f_0 * id_48[k];

        t_49[k] = -2.0 * gd_49[k]
                  + f_0 * id_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gd_50, gd_51, gd_52, gd_53, gd_54, \
                         id_50, id_51, id_52, id_53, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * gd_50[k]
                  + f_0 * id_50[k];

        t_51[k] = -2.0 * gd_51[k]
                  + f_0 * id_51[k];

        t_52[k] = -2.0 * gd_52[k]
                  + f_0 * id_52[k];

        t_53[k] = -2.0 * gd_53[k]
                  + f_0 * id_53[k];

        t_54[k] = -2.0 * gd_54[k]
                  + f_0 * id_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gd_55, gd_56, gd_57, gd_58, gd_59, \
                         id_55, id_56, id_57, id_58, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * gd_55[k]
                  + f_0 * id_55[k];

        t_56[k] = -2.0 * gd_56[k]
                  + f_0 * id_56[k];

        t_57[k] = -2.0 * gd_57[k]
                  + f_0 * id_57[k];

        t_58[k] = -2.0 * gd_58[k]
                  + f_0 * id_58[k];

        t_59[k] = -2.0 * gd_59[k]
                  + f_0 * id_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, gd_60, gd_61, gd_62, gd_63, gd_64, \
                         id_60, id_61, id_62, id_63, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -gd_60[k]
                  + f_0 * id_60[k];

        t_61[k] = -gd_61[k]
                  + f_0 * id_61[k];

        t_62[k] = -gd_62[k]
                  + f_0 * id_62[k];

        t_63[k] = -gd_63[k]
                  + f_0 * id_63[k];

        t_64[k] = -gd_64[k]
                  + f_0 * id_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gd_65, gd_66, gd_67, gd_68, gd_69, \
                         id_65, id_66, id_67, id_68, id_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -gd_65[k]
                  + f_0 * id_65[k];

        t_66[k] = -gd_66[k]
                  + f_0 * id_66[k];

        t_67[k] = -gd_67[k]
                  + f_0 * id_67[k];

        t_68[k] = -gd_68[k]
                  + f_0 * id_68[k];

        t_69[k] = -gd_69[k]
                  + f_0 * id_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gd_70, gd_71, gd_72, gd_73, gd_74, \
                         id_70, id_71, id_72, id_73, id_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -gd_70[k]
                  + f_0 * id_70[k];

        t_71[k] = -gd_71[k]
                  + f_0 * id_71[k];

        t_72[k] = -gd_72[k]
                  + f_0 * id_72[k];

        t_73[k] = -gd_73[k]
                  + f_0 * id_73[k];

        t_74[k] = -gd_74[k]
                  + f_0 * id_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gd_75, gd_76, gd_77, gd_78, gd_79, \
                         id_75, id_76, id_77, id_78, id_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -gd_75[k]
                  + f_0 * id_75[k];

        t_76[k] = -gd_76[k]
                  + f_0 * id_76[k];

        t_77[k] = -gd_77[k]
                  + f_0 * id_77[k];

        t_78[k] = -gd_78[k]
                  + f_0 * id_78[k];

        t_79[k] = -gd_79[k]
                  + f_0 * id_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gd_80, gd_81, gd_82, gd_83, gd_84, \
                         id_80, id_81, id_82, id_83, id_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -gd_80[k]
                  + f_0 * id_80[k];

        t_81[k] = -gd_81[k]
                  + f_0 * id_81[k];

        t_82[k] = -gd_82[k]
                  + f_0 * id_82[k];

        t_83[k] = -gd_83[k]
                  + f_0 * id_83[k];

        t_84[k] = -gd_84[k]
                  + f_0 * id_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gd_85, gd_86, gd_87, gd_88, gd_89, \
                         id_85, id_86, id_87, id_88, id_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -gd_85[k]
                  + f_0 * id_85[k];

        t_86[k] = -gd_86[k]
                  + f_0 * id_86[k];

        t_87[k] = -gd_87[k]
                  + f_0 * id_87[k];

        t_88[k] = -gd_88[k]
                  + f_0 * id_88[k];

        t_89[k] = -gd_89[k]
                  + f_0 * id_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, t_97, id_90, id_91, id_92, \
                         id_93, id_94, id_95, id_96, id_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * id_90[k];

        t_91[k] = f_0 * id_91[k];

        t_92[k] = f_0 * id_92[k];

        t_93[k] = f_0 * id_93[k];

        t_94[k] = f_0 * id_94[k];

        t_95[k] = f_0 * id_95[k];

        t_96[k] = f_0 * id_96[k];

        t_97[k] = f_0 * id_97[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, t_103, t_104, t_105, id_98, id_99, \
                         id_100, id_101, id_102, id_103, id_104, \
                         id_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * id_98[k];

        t_99[k] = f_0 * id_99[k];

        t_100[k] = f_0 * id_100[k];

        t_101[k] = f_0 * id_101[k];

        t_102[k] = f_0 * id_102[k];

        t_103[k] = f_0 * id_103[k];

        t_104[k] = f_0 * id_104[k];

        t_105[k] = f_0 * id_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, t_111, t_112, t_113, id_106, \
                         id_107, id_108, id_109, id_110, id_111, id_112, \
                         id_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_0 * id_106[k];

        t_107[k] = f_0 * id_107[k];

        t_108[k] = f_0 * id_108[k];

        t_109[k] = f_0 * id_109[k];

        t_110[k] = f_0 * id_110[k];

        t_111[k] = f_0 * id_111[k];

        t_112[k] = f_0 * id_112[k];

        t_113[k] = f_0 * id_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, t_119, t_120, t_121, id_114, \
                         id_115, id_116, id_117, id_118, id_119, id_120, \
                         id_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_0 * id_114[k];

        t_115[k] = f_0 * id_115[k];

        t_116[k] = f_0 * id_116[k];

        t_117[k] = f_0 * id_117[k];

        t_118[k] = f_0 * id_118[k];

        t_119[k] = f_0 * id_119[k];

        t_120[k] = f_0 * id_120[k];

        t_121[k] = f_0 * id_121[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, id_122, id_123, id_124, \
                         id_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * id_122[k];

        t_123[k] = f_0 * id_123[k];

        t_124[k] = f_0 * id_124[k];

        t_125[k] = f_0 * id_125[k];
    }
}

auto
compute_prim_geom_10_hd_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gd, const size_t id,
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

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_62 = buffer.data(gd + 62);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_85 = buffer.data(gd + 85);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_128 = buffer.data(id + 128);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
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
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, gd_0, id_6, id_7, id_8, id_9, \
                         id_10, id_11, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_6[k];

        t_1[k] = f_0 * id_7[k];

        t_2[k] = f_0 * id_8[k];

        t_3[k] = f_0 * id_9[k];

        t_4[k] = f_0 * id_10[k];

        t_5[k] = f_0 * id_11[k];

        t_6[k] = -gd_0[k]
                 + f_0 * id_18[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, gd_1, gd_2, gd_3, gd_4, gd_5, id_19, \
                         id_20, id_21, id_22, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -gd_1[k]
                 + f_0 * id_19[k];

        t_8[k] = -gd_2[k]
                 + f_0 * id_20[k];

        t_9[k] = -gd_3[k]
                 + f_0 * id_21[k];

        t_10[k] = -gd_4[k]
                  + f_0 * id_22[k];

        t_11[k] = -gd_5[k]
                  + f_0 * id_23[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, t_18, gd_6, id_24, id_25, id_26, \
                         id_27, id_28, id_29, id_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * id_24[k];

        t_13[k] = f_0 * id_25[k];

        t_14[k] = f_0 * id_26[k];

        t_15[k] = f_0 * id_27[k];

        t_16[k] = f_0 * id_28[k];

        t_17[k] = f_0 * id_29[k];

        t_18[k] = -2.0 * gd_6[k]
                  + f_0 * id_36[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, gd_7, gd_8, gd_9, gd_10, gd_11, id_37, \
                         id_38, id_39, id_40, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -2.0 * gd_7[k]
                  + f_0 * id_37[k];

        t_20[k] = -2.0 * gd_8[k]
                  + f_0 * id_38[k];

        t_21[k] = -2.0 * gd_9[k]
                  + f_0 * id_39[k];

        t_22[k] = -2.0 * gd_10[k]
                  + f_0 * id_40[k];

        t_23[k] = -2.0 * gd_11[k]
                  + f_0 * id_41[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, gd_12, gd_13, gd_14, gd_15, gd_16, \
                         id_42, id_43, id_44, id_45, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -gd_12[k]
                  + f_0 * id_42[k];

        t_25[k] = -gd_13[k]
                  + f_0 * id_43[k];

        t_26[k] = -gd_14[k]
                  + f_0 * id_44[k];

        t_27[k] = -gd_15[k]
                  + f_0 * id_45[k];

        t_28[k] = -gd_16[k]
                  + f_0 * id_46[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, t_35, gd_17, id_47, id_48, id_49, \
                         id_50, id_51, id_52, id_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -gd_17[k]
                  + f_0 * id_47[k];

        t_30[k] = f_0 * id_48[k];

        t_31[k] = f_0 * id_49[k];

        t_32[k] = f_0 * id_50[k];

        t_33[k] = f_0 * id_51[k];

        t_34[k] = f_0 * id_52[k];

        t_35[k] = f_0 * id_53[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, gd_18, gd_19, gd_20, gd_21, gd_22, \
                         id_60, id_61, id_62, id_63, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -3.0 * gd_18[k]
                  + f_0 * id_60[k];

        t_37[k] = -3.0 * gd_19[k]
                  + f_0 * id_61[k];

        t_38[k] = -3.0 * gd_20[k]
                  + f_0 * id_62[k];

        t_39[k] = -3.0 * gd_21[k]
                  + f_0 * id_63[k];

        t_40[k] = -3.0 * gd_22[k]
                  + f_0 * id_64[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, gd_23, gd_24, gd_25, gd_26, gd_27, \
                         id_65, id_66, id_67, id_68, id_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -3.0 * gd_23[k]
                  + f_0 * id_65[k];

        t_42[k] = -2.0 * gd_24[k]
                  + f_0 * id_66[k];

        t_43[k] = -2.0 * gd_25[k]
                  + f_0 * id_67[k];

        t_44[k] = -2.0 * gd_26[k]
                  + f_0 * id_68[k];

        t_45[k] = -2.0 * gd_27[k]
                  + f_0 * id_69[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, gd_28, gd_29, gd_30, gd_31, gd_32, \
                         id_70, id_71, id_72, id_73, id_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -2.0 * gd_28[k]
                  + f_0 * id_70[k];

        t_47[k] = -2.0 * gd_29[k]
                  + f_0 * id_71[k];

        t_48[k] = -gd_30[k]
                  + f_0 * id_72[k];

        t_49[k] = -gd_31[k]
                  + f_0 * id_73[k];

        t_50[k] = -gd_32[k]
                  + f_0 * id_74[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, gd_33, gd_34, gd_35, id_75, \
                         id_76, id_77, id_78, id_79, id_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -gd_33[k]
                  + f_0 * id_75[k];

        t_52[k] = -gd_34[k]
                  + f_0 * id_76[k];

        t_53[k] = -gd_35[k]
                  + f_0 * id_77[k];

        t_54[k] = f_0 * id_78[k];

        t_55[k] = f_0 * id_79[k];

        t_56[k] = f_0 * id_80[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, t_62, gd_36, gd_37, gd_38, id_81, \
                         id_82, id_83, id_90, id_91, id_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_0 * id_81[k];

        t_58[k] = f_0 * id_82[k];

        t_59[k] = f_0 * id_83[k];

        t_60[k] = -4.0 * gd_36[k]
                  + f_0 * id_90[k];

        t_61[k] = -4.0 * gd_37[k]
                  + f_0 * id_91[k];

        t_62[k] = -4.0 * gd_38[k]
                  + f_0 * id_92[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, gd_39, gd_40, gd_41, gd_42, gd_43, \
                         id_93, id_94, id_95, id_96, id_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -4.0 * gd_39[k]
                  + f_0 * id_93[k];

        t_64[k] = -4.0 * gd_40[k]
                  + f_0 * id_94[k];

        t_65[k] = -4.0 * gd_41[k]
                  + f_0 * id_95[k];

        t_66[k] = -3.0 * gd_42[k]
                  + f_0 * id_96[k];

        t_67[k] = -3.0 * gd_43[k]
                  + f_0 * id_97[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, gd_44, gd_45, gd_46, gd_47, gd_48, \
                         id_98, id_99, id_100, id_101, id_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -3.0 * gd_44[k]
                  + f_0 * id_98[k];

        t_69[k] = -3.0 * gd_45[k]
                  + f_0 * id_99[k];

        t_70[k] = -3.0 * gd_46[k]
                  + f_0 * id_100[k];

        t_71[k] = -3.0 * gd_47[k]
                  + f_0 * id_101[k];

        t_72[k] = -2.0 * gd_48[k]
                  + f_0 * id_102[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, gd_49, gd_50, gd_51, gd_52, gd_53, \
                         id_103, id_104, id_105, id_106, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = -2.0 * gd_49[k]
                  + f_0 * id_103[k];

        t_74[k] = -2.0 * gd_50[k]
                  + f_0 * id_104[k];

        t_75[k] = -2.0 * gd_51[k]
                  + f_0 * id_105[k];

        t_76[k] = -2.0 * gd_52[k]
                  + f_0 * id_106[k];

        t_77[k] = -2.0 * gd_53[k]
                  + f_0 * id_107[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, gd_54, gd_55, gd_56, gd_57, gd_58, \
                         id_108, id_109, id_110, id_111, id_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -gd_54[k]
                  + f_0 * id_108[k];

        t_79[k] = -gd_55[k]
                  + f_0 * id_109[k];

        t_80[k] = -gd_56[k]
                  + f_0 * id_110[k];

        t_81[k] = -gd_57[k]
                  + f_0 * id_111[k];

        t_82[k] = -gd_58[k]
                  + f_0 * id_112[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, t_88, t_89, gd_59, id_113, id_114, \
                         id_115, id_116, id_117, id_118, id_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = -gd_59[k]
                  + f_0 * id_113[k];

        t_84[k] = f_0 * id_114[k];

        t_85[k] = f_0 * id_115[k];

        t_86[k] = f_0 * id_116[k];

        t_87[k] = f_0 * id_117[k];

        t_88[k] = f_0 * id_118[k];

        t_89[k] = f_0 * id_119[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, gd_60, gd_61, gd_62, gd_63, gd_64, \
                         id_126, id_127, id_128, id_129, id_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = -5.0 * gd_60[k]
                  + f_0 * id_126[k];

        t_91[k] = -5.0 * gd_61[k]
                  + f_0 * id_127[k];

        t_92[k] = -5.0 * gd_62[k]
                  + f_0 * id_128[k];

        t_93[k] = -5.0 * gd_63[k]
                  + f_0 * id_129[k];

        t_94[k] = -5.0 * gd_64[k]
                  + f_0 * id_130[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, gd_65, gd_66, gd_67, gd_68, gd_69, \
                         id_131, id_132, id_133, id_134, id_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = -5.0 * gd_65[k]
                  + f_0 * id_131[k];

        t_96[k] = -4.0 * gd_66[k]
                  + f_0 * id_132[k];

        t_97[k] = -4.0 * gd_67[k]
                  + f_0 * id_133[k];

        t_98[k] = -4.0 * gd_68[k]
                  + f_0 * id_134[k];

        t_99[k] = -4.0 * gd_69[k]
                  + f_0 * id_135[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, gd_70, gd_71, gd_72, gd_73, gd_74, \
                         id_136, id_137, id_138, id_139, id_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -4.0 * gd_70[k]
                   + f_0 * id_136[k];

        t_101[k] = -4.0 * gd_71[k]
                   + f_0 * id_137[k];

        t_102[k] = -3.0 * gd_72[k]
                   + f_0 * id_138[k];

        t_103[k] = -3.0 * gd_73[k]
                   + f_0 * id_139[k];

        t_104[k] = -3.0 * gd_74[k]
                   + f_0 * id_140[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, gd_75, gd_76, gd_77, gd_78, gd_79, \
                         id_141, id_142, id_143, id_144, id_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -3.0 * gd_75[k]
                   + f_0 * id_141[k];

        t_106[k] = -3.0 * gd_76[k]
                   + f_0 * id_142[k];

        t_107[k] = -3.0 * gd_77[k]
                   + f_0 * id_143[k];

        t_108[k] = -2.0 * gd_78[k]
                   + f_0 * id_144[k];

        t_109[k] = -2.0 * gd_79[k]
                   + f_0 * id_145[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, gd_80, gd_81, gd_82, gd_83, gd_84, \
                         id_146, id_147, id_148, id_149, id_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -2.0 * gd_80[k]
                   + f_0 * id_146[k];

        t_111[k] = -2.0 * gd_81[k]
                   + f_0 * id_147[k];

        t_112[k] = -2.0 * gd_82[k]
                   + f_0 * id_148[k];

        t_113[k] = -2.0 * gd_83[k]
                   + f_0 * id_149[k];

        t_114[k] = -gd_84[k]
                   + f_0 * id_150[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, gd_85, gd_86, gd_87, gd_88, gd_89, \
                         id_151, id_152, id_153, id_154, id_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = -gd_85[k]
                   + f_0 * id_151[k];

        t_116[k] = -gd_86[k]
                   + f_0 * id_152[k];

        t_117[k] = -gd_87[k]
                   + f_0 * id_153[k];

        t_118[k] = -gd_88[k]
                   + f_0 * id_154[k];

        t_119[k] = -gd_89[k]
                   + f_0 * id_155[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, id_156, id_157, id_158, \
                         id_159, id_160, id_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_0 * id_156[k];

        t_121[k] = f_0 * id_157[k];

        t_122[k] = f_0 * id_158[k];

        t_123[k] = f_0 * id_159[k];

        t_124[k] = f_0 * id_160[k];

        t_125[k] = f_0 * id_161[k];
    }
}

auto
compute_prim_geom_10_hd_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gd, const size_t id,
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

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_1 = buffer.data(gd + 1);
    const auto *gd_2 = buffer.data(gd + 2);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_4 = buffer.data(gd + 4);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_7 = buffer.data(gd + 7);
    const auto *gd_8 = buffer.data(gd + 8);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_10 = buffer.data(gd + 10);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_13 = buffer.data(gd + 13);
    const auto *gd_14 = buffer.data(gd + 14);
    const auto *gd_15 = buffer.data(gd + 15);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_19 = buffer.data(gd + 19);
    const auto *gd_20 = buffer.data(gd + 20);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_22 = buffer.data(gd + 22);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_24 = buffer.data(gd + 24);
    const auto *gd_25 = buffer.data(gd + 25);
    const auto *gd_26 = buffer.data(gd + 26);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_31 = buffer.data(gd + 31);
    const auto *gd_32 = buffer.data(gd + 32);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_37 = buffer.data(gd + 37);
    const auto *gd_38 = buffer.data(gd + 38);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_40 = buffer.data(gd + 40);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_43 = buffer.data(gd + 43);
    const auto *gd_44 = buffer.data(gd + 44);
    const auto *gd_45 = buffer.data(gd + 45);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_49 = buffer.data(gd + 49);
    const auto *gd_50 = buffer.data(gd + 50);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_53 = buffer.data(gd + 53);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_55 = buffer.data(gd + 55);
    const auto *gd_56 = buffer.data(gd + 56);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_58 = buffer.data(gd + 58);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_62 = buffer.data(gd + 62);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_64 = buffer.data(gd + 64);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_66 = buffer.data(gd + 66);
    const auto *gd_67 = buffer.data(gd + 67);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_78 = buffer.data(gd + 78);
    const auto *gd_79 = buffer.data(gd + 79);
    const auto *gd_80 = buffer.data(gd + 80);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_85 = buffer.data(gd + 85);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_97 = buffer.data(id + 97);
    const auto *id_98 = buffer.data(id + 98);
    const auto *id_99 = buffer.data(id + 99);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_103 = buffer.data(id + 103);
    const auto *id_104 = buffer.data(id + 104);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_109 = buffer.data(id + 109);
    const auto *id_110 = buffer.data(id + 110);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_115 = buffer.data(id + 115);
    const auto *id_116 = buffer.data(id + 116);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_119 = buffer.data(id + 119);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_121 = buffer.data(id + 121);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_124 = buffer.data(id + 124);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_132 = buffer.data(id + 132);
    const auto *id_133 = buffer.data(id + 133);
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
    const auto *id_156 = buffer.data(id + 156);
    const auto *id_157 = buffer.data(id + 157);
    const auto *id_158 = buffer.data(id + 158);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_163 = buffer.data(id + 163);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, id_12, id_13, id_14, id_15, \
                         id_16, id_17, id_24, id_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_12[k];

        t_1[k] = f_0 * id_13[k];

        t_2[k] = f_0 * id_14[k];

        t_3[k] = f_0 * id_15[k];

        t_4[k] = f_0 * id_16[k];

        t_5[k] = f_0 * id_17[k];

        t_6[k] = f_0 * id_24[k];

        t_7[k] = f_0 * id_25[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, t_12, t_13, gd_0, gd_1, id_26, id_27, id_28, \
                         id_29, id_30, id_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * id_26[k];

        t_9[k] = f_0 * id_27[k];

        t_10[k] = f_0 * id_28[k];

        t_11[k] = f_0 * id_29[k];

        t_12[k] = -gd_0[k]
                  + f_0 * id_30[k];

        t_13[k] = -gd_1[k]
                  + f_0 * id_31[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, gd_2, gd_3, gd_4, gd_5, id_32, \
                         id_33, id_34, id_35, id_42, id_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -gd_2[k]
                  + f_0 * id_32[k];

        t_15[k] = -gd_3[k]
                  + f_0 * id_33[k];

        t_16[k] = -gd_4[k]
                  + f_0 * id_34[k];

        t_17[k] = -gd_5[k]
                  + f_0 * id_35[k];

        t_18[k] = f_0 * id_42[k];

        t_19[k] = f_0 * id_43[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, gd_6, gd_7, id_44, id_45, id_46, \
                         id_47, id_48, id_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * id_44[k];

        t_21[k] = f_0 * id_45[k];

        t_22[k] = f_0 * id_46[k];

        t_23[k] = f_0 * id_47[k];

        t_24[k] = -gd_6[k]
                  + f_0 * id_48[k];

        t_25[k] = -gd_7[k]
                  + f_0 * id_49[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, gd_8, gd_9, gd_10, gd_11, gd_12, id_50, \
                         id_51, id_52, id_53, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -gd_8[k]
                  + f_0 * id_50[k];

        t_27[k] = -gd_9[k]
                  + f_0 * id_51[k];

        t_28[k] = -gd_10[k]
                  + f_0 * id_52[k];

        t_29[k] = -gd_11[k]
                  + f_0 * id_53[k];

        t_30[k] = -2.0 * gd_12[k]
                  + f_0 * id_54[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, gd_13, gd_14, gd_15, gd_16, gd_17, \
                         id_55, id_56, id_57, id_58, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -2.0 * gd_13[k]
                  + f_0 * id_55[k];

        t_32[k] = -2.0 * gd_14[k]
                  + f_0 * id_56[k];

        t_33[k] = -2.0 * gd_15[k]
                  + f_0 * id_57[k];

        t_34[k] = -2.0 * gd_16[k]
                  + f_0 * id_58[k];

        t_35[k] = -2.0 * gd_17[k]
                  + f_0 * id_59[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, t_42, gd_18, id_66, id_67, id_68, \
                         id_69, id_70, id_71, id_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * id_66[k];

        t_37[k] = f_0 * id_67[k];

        t_38[k] = f_0 * id_68[k];

        t_39[k] = f_0 * id_69[k];

        t_40[k] = f_0 * id_70[k];

        t_41[k] = f_0 * id_71[k];

        t_42[k] = -gd_18[k]
                  + f_0 * id_72[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, gd_19, gd_20, gd_21, gd_22, gd_23, \
                         id_73, id_74, id_75, id_76, id_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -gd_19[k]
                  + f_0 * id_73[k];

        t_44[k] = -gd_20[k]
                  + f_0 * id_74[k];

        t_45[k] = -gd_21[k]
                  + f_0 * id_75[k];

        t_46[k] = -gd_22[k]
                  + f_0 * id_76[k];

        t_47[k] = -gd_23[k]
                  + f_0 * id_77[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, gd_24, gd_25, gd_26, gd_27, gd_28, \
                         id_78, id_79, id_80, id_81, id_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -2.0 * gd_24[k]
                  + f_0 * id_78[k];

        t_49[k] = -2.0 * gd_25[k]
                  + f_0 * id_79[k];

        t_50[k] = -2.0 * gd_26[k]
                  + f_0 * id_80[k];

        t_51[k] = -2.0 * gd_27[k]
                  + f_0 * id_81[k];

        t_52[k] = -2.0 * gd_28[k]
                  + f_0 * id_82[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, gd_29, gd_30, gd_31, gd_32, gd_33, \
                         id_83, id_84, id_85, id_86, id_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -2.0 * gd_29[k]
                  + f_0 * id_83[k];

        t_54[k] = -3.0 * gd_30[k]
                  + f_0 * id_84[k];

        t_55[k] = -3.0 * gd_31[k]
                  + f_0 * id_85[k];

        t_56[k] = -3.0 * gd_32[k]
                  + f_0 * id_86[k];

        t_57[k] = -3.0 * gd_33[k]
                  + f_0 * id_87[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, t_63, t_64, gd_34, gd_35, id_88, id_89, \
                         id_96, id_97, id_98, id_99, id_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -3.0 * gd_34[k]
                  + f_0 * id_88[k];

        t_59[k] = -3.0 * gd_35[k]
                  + f_0 * id_89[k];

        t_60[k] = f_0 * id_96[k];

        t_61[k] = f_0 * id_97[k];

        t_62[k] = f_0 * id_98[k];

        t_63[k] = f_0 * id_99[k];

        t_64[k] = f_0 * id_100[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, gd_36, gd_37, gd_38, gd_39, id_101, \
                         id_102, id_103, id_104, id_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * id_101[k];

        t_66[k] = -gd_36[k]
                  + f_0 * id_102[k];

        t_67[k] = -gd_37[k]
                  + f_0 * id_103[k];

        t_68[k] = -gd_38[k]
                  + f_0 * id_104[k];

        t_69[k] = -gd_39[k]
                  + f_0 * id_105[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, gd_40, gd_41, gd_42, gd_43, gd_44, \
                         id_106, id_107, id_108, id_109, id_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -gd_40[k]
                  + f_0 * id_106[k];

        t_71[k] = -gd_41[k]
                  + f_0 * id_107[k];

        t_72[k] = -2.0 * gd_42[k]
                  + f_0 * id_108[k];

        t_73[k] = -2.0 * gd_43[k]
                  + f_0 * id_109[k];

        t_74[k] = -2.0 * gd_44[k]
                  + f_0 * id_110[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, gd_45, gd_46, gd_47, gd_48, gd_49, \
                         id_111, id_112, id_113, id_114, id_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -2.0 * gd_45[k]
                  + f_0 * id_111[k];

        t_76[k] = -2.0 * gd_46[k]
                  + f_0 * id_112[k];

        t_77[k] = -2.0 * gd_47[k]
                  + f_0 * id_113[k];

        t_78[k] = -3.0 * gd_48[k]
                  + f_0 * id_114[k];

        t_79[k] = -3.0 * gd_49[k]
                  + f_0 * id_115[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, gd_50, gd_51, gd_52, gd_53, gd_54, \
                         id_116, id_117, id_118, id_119, id_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -3.0 * gd_50[k]
                  + f_0 * id_116[k];

        t_81[k] = -3.0 * gd_51[k]
                  + f_0 * id_117[k];

        t_82[k] = -3.0 * gd_52[k]
                  + f_0 * id_118[k];

        t_83[k] = -3.0 * gd_53[k]
                  + f_0 * id_119[k];

        t_84[k] = -4.0 * gd_54[k]
                  + f_0 * id_120[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, gd_55, gd_56, gd_57, gd_58, gd_59, \
                         id_121, id_122, id_123, id_124, id_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = -4.0 * gd_55[k]
                  + f_0 * id_121[k];

        t_86[k] = -4.0 * gd_56[k]
                  + f_0 * id_122[k];

        t_87[k] = -4.0 * gd_57[k]
                  + f_0 * id_123[k];

        t_88[k] = -4.0 * gd_58[k]
                  + f_0 * id_124[k];

        t_89[k] = -4.0 * gd_59[k]
                  + f_0 * id_125[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, t_95, t_96, gd_60, id_132, id_133, \
                         id_134, id_135, id_136, id_137, id_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_0 * id_132[k];

        t_91[k] = f_0 * id_133[k];

        t_92[k] = f_0 * id_134[k];

        t_93[k] = f_0 * id_135[k];

        t_94[k] = f_0 * id_136[k];

        t_95[k] = f_0 * id_137[k];

        t_96[k] = -gd_60[k]
                  + f_0 * id_138[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, gd_61, gd_62, gd_63, gd_64, gd_65, \
                         id_139, id_140, id_141, id_142, id_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -gd_61[k]
                  + f_0 * id_139[k];

        t_98[k] = -gd_62[k]
                  + f_0 * id_140[k];

        t_99[k] = -gd_63[k]
                  + f_0 * id_141[k];

        t_100[k] = -gd_64[k]
                   + f_0 * id_142[k];

        t_101[k] = -gd_65[k]
                   + f_0 * id_143[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, gd_66, gd_67, gd_68, gd_69, gd_70, \
                         id_144, id_145, id_146, id_147, id_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -2.0 * gd_66[k]
                   + f_0 * id_144[k];

        t_103[k] = -2.0 * gd_67[k]
                   + f_0 * id_145[k];

        t_104[k] = -2.0 * gd_68[k]
                   + f_0 * id_146[k];

        t_105[k] = -2.0 * gd_69[k]
                   + f_0 * id_147[k];

        t_106[k] = -2.0 * gd_70[k]
                   + f_0 * id_148[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, gd_71, gd_72, gd_73, gd_74, gd_75, \
                         id_149, id_150, id_151, id_152, id_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -2.0 * gd_71[k]
                   + f_0 * id_149[k];

        t_108[k] = -3.0 * gd_72[k]
                   + f_0 * id_150[k];

        t_109[k] = -3.0 * gd_73[k]
                   + f_0 * id_151[k];

        t_110[k] = -3.0 * gd_74[k]
                   + f_0 * id_152[k];

        t_111[k] = -3.0 * gd_75[k]
                   + f_0 * id_153[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, gd_76, gd_77, gd_78, gd_79, gd_80, \
                         id_154, id_155, id_156, id_157, id_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -3.0 * gd_76[k]
                   + f_0 * id_154[k];

        t_113[k] = -3.0 * gd_77[k]
                   + f_0 * id_155[k];

        t_114[k] = -4.0 * gd_78[k]
                   + f_0 * id_156[k];

        t_115[k] = -4.0 * gd_79[k]
                   + f_0 * id_157[k];

        t_116[k] = -4.0 * gd_80[k]
                   + f_0 * id_158[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, gd_81, gd_82, gd_83, gd_84, gd_85, \
                         id_159, id_160, id_161, id_162, id_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = -4.0 * gd_81[k]
                   + f_0 * id_159[k];

        t_118[k] = -4.0 * gd_82[k]
                   + f_0 * id_160[k];

        t_119[k] = -4.0 * gd_83[k]
                   + f_0 * id_161[k];

        t_120[k] = -5.0 * gd_84[k]
                   + f_0 * id_162[k];

        t_121[k] = -5.0 * gd_85[k]
                   + f_0 * id_163[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, gd_86, gd_87, gd_88, gd_89, id_164, \
                         id_165, id_166, id_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -5.0 * gd_86[k]
                   + f_0 * id_164[k];

        t_123[k] = -5.0 * gd_87[k]
                   + f_0 * id_165[k];

        t_124[k] = -5.0 * gd_88[k]
                   + f_0 * id_166[k];

        t_125[k] = -5.0 * gd_89[k]
                   + f_0 * id_167[k];
    }
}

}  // namespace simdt2ceri
