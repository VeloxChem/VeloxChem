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


#include "SimdElectronRepulsionGeom10VrrRecKP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_kp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t ip, const size_t lp,
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
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
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
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_24 = buffer.data(lp + 24);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_33 = buffer.data(lp + 33);
    const auto *lp_34 = buffer.data(lp + 34);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_36 = buffer.data(lp + 36);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_39 = buffer.data(lp + 39);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_48 = buffer.data(lp + 48);
    const auto *lp_49 = buffer.data(lp + 49);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_51 = buffer.data(lp + 51);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_54 = buffer.data(lp + 54);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_57 = buffer.data(lp + 57);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_60 = buffer.data(lp + 60);
    const auto *lp_61 = buffer.data(lp + 61);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_63 = buffer.data(lp + 63);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_65 = buffer.data(lp + 65);
    const auto *lp_66 = buffer.data(lp + 66);
    const auto *lp_67 = buffer.data(lp + 67);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_69 = buffer.data(lp + 69);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_72 = buffer.data(lp + 72);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_75 = buffer.data(lp + 75);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_78 = buffer.data(lp + 78);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_81 = buffer.data(lp + 81);
    const auto *lp_82 = buffer.data(lp + 82);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_84 = buffer.data(lp + 84);
    const auto *lp_85 = buffer.data(lp + 85);
    const auto *lp_86 = buffer.data(lp + 86);
    const auto *lp_87 = buffer.data(lp + 87);
    const auto *lp_88 = buffer.data(lp + 88);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_90 = buffer.data(lp + 90);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_93 = buffer.data(lp + 93);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_96 = buffer.data(lp + 96);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_99 = buffer.data(lp + 99);
    const auto *lp_100 = buffer.data(lp + 100);
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_102 = buffer.data(lp + 102);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_104 = buffer.data(lp + 104);
    const auto *lp_105 = buffer.data(lp + 105);
    const auto *lp_106 = buffer.data(lp + 106);
    const auto *lp_107 = buffer.data(lp + 107);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ip_0, ip_1, ip_2, ip_3, ip_4, lp_0, lp_1, \
                         lp_2, lp_3, lp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -7.0 * ip_0[k]
                 + f_0 * lp_0[k];

        t_1[k] = -7.0 * ip_1[k]
                 + f_0 * lp_1[k];

        t_2[k] = -7.0 * ip_2[k]
                 + f_0 * lp_2[k];

        t_3[k] = -6.0 * ip_3[k]
                 + f_0 * lp_3[k];

        t_4[k] = -6.0 * ip_4[k]
                 + f_0 * lp_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ip_5, ip_6, ip_7, ip_8, ip_9, lp_5, lp_6, \
                         lp_7, lp_8, lp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -6.0 * ip_5[k]
                 + f_0 * lp_5[k];

        t_6[k] = -6.0 * ip_6[k]
                 + f_0 * lp_6[k];

        t_7[k] = -6.0 * ip_7[k]
                 + f_0 * lp_7[k];

        t_8[k] = -6.0 * ip_8[k]
                 + f_0 * lp_8[k];

        t_9[k] = -5.0 * ip_9[k]
                 + f_0 * lp_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ip_10, ip_11, ip_12, ip_13, ip_14, \
                         lp_10, lp_11, lp_12, lp_13, lp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -5.0 * ip_10[k]
                  + f_0 * lp_10[k];

        t_11[k] = -5.0 * ip_11[k]
                  + f_0 * lp_11[k];

        t_12[k] = -5.0 * ip_12[k]
                  + f_0 * lp_12[k];

        t_13[k] = -5.0 * ip_13[k]
                  + f_0 * lp_13[k];

        t_14[k] = -5.0 * ip_14[k]
                  + f_0 * lp_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ip_15, ip_16, ip_17, ip_18, ip_19, \
                         lp_15, lp_16, lp_17, lp_18, lp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -5.0 * ip_15[k]
                  + f_0 * lp_15[k];

        t_16[k] = -5.0 * ip_16[k]
                  + f_0 * lp_16[k];

        t_17[k] = -5.0 * ip_17[k]
                  + f_0 * lp_17[k];

        t_18[k] = -4.0 * ip_18[k]
                  + f_0 * lp_18[k];

        t_19[k] = -4.0 * ip_19[k]
                  + f_0 * lp_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ip_20, ip_21, ip_22, ip_23, ip_24, \
                         lp_20, lp_21, lp_22, lp_23, lp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -4.0 * ip_20[k]
                  + f_0 * lp_20[k];

        t_21[k] = -4.0 * ip_21[k]
                  + f_0 * lp_21[k];

        t_22[k] = -4.0 * ip_22[k]
                  + f_0 * lp_22[k];

        t_23[k] = -4.0 * ip_23[k]
                  + f_0 * lp_23[k];

        t_24[k] = -4.0 * ip_24[k]
                  + f_0 * lp_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ip_25, ip_26, ip_27, ip_28, ip_29, \
                         lp_25, lp_26, lp_27, lp_28, lp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -4.0 * ip_25[k]
                  + f_0 * lp_25[k];

        t_26[k] = -4.0 * ip_26[k]
                  + f_0 * lp_26[k];

        t_27[k] = -4.0 * ip_27[k]
                  + f_0 * lp_27[k];

        t_28[k] = -4.0 * ip_28[k]
                  + f_0 * lp_28[k];

        t_29[k] = -4.0 * ip_29[k]
                  + f_0 * lp_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ip_30, ip_31, ip_32, ip_33, ip_34, \
                         lp_30, lp_31, lp_32, lp_33, lp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -3.0 * ip_30[k]
                  + f_0 * lp_30[k];

        t_31[k] = -3.0 * ip_31[k]
                  + f_0 * lp_31[k];

        t_32[k] = -3.0 * ip_32[k]
                  + f_0 * lp_32[k];

        t_33[k] = -3.0 * ip_33[k]
                  + f_0 * lp_33[k];

        t_34[k] = -3.0 * ip_34[k]
                  + f_0 * lp_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ip_35, ip_36, ip_37, ip_38, ip_39, \
                         lp_35, lp_36, lp_37, lp_38, lp_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -3.0 * ip_35[k]
                  + f_0 * lp_35[k];

        t_36[k] = -3.0 * ip_36[k]
                  + f_0 * lp_36[k];

        t_37[k] = -3.0 * ip_37[k]
                  + f_0 * lp_37[k];

        t_38[k] = -3.0 * ip_38[k]
                  + f_0 * lp_38[k];

        t_39[k] = -3.0 * ip_39[k]
                  + f_0 * lp_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ip_40, ip_41, ip_42, ip_43, ip_44, \
                         lp_40, lp_41, lp_42, lp_43, lp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * ip_40[k]
                  + f_0 * lp_40[k];

        t_41[k] = -3.0 * ip_41[k]
                  + f_0 * lp_41[k];

        t_42[k] = -3.0 * ip_42[k]
                  + f_0 * lp_42[k];

        t_43[k] = -3.0 * ip_43[k]
                  + f_0 * lp_43[k];

        t_44[k] = -3.0 * ip_44[k]
                  + f_0 * lp_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ip_45, ip_46, ip_47, ip_48, ip_49, \
                         lp_45, lp_46, lp_47, lp_48, lp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -2.0 * ip_45[k]
                  + f_0 * lp_45[k];

        t_46[k] = -2.0 * ip_46[k]
                  + f_0 * lp_46[k];

        t_47[k] = -2.0 * ip_47[k]
                  + f_0 * lp_47[k];

        t_48[k] = -2.0 * ip_48[k]
                  + f_0 * lp_48[k];

        t_49[k] = -2.0 * ip_49[k]
                  + f_0 * lp_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ip_50, ip_51, ip_52, ip_53, ip_54, \
                         lp_50, lp_51, lp_52, lp_53, lp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -2.0 * ip_50[k]
                  + f_0 * lp_50[k];

        t_51[k] = -2.0 * ip_51[k]
                  + f_0 * lp_51[k];

        t_52[k] = -2.0 * ip_52[k]
                  + f_0 * lp_52[k];

        t_53[k] = -2.0 * ip_53[k]
                  + f_0 * lp_53[k];

        t_54[k] = -2.0 * ip_54[k]
                  + f_0 * lp_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ip_55, ip_56, ip_57, ip_58, ip_59, \
                         lp_55, lp_56, lp_57, lp_58, lp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * ip_55[k]
                  + f_0 * lp_55[k];

        t_56[k] = -2.0 * ip_56[k]
                  + f_0 * lp_56[k];

        t_57[k] = -2.0 * ip_57[k]
                  + f_0 * lp_57[k];

        t_58[k] = -2.0 * ip_58[k]
                  + f_0 * lp_58[k];

        t_59[k] = -2.0 * ip_59[k]
                  + f_0 * lp_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ip_60, ip_61, ip_62, ip_63, ip_64, \
                         lp_60, lp_61, lp_62, lp_63, lp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -2.0 * ip_60[k]
                  + f_0 * lp_60[k];

        t_61[k] = -2.0 * ip_61[k]
                  + f_0 * lp_61[k];

        t_62[k] = -2.0 * ip_62[k]
                  + f_0 * lp_62[k];

        t_63[k] = -ip_63[k]
                  + f_0 * lp_63[k];

        t_64[k] = -ip_64[k]
                  + f_0 * lp_64[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ip_65, ip_66, ip_67, ip_68, ip_69, \
                         lp_65, lp_66, lp_67, lp_68, lp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = -ip_65[k]
                  + f_0 * lp_65[k];

        t_66[k] = -ip_66[k]
                  + f_0 * lp_66[k];

        t_67[k] = -ip_67[k]
                  + f_0 * lp_67[k];

        t_68[k] = -ip_68[k]
                  + f_0 * lp_68[k];

        t_69[k] = -ip_69[k]
                  + f_0 * lp_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ip_70, ip_71, ip_72, ip_73, ip_74, \
                         lp_70, lp_71, lp_72, lp_73, lp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = -ip_70[k]
                  + f_0 * lp_70[k];

        t_71[k] = -ip_71[k]
                  + f_0 * lp_71[k];

        t_72[k] = -ip_72[k]
                  + f_0 * lp_72[k];

        t_73[k] = -ip_73[k]
                  + f_0 * lp_73[k];

        t_74[k] = -ip_74[k]
                  + f_0 * lp_74[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ip_75, ip_76, ip_77, ip_78, ip_79, \
                         lp_75, lp_76, lp_77, lp_78, lp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -ip_75[k]
                  + f_0 * lp_75[k];

        t_76[k] = -ip_76[k]
                  + f_0 * lp_76[k];

        t_77[k] = -ip_77[k]
                  + f_0 * lp_77[k];

        t_78[k] = -ip_78[k]
                  + f_0 * lp_78[k];

        t_79[k] = -ip_79[k]
                  + f_0 * lp_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, t_85, ip_80, ip_81, ip_82, ip_83, \
                         lp_80, lp_81, lp_82, lp_83, lp_84, lp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = -ip_80[k]
                  + f_0 * lp_80[k];

        t_81[k] = -ip_81[k]
                  + f_0 * lp_81[k];

        t_82[k] = -ip_82[k]
                  + f_0 * lp_82[k];

        t_83[k] = -ip_83[k]
                  + f_0 * lp_83[k];

        t_84[k] = f_0 * lp_84[k];

        t_85[k] = f_0 * lp_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, t_90, t_91, t_92, t_93, lp_86, lp_87, lp_88, \
                         lp_89, lp_90, lp_91, lp_92, lp_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_0 * lp_86[k];

        t_87[k] = f_0 * lp_87[k];

        t_88[k] = f_0 * lp_88[k];

        t_89[k] = f_0 * lp_89[k];

        t_90[k] = f_0 * lp_90[k];

        t_91[k] = f_0 * lp_91[k];

        t_92[k] = f_0 * lp_92[k];

        t_93[k] = f_0 * lp_93[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, t_99, t_100, t_101, lp_94, lp_95, \
                         lp_96, lp_97, lp_98, lp_99, lp_100, lp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_0 * lp_94[k];

        t_95[k] = f_0 * lp_95[k];

        t_96[k] = f_0 * lp_96[k];

        t_97[k] = f_0 * lp_97[k];

        t_98[k] = f_0 * lp_98[k];

        t_99[k] = f_0 * lp_99[k];

        t_100[k] = f_0 * lp_100[k];

        t_101[k] = f_0 * lp_101[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, lp_102, lp_103, lp_104, \
                         lp_105, lp_106, lp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_0 * lp_102[k];

        t_103[k] = f_0 * lp_103[k];

        t_104[k] = f_0 * lp_104[k];

        t_105[k] = f_0 * lp_105[k];

        t_106[k] = f_0 * lp_106[k];

        t_107[k] = f_0 * lp_107[k];
    }
}

auto
compute_prim_geom_10_kp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t ip, const size_t lp,
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
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
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
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *lp_3 = buffer.data(lp + 3);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_5 = buffer.data(lp + 5);
    const auto *lp_9 = buffer.data(lp + 9);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_18 = buffer.data(lp + 18);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_24 = buffer.data(lp + 24);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_30 = buffer.data(lp + 30);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_33 = buffer.data(lp + 33);
    const auto *lp_34 = buffer.data(lp + 34);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_36 = buffer.data(lp + 36);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_39 = buffer.data(lp + 39);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_45 = buffer.data(lp + 45);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_48 = buffer.data(lp + 48);
    const auto *lp_49 = buffer.data(lp + 49);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_51 = buffer.data(lp + 51);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_54 = buffer.data(lp + 54);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_57 = buffer.data(lp + 57);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_63 = buffer.data(lp + 63);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_65 = buffer.data(lp + 65);
    const auto *lp_66 = buffer.data(lp + 66);
    const auto *lp_67 = buffer.data(lp + 67);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_69 = buffer.data(lp + 69);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_72 = buffer.data(lp + 72);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_75 = buffer.data(lp + 75);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_78 = buffer.data(lp + 78);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_84 = buffer.data(lp + 84);
    const auto *lp_85 = buffer.data(lp + 85);
    const auto *lp_86 = buffer.data(lp + 86);
    const auto *lp_87 = buffer.data(lp + 87);
    const auto *lp_88 = buffer.data(lp + 88);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_90 = buffer.data(lp + 90);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_93 = buffer.data(lp + 93);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_96 = buffer.data(lp + 96);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_99 = buffer.data(lp + 99);
    const auto *lp_100 = buffer.data(lp + 100);
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_102 = buffer.data(lp + 102);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_104 = buffer.data(lp + 104);
    const auto *lp_108 = buffer.data(lp + 108);
    const auto *lp_109 = buffer.data(lp + 109);
    const auto *lp_110 = buffer.data(lp + 110);
    const auto *lp_111 = buffer.data(lp + 111);
    const auto *lp_112 = buffer.data(lp + 112);
    const auto *lp_113 = buffer.data(lp + 113);
    const auto *lp_114 = buffer.data(lp + 114);
    const auto *lp_115 = buffer.data(lp + 115);
    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_117 = buffer.data(lp + 117);
    const auto *lp_118 = buffer.data(lp + 118);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_120 = buffer.data(lp + 120);
    const auto *lp_121 = buffer.data(lp + 121);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_123 = buffer.data(lp + 123);
    const auto *lp_124 = buffer.data(lp + 124);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_126 = buffer.data(lp + 126);
    const auto *lp_127 = buffer.data(lp + 127);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_129 = buffer.data(lp + 129);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_131 = buffer.data(lp + 131);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ip_0, ip_1, ip_2, lp_3, lp_4, lp_5, \
                         lp_9, lp_10, lp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lp_3[k];

        t_1[k] = f_0 * lp_4[k];

        t_2[k] = f_0 * lp_5[k];

        t_3[k] = -ip_0[k]
                 + f_0 * lp_9[k];

        t_4[k] = -ip_1[k]
                 + f_0 * lp_10[k];

        t_5[k] = -ip_2[k]
                 + f_0 * lp_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, ip_3, ip_4, ip_5, lp_12, lp_13, \
                         lp_14, lp_18, lp_19, lp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * lp_12[k];

        t_7[k] = f_0 * lp_13[k];

        t_8[k] = f_0 * lp_14[k];

        t_9[k] = -2.0 * ip_3[k]
                 + f_0 * lp_18[k];

        t_10[k] = -2.0 * ip_4[k]
                  + f_0 * lp_19[k];

        t_11[k] = -2.0 * ip_5[k]
                  + f_0 * lp_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, ip_6, ip_7, ip_8, lp_21, lp_22, \
                         lp_23, lp_24, lp_25, lp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -ip_6[k]
                  + f_0 * lp_21[k];

        t_13[k] = -ip_7[k]
                  + f_0 * lp_22[k];

        t_14[k] = -ip_8[k]
                  + f_0 * lp_23[k];

        t_15[k] = f_0 * lp_24[k];

        t_16[k] = f_0 * lp_25[k];

        t_17[k] = f_0 * lp_26[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ip_9, ip_10, ip_11, ip_12, ip_13, \
                         lp_30, lp_31, lp_32, lp_33, lp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -3.0 * ip_9[k]
                  + f_0 * lp_30[k];

        t_19[k] = -3.0 * ip_10[k]
                  + f_0 * lp_31[k];

        t_20[k] = -3.0 * ip_11[k]
                  + f_0 * lp_32[k];

        t_21[k] = -2.0 * ip_12[k]
                  + f_0 * lp_33[k];

        t_22[k] = -2.0 * ip_13[k]
                  + f_0 * lp_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, ip_14, ip_15, ip_16, ip_17, \
                         lp_35, lp_36, lp_37, lp_38, lp_39, lp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -2.0 * ip_14[k]
                  + f_0 * lp_35[k];

        t_24[k] = -ip_15[k]
                  + f_0 * lp_36[k];

        t_25[k] = -ip_16[k]
                  + f_0 * lp_37[k];

        t_26[k] = -ip_17[k]
                  + f_0 * lp_38[k];

        t_27[k] = f_0 * lp_39[k];

        t_28[k] = f_0 * lp_40[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ip_18, ip_19, ip_20, ip_21, lp_41, \
                         lp_45, lp_46, lp_47, lp_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * lp_41[k];

        t_30[k] = -4.0 * ip_18[k]
                  + f_0 * lp_45[k];

        t_31[k] = -4.0 * ip_19[k]
                  + f_0 * lp_46[k];

        t_32[k] = -4.0 * ip_20[k]
                  + f_0 * lp_47[k];

        t_33[k] = -3.0 * ip_21[k]
                  + f_0 * lp_48[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, ip_22, ip_23, ip_24, ip_25, ip_26, \
                         lp_49, lp_50, lp_51, lp_52, lp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -3.0 * ip_22[k]
                  + f_0 * lp_49[k];

        t_35[k] = -3.0 * ip_23[k]
                  + f_0 * lp_50[k];

        t_36[k] = -2.0 * ip_24[k]
                  + f_0 * lp_51[k];

        t_37[k] = -2.0 * ip_25[k]
                  + f_0 * lp_52[k];

        t_38[k] = -2.0 * ip_26[k]
                  + f_0 * lp_53[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, ip_27, ip_28, ip_29, lp_54, \
                         lp_55, lp_56, lp_57, lp_58, lp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -ip_27[k]
                  + f_0 * lp_54[k];

        t_40[k] = -ip_28[k]
                  + f_0 * lp_55[k];

        t_41[k] = -ip_29[k]
                  + f_0 * lp_56[k];

        t_42[k] = f_0 * lp_57[k];

        t_43[k] = f_0 * lp_58[k];

        t_44[k] = f_0 * lp_59[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ip_30, ip_31, ip_32, ip_33, ip_34, \
                         lp_63, lp_64, lp_65, lp_66, lp_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * ip_30[k]
                  + f_0 * lp_63[k];

        t_46[k] = -5.0 * ip_31[k]
                  + f_0 * lp_64[k];

        t_47[k] = -5.0 * ip_32[k]
                  + f_0 * lp_65[k];

        t_48[k] = -4.0 * ip_33[k]
                  + f_0 * lp_66[k];

        t_49[k] = -4.0 * ip_34[k]
                  + f_0 * lp_67[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ip_35, ip_36, ip_37, ip_38, ip_39, \
                         lp_68, lp_69, lp_70, lp_71, lp_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * ip_35[k]
                  + f_0 * lp_68[k];

        t_51[k] = -3.0 * ip_36[k]
                  + f_0 * lp_69[k];

        t_52[k] = -3.0 * ip_37[k]
                  + f_0 * lp_70[k];

        t_53[k] = -3.0 * ip_38[k]
                  + f_0 * lp_71[k];

        t_54[k] = -2.0 * ip_39[k]
                  + f_0 * lp_72[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ip_40, ip_41, ip_42, ip_43, ip_44, \
                         lp_73, lp_74, lp_75, lp_76, lp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * ip_40[k]
                  + f_0 * lp_73[k];

        t_56[k] = -2.0 * ip_41[k]
                  + f_0 * lp_74[k];

        t_57[k] = -ip_42[k]
                  + f_0 * lp_75[k];

        t_58[k] = -ip_43[k]
                  + f_0 * lp_76[k];

        t_59[k] = -ip_44[k]
                  + f_0 * lp_77[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, ip_45, ip_46, ip_47, lp_78, \
                         lp_79, lp_80, lp_84, lp_85, lp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * lp_78[k];

        t_61[k] = f_0 * lp_79[k];

        t_62[k] = f_0 * lp_80[k];

        t_63[k] = -6.0 * ip_45[k]
                  + f_0 * lp_84[k];

        t_64[k] = -6.0 * ip_46[k]
                  + f_0 * lp_85[k];

        t_65[k] = -6.0 * ip_47[k]
                  + f_0 * lp_86[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, ip_48, ip_49, ip_50, ip_51, ip_52, \
                         lp_87, lp_88, lp_89, lp_90, lp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -5.0 * ip_48[k]
                  + f_0 * lp_87[k];

        t_67[k] = -5.0 * ip_49[k]
                  + f_0 * lp_88[k];

        t_68[k] = -5.0 * ip_50[k]
                  + f_0 * lp_89[k];

        t_69[k] = -4.0 * ip_51[k]
                  + f_0 * lp_90[k];

        t_70[k] = -4.0 * ip_52[k]
                  + f_0 * lp_91[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, t_75, ip_53, ip_54, ip_55, ip_56, ip_57, \
                         lp_92, lp_93, lp_94, lp_95, lp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -4.0 * ip_53[k]
                  + f_0 * lp_92[k];

        t_72[k] = -3.0 * ip_54[k]
                  + f_0 * lp_93[k];

        t_73[k] = -3.0 * ip_55[k]
                  + f_0 * lp_94[k];

        t_74[k] = -3.0 * ip_56[k]
                  + f_0 * lp_95[k];

        t_75[k] = -2.0 * ip_57[k]
                  + f_0 * lp_96[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, ip_58, ip_59, ip_60, ip_61, ip_62, \
                         lp_97, lp_98, lp_99, lp_100, lp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -2.0 * ip_58[k]
                  + f_0 * lp_97[k];

        t_77[k] = -2.0 * ip_59[k]
                  + f_0 * lp_98[k];

        t_78[k] = -ip_60[k]
                  + f_0 * lp_99[k];

        t_79[k] = -ip_61[k]
                  + f_0 * lp_100[k];

        t_80[k] = -ip_62[k]
                  + f_0 * lp_101[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, t_85, t_86, ip_63, ip_64, ip_65, lp_102, \
                         lp_103, lp_104, lp_108, lp_109, lp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_0 * lp_102[k];

        t_82[k] = f_0 * lp_103[k];

        t_83[k] = f_0 * lp_104[k];

        t_84[k] = -7.0 * ip_63[k]
                  + f_0 * lp_108[k];

        t_85[k] = -7.0 * ip_64[k]
                  + f_0 * lp_109[k];

        t_86[k] = -7.0 * ip_65[k]
                  + f_0 * lp_110[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, ip_66, ip_67, ip_68, ip_69, ip_70, \
                         lp_111, lp_112, lp_113, lp_114, lp_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -6.0 * ip_66[k]
                  + f_0 * lp_111[k];

        t_88[k] = -6.0 * ip_67[k]
                  + f_0 * lp_112[k];

        t_89[k] = -6.0 * ip_68[k]
                  + f_0 * lp_113[k];

        t_90[k] = -5.0 * ip_69[k]
                  + f_0 * lp_114[k];

        t_91[k] = -5.0 * ip_70[k]
                  + f_0 * lp_115[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, ip_71, ip_72, ip_73, ip_74, ip_75, \
                         lp_116, lp_117, lp_118, lp_119, lp_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = -5.0 * ip_71[k]
                  + f_0 * lp_116[k];

        t_93[k] = -4.0 * ip_72[k]
                  + f_0 * lp_117[k];

        t_94[k] = -4.0 * ip_73[k]
                  + f_0 * lp_118[k];

        t_95[k] = -4.0 * ip_74[k]
                  + f_0 * lp_119[k];

        t_96[k] = -3.0 * ip_75[k]
                  + f_0 * lp_120[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, ip_76, ip_77, ip_78, ip_79, ip_80, \
                         lp_121, lp_122, lp_123, lp_124, lp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = -3.0 * ip_76[k]
                  + f_0 * lp_121[k];

        t_98[k] = -3.0 * ip_77[k]
                  + f_0 * lp_122[k];

        t_99[k] = -2.0 * ip_78[k]
                  + f_0 * lp_123[k];

        t_100[k] = -2.0 * ip_79[k]
                   + f_0 * lp_124[k];

        t_101[k] = -2.0 * ip_80[k]
                   + f_0 * lp_125[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, t_107, ip_81, ip_82, ip_83, \
                         lp_126, lp_127, lp_128, lp_129, lp_130, \
                         lp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = -ip_81[k]
                   + f_0 * lp_126[k];

        t_103[k] = -ip_82[k]
                   + f_0 * lp_127[k];

        t_104[k] = -ip_83[k]
                   + f_0 * lp_128[k];

        t_105[k] = f_0 * lp_129[k];

        t_106[k] = f_0 * lp_130[k];

        t_107[k] = f_0 * lp_131[k];
    }
}

auto
compute_prim_geom_10_kp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t ip, const size_t lp,
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
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
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
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *lp_6 = buffer.data(lp + 6);
    const auto *lp_7 = buffer.data(lp + 7);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_12 = buffer.data(lp + 12);
    const auto *lp_13 = buffer.data(lp + 13);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_15 = buffer.data(lp + 15);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_21 = buffer.data(lp + 21);
    const auto *lp_22 = buffer.data(lp + 22);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_24 = buffer.data(lp + 24);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_27 = buffer.data(lp + 27);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_33 = buffer.data(lp + 33);
    const auto *lp_34 = buffer.data(lp + 34);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_36 = buffer.data(lp + 36);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_39 = buffer.data(lp + 39);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_42 = buffer.data(lp + 42);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_48 = buffer.data(lp + 48);
    const auto *lp_49 = buffer.data(lp + 49);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_51 = buffer.data(lp + 51);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_54 = buffer.data(lp + 54);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_57 = buffer.data(lp + 57);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_60 = buffer.data(lp + 60);
    const auto *lp_61 = buffer.data(lp + 61);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_66 = buffer.data(lp + 66);
    const auto *lp_67 = buffer.data(lp + 67);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_69 = buffer.data(lp + 69);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_72 = buffer.data(lp + 72);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_75 = buffer.data(lp + 75);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_78 = buffer.data(lp + 78);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_81 = buffer.data(lp + 81);
    const auto *lp_82 = buffer.data(lp + 82);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_87 = buffer.data(lp + 87);
    const auto *lp_88 = buffer.data(lp + 88);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_90 = buffer.data(lp + 90);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_93 = buffer.data(lp + 93);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_96 = buffer.data(lp + 96);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_99 = buffer.data(lp + 99);
    const auto *lp_100 = buffer.data(lp + 100);
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_102 = buffer.data(lp + 102);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_104 = buffer.data(lp + 104);
    const auto *lp_105 = buffer.data(lp + 105);
    const auto *lp_106 = buffer.data(lp + 106);
    const auto *lp_107 = buffer.data(lp + 107);
    const auto *lp_111 = buffer.data(lp + 111);
    const auto *lp_112 = buffer.data(lp + 112);
    const auto *lp_113 = buffer.data(lp + 113);
    const auto *lp_114 = buffer.data(lp + 114);
    const auto *lp_115 = buffer.data(lp + 115);
    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_117 = buffer.data(lp + 117);
    const auto *lp_118 = buffer.data(lp + 118);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_120 = buffer.data(lp + 120);
    const auto *lp_121 = buffer.data(lp + 121);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_123 = buffer.data(lp + 123);
    const auto *lp_124 = buffer.data(lp + 124);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_126 = buffer.data(lp + 126);
    const auto *lp_127 = buffer.data(lp + 127);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_129 = buffer.data(lp + 129);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_131 = buffer.data(lp + 131);
    const auto *lp_132 = buffer.data(lp + 132);
    const auto *lp_133 = buffer.data(lp + 133);
    const auto *lp_134 = buffer.data(lp + 134);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, ip_0, lp_6, lp_7, lp_8, lp_12, \
                         lp_13, lp_14, lp_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lp_6[k];

        t_1[k] = f_0 * lp_7[k];

        t_2[k] = f_0 * lp_8[k];

        t_3[k] = f_0 * lp_12[k];

        t_4[k] = f_0 * lp_13[k];

        t_5[k] = f_0 * lp_14[k];

        t_6[k] = -ip_0[k]
                 + f_0 * lp_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, ip_1, ip_2, ip_3, lp_16, lp_17, \
                         lp_21, lp_22, lp_23, lp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -ip_1[k]
                 + f_0 * lp_16[k];

        t_8[k] = -ip_2[k]
                 + f_0 * lp_17[k];

        t_9[k] = f_0 * lp_21[k];

        t_10[k] = f_0 * lp_22[k];

        t_11[k] = f_0 * lp_23[k];

        t_12[k] = -ip_3[k]
                  + f_0 * lp_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ip_4, ip_5, ip_6, ip_7, ip_8, lp_25, \
                         lp_26, lp_27, lp_28, lp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -ip_4[k]
                  + f_0 * lp_25[k];

        t_14[k] = -ip_5[k]
                  + f_0 * lp_26[k];

        t_15[k] = -2.0 * ip_6[k]
                  + f_0 * lp_27[k];

        t_16[k] = -2.0 * ip_7[k]
                  + f_0 * lp_28[k];

        t_17[k] = -2.0 * ip_8[k]
                  + f_0 * lp_29[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, ip_9, ip_10, ip_11, lp_33, lp_34, \
                         lp_35, lp_36, lp_37, lp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * lp_33[k];

        t_19[k] = f_0 * lp_34[k];

        t_20[k] = f_0 * lp_35[k];

        t_21[k] = -ip_9[k]
                  + f_0 * lp_36[k];

        t_22[k] = -ip_10[k]
                  + f_0 * lp_37[k];

        t_23[k] = -ip_11[k]
                  + f_0 * lp_38[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ip_12, ip_13, ip_14, ip_15, ip_16, \
                         lp_39, lp_40, lp_41, lp_42, lp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -2.0 * ip_12[k]
                  + f_0 * lp_39[k];

        t_25[k] = -2.0 * ip_13[k]
                  + f_0 * lp_40[k];

        t_26[k] = -2.0 * ip_14[k]
                  + f_0 * lp_41[k];

        t_27[k] = -3.0 * ip_15[k]
                  + f_0 * lp_42[k];

        t_28[k] = -3.0 * ip_16[k]
                  + f_0 * lp_43[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, ip_17, ip_18, ip_19, lp_44, \
                         lp_48, lp_49, lp_50, lp_51, lp_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -3.0 * ip_17[k]
                  + f_0 * lp_44[k];

        t_30[k] = f_0 * lp_48[k];

        t_31[k] = f_0 * lp_49[k];

        t_32[k] = f_0 * lp_50[k];

        t_33[k] = -ip_18[k]
                  + f_0 * lp_51[k];

        t_34[k] = -ip_19[k]
                  + f_0 * lp_52[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ip_20, ip_21, ip_22, ip_23, ip_24, \
                         lp_53, lp_54, lp_55, lp_56, lp_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -ip_20[k]
                  + f_0 * lp_53[k];

        t_36[k] = -2.0 * ip_21[k]
                  + f_0 * lp_54[k];

        t_37[k] = -2.0 * ip_22[k]
                  + f_0 * lp_55[k];

        t_38[k] = -2.0 * ip_23[k]
                  + f_0 * lp_56[k];

        t_39[k] = -3.0 * ip_24[k]
                  + f_0 * lp_57[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ip_25, ip_26, ip_27, ip_28, ip_29, \
                         lp_58, lp_59, lp_60, lp_61, lp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * ip_25[k]
                  + f_0 * lp_58[k];

        t_41[k] = -3.0 * ip_26[k]
                  + f_0 * lp_59[k];

        t_42[k] = -4.0 * ip_27[k]
                  + f_0 * lp_60[k];

        t_43[k] = -4.0 * ip_28[k]
                  + f_0 * lp_61[k];

        t_44[k] = -4.0 * ip_29[k]
                  + f_0 * lp_62[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, ip_30, ip_31, ip_32, lp_66, \
                         lp_67, lp_68, lp_69, lp_70, lp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * lp_66[k];

        t_46[k] = f_0 * lp_67[k];

        t_47[k] = f_0 * lp_68[k];

        t_48[k] = -ip_30[k]
                  + f_0 * lp_69[k];

        t_49[k] = -ip_31[k]
                  + f_0 * lp_70[k];

        t_50[k] = -ip_32[k]
                  + f_0 * lp_71[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, ip_33, ip_34, ip_35, ip_36, ip_37, \
                         lp_72, lp_73, lp_74, lp_75, lp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -2.0 * ip_33[k]
                  + f_0 * lp_72[k];

        t_52[k] = -2.0 * ip_34[k]
                  + f_0 * lp_73[k];

        t_53[k] = -2.0 * ip_35[k]
                  + f_0 * lp_74[k];

        t_54[k] = -3.0 * ip_36[k]
                  + f_0 * lp_75[k];

        t_55[k] = -3.0 * ip_37[k]
                  + f_0 * lp_76[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, ip_38, ip_39, ip_40, ip_41, ip_42, \
                         lp_77, lp_78, lp_79, lp_80, lp_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -3.0 * ip_38[k]
                  + f_0 * lp_77[k];

        t_57[k] = -4.0 * ip_39[k]
                  + f_0 * lp_78[k];

        t_58[k] = -4.0 * ip_40[k]
                  + f_0 * lp_79[k];

        t_59[k] = -4.0 * ip_41[k]
                  + f_0 * lp_80[k];

        t_60[k] = -5.0 * ip_42[k]
                  + f_0 * lp_81[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, ip_43, ip_44, ip_45, lp_82, \
                         lp_83, lp_87, lp_88, lp_89, lp_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -5.0 * ip_43[k]
                  + f_0 * lp_82[k];

        t_62[k] = -5.0 * ip_44[k]
                  + f_0 * lp_83[k];

        t_63[k] = f_0 * lp_87[k];

        t_64[k] = f_0 * lp_88[k];

        t_65[k] = f_0 * lp_89[k];

        t_66[k] = -ip_45[k]
                  + f_0 * lp_90[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ip_46, ip_47, ip_48, ip_49, ip_50, \
                         lp_91, lp_92, lp_93, lp_94, lp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -ip_46[k]
                  + f_0 * lp_91[k];

        t_68[k] = -ip_47[k]
                  + f_0 * lp_92[k];

        t_69[k] = -2.0 * ip_48[k]
                  + f_0 * lp_93[k];

        t_70[k] = -2.0 * ip_49[k]
                  + f_0 * lp_94[k];

        t_71[k] = -2.0 * ip_50[k]
                  + f_0 * lp_95[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ip_51, ip_52, ip_53, ip_54, ip_55, \
                         lp_96, lp_97, lp_98, lp_99, lp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -3.0 * ip_51[k]
                  + f_0 * lp_96[k];

        t_73[k] = -3.0 * ip_52[k]
                  + f_0 * lp_97[k];

        t_74[k] = -3.0 * ip_53[k]
                  + f_0 * lp_98[k];

        t_75[k] = -4.0 * ip_54[k]
                  + f_0 * lp_99[k];

        t_76[k] = -4.0 * ip_55[k]
                  + f_0 * lp_100[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, ip_56, ip_57, ip_58, ip_59, ip_60, \
                         lp_101, lp_102, lp_103, lp_104, lp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -4.0 * ip_56[k]
                  + f_0 * lp_101[k];

        t_78[k] = -5.0 * ip_57[k]
                  + f_0 * lp_102[k];

        t_79[k] = -5.0 * ip_58[k]
                  + f_0 * lp_103[k];

        t_80[k] = -5.0 * ip_59[k]
                  + f_0 * lp_104[k];

        t_81[k] = -6.0 * ip_60[k]
                  + f_0 * lp_105[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, t_87, ip_61, ip_62, ip_63, lp_106, \
                         lp_107, lp_111, lp_112, lp_113, lp_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = -6.0 * ip_61[k]
                  + f_0 * lp_106[k];

        t_83[k] = -6.0 * ip_62[k]
                  + f_0 * lp_107[k];

        t_84[k] = f_0 * lp_111[k];

        t_85[k] = f_0 * lp_112[k];

        t_86[k] = f_0 * lp_113[k];

        t_87[k] = -ip_63[k]
                  + f_0 * lp_114[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, ip_64, ip_65, ip_66, ip_67, ip_68, \
                         lp_115, lp_116, lp_117, lp_118, lp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = -ip_64[k]
                  + f_0 * lp_115[k];

        t_89[k] = -ip_65[k]
                  + f_0 * lp_116[k];

        t_90[k] = -2.0 * ip_66[k]
                  + f_0 * lp_117[k];

        t_91[k] = -2.0 * ip_67[k]
                  + f_0 * lp_118[k];

        t_92[k] = -2.0 * ip_68[k]
                  + f_0 * lp_119[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, ip_69, ip_70, ip_71, ip_72, ip_73, \
                         lp_120, lp_121, lp_122, lp_123, lp_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = -3.0 * ip_69[k]
                  + f_0 * lp_120[k];

        t_94[k] = -3.0 * ip_70[k]
                  + f_0 * lp_121[k];

        t_95[k] = -3.0 * ip_71[k]
                  + f_0 * lp_122[k];

        t_96[k] = -4.0 * ip_72[k]
                  + f_0 * lp_123[k];

        t_97[k] = -4.0 * ip_73[k]
                  + f_0 * lp_124[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, ip_74, ip_75, ip_76, ip_77, ip_78, \
                         lp_125, lp_126, lp_127, lp_128, lp_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = -4.0 * ip_74[k]
                  + f_0 * lp_125[k];

        t_99[k] = -5.0 * ip_75[k]
                  + f_0 * lp_126[k];

        t_100[k] = -5.0 * ip_76[k]
                   + f_0 * lp_127[k];

        t_101[k] = -5.0 * ip_77[k]
                   + f_0 * lp_128[k];

        t_102[k] = -6.0 * ip_78[k]
                   + f_0 * lp_129[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, ip_79, ip_80, ip_81, ip_82, ip_83, \
                         lp_130, lp_131, lp_132, lp_133, lp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = -6.0 * ip_79[k]
                   + f_0 * lp_130[k];

        t_104[k] = -6.0 * ip_80[k]
                   + f_0 * lp_131[k];

        t_105[k] = -7.0 * ip_81[k]
                   + f_0 * lp_132[k];

        t_106[k] = -7.0 * ip_82[k]
                   + f_0 * lp_133[k];

        t_107[k] = -7.0 * ip_83[k]
                   + f_0 * lp_134[k];
    }
}

}  // namespace simdt2ceri
