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


#include "SimdKineticEnergyVrrRecFH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_fh_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t ph_s,
                                        const size_t ph, const size_t dg, const size_t dh,
                                        const size_t ff_s, const size_t fh_s, const size_t ff,
                                        const size_t fg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_36 = buffer.data(ph_s + 36);
    const auto *ph_s_62 = buffer.data(ph_s + 62);

    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_62 = buffer.data(ph + 62);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_80 = buffer.data(dg + 80);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_7 = buffer.data(dh + 7);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_110 = buffer.data(dh + 110);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_1 = buffer.data(ff_s + 1);
    const auto *ff_s_2 = buffer.data(ff_s + 2);
    const auto *ff_s_6 = buffer.data(ff_s + 6);
    const auto *ff_s_8 = buffer.data(ff_s + 8);
    const auto *ff_s_9 = buffer.data(ff_s + 9);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_17 = buffer.data(ff_s + 17);
    const auto *ff_s_27 = buffer.data(ff_s + 27);
    const auto *ff_s_28 = buffer.data(ff_s + 28);
    const auto *ff_s_29 = buffer.data(ff_s + 29);

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_1 = buffer.data(fh_s + 1);
    const auto *fh_s_2 = buffer.data(fh_s + 2);
    const auto *fh_s_3 = buffer.data(fh_s + 3);
    const auto *fh_s_4 = buffer.data(fh_s + 4);
    const auto *fh_s_5 = buffer.data(fh_s + 5);
    const auto *fh_s_6 = buffer.data(fh_s + 6);
    const auto *fh_s_7 = buffer.data(fh_s + 7);
    const auto *fh_s_8 = buffer.data(fh_s + 8);
    const auto *fh_s_9 = buffer.data(fh_s + 9);
    const auto *fh_s_10 = buffer.data(fh_s + 10);
    const auto *fh_s_11 = buffer.data(fh_s + 11);
    const auto *fh_s_12 = buffer.data(fh_s + 12);
    const auto *fh_s_13 = buffer.data(fh_s + 13);
    const auto *fh_s_14 = buffer.data(fh_s + 14);
    const auto *fh_s_15 = buffer.data(fh_s + 15);
    const auto *fh_s_16 = buffer.data(fh_s + 16);
    const auto *fh_s_17 = buffer.data(fh_s + 17);
    const auto *fh_s_18 = buffer.data(fh_s + 18);
    const auto *fh_s_19 = buffer.data(fh_s + 19);
    const auto *fh_s_20 = buffer.data(fh_s + 20);
    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_22 = buffer.data(fh_s + 22);
    const auto *fh_s_23 = buffer.data(fh_s + 23);
    const auto *fh_s_24 = buffer.data(fh_s + 24);
    const auto *fh_s_25 = buffer.data(fh_s + 25);
    const auto *fh_s_26 = buffer.data(fh_s + 26);
    const auto *fh_s_27 = buffer.data(fh_s + 27);
    const auto *fh_s_28 = buffer.data(fh_s + 28);
    const auto *fh_s_29 = buffer.data(fh_s + 29);
    const auto *fh_s_30 = buffer.data(fh_s + 30);
    const auto *fh_s_31 = buffer.data(fh_s + 31);
    const auto *fh_s_32 = buffer.data(fh_s + 32);
    const auto *fh_s_33 = buffer.data(fh_s + 33);
    const auto *fh_s_34 = buffer.data(fh_s + 34);
    const auto *fh_s_35 = buffer.data(fh_s + 35);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_37 = buffer.data(fh_s + 37);
    const auto *fh_s_38 = buffer.data(fh_s + 38);
    const auto *fh_s_39 = buffer.data(fh_s + 39);
    const auto *fh_s_40 = buffer.data(fh_s + 40);
    const auto *fh_s_41 = buffer.data(fh_s + 41);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_43 = buffer.data(fh_s + 43);
    const auto *fh_s_44 = buffer.data(fh_s + 44);
    const auto *fh_s_45 = buffer.data(fh_s + 45);
    const auto *fh_s_46 = buffer.data(fh_s + 46);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_48 = buffer.data(fh_s + 48);
    const auto *fh_s_49 = buffer.data(fh_s + 49);
    const auto *fh_s_50 = buffer.data(fh_s + 50);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_52 = buffer.data(fh_s + 52);
    const auto *fh_s_53 = buffer.data(fh_s + 53);
    const auto *fh_s_54 = buffer.data(fh_s + 54);
    const auto *fh_s_55 = buffer.data(fh_s + 55);
    const auto *fh_s_56 = buffer.data(fh_s + 56);
    const auto *fh_s_57 = buffer.data(fh_s + 57);
    const auto *fh_s_58 = buffer.data(fh_s + 58);
    const auto *fh_s_59 = buffer.data(fh_s + 59);
    const auto *fh_s_60 = buffer.data(fh_s + 60);
    const auto *fh_s_61 = buffer.data(fh_s + 61);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_63 = buffer.data(fh_s + 63);
    const auto *fh_s_64 = buffer.data(fh_s + 64);
    const auto *fh_s_65 = buffer.data(fh_s + 65);
    const auto *fh_s_66 = buffer.data(fh_s + 66);
    const auto *fh_s_67 = buffer.data(fh_s + 67);
    const auto *fh_s_68 = buffer.data(fh_s + 68);
    const auto *fh_s_69 = buffer.data(fh_s + 69);
    const auto *fh_s_70 = buffer.data(fh_s + 70);
    const auto *fh_s_71 = buffer.data(fh_s + 71);
    const auto *fh_s_72 = buffer.data(fh_s + 72);
    const auto *fh_s_73 = buffer.data(fh_s + 73);
    const auto *fh_s_74 = buffer.data(fh_s + 74);
    const auto *fh_s_75 = buffer.data(fh_s + 75);
    const auto *fh_s_76 = buffer.data(fh_s + 76);
    const auto *fh_s_77 = buffer.data(fh_s + 77);
    const auto *fh_s_78 = buffer.data(fh_s + 78);
    const auto *fh_s_79 = buffer.data(fh_s + 79);
    const auto *fh_s_80 = buffer.data(fh_s + 80);
    const auto *fh_s_81 = buffer.data(fh_s + 81);
    const auto *fh_s_82 = buffer.data(fh_s + 82);
    const auto *fh_s_83 = buffer.data(fh_s + 83);
    const auto *fh_s_84 = buffer.data(fh_s + 84);
    const auto *fh_s_85 = buffer.data(fh_s + 85);
    const auto *fh_s_86 = buffer.data(fh_s + 86);
    const auto *fh_s_87 = buffer.data(fh_s + 87);
    const auto *fh_s_88 = buffer.data(fh_s + 88);
    const auto *fh_s_89 = buffer.data(fh_s + 89);
    const auto *fh_s_90 = buffer.data(fh_s + 90);
    const auto *fh_s_91 = buffer.data(fh_s + 91);
    const auto *fh_s_92 = buffer.data(fh_s + 92);
    const auto *fh_s_93 = buffer.data(fh_s + 93);
    const auto *fh_s_94 = buffer.data(fh_s + 94);
    const auto *fh_s_95 = buffer.data(fh_s + 95);
    const auto *fh_s_96 = buffer.data(fh_s + 96);
    const auto *fh_s_97 = buffer.data(fh_s + 97);
    const auto *fh_s_98 = buffer.data(fh_s + 98);
    const auto *fh_s_99 = buffer.data(fh_s + 99);
    const auto *fh_s_100 = buffer.data(fh_s + 100);
    const auto *fh_s_101 = buffer.data(fh_s + 101);
    const auto *fh_s_102 = buffer.data(fh_s + 102);
    const auto *fh_s_103 = buffer.data(fh_s + 103);
    const auto *fh_s_104 = buffer.data(fh_s + 104);
    const auto *fh_s_105 = buffer.data(fh_s + 105);
    const auto *fh_s_106 = buffer.data(fh_s + 106);
    const auto *fh_s_107 = buffer.data(fh_s + 107);
    const auto *fh_s_108 = buffer.data(fh_s + 108);
    const auto *fh_s_109 = buffer.data(fh_s + 109);
    const auto *fh_s_110 = buffer.data(fh_s + 110);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_8 = buffer.data(ff + 8);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_17 = buffer.data(ff + 17);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_28 = buffer.data(ff + 28);
    const auto *ff_29 = buffer.data(ff + 29);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_13 = buffer.data(fg + 13);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_26 = buffer.data(fg + 26);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_43 = buffer.data(fg + 43);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_77 = buffer.data(fg + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, dg_0, ff_s_0, fh_s_0, fh_s_1, \
                         fh_s_2, fh_s_3, ff_0, fg_0, fg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dg_0[k]
                 - f_1 * ff_s_0[k]
                 + f_2 * fh_s_0[k]
                 + f_3 * ff_0[k]
                 + pb_x[k] * fg_0[k];

        t_1[k] = f_2 * fh_s_1[k]
                 + pb_y[k] * fg_0[k];

        t_2[k] = f_2 * fh_s_2[k]
                 + pb_z[k] * fg_0[k];

        t_3[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_3[k]
                 + f_5 * ff_0[k]
                 + pb_y[k] * fg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, ff_s_0, ff_s_1, fh_s_4, fh_s_5, \
                         fh_s_6, fh_s_7, ff_0, ff_1, fg_2, fg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * fh_s_4[k]
                 + pb_y[k] * fg_2[k];

        t_5[k] = -f_4 * ff_s_0[k]
                 + f_2 * fh_s_5[k]
                 + f_5 * ff_0[k]
                 + pb_z[k] * fg_2[k];

        t_6[k] = -f_6 * ff_s_1[k]
                 + f_2 * fh_s_6[k]
                 + f_7 * ff_1[k]
                 + pb_y[k] * fg_3[k];

        t_7[k] = f_2 * fh_s_7[k]
                 + pb_z[k] * fg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, dg_10, ff_s_2, fh_s_8, fh_s_9, \
                         fh_s_10, ff_2, fg_5, fg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * fh_s_8[k]
                 + pb_y[k] * fg_5[k];

        t_9[k] = -f_6 * ff_s_2[k]
                 + f_2 * fh_s_9[k]
                 + f_7 * ff_2[k]
                 + pb_z[k] * fg_5[k];

        t_10[k] = f_0 * dg_10[k]
                  + f_2 * fh_s_10[k]
                  + pb_x[k] * fg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, dg_12, fh_s_11, fh_s_12, fh_s_13, \
                         fg_6, fg_9, fg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * fh_s_11[k]
                  + pb_z[k] * fg_6[k];

        t_12[k] = f_0 * dg_12[k]
                  + f_2 * fh_s_12[k]
                  + pb_x[k] * fg_12[k];

        t_13[k] = f_2 * fh_s_13[k]
                  + pb_y[k] * fg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, dg_14, ff_s_6, fh_s_14, fh_s_15, \
                         fh_s_16, ff_6, fg_10, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * dg_14[k]
                  + f_2 * fh_s_14[k]
                  + pb_x[k] * fg_14[k];

        t_15[k] = -f_1 * ff_s_6[k]
                  + f_2 * fh_s_15[k]
                  + f_3 * ff_6[k]
                  + pb_y[k] * fg_10[k];

        t_16[k] = f_2 * fh_s_16[k]
                  + pb_z[k] * fg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, ff_s_8, ff_s_9, fh_s_17, fh_s_18, fh_s_19, \
                         ff_8, ff_9, fg_12, fg_13, fg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_6 * ff_s_8[k]
                  + f_2 * fh_s_17[k]
                  + f_7 * ff_8[k]
                  + pb_y[k] * fg_12[k];

        t_18[k] = -f_4 * ff_s_9[k]
                  + f_2 * fh_s_18[k]
                  + f_5 * ff_9[k]
                  + pb_y[k] * fg_13[k];

        t_19[k] = f_2 * fh_s_19[k]
                  + pb_y[k] * fg_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, dg_0, dh_0, ff_s_9, fh_s_20, \
                         fh_s_21, fh_s_22, ff_9, fg_14, fg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * ff_s_9[k]
                  + f_2 * fh_s_20[k]
                  + f_3 * ff_9[k]
                  + pb_z[k] * fg_14[k];

        t_21[k] = pa_y[k] * dh_0[k]
                  + f_2 * fh_s_21[k];

        t_22[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_22[k]
                  + pb_y[k] * fg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_z, dg_1, dh_3, dh_5, fh_s_23, \
                         fh_s_24, fh_s_25, fh_s_26, fg_15, fg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * fh_s_23[k]
                  + pb_z[k] * fg_15[k];

        t_24[k] = f_7 * dg_1[k]
                  + pa_y[k] * dh_3[k]
                  + f_2 * fh_s_24[k];

        t_25[k] = f_2 * fh_s_25[k]
                  + pb_z[k] * fg_16[k];

        t_26[k] = pa_y[k] * dh_5[k]
                  + f_2 * fh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, dg_3, dg_5, dh_6, fh_s_27, \
                         fh_s_28, fh_s_29, fg_18, fg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * dg_3[k]
                  + pa_y[k] * dh_6[k]
                  + f_2 * fh_s_27[k];

        t_28[k] = f_2 * fh_s_28[k]
                  + pb_z[k] * fg_18[k];

        t_29[k] = f_5 * dg_5[k]
                  + f_2 * fh_s_29[k]
                  + pb_y[k] * fg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, dg_25, dh_9, fh_s_30, fh_s_31, \
                         fh_s_32, fg_21, fg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * dh_9[k]
                  + f_2 * fh_s_30[k];

        t_31[k] = f_7 * dg_25[k]
                  + f_2 * fh_s_31[k]
                  + pb_x[k] * fg_25[k];

        t_32[k] = f_2 * fh_s_32[k]
                  + pb_z[k] * fg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, dg_27, dg_28, dh_14, fh_s_33, fh_s_34, \
                         fh_s_35, fg_27, fg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * dg_27[k]
                  + f_2 * fh_s_33[k]
                  + pb_x[k] * fg_27[k];

        t_34[k] = f_7 * dg_28[k]
                  + f_2 * fh_s_34[k]
                  + pb_x[k] * fg_28[k];

        t_35[k] = pa_y[k] * dh_14[k]
                  + f_2 * fh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_z, ph_s_36, ph_36, dh_36, ff_s_16, \
                         fh_s_36, fh_s_37, fh_s_38, ff_16, fg_25, \
                         fg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * ph_s_36[k]
                  + f_5 * ph_36[k]
                  + pa_x[k] * dh_36[k]
                  + f_2 * fh_s_36[k];

        t_37[k] = f_2 * fh_s_37[k]
                  + pb_z[k] * fg_25[k];

        t_38[k] = -f_4 * ff_s_16[k]
                  + f_2 * fh_s_38[k]
                  + f_5 * ff_16[k]
                  + pb_z[k] * fg_26[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, pb_z, dg_14, dh_20, ff_s_17, fh_s_39, \
                         fh_s_40, fh_s_41, ff_17, fg_27, fg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_6 * ff_s_17[k]
                  + f_2 * fh_s_39[k]
                  + f_7 * ff_17[k]
                  + pb_z[k] * fg_27[k];

        t_40[k] = f_5 * dg_14[k]
                  + f_2 * fh_s_40[k]
                  + pb_y[k] * fg_29[k];

        t_41[k] = pa_y[k] * dh_20[k]
                  + f_2 * fh_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, dg_0, dh_0, dh_3, fh_s_42, \
                         fh_s_43, fh_s_44, fh_s_45, fg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * dh_0[k]
                  + f_2 * fh_s_42[k];

        t_43[k] = f_2 * fh_s_43[k]
                  + pb_y[k] * fg_30[k];

        t_44[k] = f_5 * dg_0[k]
                  + f_2 * fh_s_44[k]
                  + pb_z[k] * fg_30[k];

        t_45[k] = pa_z[k] * dh_3[k]
                  + f_2 * fh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_y, dg_2, dg_3, dh_5, dh_6, dh_7, \
                         fh_s_46, fh_s_47, fh_s_48, fh_s_49, fg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * fh_s_46[k]
                  + pb_y[k] * fg_32[k];

        t_47[k] = f_7 * dg_2[k]
                  + pa_z[k] * dh_5[k]
                  + f_2 * fh_s_47[k];

        t_48[k] = pa_z[k] * dh_6[k]
                  + f_2 * fh_s_48[k];

        t_49[k] = f_5 * dg_3[k]
                  + pa_z[k] * dh_7[k]
                  + f_2 * fh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_z, pb_y, dg_5, dh_9, dh_10, fh_s_50, fh_s_51, \
                         fh_s_52, fg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * fh_s_50[k]
                  + pb_y[k] * fg_35[k];

        t_51[k] = f_0 * dg_5[k]
                  + pa_z[k] * dh_9[k]
                  + f_2 * fh_s_51[k];

        t_52[k] = pa_z[k] * dh_10[k]
                  + f_2 * fh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, dg_41, dg_42, fh_s_53, fh_s_54, \
                         fh_s_55, fg_39, fg_41, fg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * dg_41[k]
                  + f_2 * fh_s_53[k]
                  + pb_x[k] * fg_41[k];

        t_54[k] = f_7 * dg_42[k]
                  + f_2 * fh_s_54[k]
                  + pb_x[k] * fg_42[k];

        t_55[k] = f_2 * fh_s_55[k]
                  + pb_y[k] * fg_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_z, pb_x, pb_y, dg_44, dh_15, ff_s_27, fh_s_56, \
                         fh_s_57, fh_s_58, ff_27, fg_41, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_7 * dg_44[k]
                  + f_2 * fh_s_56[k]
                  + pb_x[k] * fg_44[k];

        t_57[k] = pa_z[k] * dh_15[k]
                  + f_2 * fh_s_57[k];

        t_58[k] = -f_9 * ff_s_27[k]
                  + f_2 * fh_s_58[k]
                  + f_0 * ff_27[k]
                  + pb_y[k] * fg_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, ff_s_28, ff_s_29, fh_s_59, fh_s_60, fh_s_61, \
                         ff_28, ff_29, fg_42, fg_43, fg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_6 * ff_s_28[k]
                  + f_2 * fh_s_59[k]
                  + f_7 * ff_28[k]
                  + pb_y[k] * fg_42[k];

        t_60[k] = -f_4 * ff_s_29[k]
                  + f_2 * fh_s_60[k]
                  + f_5 * ff_29[k]
                  + pb_y[k] * fg_43[k];

        t_61[k] = f_2 * fh_s_61[k]
                  + pb_y[k] * fg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_x, pb_y, ph_s_62, ph_62, dg_15, dg_45, dh_62, \
                         dh_63, fh_s_62, fh_s_63, fh_s_64, fg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -f_8 * ph_s_62[k]
                  + f_5 * ph_62[k]
                  + pa_x[k] * dh_62[k]
                  + f_2 * fh_s_62[k];

        t_63[k] = f_10 * dg_45[k]
                  + pa_x[k] * dh_63[k]
                  + f_2 * fh_s_63[k];

        t_64[k] = f_7 * dg_15[k]
                  + f_2 * fh_s_64[k]
                  + pb_y[k] * fg_45[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_z, dg_48, dg_50, dh_66, dh_68, \
                         fh_s_65, fh_s_66, fh_s_67, fh_s_68, fg_45, \
                         fg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * fh_s_65[k]
                  + pb_z[k] * fg_45[k];

        t_66[k] = f_0 * dg_48[k]
                  + pa_x[k] * dh_66[k]
                  + f_2 * fh_s_66[k];

        t_67[k] = f_2 * fh_s_67[k]
                  + pb_z[k] * fg_46[k];

        t_68[k] = f_0 * dg_50[k]
                  + pa_x[k] * dh_68[k]
                  + f_2 * fh_s_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_x, pb_y, pb_z, dg_20, dg_51, dh_69, fh_s_69, \
                         fh_s_70, fh_s_71, fg_48, fg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_7 * dg_51[k]
                  + pa_x[k] * dh_69[k]
                  + f_2 * fh_s_69[k];

        t_70[k] = f_2 * fh_s_70[k]
                  + pb_z[k] * fg_48[k];

        t_71[k] = f_7 * dg_20[k]
                  + f_2 * fh_s_71[k]
                  + pb_y[k] * fg_50[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_x, pb_x, pb_z, dg_54, dg_55, dh_72, fh_s_72, \
                         fh_s_73, fh_s_74, fg_51, fg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * dg_54[k]
                  + pa_x[k] * dh_72[k]
                  + f_2 * fh_s_72[k];

        t_73[k] = f_5 * dg_55[k]
                  + f_2 * fh_s_73[k]
                  + pb_x[k] * fg_55[k];

        t_74[k] = f_2 * fh_s_74[k]
                  + pb_z[k] * fg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_x, dg_57, dg_58, dg_59, fh_s_75, fh_s_76, \
                         fh_s_77, fg_57, fg_58, fg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_5 * dg_57[k]
                  + f_2 * fh_s_75[k]
                  + pb_x[k] * fg_57[k];

        t_76[k] = f_5 * dg_58[k]
                  + f_2 * fh_s_76[k]
                  + pb_x[k] * fg_58[k];

        t_77[k] = f_5 * dg_59[k]
                  + f_2 * fh_s_77[k]
                  + pb_x[k] * fg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, dh_78, dh_80, dh_81, fh_s_78, \
                         fh_s_79, fh_s_80, fh_s_81, fg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * dh_78[k]
                  + f_2 * fh_s_78[k];

        t_79[k] = f_2 * fh_s_79[k]
                  + pb_z[k] * fg_55[k];

        t_80[k] = pa_x[k] * dh_80[k]
                  + f_2 * fh_s_80[k];

        t_81[k] = pa_x[k] * dh_81[k]
                  + f_2 * fh_s_81[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, pa_z, dh_22, dh_42, dh_82, dh_83, \
                         fh_s_82, fh_s_83, fh_s_84, fh_s_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_x[k] * dh_82[k]
                  + f_2 * fh_s_82[k];

        t_83[k] = pa_x[k] * dh_83[k]
                  + f_2 * fh_s_83[k];

        t_84[k] = pa_y[k] * dh_42[k]
                  + f_2 * fh_s_84[k];

        t_85[k] = pa_z[k] * dh_22[k]
                  + f_2 * fh_s_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pa_z, pb_y, dg_32, dh_24, dh_44, dh_47, \
                         fh_s_86, fh_s_87, fh_s_88, fh_s_89, fg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_y[k] * dh_44[k]
                  + f_2 * fh_s_86[k];

        t_87[k] = pa_z[k] * dh_24[k]
                  + f_2 * fh_s_87[k];

        t_88[k] = f_5 * dg_32[k]
                  + f_2 * fh_s_88[k]
                  + pb_y[k] * fg_62[k];

        t_89[k] = pa_y[k] * dh_47[k]
                  + f_2 * fh_s_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pa_z, pb_y, pb_z, dg_18, dg_35, dh_27, fh_s_90, \
                         fh_s_91, fh_s_92, fg_63, fg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_z[k] * dh_27[k]
                  + f_2 * fh_s_90[k];

        t_91[k] = f_5 * dg_18[k]
                  + f_2 * fh_s_91[k]
                  + pb_z[k] * fg_63[k];

        t_92[k] = f_5 * dg_35[k]
                  + f_2 * fh_s_92[k]
                  + pb_y[k] * fg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_y, pa_z, pb_x, dg_71, dh_31, dh_51, fh_s_93, \
                         fh_s_94, fh_s_95, fg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * dh_51[k]
                  + f_2 * fh_s_93[k];

        t_94[k] = pa_z[k] * dh_31[k]
                  + f_2 * fh_s_94[k];

        t_95[k] = f_5 * dg_71[k]
                  + f_2 * fh_s_95[k]
                  + pb_x[k] * fg_71[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_y, pb_x, dg_72, dg_73, dh_56, fh_s_96, fh_s_97, \
                         fh_s_98, fg_72, fg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * dg_72[k]
                  + f_2 * fh_s_96[k]
                  + pb_x[k] * fg_72[k];

        t_97[k] = f_5 * dg_73[k]
                  + f_2 * fh_s_97[k]
                  + pb_x[k] * fg_73[k];

        t_98[k] = pa_y[k] * dh_56[k]
                  + f_2 * fh_s_98[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_x, dh_99, dh_100, dh_101, \
                         dh_102, dh_103, fh_s_99, fh_s_100, fh_s_101, fh_s_102, \
                         fh_s_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * dh_99[k]
                  + f_2 * fh_s_99[k];

        t_100[k] = pa_x[k] * dh_100[k]
                   + f_2 * fh_s_100[k];

        t_101[k] = pa_x[k] * dh_101[k]
                   + f_2 * fh_s_101[k];

        t_102[k] = pa_x[k] * dh_102[k]
                   + f_2 * fh_s_102[k];

        t_103[k] = pa_x[k] * dh_103[k]
                   + f_2 * fh_s_103[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_y, pb_z, dg_30, dg_75, dh_104, \
                         dh_105, fh_s_104, fh_s_105, fh_s_106, fh_s_107, \
                         fg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_x[k] * dh_104[k]
                   + f_2 * fh_s_104[k];

        t_105[k] = f_10 * dg_75[k]
                   + pa_x[k] * dh_105[k]
                   + f_2 * fh_s_105[k];

        t_106[k] = f_2 * fh_s_106[k]
                   + pb_y[k] * fg_75[k];

        t_107[k] = f_7 * dg_30[k]
                   + f_2 * fh_s_107[k]
                   + pb_z[k] * fg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_x, pb_y, dg_78, dg_80, dh_108, dh_110, \
                         fh_s_108, fh_s_109, fh_s_110, fg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * dg_78[k]
                   + pa_x[k] * dh_108[k]
                   + f_2 * fh_s_108[k];

        t_109[k] = f_2 * fh_s_109[k]
                   + pb_y[k] * fg_77[k];

        t_110[k] = f_0 * dg_80[k]
                   + pa_x[k] * dh_110[k]
                   + f_2 * fh_s_110[k];
    }
}

static auto
compute_prim_fh_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t ph_s,
                                        const size_t ph, const size_t dg, const size_t dh,
                                        const size_t ff_s, const size_t fh_s, const size_t ff,
                                        const size_t fg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_s_62 = buffer.data(ph_s + 62);

    const auto *ph_62 = buffer.data(ph + 62);

    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_76 = buffer.data(dg + 76);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_79 = buffer.data(dg + 79);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_67 = buffer.data(dh + 67);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_70 = buffer.data(dh + 70);
    const auto *dh_71 = buffer.data(dh + 71);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_106 = buffer.data(dh + 106);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_109 = buffer.data(dh + 109);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_112 = buffer.data(dh + 112);
    const auto *dh_113 = buffer.data(dh + 113);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *ff_s_60 = buffer.data(ff_s + 60);
    const auto *ff_s_61 = buffer.data(ff_s + 61);
    const auto *ff_s_63 = buffer.data(ff_s + 63);
    const auto *ff_s_65 = buffer.data(ff_s + 65);
    const auto *ff_s_66 = buffer.data(ff_s + 66);
    const auto *ff_s_67 = buffer.data(ff_s + 67);
    const auto *ff_s_68 = buffer.data(ff_s + 68);
    const auto *ff_s_69 = buffer.data(ff_s + 69);
    const auto *ff_s_72 = buffer.data(ff_s + 72);
    const auto *ff_s_75 = buffer.data(ff_s + 75);
    const auto *ff_s_79 = buffer.data(ff_s + 79);
    const auto *ff_s_90 = buffer.data(ff_s + 90);
    const auto *ff_s_92 = buffer.data(ff_s + 92);
    const auto *ff_s_93 = buffer.data(ff_s + 93);
    const auto *ff_s_95 = buffer.data(ff_s + 95);
    const auto *ff_s_96 = buffer.data(ff_s + 96);
    const auto *ff_s_97 = buffer.data(ff_s + 97);
    const auto *ff_s_98 = buffer.data(ff_s + 98);
    const auto *ff_s_99 = buffer.data(ff_s + 99);

    const auto *fh_s_111 = buffer.data(fh_s + 111);
    const auto *fh_s_112 = buffer.data(fh_s + 112);
    const auto *fh_s_113 = buffer.data(fh_s + 113);
    const auto *fh_s_114 = buffer.data(fh_s + 114);
    const auto *fh_s_115 = buffer.data(fh_s + 115);
    const auto *fh_s_116 = buffer.data(fh_s + 116);
    const auto *fh_s_117 = buffer.data(fh_s + 117);
    const auto *fh_s_118 = buffer.data(fh_s + 118);
    const auto *fh_s_119 = buffer.data(fh_s + 119);
    const auto *fh_s_120 = buffer.data(fh_s + 120);
    const auto *fh_s_121 = buffer.data(fh_s + 121);
    const auto *fh_s_122 = buffer.data(fh_s + 122);
    const auto *fh_s_123 = buffer.data(fh_s + 123);
    const auto *fh_s_124 = buffer.data(fh_s + 124);
    const auto *fh_s_125 = buffer.data(fh_s + 125);
    const auto *fh_s_126 = buffer.data(fh_s + 126);
    const auto *fh_s_127 = buffer.data(fh_s + 127);
    const auto *fh_s_128 = buffer.data(fh_s + 128);
    const auto *fh_s_129 = buffer.data(fh_s + 129);
    const auto *fh_s_130 = buffer.data(fh_s + 130);
    const auto *fh_s_131 = buffer.data(fh_s + 131);
    const auto *fh_s_132 = buffer.data(fh_s + 132);
    const auto *fh_s_133 = buffer.data(fh_s + 133);
    const auto *fh_s_134 = buffer.data(fh_s + 134);
    const auto *fh_s_135 = buffer.data(fh_s + 135);
    const auto *fh_s_136 = buffer.data(fh_s + 136);
    const auto *fh_s_137 = buffer.data(fh_s + 137);
    const auto *fh_s_138 = buffer.data(fh_s + 138);
    const auto *fh_s_139 = buffer.data(fh_s + 139);
    const auto *fh_s_140 = buffer.data(fh_s + 140);
    const auto *fh_s_141 = buffer.data(fh_s + 141);
    const auto *fh_s_142 = buffer.data(fh_s + 142);
    const auto *fh_s_143 = buffer.data(fh_s + 143);
    const auto *fh_s_144 = buffer.data(fh_s + 144);
    const auto *fh_s_145 = buffer.data(fh_s + 145);
    const auto *fh_s_146 = buffer.data(fh_s + 146);
    const auto *fh_s_147 = buffer.data(fh_s + 147);
    const auto *fh_s_148 = buffer.data(fh_s + 148);
    const auto *fh_s_149 = buffer.data(fh_s + 149);
    const auto *fh_s_150 = buffer.data(fh_s + 150);
    const auto *fh_s_151 = buffer.data(fh_s + 151);
    const auto *fh_s_152 = buffer.data(fh_s + 152);
    const auto *fh_s_153 = buffer.data(fh_s + 153);
    const auto *fh_s_154 = buffer.data(fh_s + 154);
    const auto *fh_s_155 = buffer.data(fh_s + 155);
    const auto *fh_s_156 = buffer.data(fh_s + 156);
    const auto *fh_s_157 = buffer.data(fh_s + 157);
    const auto *fh_s_158 = buffer.data(fh_s + 158);
    const auto *fh_s_159 = buffer.data(fh_s + 159);
    const auto *fh_s_160 = buffer.data(fh_s + 160);
    const auto *fh_s_161 = buffer.data(fh_s + 161);
    const auto *fh_s_162 = buffer.data(fh_s + 162);
    const auto *fh_s_163 = buffer.data(fh_s + 163);
    const auto *fh_s_164 = buffer.data(fh_s + 164);
    const auto *fh_s_165 = buffer.data(fh_s + 165);
    const auto *fh_s_166 = buffer.data(fh_s + 166);
    const auto *fh_s_167 = buffer.data(fh_s + 167);
    const auto *fh_s_168 = buffer.data(fh_s + 168);
    const auto *fh_s_169 = buffer.data(fh_s + 169);
    const auto *fh_s_170 = buffer.data(fh_s + 170);
    const auto *fh_s_171 = buffer.data(fh_s + 171);
    const auto *fh_s_172 = buffer.data(fh_s + 172);
    const auto *fh_s_173 = buffer.data(fh_s + 173);
    const auto *fh_s_174 = buffer.data(fh_s + 174);
    const auto *fh_s_175 = buffer.data(fh_s + 175);
    const auto *fh_s_176 = buffer.data(fh_s + 176);
    const auto *fh_s_177 = buffer.data(fh_s + 177);
    const auto *fh_s_178 = buffer.data(fh_s + 178);
    const auto *fh_s_179 = buffer.data(fh_s + 179);
    const auto *fh_s_180 = buffer.data(fh_s + 180);
    const auto *fh_s_181 = buffer.data(fh_s + 181);
    const auto *fh_s_182 = buffer.data(fh_s + 182);
    const auto *fh_s_183 = buffer.data(fh_s + 183);
    const auto *fh_s_184 = buffer.data(fh_s + 184);
    const auto *fh_s_185 = buffer.data(fh_s + 185);
    const auto *fh_s_186 = buffer.data(fh_s + 186);
    const auto *fh_s_187 = buffer.data(fh_s + 187);
    const auto *fh_s_188 = buffer.data(fh_s + 188);
    const auto *fh_s_189 = buffer.data(fh_s + 189);
    const auto *fh_s_190 = buffer.data(fh_s + 190);
    const auto *fh_s_191 = buffer.data(fh_s + 191);
    const auto *fh_s_192 = buffer.data(fh_s + 192);
    const auto *fh_s_193 = buffer.data(fh_s + 193);
    const auto *fh_s_194 = buffer.data(fh_s + 194);
    const auto *fh_s_195 = buffer.data(fh_s + 195);
    const auto *fh_s_196 = buffer.data(fh_s + 196);
    const auto *fh_s_197 = buffer.data(fh_s + 197);
    const auto *fh_s_198 = buffer.data(fh_s + 198);
    const auto *fh_s_199 = buffer.data(fh_s + 199);
    const auto *fh_s_200 = buffer.data(fh_s + 200);
    const auto *fh_s_201 = buffer.data(fh_s + 201);
    const auto *fh_s_202 = buffer.data(fh_s + 202);
    const auto *fh_s_203 = buffer.data(fh_s + 203);
    const auto *fh_s_204 = buffer.data(fh_s + 204);
    const auto *fh_s_205 = buffer.data(fh_s + 205);
    const auto *fh_s_206 = buffer.data(fh_s + 206);
    const auto *fh_s_207 = buffer.data(fh_s + 207);
    const auto *fh_s_208 = buffer.data(fh_s + 208);
    const auto *fh_s_209 = buffer.data(fh_s + 209);

    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_72 = buffer.data(ff + 72);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_89 = buffer.data(fg + 89);
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
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
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

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pb_y, dg_81, dg_82, dh_111, dh_112, \
                         fh_s_111, fh_s_112, fh_s_113, fg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_7 * dg_81[k]
                   + pa_x[k] * dh_111[k]
                   + f_2 * fh_s_111[k];

        t_112[k] = f_7 * dg_82[k]
                   + pa_x[k] * dh_112[k]
                   + f_2 * fh_s_112[k];

        t_113[k] = f_2 * fh_s_113[k]
                   + pb_y[k] * fg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_x, pb_x, dg_84, dg_85, dg_86, dh_114, \
                         fh_s_114, fh_s_115, fh_s_116, fg_85, fg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * dg_84[k]
                   + pa_x[k] * dh_114[k]
                   + f_2 * fh_s_114[k];

        t_115[k] = f_5 * dg_85[k]
                   + f_2 * fh_s_115[k]
                   + pb_x[k] * fg_85[k];

        t_116[k] = f_5 * dg_86[k]
                   + f_2 * fh_s_116[k]
                   + pb_x[k] * fg_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_x, pb_y, dg_87, dg_89, fh_s_117, fh_s_118, \
                         fh_s_119, fg_84, fg_87, fg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * dg_87[k]
                   + f_2 * fh_s_117[k]
                   + pb_x[k] * fg_87[k];

        t_118[k] = f_2 * fh_s_118[k]
                   + pb_y[k] * fg_84[k];

        t_119[k] = f_5 * dg_89[k]
                   + f_2 * fh_s_119[k]
                   + pb_x[k] * fg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, dh_120, dh_121, dh_122, dh_123, \
                         fh_s_120, fh_s_121, fh_s_122, fh_s_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_x[k] * dh_120[k]
                   + f_2 * fh_s_120[k];

        t_121[k] = pa_x[k] * dh_121[k]
                   + f_2 * fh_s_121[k];

        t_122[k] = pa_x[k] * dh_122[k]
                   + f_2 * fh_s_122[k];

        t_123[k] = pa_x[k] * dh_123[k]
                   + f_2 * fh_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pb_x, pb_y, dh_125, ff_s_60, fh_s_124, \
                         fh_s_125, fh_s_126, ff_60, fg_89, fg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * fh_s_124[k]
                   + pb_y[k] * fg_89[k];

        t_125[k] = pa_x[k] * dh_125[k]
                   + f_2 * fh_s_125[k];

        t_126[k] = -f_1 * ff_s_60[k]
                   + f_2 * fh_s_126[k]
                   + f_3 * ff_60[k]
                   + pb_x[k] * fg_90[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pb_x, pb_z, ff_s_61, ff_s_63, fh_s_127, \
                         fh_s_128, fh_s_129, ff_61, ff_63, fg_90, fg_91, \
                         fg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = -f_9 * ff_s_61[k]
                   + f_2 * fh_s_127[k]
                   + f_0 * ff_61[k]
                   + pb_x[k] * fg_91[k];

        t_128[k] = f_2 * fh_s_128[k]
                   + pb_z[k] * fg_90[k];

        t_129[k] = -f_6 * ff_s_63[k]
                   + f_2 * fh_s_129[k]
                   + f_7 * ff_63[k]
                   + pb_x[k] * fg_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pb_x, pb_z, ff_s_65, ff_s_66, fh_s_130, \
                         fh_s_131, fh_s_132, ff_65, ff_66, fg_91, fg_95, \
                         fg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * fh_s_130[k]
                   + pb_z[k] * fg_91[k];

        t_131[k] = -f_6 * ff_s_65[k]
                   + f_2 * fh_s_131[k]
                   + f_7 * ff_65[k]
                   + pb_x[k] * fg_95[k];

        t_132[k] = -f_4 * ff_s_66[k]
                   + f_2 * fh_s_132[k]
                   + f_5 * ff_66[k]
                   + pb_x[k] * fg_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_x, pb_z, ff_s_68, ff_s_69, fh_s_133, \
                         fh_s_134, fh_s_135, ff_68, ff_69, fg_93, fg_98, \
                         fg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_2 * fh_s_133[k]
                   + pb_z[k] * fg_93[k];

        t_134[k] = -f_4 * ff_s_68[k]
                   + f_2 * fh_s_134[k]
                   + f_5 * ff_68[k]
                   + pb_x[k] * fg_98[k];

        t_135[k] = -f_4 * ff_s_69[k]
                   + f_2 * fh_s_135[k]
                   + f_5 * ff_69[k]
                   + pb_x[k] * fg_99[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, fh_s_136, fh_s_137, \
                         fh_s_138, fh_s_139, fh_s_140, fg_100, fg_101, fg_102, fg_103, \
                         fg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * fh_s_136[k]
                   + pb_x[k] * fg_100[k];

        t_137[k] = f_2 * fh_s_137[k]
                   + pb_x[k] * fg_101[k];

        t_138[k] = f_2 * fh_s_138[k]
                   + pb_x[k] * fg_102[k];

        t_139[k] = f_2 * fh_s_139[k]
                   + pb_x[k] * fg_103[k];

        t_140[k] = f_2 * fh_s_140[k]
                   + pb_x[k] * fg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pb_y, pb_z, dg_55, ff_s_66, fh_s_141, fh_s_142, \
                         fh_s_143, ff_66, fg_100, fg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_0 * dg_55[k]
                   - f_1 * ff_s_66[k]
                   + f_2 * fh_s_141[k]
                   + f_3 * ff_66[k]
                   + pb_y[k] * fg_100[k];

        t_142[k] = f_2 * fh_s_142[k]
                   + pb_z[k] * fg_100[k];

        t_143[k] = -f_4 * ff_s_66[k]
                   + f_2 * fh_s_143[k]
                   + f_5 * ff_66[k]
                   + pb_z[k] * fg_101[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pb_y, pb_z, dg_59, ff_s_67, ff_s_69, fh_s_144, \
                         fh_s_145, fh_s_146, ff_67, ff_69, fg_102, \
                         fg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = -f_6 * ff_s_67[k]
                   + f_2 * fh_s_144[k]
                   + f_7 * ff_67[k]
                   + pb_z[k] * fg_102[k];

        t_145[k] = f_0 * dg_59[k]
                   + f_2 * fh_s_145[k]
                   + pb_y[k] * fg_104[k];

        t_146[k] = -f_1 * ff_s_69[k]
                   + f_2 * fh_s_146[k]
                   + f_3 * ff_69[k]
                   + pb_z[k] * fg_104[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pb_x, dh_63, dh_64, dh_66, ff_s_72, \
                         fh_s_147, fh_s_148, fh_s_149, fh_s_150, ff_72, \
                         fg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pa_z[k] * dh_63[k]
                   + f_2 * fh_s_147[k];

        t_148[k] = pa_z[k] * dh_64[k]
                   + f_2 * fh_s_148[k];

        t_149[k] = -f_9 * ff_s_72[k]
                   + f_2 * fh_s_149[k]
                   + f_0 * ff_72[k]
                   + pb_x[k] * fg_107[k];

        t_150[k] = pa_z[k] * dh_66[k]
                   + f_2 * fh_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_z, pb_x, dg_46, dh_67, dh_69, ff_s_75, \
                         fh_s_151, fh_s_152, fh_s_153, ff_75, fg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_5 * dg_46[k]
                   + pa_z[k] * dh_67[k]
                   + f_2 * fh_s_151[k];

        t_152[k] = -f_6 * ff_s_75[k]
                   + f_2 * fh_s_152[k]
                   + f_7 * ff_75[k]
                   + pb_x[k] * fg_110[k];

        t_153[k] = pa_z[k] * dh_69[k]
                   + f_2 * fh_s_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_z, pb_x, dg_48, dg_49, dh_70, dh_71, ff_s_79, \
                         fh_s_154, fh_s_155, fh_s_156, ff_79, fg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_5 * dg_48[k]
                   + pa_z[k] * dh_70[k]
                   + f_2 * fh_s_154[k];

        t_155[k] = f_7 * dg_49[k]
                   + pa_z[k] * dh_71[k]
                   + f_2 * fh_s_155[k];

        t_156[k] = -f_4 * ff_s_79[k]
                   + f_2 * fh_s_156[k]
                   + f_5 * ff_79[k]
                   + pb_x[k] * fg_114[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pb_x, fh_s_157, fh_s_158, \
                         fh_s_159, fh_s_160, fh_s_161, fg_115, fg_116, fg_117, fg_118, \
                         fg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_2 * fh_s_157[k]
                   + pb_x[k] * fg_115[k];

        t_158[k] = f_2 * fh_s_158[k]
                   + pb_x[k] * fg_116[k];

        t_159[k] = f_2 * fh_s_159[k]
                   + pb_x[k] * fg_117[k];

        t_160[k] = f_2 * fh_s_160[k]
                   + pb_x[k] * fg_118[k];

        t_161[k] = f_2 * fh_s_161[k]
                   + pb_x[k] * fg_119[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pa_z, pb_z, dg_55, dg_56, dh_78, dh_80, \
                         fh_s_162, fh_s_163, fh_s_164, fg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_z[k] * dh_78[k]
                   + f_2 * fh_s_162[k];

        t_163[k] = f_5 * dg_55[k]
                   + f_2 * fh_s_163[k]
                   + pb_z[k] * fg_115[k];

        t_164[k] = f_7 * dg_56[k]
                   + pa_z[k] * dh_80[k]
                   + f_2 * fh_s_164[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_y, pa_z, pb_y, ph_s_62, ph_62, dg_57, dg_74, \
                         dh_81, dh_104, fh_s_165, fh_s_166, fh_s_167, \
                         fg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_0 * dg_57[k]
                   + pa_z[k] * dh_81[k]
                   + f_2 * fh_s_165[k];

        t_166[k] = f_7 * dg_74[k]
                   + f_2 * fh_s_166[k]
                   + pb_y[k] * fg_119[k];

        t_167[k] = -f_8 * ph_s_62[k]
                   + f_5 * ph_62[k]
                   + pa_y[k] * dh_104[k]
                   + f_2 * fh_s_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, dg_75, dg_76, dh_105, dh_106, \
                         dh_107, dh_108, fh_s_168, fh_s_169, fh_s_170, \
                         fh_s_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * dh_105[k]
                   + f_2 * fh_s_168[k];

        t_169[k] = f_5 * dg_75[k]
                   + pa_y[k] * dh_106[k]
                   + f_2 * fh_s_169[k];

        t_170[k] = pa_y[k] * dh_107[k]
                   + f_2 * fh_s_170[k];

        t_171[k] = f_7 * dg_76[k]
                   + pa_y[k] * dh_108[k]
                   + f_2 * fh_s_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, dg_77, dg_78, dg_79, dh_109, \
                         dh_110, dh_111, dh_112, fh_s_172, fh_s_173, fh_s_174, \
                         fh_s_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * dg_77[k]
                   + pa_y[k] * dh_109[k]
                   + f_2 * fh_s_172[k];

        t_173[k] = pa_y[k] * dh_110[k]
                   + f_2 * fh_s_173[k];

        t_174[k] = f_0 * dg_78[k]
                   + pa_y[k] * dh_111[k]
                   + f_2 * fh_s_174[k];

        t_175[k] = f_7 * dg_79[k]
                   + pa_y[k] * dh_112[k]
                   + f_2 * fh_s_175[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, dg_80, dh_113, dh_114, \
                         fh_s_176, fh_s_177, fh_s_178, fh_s_179, fg_130, \
                         fg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_5 * dg_80[k]
                   + pa_y[k] * dh_113[k]
                   + f_2 * fh_s_176[k];

        t_177[k] = pa_y[k] * dh_114[k]
                   + f_2 * fh_s_177[k];

        t_178[k] = f_2 * fh_s_178[k]
                   + pb_x[k] * fg_130[k];

        t_179[k] = f_2 * fh_s_179[k]
                   + pb_x[k] * fg_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pb_x, dg_85, dh_120, fh_s_180, \
                         fh_s_181, fh_s_182, fh_s_183, fg_132, fg_133, \
                         fg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_2 * fh_s_180[k]
                   + pb_x[k] * fg_132[k];

        t_181[k] = f_2 * fh_s_181[k]
                   + pb_x[k] * fg_133[k];

        t_182[k] = f_2 * fh_s_182[k]
                   + pb_x[k] * fg_134[k];

        t_183[k] = f_10 * dg_85[k]
                   + pa_y[k] * dh_120[k]
                   + f_2 * fh_s_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_y, pb_z, dg_70, dg_87, dg_88, dh_122, dh_123, \
                         fh_s_184, fh_s_185, fh_s_186, fg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_7 * dg_70[k]
                   + f_2 * fh_s_184[k]
                   + pb_z[k] * fg_130[k];

        t_185[k] = f_0 * dg_87[k]
                   + pa_y[k] * dh_122[k]
                   + f_2 * fh_s_185[k];

        t_186[k] = f_7 * dg_88[k]
                   + pa_y[k] * dh_123[k]
                   + f_2 * fh_s_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pb_x, pb_y, dg_89, dh_125, ff_s_90, \
                         fh_s_187, fh_s_188, fh_s_189, ff_90, fg_134, \
                         fg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_5 * dg_89[k]
                   + f_2 * fh_s_187[k]
                   + pb_y[k] * fg_134[k];

        t_188[k] = pa_y[k] * dh_125[k]
                   + f_2 * fh_s_188[k];

        t_189[k] = -f_1 * ff_s_90[k]
                   + f_2 * fh_s_189[k]
                   + f_3 * ff_90[k]
                   + pb_x[k] * fg_135[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pb_x, pb_y, ff_s_92, ff_s_93, fh_s_190, \
                         fh_s_191, fh_s_192, ff_92, ff_93, fg_135, fg_137, \
                         fg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_2 * fh_s_190[k]
                   + pb_y[k] * fg_135[k];

        t_191[k] = -f_9 * ff_s_92[k]
                   + f_2 * fh_s_191[k]
                   + f_0 * ff_92[k]
                   + pb_x[k] * fg_137[k];

        t_192[k] = -f_6 * ff_s_93[k]
                   + f_2 * fh_s_192[k]
                   + f_7 * ff_93[k]
                   + pb_x[k] * fg_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pb_x, pb_y, ff_s_95, ff_s_96, fh_s_193, \
                         fh_s_194, fh_s_195, ff_95, ff_96, fg_137, fg_140, \
                         fg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_2 * fh_s_193[k]
                   + pb_y[k] * fg_137[k];

        t_194[k] = -f_6 * ff_s_95[k]
                   + f_2 * fh_s_194[k]
                   + f_7 * ff_95[k]
                   + pb_x[k] * fg_140[k];

        t_195[k] = -f_4 * ff_s_96[k]
                   + f_2 * fh_s_195[k]
                   + f_5 * ff_96[k]
                   + pb_x[k] * fg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pb_x, pb_y, ff_s_97, ff_s_99, fh_s_196, \
                         fh_s_197, fh_s_198, ff_97, ff_99, fg_140, fg_142, \
                         fg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -f_4 * ff_s_97[k]
                   + f_2 * fh_s_196[k]
                   + f_5 * ff_97[k]
                   + pb_x[k] * fg_142[k];

        t_197[k] = f_2 * fh_s_197[k]
                   + pb_y[k] * fg_140[k];

        t_198[k] = -f_4 * ff_s_99[k]
                   + f_2 * fh_s_198[k]
                   + f_5 * ff_99[k]
                   + pb_x[k] * fg_144[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pb_x, fh_s_199, fh_s_200, \
                         fh_s_201, fh_s_202, fh_s_203, fg_145, fg_146, fg_147, fg_148, \
                         fg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_2 * fh_s_199[k]
                   + pb_x[k] * fg_145[k];

        t_200[k] = f_2 * fh_s_200[k]
                   + pb_x[k] * fg_146[k];

        t_201[k] = f_2 * fh_s_201[k]
                   + pb_x[k] * fg_147[k];

        t_202[k] = f_2 * fh_s_202[k]
                   + pb_x[k] * fg_148[k];

        t_203[k] = f_2 * fh_s_203[k]
                   + pb_x[k] * fg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pb_y, ff_s_96, ff_s_97, ff_s_98, fh_s_204, \
                         fh_s_205, fh_s_206, ff_96, ff_97, ff_98, fg_145, fg_146, \
                         fg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -f_1 * ff_s_96[k]
                   + f_2 * fh_s_204[k]
                   + f_3 * ff_96[k]
                   + pb_y[k] * fg_145[k];

        t_205[k] = -f_9 * ff_s_97[k]
                   + f_2 * fh_s_205[k]
                   + f_0 * ff_97[k]
                   + pb_y[k] * fg_146[k];

        t_206[k] = -f_6 * ff_s_98[k]
                   + f_2 * fh_s_206[k]
                   + f_7 * ff_98[k]
                   + pb_y[k] * fg_147[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pb_y, pb_z, dg_89, ff_s_99, fh_s_207, fh_s_208, \
                         fh_s_209, ff_99, fg_148, fg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -f_4 * ff_s_99[k]
                   + f_2 * fh_s_207[k]
                   + f_5 * ff_99[k]
                   + pb_y[k] * fg_148[k];

        t_208[k] = f_2 * fh_s_208[k]
                   + pb_y[k] * fg_149[k];

        t_209[k] = f_0 * dg_89[k]
                   - f_1 * ff_s_99[k]
                   + f_2 * fh_s_209[k]
                   + f_3 * ff_99[k]
                   + pb_z[k] * fg_149[k];
    }
}

auto
compute_prim_fh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ph_s, const size_t ph,
                                 const size_t dg, const size_t dh, const size_t ff_s,
                                 const size_t fh_s, const size_t ff, const size_t fg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_fh_kinetic_energy_0_piece0(buffer, target, pa, pb, ph_s, ph, dg, dh, ff_s,
                                            fh_s, ff, fg, ncols, alpha, beta, p);

    compute_prim_fh_kinetic_energy_0_piece1(buffer, target, pa, pb, ph_s, ph, dg, dh, ff_s,
                                            fh_s, ff, fg, ncols, alpha, beta, p);
}

}  // namespace simdkin
