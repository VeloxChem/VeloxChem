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


#include "SimdKineticEnergyVrrRecDH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_dh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);
    const auto *pg_16 = buffer.data(pg + 16);
    const auto *pg_17 = buffer.data(pg + 17);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_19 = buffer.data(pg + 19);

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

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_3 = buffer.data(df_s + 3);
    const auto *df_s_4 = buffer.data(df_s + 4);
    const auto *df_s_5 = buffer.data(df_s + 5);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_13 = buffer.data(df_s + 13);
    const auto *df_s_14 = buffer.data(df_s + 14);
    const auto *df_s_15 = buffer.data(df_s + 15);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_17 = buffer.data(df_s + 17);
    const auto *df_s_18 = buffer.data(df_s + 18);
    const auto *df_s_24 = buffer.data(df_s + 24);
    const auto *df_s_25 = buffer.data(df_s + 25);
    const auto *df_s_26 = buffer.data(df_s + 26);
    const auto *df_s_27 = buffer.data(df_s + 27);
    const auto *df_s_28 = buffer.data(df_s + 28);
    const auto *df_s_29 = buffer.data(df_s + 29);
    const auto *df_s_30 = buffer.data(df_s + 30);
    const auto *df_s_31 = buffer.data(df_s + 31);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);
    const auto *dh_s_25 = buffer.data(dh_s + 25);
    const auto *dh_s_26 = buffer.data(dh_s + 26);
    const auto *dh_s_27 = buffer.data(dh_s + 27);
    const auto *dh_s_28 = buffer.data(dh_s + 28);
    const auto *dh_s_29 = buffer.data(dh_s + 29);
    const auto *dh_s_30 = buffer.data(dh_s + 30);
    const auto *dh_s_31 = buffer.data(dh_s + 31);
    const auto *dh_s_32 = buffer.data(dh_s + 32);
    const auto *dh_s_33 = buffer.data(dh_s + 33);
    const auto *dh_s_34 = buffer.data(dh_s + 34);
    const auto *dh_s_35 = buffer.data(dh_s + 35);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_37 = buffer.data(dh_s + 37);
    const auto *dh_s_38 = buffer.data(dh_s + 38);
    const auto *dh_s_39 = buffer.data(dh_s + 39);
    const auto *dh_s_40 = buffer.data(dh_s + 40);
    const auto *dh_s_41 = buffer.data(dh_s + 41);
    const auto *dh_s_42 = buffer.data(dh_s + 42);
    const auto *dh_s_43 = buffer.data(dh_s + 43);
    const auto *dh_s_44 = buffer.data(dh_s + 44);
    const auto *dh_s_45 = buffer.data(dh_s + 45);
    const auto *dh_s_46 = buffer.data(dh_s + 46);
    const auto *dh_s_47 = buffer.data(dh_s + 47);
    const auto *dh_s_48 = buffer.data(dh_s + 48);
    const auto *dh_s_49 = buffer.data(dh_s + 49);
    const auto *dh_s_50 = buffer.data(dh_s + 50);
    const auto *dh_s_51 = buffer.data(dh_s + 51);
    const auto *dh_s_52 = buffer.data(dh_s + 52);
    const auto *dh_s_53 = buffer.data(dh_s + 53);
    const auto *dh_s_54 = buffer.data(dh_s + 54);
    const auto *dh_s_55 = buffer.data(dh_s + 55);
    const auto *dh_s_56 = buffer.data(dh_s + 56);
    const auto *dh_s_57 = buffer.data(dh_s + 57);
    const auto *dh_s_58 = buffer.data(dh_s + 58);
    const auto *dh_s_59 = buffer.data(dh_s + 59);
    const auto *dh_s_60 = buffer.data(dh_s + 60);
    const auto *dh_s_61 = buffer.data(dh_s + 61);
    const auto *dh_s_62 = buffer.data(dh_s + 62);
    const auto *dh_s_63 = buffer.data(dh_s + 63);
    const auto *dh_s_64 = buffer.data(dh_s + 64);
    const auto *dh_s_65 = buffer.data(dh_s + 65);
    const auto *dh_s_66 = buffer.data(dh_s + 66);
    const auto *dh_s_67 = buffer.data(dh_s + 67);
    const auto *dh_s_68 = buffer.data(dh_s + 68);
    const auto *dh_s_69 = buffer.data(dh_s + 69);
    const auto *dh_s_70 = buffer.data(dh_s + 70);
    const auto *dh_s_71 = buffer.data(dh_s + 71);
    const auto *dh_s_72 = buffer.data(dh_s + 72);
    const auto *dh_s_73 = buffer.data(dh_s + 73);
    const auto *dh_s_74 = buffer.data(dh_s + 74);
    const auto *dh_s_75 = buffer.data(dh_s + 75);
    const auto *dh_s_76 = buffer.data(dh_s + 76);
    const auto *dh_s_77 = buffer.data(dh_s + 77);
    const auto *dh_s_78 = buffer.data(dh_s + 78);
    const auto *dh_s_79 = buffer.data(dh_s + 79);
    const auto *dh_s_80 = buffer.data(dh_s + 80);
    const auto *dh_s_81 = buffer.data(dh_s + 81);
    const auto *dh_s_82 = buffer.data(dh_s + 82);
    const auto *dh_s_83 = buffer.data(dh_s + 83);
    const auto *dh_s_84 = buffer.data(dh_s + 84);
    const auto *dh_s_85 = buffer.data(dh_s + 85);
    const auto *dh_s_86 = buffer.data(dh_s + 86);
    const auto *dh_s_87 = buffer.data(dh_s + 87);
    const auto *dh_s_88 = buffer.data(dh_s + 88);
    const auto *dh_s_89 = buffer.data(dh_s + 89);
    const auto *dh_s_90 = buffer.data(dh_s + 90);
    const auto *dh_s_91 = buffer.data(dh_s + 91);
    const auto *dh_s_92 = buffer.data(dh_s + 92);
    const auto *dh_s_93 = buffer.data(dh_s + 93);
    const auto *dh_s_94 = buffer.data(dh_s + 94);
    const auto *dh_s_95 = buffer.data(dh_s + 95);
    const auto *dh_s_96 = buffer.data(dh_s + 96);
    const auto *dh_s_97 = buffer.data(dh_s + 97);
    const auto *dh_s_98 = buffer.data(dh_s + 98);
    const auto *dh_s_99 = buffer.data(dh_s + 99);
    const auto *dh_s_100 = buffer.data(dh_s + 100);
    const auto *dh_s_101 = buffer.data(dh_s + 101);
    const auto *dh_s_102 = buffer.data(dh_s + 102);
    const auto *dh_s_103 = buffer.data(dh_s + 103);
    const auto *dh_s_104 = buffer.data(dh_s + 104);
    const auto *dh_s_105 = buffer.data(dh_s + 105);
    const auto *dh_s_106 = buffer.data(dh_s + 106);
    const auto *dh_s_107 = buffer.data(dh_s + 107);
    const auto *dh_s_108 = buffer.data(dh_s + 108);
    const auto *dh_s_109 = buffer.data(dh_s + 109);
    const auto *dh_s_110 = buffer.data(dh_s + 110);
    const auto *dh_s_111 = buffer.data(dh_s + 111);
    const auto *dh_s_112 = buffer.data(dh_s + 112);
    const auto *dh_s_113 = buffer.data(dh_s + 113);
    const auto *dh_s_114 = buffer.data(dh_s + 114);
    const auto *dh_s_115 = buffer.data(dh_s + 115);
    const auto *dh_s_116 = buffer.data(dh_s + 116);
    const auto *dh_s_117 = buffer.data(dh_s + 117);
    const auto *dh_s_118 = buffer.data(dh_s + 118);
    const auto *dh_s_119 = buffer.data(dh_s + 119);
    const auto *dh_s_120 = buffer.data(dh_s + 120);
    const auto *dh_s_121 = buffer.data(dh_s + 121);
    const auto *dh_s_122 = buffer.data(dh_s + 122);
    const auto *dh_s_123 = buffer.data(dh_s + 123);
    const auto *dh_s_124 = buffer.data(dh_s + 124);
    const auto *dh_s_125 = buffer.data(dh_s + 125);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_52 = buffer.data(dg + 52);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, dh_s_3, df_0, dg_0, dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_2 * dh_s_1[k]
                 + pb_y[k] * dg_0[k];

        t_2[k] = f_2 * dh_s_2[k]
                 + pb_z[k] * dg_0[k];

        t_3[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_3[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, df_s_0, df_s_1, dh_s_4, dh_s_5, \
                         dh_s_6, dh_s_7, df_0, df_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * dh_s_4[k]
                 + pb_y[k] * dg_2[k];

        t_5[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_5[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];

        t_6[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_6[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];

        t_7[k] = f_2 * dh_s_7[k]
                 + pb_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, pg_2, df_s_2, dh_s_8, dh_s_9, \
                         dh_s_10, df_2, dg_4, dg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * dh_s_8[k]
                 + pb_y[k] * dg_4[k];

        t_9[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_9[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];

        t_10[k] = f_0 * pg_2[k]
                  + f_2 * dh_s_10[k]
                  + pb_x[k] * dg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, pg_3, dh_s_11, dh_s_12, dh_s_13, \
                         dg_5, dg_6, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * dh_s_11[k]
                  + pb_z[k] * dg_5[k];

        t_12[k] = f_0 * pg_3[k]
                  + f_2 * dh_s_12[k]
                  + pb_x[k] * dg_8[k];

        t_13[k] = f_2 * dh_s_13[k]
                  + pb_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, pg_4, df_s_3, dh_s_14, dh_s_15, \
                         dh_s_16, df_3, dg_7, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * pg_4[k]
                  + f_2 * dh_s_14[k]
                  + pb_x[k] * dg_10[k];

        t_15[k] = -f_1 * df_s_3[k]
                  + f_2 * dh_s_15[k]
                  + f_3 * df_3[k]
                  + pb_y[k] * dg_7[k];

        t_16[k] = f_2 * dh_s_16[k]
                  + pb_z[k] * dg_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, df_s_4, df_s_5, dh_s_17, dh_s_18, dh_s_19, \
                         df_4, df_5, dg_8, dg_9, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_6 * df_s_4[k]
                  + f_2 * dh_s_17[k]
                  + f_0 * df_4[k]
                  + pb_y[k] * dg_8[k];

        t_18[k] = -f_4 * df_s_5[k]
                  + f_2 * dh_s_18[k]
                  + f_5 * df_5[k]
                  + pb_y[k] * dg_9[k];

        t_19[k] = f_2 * dh_s_19[k]
                  + pb_y[k] * dg_10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, pg_0, ph_0, df_s_5, dh_s_20, \
                         dh_s_21, dh_s_22, df_5, dg_10, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * df_s_5[k]
                  + f_2 * dh_s_20[k]
                  + f_3 * df_5[k]
                  + pb_z[k] * dg_10[k];

        t_21[k] = pa_y[k] * ph_0[k]
                  + f_2 * dh_s_21[k];

        t_22[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_22[k]
                  + pb_y[k] * dg_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, pb_z, pg_5, ph_2, ph_8, dh_s_23, \
                         dh_s_24, dh_s_25, dh_s_26, dg_11, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * dh_s_23[k]
                  + pb_z[k] * dg_11[k];

        t_24[k] = f_7 * pg_5[k]
                  + pa_x[k] * ph_8[k]
                  + f_2 * dh_s_24[k];

        t_25[k] = f_2 * dh_s_25[k]
                  + pb_z[k] * dg_12[k];

        t_26[k] = pa_y[k] * ph_2[k]
                  + f_2 * dh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_y, pb_z, pg_1, pg_6, ph_9, dh_s_27, \
                         dh_s_28, dh_s_29, dg_13, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pg_6[k]
                  + pa_x[k] * ph_9[k]
                  + f_2 * dh_s_27[k];

        t_28[k] = f_2 * dh_s_28[k]
                  + pb_z[k] * dg_13[k];

        t_29[k] = f_5 * pg_1[k]
                  + f_2 * dh_s_29[k]
                  + pb_y[k] * dg_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, pg_7, ph_4, dh_s_30, dh_s_31, \
                         dh_s_32, dg_15, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * ph_4[k]
                  + f_2 * dh_s_30[k];

        t_31[k] = f_5 * pg_7[k]
                  + f_2 * dh_s_31[k]
                  + pb_x[k] * dg_16[k];

        t_32[k] = f_2 * dh_s_32[k]
                  + pb_z[k] * dg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, pg_8, pg_9, ph_6, dh_s_33, dh_s_34, \
                         dh_s_35, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_5 * pg_8[k]
                  + f_2 * dh_s_33[k]
                  + pb_x[k] * dg_17[k];

        t_34[k] = f_5 * pg_9[k]
                  + f_2 * dh_s_34[k]
                  + pb_x[k] * dg_18[k];

        t_35[k] = pa_y[k] * ph_6[k]
                  + f_2 * dh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_z, ph_10, ph_11, ph_12, dh_s_36, \
                         dh_s_37, dh_s_38, dh_s_39, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * ph_10[k]
                  + f_2 * dh_s_36[k];

        t_37[k] = f_2 * dh_s_37[k]
                  + pb_z[k] * dg_16[k];

        t_38[k] = pa_x[k] * ph_11[k]
                  + f_2 * dh_s_38[k];

        t_39[k] = pa_x[k] * ph_12[k]
                  + f_2 * dh_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_z, pb_y, ph_0, ph_13, ph_14, \
                         dh_s_40, dh_s_41, dh_s_42, dh_s_43, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * ph_13[k]
                  + f_2 * dh_s_40[k];

        t_41[k] = pa_x[k] * ph_14[k]
                  + f_2 * dh_s_41[k];

        t_42[k] = pa_z[k] * ph_0[k]
                  + f_2 * dh_s_42[k];

        t_43[k] = f_2 * dh_s_43[k]
                  + pb_y[k] * dg_19[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_y, pb_z, pg_0, ph_1, dh_s_44, dh_s_45, \
                         dh_s_46, dg_19, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_44[k]
                  + pb_z[k] * dg_19[k];

        t_45[k] = pa_z[k] * ph_1[k]
                  + f_2 * dh_s_45[k];

        t_46[k] = f_2 * dh_s_46[k]
                  + pb_y[k] * dg_20[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_z, pg_13, pg_14, ph_3, ph_18, ph_19, \
                         dh_s_47, dh_s_48, dh_s_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_7 * pg_13[k]
                  + pa_x[k] * ph_18[k]
                  + f_2 * dh_s_47[k];

        t_48[k] = pa_z[k] * ph_3[k]
                  + f_2 * dh_s_48[k];

        t_49[k] = f_0 * pg_14[k]
                  + pa_x[k] * ph_19[k]
                  + f_2 * dh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pa_z, pb_y, pg_15, ph_5, ph_21, dh_s_50, \
                         dh_s_51, dh_s_52, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * dh_s_50[k]
                  + pb_y[k] * dg_21[k];

        t_51[k] = f_0 * pg_15[k]
                  + pa_x[k] * ph_21[k]
                  + f_2 * dh_s_51[k];

        t_52[k] = pa_z[k] * ph_5[k]
                  + f_2 * dh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, pg_16, pg_17, dh_s_53, dh_s_54, \
                         dh_s_55, dg_22, dg_23, dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_5 * pg_16[k]
                  + f_2 * dh_s_53[k]
                  + pb_x[k] * dg_23[k];

        t_54[k] = f_5 * pg_17[k]
                  + f_2 * dh_s_54[k]
                  + pb_x[k] * dg_24[k];

        t_55[k] = f_2 * dh_s_55[k]
                  + pb_y[k] * dg_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pg_19, ph_22, ph_23, ph_24, \
                         dh_s_56, dh_s_57, dh_s_58, dh_s_59, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_5 * pg_19[k]
                  + f_2 * dh_s_56[k]
                  + pb_x[k] * dg_25[k];

        t_57[k] = pa_x[k] * ph_22[k]
                  + f_2 * dh_s_57[k];

        t_58[k] = pa_x[k] * ph_23[k]
                  + f_2 * dh_s_58[k];

        t_59[k] = pa_x[k] * ph_24[k]
                  + f_2 * dh_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pb_y, ph_25, ph_26, dh_s_60, dh_s_61, \
                         dh_s_62, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * ph_25[k]
                  + f_2 * dh_s_60[k];

        t_61[k] = f_2 * dh_s_61[k]
                  + pb_y[k] * dg_25[k];

        t_62[k] = pa_x[k] * ph_26[k]
                  + f_2 * dh_s_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, pb_z, df_s_11, df_s_12, dh_s_63, dh_s_64, \
                         dh_s_65, df_10, df_11, dg_26, dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -f_1 * df_s_11[k]
                  + f_2 * dh_s_63[k]
                  + f_3 * df_10[k]
                  + pb_x[k] * dg_26[k];

        t_64[k] = -f_8 * df_s_12[k]
                  + f_2 * dh_s_64[k]
                  + f_7 * df_11[k]
                  + pb_x[k] * dg_27[k];

        t_65[k] = f_2 * dh_s_65[k]
                  + pb_z[k] * dg_26[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, df_s_13, df_s_14, dh_s_66, dh_s_67, \
                         dh_s_68, df_12, df_13, dg_27, dg_28, dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -f_6 * df_s_13[k]
                  + f_2 * dh_s_66[k]
                  + f_0 * df_12[k]
                  + pb_x[k] * dg_28[k];

        t_67[k] = f_2 * dh_s_67[k]
                  + pb_z[k] * dg_27[k];

        t_68[k] = -f_6 * df_s_14[k]
                  + f_2 * dh_s_68[k]
                  + f_0 * df_13[k]
                  + pb_x[k] * dg_29[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_z, df_s_15, df_s_17, dh_s_69, dh_s_70, \
                         dh_s_71, df_14, df_16, dg_28, dg_30, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = -f_4 * df_s_15[k]
                  + f_2 * dh_s_69[k]
                  + f_5 * df_14[k]
                  + pb_x[k] * dg_30[k];

        t_70[k] = f_2 * dh_s_70[k]
                  + pb_z[k] * dg_28[k];

        t_71[k] = -f_4 * df_s_17[k]
                  + f_2 * dh_s_71[k]
                  + f_5 * df_16[k]
                  + pb_x[k] * dg_31[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, df_s_18, dh_s_72, dh_s_73, dh_s_74, \
                         dh_s_75, df_17, dg_32, dg_33, dg_34, dg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -f_4 * df_s_18[k]
                  + f_2 * dh_s_72[k]
                  + f_5 * df_17[k]
                  + pb_x[k] * dg_32[k];

        t_73[k] = f_2 * dh_s_73[k]
                  + pb_x[k] * dg_33[k];

        t_74[k] = f_2 * dh_s_74[k]
                  + pb_x[k] * dg_34[k];

        t_75[k] = f_2 * dh_s_75[k]
                  + pb_x[k] * dg_35[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_x, pb_y, pg_7, df_s_15, dh_s_76, dh_s_77, \
                         dh_s_78, df_14, dg_33, dg_36, dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * dh_s_76[k]
                  + pb_x[k] * dg_36[k];

        t_77[k] = f_2 * dh_s_77[k]
                  + pb_x[k] * dg_37[k];

        t_78[k] = f_0 * pg_7[k]
                  - f_1 * df_s_15[k]
                  + f_2 * dh_s_78[k]
                  + f_3 * df_14[k]
                  + pb_y[k] * dg_33[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_z, df_s_15, df_s_16, dh_s_79, dh_s_80, dh_s_81, \
                         df_14, df_15, dg_33, dg_34, dg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_2 * dh_s_79[k]
                  + pb_z[k] * dg_33[k];

        t_80[k] = -f_4 * df_s_15[k]
                  + f_2 * dh_s_80[k]
                  + f_5 * df_14[k]
                  + pb_z[k] * dg_34[k];

        t_81[k] = -f_6 * df_s_16[k]
                  + f_2 * dh_s_81[k]
                  + f_0 * df_15[k]
                  + pb_z[k] * dg_35[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pb_y, pb_z, pg_10, ph_15, df_s_18, dh_s_82, \
                         dh_s_83, dh_s_84, df_17, dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * pg_10[k]
                  + f_2 * dh_s_82[k]
                  + pb_y[k] * dg_37[k];

        t_83[k] = -f_1 * df_s_18[k]
                  + f_2 * dh_s_83[k]
                  + f_3 * df_17[k]
                  + pb_z[k] * dg_37[k];

        t_84[k] = pa_y[k] * ph_15[k]
                  + f_2 * dh_s_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pa_z, pg_11, ph_7, ph_8, ph_16, ph_17, \
                         dh_s_85, dh_s_86, dh_s_87, dh_s_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_z[k] * ph_7[k]
                  + f_2 * dh_s_85[k];

        t_86[k] = pa_y[k] * ph_16[k]
                  + f_2 * dh_s_86[k];

        t_87[k] = pa_z[k] * ph_8[k]
                  + f_2 * dh_s_87[k];

        t_88[k] = f_5 * pg_11[k]
                  + pa_y[k] * ph_17[k]
                  + f_2 * dh_s_88[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pa_z, pg_12, pg_13, ph_9, ph_18, ph_19, \
                         ph_20, dh_s_89, dh_s_90, dh_s_91, dh_s_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * ph_18[k]
                  + f_2 * dh_s_89[k];

        t_90[k] = pa_z[k] * ph_9[k]
                  + f_2 * dh_s_90[k];

        t_91[k] = f_0 * pg_12[k]
                  + pa_y[k] * ph_19[k]
                  + f_2 * dh_s_91[k];

        t_92[k] = f_5 * pg_13[k]
                  + pa_y[k] * ph_20[k]
                  + f_2 * dh_s_92[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_x, ph_21, dh_s_93, dh_s_94, dh_s_95, \
                         dh_s_96, dg_38, dg_39, dg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * ph_21[k]
                  + f_2 * dh_s_93[k];

        t_94[k] = f_2 * dh_s_94[k]
                  + pb_x[k] * dg_38[k];

        t_95[k] = f_2 * dh_s_95[k]
                  + pb_x[k] * dg_39[k];

        t_96[k] = f_2 * dh_s_96[k]
                  + pb_x[k] * dg_40[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_z, pb_x, pb_z, pg_7, ph_10, dh_s_97, \
                         dh_s_98, dh_s_99, dh_s_100, dg_38, dg_41, \
                         dg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * dh_s_97[k]
                  + pb_x[k] * dg_41[k];

        t_98[k] = f_2 * dh_s_98[k]
                  + pb_x[k] * dg_42[k];

        t_99[k] = pa_z[k] * ph_10[k]
                  + f_2 * dh_s_99[k];

        t_100[k] = f_5 * pg_7[k]
                   + f_2 * dh_s_100[k]
                   + pb_z[k] * dg_38[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_y, pb_y, pg_17, pg_18, pg_19, ph_24, ph_25, \
                         dh_s_101, dh_s_102, dh_s_103, dg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_7 * pg_17[k]
                   + pa_y[k] * ph_24[k]
                   + f_2 * dh_s_101[k];

        t_102[k] = f_0 * pg_18[k]
                   + pa_y[k] * ph_25[k]
                   + f_2 * dh_s_102[k];

        t_103[k] = f_5 * pg_19[k]
                   + f_2 * dh_s_103[k]
                   + pb_y[k] * dg_42[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pb_x, pb_y, ph_26, df_s_24, dh_s_104, \
                         dh_s_105, dh_s_106, df_20, dg_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * ph_26[k]
                   + f_2 * dh_s_104[k];

        t_105[k] = -f_1 * df_s_24[k]
                   + f_2 * dh_s_105[k]
                   + f_3 * df_20[k]
                   + pb_x[k] * dg_43[k];

        t_106[k] = f_2 * dh_s_106[k]
                   + pb_y[k] * dg_43[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, pb_y, df_s_25, df_s_26, dh_s_107, \
                         dh_s_108, dh_s_109, df_21, df_22, dg_44, \
                         dg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_8 * df_s_25[k]
                   + f_2 * dh_s_107[k]
                   + f_7 * df_21[k]
                   + pb_x[k] * dg_44[k];

        t_108[k] = -f_6 * df_s_26[k]
                   + f_2 * dh_s_108[k]
                   + f_0 * df_22[k]
                   + pb_x[k] * dg_45[k];

        t_109[k] = f_2 * dh_s_109[k]
                   + pb_y[k] * dg_44[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, df_s_27, df_s_28, df_s_29, dh_s_110, \
                         dh_s_111, dh_s_112, df_23, df_24, df_25, dg_46, dg_47, \
                         dg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -f_6 * df_s_27[k]
                   + f_2 * dh_s_110[k]
                   + f_0 * df_23[k]
                   + pb_x[k] * dg_46[k];

        t_111[k] = -f_4 * df_s_28[k]
                   + f_2 * dh_s_111[k]
                   + f_5 * df_24[k]
                   + pb_x[k] * dg_47[k];

        t_112[k] = -f_4 * df_s_29[k]
                   + f_2 * dh_s_112[k]
                   + f_5 * df_25[k]
                   + pb_x[k] * dg_48[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pb_y, df_s_31, dh_s_113, dh_s_114, \
                         dh_s_115, dh_s_116, df_27, dg_46, dg_49, dg_50, \
                         dg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_2 * dh_s_113[k]
                   + pb_y[k] * dg_46[k];

        t_114[k] = -f_4 * df_s_31[k]
                   + f_2 * dh_s_114[k]
                   + f_5 * df_27[k]
                   + pb_x[k] * dg_49[k];

        t_115[k] = f_2 * dh_s_115[k]
                   + pb_x[k] * dg_50[k];

        t_116[k] = f_2 * dh_s_116[k]
                   + pb_x[k] * dg_51[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, df_s_28, dh_s_117, dh_s_118, \
                         dh_s_119, dh_s_120, df_24, dg_50, dg_52, dg_53, \
                         dg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_2 * dh_s_117[k]
                   + pb_x[k] * dg_52[k];

        t_118[k] = f_2 * dh_s_118[k]
                   + pb_x[k] * dg_53[k];

        t_119[k] = f_2 * dh_s_119[k]
                   + pb_x[k] * dg_54[k];

        t_120[k] = -f_1 * df_s_28[k]
                   + f_2 * dh_s_120[k]
                   + f_3 * df_24[k]
                   + pb_y[k] * dg_50[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_y, df_s_29, df_s_30, df_s_31, dh_s_121, \
                         dh_s_122, dh_s_123, df_25, df_26, df_27, dg_51, dg_52, \
                         dg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_8 * df_s_29[k]
                   + f_2 * dh_s_121[k]
                   + f_7 * df_25[k]
                   + pb_y[k] * dg_51[k];

        t_122[k] = -f_6 * df_s_30[k]
                   + f_2 * dh_s_122[k]
                   + f_0 * df_26[k]
                   + pb_y[k] * dg_52[k];

        t_123[k] = -f_4 * df_s_31[k]
                   + f_2 * dh_s_123[k]
                   + f_5 * df_27[k]
                   + pb_y[k] * dg_53[k];
    }

#pragma omp simd aligned(t_124, t_125, pb_y, pb_z, pg_19, df_s_31, dh_s_124, dh_s_125, df_27, \
                         dg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * dh_s_124[k]
                   + pb_y[k] * dg_54[k];

        t_125[k] = f_0 * pg_19[k]
                   - f_1 * df_s_31[k]
                   + f_2 * dh_s_125[k]
                   + f_3 * df_27[k]
                   + pb_z[k] * dg_54[k];
    }
}

auto
compute_prim_dh_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_12 = buffer.data(ph + 12);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_25 = buffer.data(ph + 25);
    const auto *ph_26 = buffer.data(ph + 26);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_29 = buffer.data(ph + 29);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_3 = buffer.data(df_s + 3);
    const auto *df_s_4 = buffer.data(df_s + 4);
    const auto *df_s_5 = buffer.data(df_s + 5);
    const auto *df_s_8 = buffer.data(df_s + 8);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_10 = buffer.data(df_s + 10);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_13 = buffer.data(df_s + 13);
    const auto *df_s_14 = buffer.data(df_s + 14);
    const auto *df_s_15 = buffer.data(df_s + 15);
    const auto *df_s_21 = buffer.data(df_s + 21);
    const auto *df_s_22 = buffer.data(df_s + 22);
    const auto *df_s_23 = buffer.data(df_s + 23);
    const auto *df_s_24 = buffer.data(df_s + 24);
    const auto *df_s_25 = buffer.data(df_s + 25);
    const auto *df_s_26 = buffer.data(df_s + 26);
    const auto *df_s_27 = buffer.data(df_s + 27);
    const auto *df_s_28 = buffer.data(df_s + 28);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);
    const auto *dh_s_25 = buffer.data(dh_s + 25);
    const auto *dh_s_26 = buffer.data(dh_s + 26);
    const auto *dh_s_27 = buffer.data(dh_s + 27);
    const auto *dh_s_28 = buffer.data(dh_s + 28);
    const auto *dh_s_30 = buffer.data(dh_s + 30);
    const auto *dh_s_33 = buffer.data(dh_s + 33);
    const auto *dh_s_34 = buffer.data(dh_s + 34);
    const auto *dh_s_35 = buffer.data(dh_s + 35);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_37 = buffer.data(dh_s + 37);
    const auto *dh_s_38 = buffer.data(dh_s + 38);
    const auto *dh_s_39 = buffer.data(dh_s + 39);
    const auto *dh_s_40 = buffer.data(dh_s + 40);
    const auto *dh_s_41 = buffer.data(dh_s + 41);
    const auto *dh_s_42 = buffer.data(dh_s + 42);
    const auto *dh_s_43 = buffer.data(dh_s + 43);
    const auto *dh_s_44 = buffer.data(dh_s + 44);
    const auto *dh_s_45 = buffer.data(dh_s + 45);
    const auto *dh_s_46 = buffer.data(dh_s + 46);
    const auto *dh_s_47 = buffer.data(dh_s + 47);
    const auto *dh_s_48 = buffer.data(dh_s + 48);
    const auto *dh_s_49 = buffer.data(dh_s + 49);
    const auto *dh_s_50 = buffer.data(dh_s + 50);
    const auto *dh_s_51 = buffer.data(dh_s + 51);
    const auto *dh_s_52 = buffer.data(dh_s + 52);
    const auto *dh_s_53 = buffer.data(dh_s + 53);
    const auto *dh_s_54 = buffer.data(dh_s + 54);
    const auto *dh_s_55 = buffer.data(dh_s + 55);
    const auto *dh_s_62 = buffer.data(dh_s + 62);
    const auto *dh_s_63 = buffer.data(dh_s + 63);
    const auto *dh_s_64 = buffer.data(dh_s + 64);
    const auto *dh_s_65 = buffer.data(dh_s + 65);
    const auto *dh_s_66 = buffer.data(dh_s + 66);
    const auto *dh_s_67 = buffer.data(dh_s + 67);
    const auto *dh_s_68 = buffer.data(dh_s + 68);
    const auto *dh_s_69 = buffer.data(dh_s + 69);
    const auto *dh_s_70 = buffer.data(dh_s + 70);
    const auto *dh_s_71 = buffer.data(dh_s + 71);
    const auto *dh_s_72 = buffer.data(dh_s + 72);
    const auto *dh_s_73 = buffer.data(dh_s + 73);
    const auto *dh_s_74 = buffer.data(dh_s + 74);
    const auto *dh_s_75 = buffer.data(dh_s + 75);
    const auto *dh_s_76 = buffer.data(dh_s + 76);
    const auto *dh_s_77 = buffer.data(dh_s + 77);
    const auto *dh_s_78 = buffer.data(dh_s + 78);
    const auto *dh_s_79 = buffer.data(dh_s + 79);
    const auto *dh_s_80 = buffer.data(dh_s + 80);
    const auto *dh_s_81 = buffer.data(dh_s + 81);
    const auto *dh_s_82 = buffer.data(dh_s + 82);
    const auto *dh_s_83 = buffer.data(dh_s + 83);
    const auto *dh_s_84 = buffer.data(dh_s + 84);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, dh_s_3, df_0, dg_0, dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_2 * dh_s_1[k]
                 + pb_y[k] * dg_0[k];

        t_2[k] = f_2 * dh_s_2[k]
                 + pb_z[k] * dg_0[k];

        t_3[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_3[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, df_s_0, df_s_1, dh_s_4, dh_s_5, dh_s_6, \
                         df_0, df_1, dg_2, dg_3, dg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_4[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];

        t_5[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_5[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = f_2 * dh_s_6[k]
                 + pb_y[k] * dg_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_z, pg_1, pg_2, df_s_2, dh_s_7, dh_s_8, \
                         dh_s_9, df_2, dg_4, dg_5, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_7[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];

        t_8[k] = f_0 * pg_1[k]
                 + f_2 * dh_s_8[k]
                 + pb_x[k] * dg_5[k];

        t_9[k] = f_0 * pg_2[k]
                 + f_2 * dh_s_9[k]
                 + pb_x[k] * dg_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, df_s_3, df_s_4, df_s_5, dh_s_10, dh_s_11, \
                         dh_s_12, df_3, df_4, df_5, dg_5, dg_6, dg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * df_s_3[k]
                  + f_2 * dh_s_10[k]
                  + f_3 * df_3[k]
                  + pb_y[k] * dg_5[k];

        t_11[k] = -f_6 * df_s_4[k]
                  + f_2 * dh_s_11[k]
                  + f_0 * df_4[k]
                  + pb_y[k] * dg_6[k];

        t_12[k] = -f_4 * df_s_5[k]
                  + f_2 * dh_s_12[k]
                  + f_5 * df_5[k]
                  + pb_y[k] * dg_7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, ph_0, df_s_5, dh_s_13, dh_s_14, \
                         dh_s_15, df_5, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * dh_s_13[k]
                  + pb_y[k] * dg_8[k];

        t_14[k] = -f_1 * df_s_5[k]
                  + f_2 * dh_s_14[k]
                  + f_3 * df_5[k]
                  + pb_z[k] * dg_8[k];

        t_15[k] = pa_y[k] * ph_0[k]
                  + f_2 * dh_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pa_y, pb_y, pg_0, pg_3, ph_1, ph_7, dh_s_16, \
                         dh_s_17, dh_s_18, dg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_16[k]
                  + pb_y[k] * dg_9[k];

        t_17[k] = f_7 * pg_3[k]
                  + pa_x[k] * ph_7[k]
                  + f_2 * dh_s_17[k];

        t_18[k] = pa_y[k] * ph_1[k]
                  + f_2 * dh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pa_y, pb_x, pg_4, pg_5, ph_2, ph_8, dh_s_19, \
                         dh_s_21, dh_s_22, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * pg_4[k]
                  + pa_x[k] * ph_8[k]
                  + f_2 * dh_s_19[k];

        t_20[k] = pa_y[k] * ph_2[k]
                  + f_2 * dh_s_21[k];

        t_21[k] = f_5 * pg_5[k]
                  + f_2 * dh_s_22[k]
                  + pb_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, ph_11, ph_12, ph_13, ph_14, dh_s_23, \
                         dh_s_24, dh_s_25, dh_s_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pa_x[k] * ph_11[k]
                  + f_2 * dh_s_23[k];

        t_23[k] = pa_x[k] * ph_12[k]
                  + f_2 * dh_s_24[k];

        t_24[k] = pa_x[k] * ph_13[k]
                  + f_2 * dh_s_25[k];

        t_25[k] = pa_x[k] * ph_14[k]
                  + f_2 * dh_s_26[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_x, pa_z, pb_z, pg_0, pg_9, ph_0, ph_18, dh_s_27, \
                         dh_s_28, dh_s_30, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * ph_0[k]
                  + f_2 * dh_s_27[k];

        t_27[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_28[k]
                  + pb_z[k] * dg_11[k];

        t_28[k] = f_7 * pg_9[k]
                  + pa_x[k] * ph_18[k]
                  + f_2 * dh_s_30[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pb_x, pg_10, pg_14, ph_21, ph_25, \
                         ph_26, dh_s_33, dh_s_34, dh_s_35, dh_s_36, \
                         dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * pg_10[k]
                  + pa_x[k] * ph_21[k]
                  + f_2 * dh_s_33[k];

        t_30[k] = f_5 * pg_14[k]
                  + f_2 * dh_s_34[k]
                  + pb_x[k] * dg_12[k];

        t_31[k] = pa_x[k] * ph_25[k]
                  + f_2 * dh_s_35[k];

        t_32[k] = pa_x[k] * ph_26[k]
                  + f_2 * dh_s_36[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_x, pb_x, ph_27, ph_29, df_s_8, dh_s_37, dh_s_38, \
                         dh_s_39, df_6, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_x[k] * ph_27[k]
                  + f_2 * dh_s_37[k];

        t_34[k] = pa_x[k] * ph_29[k]
                  + f_2 * dh_s_38[k];

        t_35[k] = -f_1 * df_s_8[k]
                  + f_2 * dh_s_39[k]
                  + f_3 * df_6[k]
                  + pb_x[k] * dg_13[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_x, df_s_9, df_s_10, df_s_11, dh_s_40, dh_s_41, \
                         dh_s_42, df_7, df_8, df_9, dg_14, dg_15, \
                         dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * df_s_9[k]
                  + f_2 * dh_s_40[k]
                  + f_7 * df_7[k]
                  + pb_x[k] * dg_14[k];

        t_37[k] = -f_6 * df_s_10[k]
                  + f_2 * dh_s_41[k]
                  + f_0 * df_8[k]
                  + pb_x[k] * dg_15[k];

        t_38[k] = -f_6 * df_s_11[k]
                  + f_2 * dh_s_42[k]
                  + f_0 * df_9[k]
                  + pb_x[k] * dg_16[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_x, df_s_12, df_s_14, df_s_15, dh_s_43, dh_s_44, \
                         dh_s_45, df_10, df_12, df_13, dg_17, dg_18, \
                         dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_4 * df_s_12[k]
                  + f_2 * dh_s_43[k]
                  + f_5 * df_10[k]
                  + pb_x[k] * dg_17[k];

        t_40[k] = -f_4 * df_s_14[k]
                  + f_2 * dh_s_44[k]
                  + f_5 * df_12[k]
                  + pb_x[k] * dg_18[k];

        t_41[k] = -f_4 * df_s_15[k]
                  + f_2 * dh_s_45[k]
                  + f_5 * df_13[k]
                  + pb_x[k] * dg_19[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_x, dh_s_46, dh_s_47, dh_s_48, dh_s_49, \
                         dg_20, dg_22, dg_23, dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * dh_s_46[k]
                  + pb_x[k] * dg_20[k];

        t_43[k] = f_2 * dh_s_47[k]
                  + pb_x[k] * dg_22[k];

        t_44[k] = f_2 * dh_s_48[k]
                  + pb_x[k] * dg_23[k];

        t_45[k] = f_2 * dh_s_49[k]
                  + pb_x[k] * dg_24[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pb_z, pg_5, df_s_12, dh_s_50, dh_s_51, \
                         dh_s_52, df_10, dg_20, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_0 * pg_5[k]
                  - f_1 * df_s_12[k]
                  + f_2 * dh_s_50[k]
                  + f_3 * df_10[k]
                  + pb_y[k] * dg_20[k];

        t_47[k] = f_2 * dh_s_51[k]
                  + pb_z[k] * dg_20[k];

        t_48[k] = -f_4 * df_s_12[k]
                  + f_2 * dh_s_52[k]
                  + f_5 * df_10[k]
                  + pb_z[k] * dg_21[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_y, pb_z, pg_7, df_s_13, df_s_15, dh_s_53, \
                         dh_s_54, dh_s_55, df_11, df_13, dg_22, dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_6 * df_s_13[k]
                  + f_2 * dh_s_53[k]
                  + f_0 * df_11[k]
                  + pb_z[k] * dg_22[k];

        t_50[k] = f_0 * pg_7[k]
                  + f_2 * dh_s_54[k]
                  + pb_y[k] * dg_24[k];

        t_51[k] = -f_1 * df_s_15[k]
                  + f_2 * dh_s_55[k]
                  + f_3 * df_13[k]
                  + pb_z[k] * dg_24[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pa_z, pb_z, pg_5, pg_12, ph_11, ph_26, \
                         dh_s_62, dh_s_63, dh_s_64, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_z[k] * ph_11[k]
                  + f_2 * dh_s_62[k];

        t_53[k] = f_5 * pg_5[k]
                  + f_2 * dh_s_63[k]
                  + pb_z[k] * dg_25[k];

        t_54[k] = f_7 * pg_12[k]
                  + pa_y[k] * ph_26[k]
                  + f_2 * dh_s_64[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pb_y, pg_13, pg_14, ph_27, ph_29, dh_s_65, \
                         dh_s_66, dh_s_67, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * pg_13[k]
                  + pa_y[k] * ph_27[k]
                  + f_2 * dh_s_65[k];

        t_56[k] = f_5 * pg_14[k]
                  + f_2 * dh_s_66[k]
                  + pb_y[k] * dg_26[k];

        t_57[k] = pa_y[k] * ph_29[k]
                  + f_2 * dh_s_67[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, df_s_21, df_s_22, df_s_23, dh_s_68, dh_s_69, \
                         dh_s_70, df_14, df_15, df_16, dg_27, dg_28, \
                         dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -f_1 * df_s_21[k]
                  + f_2 * dh_s_68[k]
                  + f_3 * df_14[k]
                  + pb_x[k] * dg_27[k];

        t_59[k] = -f_8 * df_s_22[k]
                  + f_2 * dh_s_69[k]
                  + f_7 * df_15[k]
                  + pb_x[k] * dg_28[k];

        t_60[k] = -f_6 * df_s_23[k]
                  + f_2 * dh_s_70[k]
                  + f_0 * df_16[k]
                  + pb_x[k] * dg_29[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, df_s_24, df_s_25, df_s_26, dh_s_71, dh_s_72, \
                         dh_s_73, df_17, df_18, df_19, dg_30, dg_31, \
                         dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -f_6 * df_s_24[k]
                  + f_2 * dh_s_71[k]
                  + f_0 * df_17[k]
                  + pb_x[k] * dg_30[k];

        t_62[k] = -f_4 * df_s_25[k]
                  + f_2 * dh_s_72[k]
                  + f_5 * df_18[k]
                  + pb_x[k] * dg_31[k];

        t_63[k] = -f_4 * df_s_26[k]
                  + f_2 * dh_s_73[k]
                  + f_5 * df_19[k]
                  + pb_x[k] * dg_32[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_x, df_s_28, dh_s_74, dh_s_75, dh_s_76, \
                         dh_s_77, df_21, dg_33, dg_34, dg_35, dg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = -f_4 * df_s_28[k]
                  + f_2 * dh_s_74[k]
                  + f_5 * df_21[k]
                  + pb_x[k] * dg_33[k];

        t_65[k] = f_2 * dh_s_75[k]
                  + pb_x[k] * dg_34[k];

        t_66[k] = f_2 * dh_s_76[k]
                  + pb_x[k] * dg_35[k];

        t_67[k] = f_2 * dh_s_77[k]
                  + pb_x[k] * dg_36[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_y, df_s_25, df_s_26, dh_s_78, dh_s_79, \
                         dh_s_80, df_18, df_19, dg_34, dg_35, dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * dh_s_78[k]
                  + pb_x[k] * dg_38[k];

        t_69[k] = -f_1 * df_s_25[k]
                  + f_2 * dh_s_79[k]
                  + f_3 * df_18[k]
                  + pb_y[k] * dg_34[k];

        t_70[k] = -f_8 * df_s_26[k]
                  + f_2 * dh_s_80[k]
                  + f_7 * df_19[k]
                  + pb_y[k] * dg_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pb_y, df_s_27, df_s_28, dh_s_81, dh_s_82, dh_s_83, \
                         df_20, df_21, dg_36, dg_37, dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -f_6 * df_s_27[k]
                  + f_2 * dh_s_81[k]
                  + f_0 * df_20[k]
                  + pb_y[k] * dg_36[k];

        t_72[k] = -f_4 * df_s_28[k]
                  + f_2 * dh_s_82[k]
                  + f_5 * df_21[k]
                  + pb_y[k] * dg_37[k];

        t_73[k] = f_2 * dh_s_83[k]
                  + pb_y[k] * dg_38[k];
    }

#pragma omp simd aligned(t_74, pb_z, pg_14, df_s_28, dh_s_84, df_21, \
                         dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * pg_14[k]
                  - f_1 * df_s_28[k]
                  + f_2 * dh_s_84[k]
                  + f_3 * df_21[k]
                  + pb_z[k] * dg_38[k];
    }
}

auto
compute_prim_dh_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_9 = buffer.data(pg + 9);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_11 = buffer.data(pg + 11);
    const auto *pg_13 = buffer.data(pg + 13);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_15 = buffer.data(pg + 15);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_8 = buffer.data(ph + 8);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_3 = buffer.data(df_s + 3);
    const auto *df_s_5 = buffer.data(df_s + 5);
    const auto *df_s_10 = buffer.data(df_s + 10);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_13 = buffer.data(df_s + 13);
    const auto *df_s_14 = buffer.data(df_s + 14);
    const auto *df_s_15 = buffer.data(df_s + 15);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_17 = buffer.data(df_s + 17);
    const auto *df_s_20 = buffer.data(df_s + 20);
    const auto *df_s_21 = buffer.data(df_s + 21);
    const auto *df_s_22 = buffer.data(df_s + 22);
    const auto *df_s_23 = buffer.data(df_s + 23);
    const auto *df_s_24 = buffer.data(df_s + 24);
    const auto *df_s_25 = buffer.data(df_s + 25);
    const auto *df_s_26 = buffer.data(df_s + 26);
    const auto *df_s_27 = buffer.data(df_s + 27);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);
    const auto *dh_s_25 = buffer.data(dh_s + 25);
    const auto *dh_s_26 = buffer.data(dh_s + 26);
    const auto *dh_s_27 = buffer.data(dh_s + 27);
    const auto *dh_s_28 = buffer.data(dh_s + 28);
    const auto *dh_s_29 = buffer.data(dh_s + 29);
    const auto *dh_s_30 = buffer.data(dh_s + 30);
    const auto *dh_s_31 = buffer.data(dh_s + 31);
    const auto *dh_s_32 = buffer.data(dh_s + 32);
    const auto *dh_s_33 = buffer.data(dh_s + 33);
    const auto *dh_s_34 = buffer.data(dh_s + 34);
    const auto *dh_s_35 = buffer.data(dh_s + 35);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_37 = buffer.data(dh_s + 37);
    const auto *dh_s_38 = buffer.data(dh_s + 38);
    const auto *dh_s_39 = buffer.data(dh_s + 39);
    const auto *dh_s_40 = buffer.data(dh_s + 40);
    const auto *dh_s_41 = buffer.data(dh_s + 41);
    const auto *dh_s_42 = buffer.data(dh_s + 42);
    const auto *dh_s_43 = buffer.data(dh_s + 43);
    const auto *dh_s_44 = buffer.data(dh_s + 44);
    const auto *dh_s_45 = buffer.data(dh_s + 45);
    const auto *dh_s_46 = buffer.data(dh_s + 46);
    const auto *dh_s_47 = buffer.data(dh_s + 47);
    const auto *dh_s_48 = buffer.data(dh_s + 48);
    const auto *dh_s_49 = buffer.data(dh_s + 49);
    const auto *dh_s_50 = buffer.data(dh_s + 50);
    const auto *dh_s_51 = buffer.data(dh_s + 51);
    const auto *dh_s_52 = buffer.data(dh_s + 52);
    const auto *dh_s_53 = buffer.data(dh_s + 53);
    const auto *dh_s_54 = buffer.data(dh_s + 54);
    const auto *dh_s_55 = buffer.data(dh_s + 55);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);
    const auto *df_25 = buffer.data(df + 25);
    const auto *df_26 = buffer.data(df + 26);
    const auto *df_27 = buffer.data(df + 27);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_40 = buffer.data(dg + 40);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_43 = buffer.data(dg + 43);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_47 = buffer.data(dg + 47);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_49 = buffer.data(dg + 49);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, df_0, dg_0, dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_1[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];

        t_2[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_2[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, df_s_1, df_s_2, dh_s_3, dh_s_4, dh_s_5, \
                         df_1, df_2, dg_3, dg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_3[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];

        t_4[k] = f_2 * dh_s_4[k]
                 + pb_z[k] * dg_3[k];

        t_5[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_5[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pg_2, pg_3, df_s_3, dh_s_6, dh_s_7, \
                         dh_s_8, df_3, dg_5, dg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * pg_2[k]
                 + f_2 * dh_s_6[k]
                 + pb_x[k] * dg_5[k];

        t_7[k] = f_0 * pg_3[k]
                 + f_2 * dh_s_7[k]
                 + pb_x[k] * dg_7[k];

        t_8[k] = -f_1 * df_s_3[k]
                 + f_2 * dh_s_8[k]
                 + f_3 * df_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_y, pb_z, pg_0, pg_4, ph_1, df_s_5, dh_s_9, \
                         dh_s_10, dh_s_11, df_5, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_1 * df_s_5[k]
                 + f_2 * dh_s_9[k]
                 + f_3 * df_5[k]
                 + pb_z[k] * dg_7[k];

        t_10[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_10[k]
                  + pb_y[k] * dg_8[k];

        t_11[k] = f_7 * pg_4[k]
                  + pa_x[k] * ph_1[k]
                  + f_2 * dh_s_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_x, pg_5, pg_6, ph_2, ph_3, dh_s_12, \
                         dh_s_13, dh_s_14, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * pg_5[k]
                  + pa_x[k] * ph_2[k]
                  + f_2 * dh_s_12[k];

        t_13[k] = f_5 * pg_6[k]
                  + f_2 * dh_s_13[k]
                  + pb_x[k] * dg_11[k];

        t_14[k] = pa_x[k] * ph_3[k]
                  + f_2 * dh_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pa_z, pb_z, pg_0, pg_10, ph_0, ph_4, dh_s_15, \
                         dh_s_16, dh_s_17, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_z[k] * ph_0[k]
                  + f_2 * dh_s_15[k];

        t_16[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_16[k]
                  + pb_z[k] * dg_14[k];

        t_17[k] = f_7 * pg_10[k]
                  + pa_x[k] * ph_4[k]
                  + f_2 * dh_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pb_x, pg_11, pg_15, ph_5, ph_8, dh_s_18, \
                         dh_s_19, dh_s_20, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * pg_11[k]
                  + pa_x[k] * ph_5[k]
                  + f_2 * dh_s_18[k];

        t_19[k] = f_5 * pg_15[k]
                  + f_2 * dh_s_19[k]
                  + pb_x[k] * dg_19[k];

        t_20[k] = pa_x[k] * ph_8[k]
                  + f_2 * dh_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_x, df_s_10, df_s_11, df_s_12, dh_s_21, dh_s_22, \
                         dh_s_23, df_10, df_11, df_12, dg_20, dg_21, \
                         dg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * df_s_10[k]
                  + f_2 * dh_s_21[k]
                  + f_3 * df_10[k]
                  + pb_x[k] * dg_20[k];

        t_22[k] = -f_8 * df_s_11[k]
                  + f_2 * dh_s_22[k]
                  + f_7 * df_11[k]
                  + pb_x[k] * dg_21[k];

        t_23[k] = -f_6 * df_s_12[k]
                  + f_2 * dh_s_23[k]
                  + f_0 * df_12[k]
                  + pb_x[k] * dg_22[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_z, df_s_13, df_s_14, dh_s_24, dh_s_25, \
                         dh_s_26, df_13, df_14, dg_21, dg_24, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_2 * dh_s_24[k]
                  + pb_z[k] * dg_21[k];

        t_25[k] = -f_6 * df_s_13[k]
                  + f_2 * dh_s_25[k]
                  + f_0 * df_13[k]
                  + pb_x[k] * dg_24[k];

        t_26[k] = -f_4 * df_s_14[k]
                  + f_2 * dh_s_26[k]
                  + f_5 * df_14[k]
                  + pb_x[k] * dg_25[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pb_z, df_s_16, df_s_17, dh_s_27, dh_s_28, \
                         dh_s_29, df_16, df_17, dg_22, dg_26, dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * dh_s_27[k]
                  + pb_z[k] * dg_22[k];

        t_28[k] = -f_4 * df_s_16[k]
                  + f_2 * dh_s_28[k]
                  + f_5 * df_16[k]
                  + pb_x[k] * dg_26[k];

        t_29[k] = -f_4 * df_s_17[k]
                  + f_2 * dh_s_29[k]
                  + f_5 * df_17[k]
                  + pb_x[k] * dg_27[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, pb_z, pg_6, df_s_14, df_s_15, dh_s_30, \
                         dh_s_31, dh_s_32, df_14, df_15, dg_28, dg_29, \
                         dg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * pg_6[k]
                  - f_1 * df_s_14[k]
                  + f_2 * dh_s_30[k]
                  + f_3 * df_14[k]
                  + pb_y[k] * dg_28[k];

        t_31[k] = -f_4 * df_s_14[k]
                  + f_2 * dh_s_31[k]
                  + f_5 * df_14[k]
                  + pb_z[k] * dg_29[k];

        t_32[k] = -f_6 * df_s_15[k]
                  + f_2 * dh_s_32[k]
                  + f_0 * df_15[k]
                  + pb_z[k] * dg_30[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_z, pb_y, pb_z, pg_9, ph_3, df_s_17, dh_s_33, \
                         dh_s_34, dh_s_35, df_17, dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * pg_9[k]
                  + f_2 * dh_s_33[k]
                  + pb_y[k] * dg_32[k];

        t_34[k] = -f_1 * df_s_17[k]
                  + f_2 * dh_s_34[k]
                  + f_3 * df_17[k]
                  + pb_z[k] * dg_32[k];

        t_35[k] = pa_z[k] * ph_3[k]
                  + f_2 * dh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_y, pb_z, pg_6, pg_13, pg_14, ph_6, ph_7, \
                         dh_s_36, dh_s_37, dh_s_38, dg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_5 * pg_6[k]
                  + f_2 * dh_s_36[k]
                  + pb_z[k] * dg_33[k];

        t_37[k] = f_7 * pg_13[k]
                  + pa_y[k] * ph_6[k]
                  + f_2 * dh_s_37[k];

        t_38[k] = f_0 * pg_14[k]
                  + pa_y[k] * ph_7[k]
                  + f_2 * dh_s_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_y, pg_15, ph_8, df_s_20, dh_s_39, \
                         dh_s_40, dh_s_41, df_20, dg_37, dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_5 * pg_15[k]
                  + f_2 * dh_s_39[k]
                  + pb_y[k] * dg_37[k];

        t_40[k] = pa_y[k] * ph_8[k]
                  + f_2 * dh_s_40[k];

        t_41[k] = -f_1 * df_s_20[k]
                  + f_2 * dh_s_41[k]
                  + f_3 * df_20[k]
                  + pb_x[k] * dg_38[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_y, df_s_21, df_s_22, dh_s_42, dh_s_43, \
                         dh_s_44, df_21, df_22, dg_38, dg_40, dg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * dh_s_42[k]
                  + pb_y[k] * dg_38[k];

        t_43[k] = -f_8 * df_s_21[k]
                  + f_2 * dh_s_43[k]
                  + f_7 * df_21[k]
                  + pb_x[k] * dg_40[k];

        t_44[k] = -f_6 * df_s_22[k]
                  + f_2 * dh_s_44[k]
                  + f_0 * df_22[k]
                  + pb_x[k] * dg_41[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pb_x, pb_y, df_s_23, df_s_24, dh_s_45, dh_s_46, \
                         dh_s_47, df_23, df_24, dg_40, dg_43, dg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_2 * dh_s_45[k]
                  + pb_y[k] * dg_40[k];

        t_46[k] = -f_6 * df_s_23[k]
                  + f_2 * dh_s_46[k]
                  + f_0 * df_23[k]
                  + pb_x[k] * dg_43[k];

        t_47[k] = -f_4 * df_s_24[k]
                  + f_2 * dh_s_47[k]
                  + f_5 * df_24[k]
                  + pb_x[k] * dg_44[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_x, pb_y, df_s_25, df_s_27, dh_s_48, dh_s_49, \
                         dh_s_50, df_25, df_27, dg_43, dg_45, dg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -f_4 * df_s_25[k]
                  + f_2 * dh_s_48[k]
                  + f_5 * df_25[k]
                  + pb_x[k] * dg_45[k];

        t_49[k] = f_2 * dh_s_49[k]
                  + pb_y[k] * dg_43[k];

        t_50[k] = -f_4 * df_s_27[k]
                  + f_2 * dh_s_50[k]
                  + f_5 * df_27[k]
                  + pb_x[k] * dg_46[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, df_s_24, df_s_25, df_s_26, dh_s_51, dh_s_52, \
                         dh_s_53, df_24, df_25, df_26, dg_47, dg_48, \
                         dg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -f_1 * df_s_24[k]
                  + f_2 * dh_s_51[k]
                  + f_3 * df_24[k]
                  + pb_y[k] * dg_47[k];

        t_52[k] = -f_8 * df_s_25[k]
                  + f_2 * dh_s_52[k]
                  + f_7 * df_25[k]
                  + pb_y[k] * dg_48[k];

        t_53[k] = -f_6 * df_s_26[k]
                  + f_2 * dh_s_53[k]
                  + f_0 * df_26[k]
                  + pb_y[k] * dg_49[k];
    }

#pragma omp simd aligned(t_54, t_55, pb_y, pb_z, pg_15, df_s_27, dh_s_54, dh_s_55, df_27, \
                         dg_50, dg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = -f_4 * df_s_27[k]
                  + f_2 * dh_s_54[k]
                  + f_5 * df_27[k]
                  + pb_y[k] * dg_50[k];

        t_55[k] = f_0 * pg_15[k]
                  - f_1 * df_s_27[k]
                  + f_2 * dh_s_55[k]
                  + f_3 * df_27[k]
                  + pb_z[k] * dg_51[k];
    }
}

auto
compute_prim_dh_kinetic_energy_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_8 = buffer.data(ph + 8);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_11 = buffer.data(ph + 11);
    const auto *ph_13 = buffer.data(ph + 13);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_19 = buffer.data(ph + 19);
    const auto *ph_21 = buffer.data(ph + 21);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_3 = buffer.data(df_s + 3);
    const auto *df_s_4 = buffer.data(df_s + 4);
    const auto *df_s_5 = buffer.data(df_s + 5);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_10 = buffer.data(df_s + 10);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_13 = buffer.data(df_s + 13);
    const auto *df_s_14 = buffer.data(df_s + 14);
    const auto *df_s_15 = buffer.data(df_s + 15);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_19 = buffer.data(df_s + 19);
    const auto *df_s_20 = buffer.data(df_s + 20);
    const auto *df_s_21 = buffer.data(df_s + 21);
    const auto *df_s_22 = buffer.data(df_s + 22);
    const auto *df_s_23 = buffer.data(df_s + 23);
    const auto *df_s_24 = buffer.data(df_s + 24);
    const auto *df_s_25 = buffer.data(df_s + 25);
    const auto *df_s_26 = buffer.data(df_s + 26);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);
    const auto *dh_s_25 = buffer.data(dh_s + 25);
    const auto *dh_s_26 = buffer.data(dh_s + 26);
    const auto *dh_s_27 = buffer.data(dh_s + 27);
    const auto *dh_s_28 = buffer.data(dh_s + 28);
    const auto *dh_s_29 = buffer.data(dh_s + 29);
    const auto *dh_s_30 = buffer.data(dh_s + 30);
    const auto *dh_s_31 = buffer.data(dh_s + 31);
    const auto *dh_s_32 = buffer.data(dh_s + 32);
    const auto *dh_s_33 = buffer.data(dh_s + 33);
    const auto *dh_s_34 = buffer.data(dh_s + 34);
    const auto *dh_s_35 = buffer.data(dh_s + 35);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_37 = buffer.data(dh_s + 37);
    const auto *dh_s_38 = buffer.data(dh_s + 38);
    const auto *dh_s_39 = buffer.data(dh_s + 39);
    const auto *dh_s_40 = buffer.data(dh_s + 40);
    const auto *dh_s_41 = buffer.data(dh_s + 41);
    const auto *dh_s_42 = buffer.data(dh_s + 42);
    const auto *dh_s_43 = buffer.data(dh_s + 43);
    const auto *dh_s_44 = buffer.data(dh_s + 44);
    const auto *dh_s_45 = buffer.data(dh_s + 45);
    const auto *dh_s_46 = buffer.data(dh_s + 46);
    const auto *dh_s_47 = buffer.data(dh_s + 47);
    const auto *dh_s_48 = buffer.data(dh_s + 48);
    const auto *dh_s_49 = buffer.data(dh_s + 49);
    const auto *dh_s_50 = buffer.data(dh_s + 50);
    const auto *dh_s_51 = buffer.data(dh_s + 51);
    const auto *dh_s_52 = buffer.data(dh_s + 52);
    const auto *dh_s_53 = buffer.data(dh_s + 53);
    const auto *dh_s_54 = buffer.data(dh_s + 54);
    const auto *dh_s_55 = buffer.data(dh_s + 55);
    const auto *dh_s_56 = buffer.data(dh_s + 56);
    const auto *dh_s_57 = buffer.data(dh_s + 57);
    const auto *dh_s_58 = buffer.data(dh_s + 58);
    const auto *dh_s_59 = buffer.data(dh_s + 59);
    const auto *dh_s_60 = buffer.data(dh_s + 60);
    const auto *dh_s_61 = buffer.data(dh_s + 61);
    const auto *dh_s_62 = buffer.data(dh_s + 62);
    const auto *dh_s_63 = buffer.data(dh_s + 63);
    const auto *dh_s_64 = buffer.data(dh_s + 64);
    const auto *dh_s_65 = buffer.data(dh_s + 65);
    const auto *dh_s_66 = buffer.data(dh_s + 66);
    const auto *dh_s_67 = buffer.data(dh_s + 67);
    const auto *dh_s_68 = buffer.data(dh_s + 68);
    const auto *dh_s_69 = buffer.data(dh_s + 69);
    const auto *dh_s_70 = buffer.data(dh_s + 70);
    const auto *dh_s_71 = buffer.data(dh_s + 71);
    const auto *dh_s_72 = buffer.data(dh_s + 72);
    const auto *dh_s_73 = buffer.data(dh_s + 73);
    const auto *dh_s_74 = buffer.data(dh_s + 74);
    const auto *dh_s_75 = buffer.data(dh_s + 75);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);
    const auto *df_24 = buffer.data(df + 24);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);
    const auto *dg_39 = buffer.data(dg + 39);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, dh_s_3, df_0, dg_0, dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_2 * dh_s_1[k]
                 + pb_y[k] * dg_0[k];

        t_2[k] = f_2 * dh_s_2[k]
                 + pb_z[k] * dg_0[k];

        t_3[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_3[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, df_s_0, df_s_1, dh_s_4, dh_s_5, dh_s_6, \
                         df_0, df_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_4[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];

        t_5[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_5[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = f_2 * dh_s_6[k]
                 + pb_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, df_s_2, df_s_3, dh_s_7, dh_s_8, dh_s_9, \
                         df_2, df_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * dh_s_7[k]
                 + pb_y[k] * dg_4[k];

        t_8[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_8[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];

        t_9[k] = -f_1 * df_s_3[k]
                 + f_2 * dh_s_9[k]
                 + f_3 * df_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, df_s_4, df_s_5, dh_s_10, dh_s_11, \
                         dh_s_12, df_4, df_5, dg_6, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * df_s_4[k]
                  + f_2 * dh_s_10[k]
                  + f_0 * df_4[k]
                  + pb_y[k] * dg_6[k];

        t_11[k] = -f_4 * df_s_5[k]
                  + f_2 * dh_s_11[k]
                  + f_5 * df_5[k]
                  + pb_y[k] * dg_7[k];

        t_12[k] = -f_1 * df_s_5[k]
                  + f_2 * dh_s_12[k]
                  + f_3 * df_5[k]
                  + pb_z[k] * dg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pa_x, pa_y, pg_1, pg_2, ph_0, ph_1, ph_5, \
                         ph_6, dh_s_13, dh_s_14, dh_s_15, dh_s_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * ph_0[k]
                  + f_2 * dh_s_13[k];

        t_14[k] = f_7 * pg_1[k]
                  + pa_x[k] * ph_5[k]
                  + f_2 * dh_s_14[k];

        t_15[k] = pa_y[k] * ph_1[k]
                  + f_2 * dh_s_15[k];

        t_16[k] = f_0 * pg_2[k]
                  + pa_x[k] * ph_6[k]
                  + f_2 * dh_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, ph_2, ph_8, ph_9, ph_10, dh_s_17, \
                         dh_s_18, dh_s_19, dh_s_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_y[k] * ph_2[k]
                  + f_2 * dh_s_17[k];

        t_18[k] = pa_x[k] * ph_8[k]
                  + f_2 * dh_s_18[k];

        t_19[k] = pa_x[k] * ph_9[k]
                  + f_2 * dh_s_19[k];

        t_20[k] = pa_x[k] * ph_10[k]
                  + f_2 * dh_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_z, pb_z, pg_0, ph_0, ph_11, dh_s_21, \
                         dh_s_22, dh_s_23, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_x[k] * ph_11[k]
                  + f_2 * dh_s_21[k];

        t_22[k] = pa_z[k] * ph_0[k]
                  + f_2 * dh_s_22[k];

        t_23[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_23[k]
                  + pb_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pg_5, pg_6, ph_13, ph_15, ph_17, ph_18, \
                         dh_s_24, dh_s_25, dh_s_26, dh_s_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * pg_5[k]
                  + pa_x[k] * ph_13[k]
                  + f_2 * dh_s_24[k];

        t_25[k] = f_0 * pg_6[k]
                  + pa_x[k] * ph_15[k]
                  + f_2 * dh_s_25[k];

        t_26[k] = pa_x[k] * ph_17[k]
                  + f_2 * dh_s_26[k];

        t_27[k] = pa_x[k] * ph_18[k]
                  + f_2 * dh_s_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pa_x, pb_x, ph_19, ph_21, df_s_9, dh_s_28, dh_s_29, \
                         dh_s_30, df_8, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_x[k] * ph_19[k]
                  + f_2 * dh_s_28[k];

        t_29[k] = pa_x[k] * ph_21[k]
                  + f_2 * dh_s_29[k];

        t_30[k] = -f_1 * df_s_9[k]
                  + f_2 * dh_s_30[k]
                  + f_3 * df_8[k]
                  + pb_x[k] * dg_13[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_x, pb_z, df_s_10, df_s_11, dh_s_31, dh_s_32, \
                         dh_s_33, df_9, df_10, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_8 * df_s_10[k]
                  + f_2 * dh_s_31[k]
                  + f_7 * df_9[k]
                  + pb_x[k] * dg_14[k];

        t_32[k] = -f_6 * df_s_11[k]
                  + f_2 * dh_s_32[k]
                  + f_0 * df_10[k]
                  + pb_x[k] * dg_15[k];

        t_33[k] = f_2 * dh_s_33[k]
                  + pb_z[k] * dg_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pb_x, pb_z, df_s_12, df_s_13, dh_s_34, dh_s_35, \
                         dh_s_36, df_11, df_12, dg_15, dg_16, dg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -f_6 * df_s_12[k]
                  + f_2 * dh_s_34[k]
                  + f_0 * df_11[k]
                  + pb_x[k] * dg_16[k];

        t_35[k] = -f_4 * df_s_13[k]
                  + f_2 * dh_s_35[k]
                  + f_5 * df_12[k]
                  + pb_x[k] * dg_17[k];

        t_36[k] = f_2 * dh_s_36[k]
                  + pb_z[k] * dg_15[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, df_s_15, df_s_16, dh_s_37, dh_s_38, dh_s_39, \
                         df_14, df_15, dg_18, dg_19, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_4 * df_s_15[k]
                  + f_2 * dh_s_37[k]
                  + f_5 * df_14[k]
                  + pb_x[k] * dg_18[k];

        t_38[k] = -f_4 * df_s_16[k]
                  + f_2 * dh_s_38[k]
                  + f_5 * df_15[k]
                  + pb_x[k] * dg_19[k];

        t_39[k] = f_2 * dh_s_39[k]
                  + pb_x[k] * dg_20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, dh_s_40, dh_s_41, dh_s_42, dg_22, dg_23, \
                         dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * dh_s_40[k]
                  + pb_x[k] * dg_22[k];

        t_41[k] = f_2 * dh_s_41[k]
                  + pb_x[k] * dg_23[k];

        t_42[k] = f_2 * dh_s_42[k]
                  + pb_x[k] * dg_24[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_y, pb_z, pg_3, df_s_13, dh_s_43, dh_s_44, \
                         dh_s_45, df_12, dg_20, dg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * pg_3[k]
                  - f_1 * df_s_13[k]
                  + f_2 * dh_s_43[k]
                  + f_3 * df_12[k]
                  + pb_y[k] * dg_20[k];

        t_44[k] = f_2 * dh_s_44[k]
                  + pb_z[k] * dg_20[k];

        t_45[k] = -f_4 * df_s_13[k]
                  + f_2 * dh_s_45[k]
                  + f_5 * df_12[k]
                  + pb_z[k] * dg_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pb_z, pg_4, df_s_14, df_s_16, dh_s_46, \
                         dh_s_47, dh_s_48, df_13, df_15, dg_22, dg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -f_6 * df_s_14[k]
                  + f_2 * dh_s_46[k]
                  + f_0 * df_13[k]
                  + pb_z[k] * dg_22[k];

        t_47[k] = f_0 * pg_4[k]
                  + f_2 * dh_s_47[k]
                  + pb_y[k] * dg_24[k];

        t_48[k] = -f_1 * df_s_16[k]
                  + f_2 * dh_s_48[k]
                  + f_3 * df_15[k]
                  + pb_z[k] * dg_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_z, pb_x, pb_z, pg_3, ph_8, dh_s_49, dh_s_50, \
                         dh_s_51, dg_25, dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * dh_s_49[k]
                  + pb_x[k] * dg_26[k];

        t_50[k] = pa_z[k] * ph_8[k]
                  + f_2 * dh_s_50[k];

        t_51[k] = f_5 * pg_3[k]
                  + f_2 * dh_s_51[k]
                  + pb_z[k] * dg_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pb_y, pg_7, pg_8, pg_9, ph_18, ph_19, \
                         dh_s_52, dh_s_53, dh_s_54, dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_7 * pg_7[k]
                  + pa_y[k] * ph_18[k]
                  + f_2 * dh_s_52[k];

        t_53[k] = f_0 * pg_8[k]
                  + pa_y[k] * ph_19[k]
                  + f_2 * dh_s_53[k];

        t_54[k] = f_5 * pg_9[k]
                  + f_2 * dh_s_54[k]
                  + pb_y[k] * dg_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pb_x, pb_y, ph_21, df_s_19, dh_s_55, dh_s_56, \
                         dh_s_57, df_17, dg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * ph_21[k]
                  + f_2 * dh_s_55[k];

        t_56[k] = -f_1 * df_s_19[k]
                  + f_2 * dh_s_56[k]
                  + f_3 * df_17[k]
                  + pb_x[k] * dg_28[k];

        t_57[k] = f_2 * dh_s_57[k]
                  + pb_y[k] * dg_28[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, df_s_20, df_s_21, dh_s_58, dh_s_59, \
                         dh_s_60, df_18, df_19, dg_29, dg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = -f_8 * df_s_20[k]
                  + f_2 * dh_s_58[k]
                  + f_7 * df_18[k]
                  + pb_x[k] * dg_29[k];

        t_59[k] = -f_6 * df_s_21[k]
                  + f_2 * dh_s_59[k]
                  + f_0 * df_19[k]
                  + pb_x[k] * dg_30[k];

        t_60[k] = f_2 * dh_s_60[k]
                  + pb_y[k] * dg_29[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_x, df_s_22, df_s_23, df_s_24, dh_s_61, dh_s_62, \
                         dh_s_63, df_20, df_21, df_22, dg_31, dg_32, \
                         dg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -f_6 * df_s_22[k]
                  + f_2 * dh_s_61[k]
                  + f_0 * df_20[k]
                  + pb_x[k] * dg_31[k];

        t_62[k] = -f_4 * df_s_23[k]
                  + f_2 * dh_s_62[k]
                  + f_5 * df_21[k]
                  + pb_x[k] * dg_32[k];

        t_63[k] = -f_4 * df_s_24[k]
                  + f_2 * dh_s_63[k]
                  + f_5 * df_22[k]
                  + pb_x[k] * dg_33[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_x, pb_y, df_s_26, dh_s_64, dh_s_65, \
                         dh_s_66, dh_s_67, df_24, dg_31, dg_34, dg_35, \
                         dg_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * dh_s_64[k]
                  + pb_y[k] * dg_31[k];

        t_65[k] = -f_4 * df_s_26[k]
                  + f_2 * dh_s_65[k]
                  + f_5 * df_24[k]
                  + pb_x[k] * dg_34[k];

        t_66[k] = f_2 * dh_s_66[k]
                  + pb_x[k] * dg_35[k];

        t_67[k] = f_2 * dh_s_67[k]
                  + pb_x[k] * dg_36[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_y, df_s_23, dh_s_68, dh_s_69, dh_s_70, \
                         df_21, dg_35, dg_37, dg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * dh_s_68[k]
                  + pb_x[k] * dg_37[k];

        t_69[k] = f_2 * dh_s_69[k]
                  + pb_x[k] * dg_39[k];

        t_70[k] = -f_1 * df_s_23[k]
                  + f_2 * dh_s_70[k]
                  + f_3 * df_21[k]
                  + pb_y[k] * dg_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pb_y, df_s_24, df_s_25, df_s_26, dh_s_71, dh_s_72, \
                         dh_s_73, df_22, df_23, df_24, dg_36, dg_37, \
                         dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = -f_8 * df_s_24[k]
                  + f_2 * dh_s_71[k]
                  + f_7 * df_22[k]
                  + pb_y[k] * dg_36[k];

        t_72[k] = -f_6 * df_s_25[k]
                  + f_2 * dh_s_72[k]
                  + f_0 * df_23[k]
                  + pb_y[k] * dg_37[k];

        t_73[k] = -f_4 * df_s_26[k]
                  + f_2 * dh_s_73[k]
                  + f_5 * df_24[k]
                  + pb_y[k] * dg_38[k];
    }

#pragma omp simd aligned(t_74, t_75, pb_y, pb_z, pg_9, df_s_26, dh_s_74, dh_s_75, df_24, \
                         dg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_2 * dh_s_74[k]
                  + pb_y[k] * dg_39[k];

        t_75[k] = f_0 * pg_9[k]
                  - f_1 * df_s_26[k]
                  + f_2 * dh_s_75[k]
                  + f_3 * df_24[k]
                  + pb_z[k] * dg_39[k];
    }
}

auto
compute_prim_dh_kinetic_energy_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_8 = buffer.data(df_s + 8);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_18 = buffer.data(df_s + 18);
    const auto *df_s_19 = buffer.data(df_s + 19);
    const auto *df_s_20 = buffer.data(df_s + 20);
    const auto *df_s_22 = buffer.data(df_s + 22);
    const auto *df_s_23 = buffer.data(df_s + 23);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_22 = buffer.data(df + 22);
    const auto *df_23 = buffer.data(df + 23);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_33 = buffer.data(dg + 33);
    const auto *dg_34 = buffer.data(dg + 34);
    const auto *dg_36 = buffer.data(dg + 36);
    const auto *dg_37 = buffer.data(dg + 37);
    const auto *dg_38 = buffer.data(dg + 38);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, df_0, dg_0, dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_1[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];

        t_2[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_2[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, pb_z, ph_0, df_s_1, df_s_2, dh_s_3, \
                         dh_s_4, dh_s_5, df_1, df_2, dg_3, dg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_3[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];

        t_4[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_4[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];

        t_5[k] = pa_x[k] * ph_0[k]
                 + f_2 * dh_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pb_x, ph_3, df_s_8, df_s_9, dh_s_6, dh_s_7, \
                         dh_s_8, df_8, df_9, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_x[k] * ph_3[k]
                 + f_2 * dh_s_6[k];

        t_7[k] = -f_1 * df_s_8[k]
                 + f_2 * dh_s_7[k]
                 + f_3 * df_8[k]
                 + pb_x[k] * dg_14[k];

        t_8[k] = -f_6 * df_s_9[k]
                 + f_2 * dh_s_8[k]
                 + f_0 * df_9[k]
                 + pb_x[k] * dg_15[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, pb_z, pg_2, df_s_11, dh_s_9, dh_s_10, \
                         dh_s_11, df_11, dg_17, dg_18, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_4 * df_s_11[k]
                 + f_2 * dh_s_9[k]
                 + f_5 * df_11[k]
                 + pb_x[k] * dg_17[k];

        t_10[k] = f_0 * pg_2[k]
                  - f_1 * df_s_11[k]
                  + f_2 * dh_s_10[k]
                  + f_3 * df_11[k]
                  + pb_y[k] * dg_18[k];

        t_11[k] = -f_4 * df_s_11[k]
                  + f_2 * dh_s_11[k]
                  + f_5 * df_11[k]
                  + pb_z[k] * dg_19[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_z, pg_4, pg_5, ph_1, ph_2, df_s_12, \
                         dh_s_12, dh_s_13, dh_s_14, df_12, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_6 * df_s_12[k]
                  + f_2 * dh_s_12[k]
                  + f_0 * df_12[k]
                  + pb_z[k] * dg_20[k];

        t_13[k] = f_7 * pg_4[k]
                  + pa_y[k] * ph_1[k]
                  + f_2 * dh_s_13[k];

        t_14[k] = f_0 * pg_5[k]
                  + pa_y[k] * ph_2[k]
                  + f_2 * dh_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, pb_x, ph_3, df_s_16, df_s_18, dh_s_15, \
                         dh_s_16, dh_s_17, df_16, df_18, dg_28, dg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * ph_3[k]
                  + f_2 * dh_s_15[k];

        t_16[k] = -f_1 * df_s_16[k]
                  + f_2 * dh_s_16[k]
                  + f_3 * df_16[k]
                  + pb_x[k] * dg_28[k];

        t_17[k] = -f_6 * df_s_18[k]
                  + f_2 * dh_s_17[k]
                  + f_0 * df_18[k]
                  + pb_x[k] * dg_30[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, df_s_19, df_s_20, df_s_23, dh_s_18, dh_s_19, \
                         dh_s_20, df_19, df_20, df_23, dg_31, dg_32, \
                         dg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -f_6 * df_s_19[k]
                  + f_2 * dh_s_18[k]
                  + f_0 * df_19[k]
                  + pb_x[k] * dg_31[k];

        t_19[k] = -f_4 * df_s_20[k]
                  + f_2 * dh_s_19[k]
                  + f_5 * df_20[k]
                  + pb_x[k] * dg_32[k];

        t_20[k] = -f_4 * df_s_23[k]
                  + f_2 * dh_s_20[k]
                  + f_5 * df_23[k]
                  + pb_x[k] * dg_33[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, df_s_20, df_s_22, df_s_23, dh_s_21, dh_s_22, \
                         dh_s_23, df_20, df_22, df_23, dg_34, dg_36, \
                         dg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * df_s_20[k]
                  + f_2 * dh_s_21[k]
                  + f_3 * df_20[k]
                  + pb_y[k] * dg_34[k];

        t_22[k] = -f_6 * df_s_22[k]
                  + f_2 * dh_s_22[k]
                  + f_0 * df_22[k]
                  + pb_y[k] * dg_36[k];

        t_23[k] = -f_4 * df_s_23[k]
                  + f_2 * dh_s_23[k]
                  + f_5 * df_23[k]
                  + pb_y[k] * dg_37[k];
    }

#pragma omp simd aligned(t_24, pb_z, pg_6, df_s_23, dh_s_24, df_23, \
                         dg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_0 * pg_6[k]
                  - f_1 * df_s_23[k]
                  + f_2 * dh_s_24[k]
                  + f_3 * df_23[k]
                  + pb_z[k] * dg_38[k];
    }
}

auto
compute_prim_dh_kinetic_energy_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 3.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_6 = buffer.data(pg + 6);
    const auto *pg_7 = buffer.data(pg + 7);
    const auto *pg_8 = buffer.data(pg + 8);
    const auto *pg_9 = buffer.data(pg + 9);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_4 = buffer.data(ph + 4);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_7 = buffer.data(ph + 7);
    const auto *ph_9 = buffer.data(ph + 9);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_3 = buffer.data(df_s + 3);
    const auto *df_s_4 = buffer.data(df_s + 4);
    const auto *df_s_7 = buffer.data(df_s + 7);
    const auto *df_s_8 = buffer.data(df_s + 8);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_10 = buffer.data(df_s + 10);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_14 = buffer.data(df_s + 14);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_17 = buffer.data(df_s + 17);
    const auto *df_s_18 = buffer.data(df_s + 18);
    const auto *df_s_19 = buffer.data(df_s + 19);
    const auto *df_s_20 = buffer.data(df_s + 20);
    const auto *df_s_21 = buffer.data(df_s + 21);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);
    const auto *dh_s_25 = buffer.data(dh_s + 25);
    const auto *dh_s_26 = buffer.data(dh_s + 26);
    const auto *dh_s_27 = buffer.data(dh_s + 27);
    const auto *dh_s_28 = buffer.data(dh_s + 28);
    const auto *dh_s_29 = buffer.data(dh_s + 29);
    const auto *dh_s_30 = buffer.data(dh_s + 30);
    const auto *dh_s_31 = buffer.data(dh_s + 31);
    const auto *dh_s_32 = buffer.data(dh_s + 32);
    const auto *dh_s_33 = buffer.data(dh_s + 33);
    const auto *dh_s_34 = buffer.data(dh_s + 34);
    const auto *dh_s_35 = buffer.data(dh_s + 35);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_37 = buffer.data(dh_s + 37);
    const auto *dh_s_38 = buffer.data(dh_s + 38);
    const auto *dh_s_39 = buffer.data(dh_s + 39);
    const auto *dh_s_40 = buffer.data(dh_s + 40);
    const auto *dh_s_41 = buffer.data(dh_s + 41);
    const auto *dh_s_42 = buffer.data(dh_s + 42);
    const auto *dh_s_43 = buffer.data(dh_s + 43);
    const auto *dh_s_44 = buffer.data(dh_s + 44);
    const auto *dh_s_45 = buffer.data(dh_s + 45);
    const auto *dh_s_46 = buffer.data(dh_s + 46);
    const auto *dh_s_47 = buffer.data(dh_s + 47);
    const auto *dh_s_48 = buffer.data(dh_s + 48);
    const auto *dh_s_49 = buffer.data(dh_s + 49);
    const auto *dh_s_50 = buffer.data(dh_s + 50);
    const auto *dh_s_51 = buffer.data(dh_s + 51);
    const auto *dh_s_52 = buffer.data(dh_s + 52);
    const auto *dh_s_53 = buffer.data(dh_s + 53);
    const auto *dh_s_54 = buffer.data(dh_s + 54);
    const auto *dh_s_55 = buffer.data(dh_s + 55);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_3 = buffer.data(df + 3);
    const auto *df_4 = buffer.data(df + 4);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);
    const auto *df_21 = buffer.data(df + 21);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_29 = buffer.data(dg + 29);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_31 = buffer.data(dg + 31);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, dh_s_3, df_0, dg_0, dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_2 * dh_s_1[k]
                 + pb_y[k] * dg_0[k];

        t_2[k] = f_2 * dh_s_2[k]
                 + pb_z[k] * dg_0[k];

        t_3[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_3[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_y, pb_z, df_s_0, df_s_1, dh_s_4, dh_s_5, dh_s_6, \
                         df_0, df_1, dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_4[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];

        t_5[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_5[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];

        t_6[k] = f_2 * dh_s_6[k]
                 + pb_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, pb_z, df_s_2, df_s_3, dh_s_7, dh_s_8, dh_s_9, \
                         df_2, df_3, dg_4, dg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * dh_s_7[k]
                 + pb_y[k] * dg_4[k];

        t_8[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_8[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];

        t_9[k] = -f_1 * df_s_3[k]
                 + f_2 * dh_s_9[k]
                 + f_3 * df_3[k]
                 + pb_y[k] * dg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pa_y, pb_z, pg_1, ph_0, ph_1, df_s_4, \
                         dh_s_10, dh_s_11, dh_s_12, df_4, dg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * df_s_4[k]
                  + f_2 * dh_s_10[k]
                  + f_3 * df_4[k]
                  + pb_z[k] * dg_6[k];

        t_11[k] = pa_y[k] * ph_0[k]
                  + f_2 * dh_s_11[k];

        t_12[k] = f_7 * pg_1[k]
                  + pa_x[k] * ph_1[k]
                  + f_2 * dh_s_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pa_z, pg_2, ph_0, ph_2, ph_3, dh_s_13, \
                         dh_s_14, dh_s_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * pg_2[k]
                  + pa_x[k] * ph_2[k]
                  + f_2 * dh_s_13[k];

        t_14[k] = pa_x[k] * ph_3[k]
                  + f_2 * dh_s_14[k];

        t_15[k] = pa_z[k] * ph_0[k]
                  + f_2 * dh_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_z, pg_0, pg_5, pg_6, ph_4, ph_5, dh_s_16, \
                         dh_s_17, dh_s_18, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_16[k]
                  + pb_z[k] * dg_8[k];

        t_17[k] = f_7 * pg_5[k]
                  + pa_x[k] * ph_4[k]
                  + f_2 * dh_s_17[k];

        t_18[k] = f_0 * pg_6[k]
                  + pa_x[k] * ph_5[k]
                  + f_2 * dh_s_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_x, pb_x, ph_9, df_s_7, df_s_8, dh_s_19, dh_s_20, \
                         dh_s_21, df_7, df_8, dg_10, dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_x[k] * ph_9[k]
                  + f_2 * dh_s_19[k];

        t_20[k] = -f_1 * df_s_7[k]
                  + f_2 * dh_s_20[k]
                  + f_3 * df_7[k]
                  + pb_x[k] * dg_10[k];

        t_21[k] = -f_6 * df_s_8[k]
                  + f_2 * dh_s_21[k]
                  + f_0 * df_8[k]
                  + pb_x[k] * dg_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pb_x, pb_z, df_s_9, df_s_10, dh_s_22, dh_s_23, \
                         dh_s_24, df_9, df_10, dg_11, dg_12, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = -f_6 * df_s_9[k]
                  + f_2 * dh_s_22[k]
                  + f_0 * df_9[k]
                  + pb_x[k] * dg_12[k];

        t_23[k] = -f_4 * df_s_10[k]
                  + f_2 * dh_s_23[k]
                  + f_5 * df_10[k]
                  + pb_x[k] * dg_13[k];

        t_24[k] = f_2 * dh_s_24[k]
                  + pb_z[k] * dg_11[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_x, df_s_12, dh_s_25, dh_s_26, dh_s_27, df_12, \
                         dg_14, dg_15, dg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_4 * df_s_12[k]
                  + f_2 * dh_s_25[k]
                  + f_5 * df_12[k]
                  + pb_x[k] * dg_14[k];

        t_26[k] = f_2 * dh_s_26[k]
                  + pb_x[k] * dg_15[k];

        t_27[k] = f_2 * dh_s_27[k]
                  + pb_x[k] * dg_17[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, pb_y, pb_z, pg_3, df_s_10, dh_s_28, dh_s_29, \
                         dh_s_30, df_10, dg_15, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * pg_3[k]
                  - f_1 * df_s_10[k]
                  + f_2 * dh_s_28[k]
                  + f_3 * df_10[k]
                  + pb_y[k] * dg_15[k];

        t_29[k] = f_2 * dh_s_29[k]
                  + pb_z[k] * dg_15[k];

        t_30[k] = -f_4 * df_s_10[k]
                  + f_2 * dh_s_30[k]
                  + f_5 * df_10[k]
                  + pb_z[k] * dg_16[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pb_y, pb_z, pg_4, df_s_11, df_s_12, dh_s_31, \
                         dh_s_32, dh_s_33, df_11, df_12, dg_17, dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = -f_6 * df_s_11[k]
                  + f_2 * dh_s_31[k]
                  + f_0 * df_11[k]
                  + pb_z[k] * dg_17[k];

        t_32[k] = f_0 * pg_4[k]
                  + f_2 * dh_s_32[k]
                  + pb_y[k] * dg_18[k];

        t_33[k] = -f_1 * df_s_12[k]
                  + f_2 * dh_s_33[k]
                  + f_3 * df_12[k]
                  + pb_z[k] * dg_18[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pa_z, pb_z, pg_3, pg_7, ph_3, ph_6, dh_s_34, \
                         dh_s_35, dh_s_36, dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * ph_3[k]
                  + f_2 * dh_s_34[k];

        t_35[k] = f_5 * pg_3[k]
                  + f_2 * dh_s_35[k]
                  + pb_z[k] * dg_19[k];

        t_36[k] = f_7 * pg_7[k]
                  + pa_y[k] * ph_6[k]
                  + f_2 * dh_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_y, pb_y, pg_8, pg_9, ph_7, ph_9, dh_s_37, \
                         dh_s_38, dh_s_39, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_0 * pg_8[k]
                  + pa_y[k] * ph_7[k]
                  + f_2 * dh_s_37[k];

        t_38[k] = f_5 * pg_9[k]
                  + f_2 * dh_s_38[k]
                  + pb_y[k] * dg_20[k];

        t_39[k] = pa_y[k] * ph_9[k]
                  + f_2 * dh_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, pb_y, df_s_14, df_s_16, dh_s_40, dh_s_41, \
                         dh_s_42, df_14, df_16, dg_21, dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -f_1 * df_s_14[k]
                  + f_2 * dh_s_40[k]
                  + f_3 * df_14[k]
                  + pb_x[k] * dg_21[k];

        t_41[k] = f_2 * dh_s_41[k]
                  + pb_y[k] * dg_21[k];

        t_42[k] = -f_6 * df_s_16[k]
                  + f_2 * dh_s_42[k]
                  + f_0 * df_16[k]
                  + pb_x[k] * dg_23[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_y, df_s_17, df_s_18, dh_s_43, dh_s_44, \
                         dh_s_45, df_17, df_18, dg_24, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_6 * df_s_17[k]
                  + f_2 * dh_s_43[k]
                  + f_0 * df_17[k]
                  + pb_x[k] * dg_24[k];

        t_44[k] = -f_4 * df_s_18[k]
                  + f_2 * dh_s_44[k]
                  + f_5 * df_18[k]
                  + pb_x[k] * dg_25[k];

        t_45[k] = f_2 * dh_s_45[k]
                  + pb_y[k] * dg_24[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, df_s_21, dh_s_46, dh_s_47, dh_s_48, \
                         dh_s_49, df_21, dg_26, dg_27, dg_29, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = -f_4 * df_s_21[k]
                  + f_2 * dh_s_46[k]
                  + f_5 * df_21[k]
                  + pb_x[k] * dg_26[k];

        t_47[k] = f_2 * dh_s_47[k]
                  + pb_x[k] * dg_27[k];

        t_48[k] = f_2 * dh_s_48[k]
                  + pb_x[k] * dg_29[k];

        t_49[k] = f_2 * dh_s_49[k]
                  + pb_x[k] * dg_31[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_y, df_s_18, df_s_19, df_s_20, dh_s_50, dh_s_51, \
                         dh_s_52, df_18, df_19, df_20, dg_27, dg_28, \
                         dg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -f_1 * df_s_18[k]
                  + f_2 * dh_s_50[k]
                  + f_3 * df_18[k]
                  + pb_y[k] * dg_27[k];

        t_51[k] = -f_8 * df_s_19[k]
                  + f_2 * dh_s_51[k]
                  + f_7 * df_19[k]
                  + pb_y[k] * dg_28[k];

        t_52[k] = -f_6 * df_s_20[k]
                  + f_2 * dh_s_52[k]
                  + f_0 * df_20[k]
                  + pb_y[k] * dg_29[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pb_z, pg_9, df_s_21, dh_s_53, dh_s_54, \
                         dh_s_55, df_21, dg_30, dg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_4 * df_s_21[k]
                  + f_2 * dh_s_53[k]
                  + f_5 * df_21[k]
                  + pb_y[k] * dg_30[k];

        t_54[k] = f_2 * dh_s_54[k]
                  + pb_y[k] * dg_31[k];

        t_55[k] = f_0 * pg_9[k]
                  - f_1 * df_s_21[k]
                  + f_2 * dh_s_55[k]
                  + f_3 * df_21[k]
                  + pb_z[k] * dg_31[k];
    }
}

auto
compute_prim_dh_kinetic_energy_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_10 = buffer.data(df_s + 10);
    const auto *df_s_18 = buffer.data(df_s + 18);
    const auto *df_s_19 = buffer.data(df_s + 19);
    const auto *df_s_20 = buffer.data(df_s + 20);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_18 = buffer.data(df + 18);
    const auto *df_19 = buffer.data(df + 19);
    const auto *df_20 = buffer.data(df + 20);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_19 = buffer.data(dg + 19);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, pg_0, df_s_0, df_s_2, dh_s_0, dh_s_1, \
                         dh_s_2, df_0, df_2, dg_0, dg_1, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_1[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_1[k];

        t_2[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_2[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_3[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, pg_1, ph_0, ph_3, df_s_9, dh_s_3, dh_s_4, \
                         dh_s_5, df_9, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * ph_0[k]
                 + f_2 * dh_s_3[k];

        t_4[k] = pa_x[k] * ph_3[k]
                 + f_2 * dh_s_4[k];

        t_5[k] = f_0 * pg_1[k]
                 - f_1 * df_s_9[k]
                 + f_2 * dh_s_5[k]
                 + f_3 * df_9[k]
                 + pb_y[k] * dg_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_z, pg_2, ph_1, df_s_9, df_s_10, dh_s_6, \
                         dh_s_7, dh_s_8, df_9, df_10, dg_9, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * df_s_9[k]
                 + f_2 * dh_s_6[k]
                 + f_5 * df_9[k]
                 + pb_z[k] * dg_9[k];

        t_7[k] = -f_6 * df_s_10[k]
                 + f_2 * dh_s_7[k]
                 + f_0 * df_10[k]
                 + pb_z[k] * dg_10[k];

        t_8[k] = f_7 * pg_2[k]
                 + pa_y[k] * ph_1[k]
                 + f_2 * dh_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_y, pg_3, ph_2, ph_3, df_s_18, dh_s_9, \
                         dh_s_10, dh_s_11, df_18, dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pg_3[k]
                 + pa_y[k] * ph_2[k]
                 + f_2 * dh_s_9[k];

        t_10[k] = pa_y[k] * ph_3[k]
                  + f_2 * dh_s_10[k];

        t_11[k] = -f_1 * df_s_18[k]
                  + f_2 * dh_s_11[k]
                  + f_3 * df_18[k]
                  + pb_y[k] * dg_16[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, pg_4, df_s_19, df_s_20, dh_s_12, \
                         dh_s_13, dh_s_14, df_19, df_20, dg_17, dg_18, \
                         dg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_6 * df_s_19[k]
                  + f_2 * dh_s_12[k]
                  + f_0 * df_19[k]
                  + pb_y[k] * dg_17[k];

        t_13[k] = -f_4 * df_s_20[k]
                  + f_2 * dh_s_13[k]
                  + f_5 * df_20[k]
                  + pb_y[k] * dg_18[k];

        t_14[k] = f_0 * pg_4[k]
                  - f_1 * df_s_20[k]
                  + f_2 * dh_s_14[k]
                  + f_3 * df_20[k]
                  + pb_z[k] * dg_19[k];
    }
}

auto
compute_prim_dh_kinetic_energy_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_5 = buffer.data(ph + 5);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_6 = buffer.data(df_s + 6);
    const auto *df_s_7 = buffer.data(df_s + 7);
    const auto *df_s_8 = buffer.data(df_s + 8);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_14 = buffer.data(df_s + 14);
    const auto *df_s_15 = buffer.data(df_s + 15);
    const auto *df_s_16 = buffer.data(df_s + 16);
    const auto *df_s_17 = buffer.data(df_s + 17);
    const auto *df_s_18 = buffer.data(df_s + 18);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);
    const auto *dh_s_22 = buffer.data(dh_s + 22);
    const auto *dh_s_23 = buffer.data(dh_s + 23);
    const auto *dh_s_24 = buffer.data(dh_s + 24);
    const auto *dh_s_25 = buffer.data(dh_s + 25);
    const auto *dh_s_26 = buffer.data(dh_s + 26);
    const auto *dh_s_27 = buffer.data(dh_s + 27);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_2 = buffer.data(df + 2);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_14 = buffer.data(df + 14);
    const auto *df_15 = buffer.data(df + 15);
    const auto *df_16 = buffer.data(df + 16);
    const auto *df_17 = buffer.data(df + 17);
    const auto *df_18 = buffer.data(df + 18);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_4 = buffer.data(dg + 4);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_11 = buffer.data(dg + 11);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_22 = buffer.data(dg + 22);
    const auto *dg_23 = buffer.data(dg + 23);
    const auto *dg_24 = buffer.data(dg + 24);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_26 = buffer.data(dg + 26);
    const auto *dg_27 = buffer.data(dg + 27);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, \
                         dh_s_2, df_0, dg_0, dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_2 * dh_s_1[k]
                 + pb_z[k] * dg_0[k];

        t_2[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_2[k]
                 + f_5 * df_0[k]
                 + pb_y[k] * dg_1[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_y, pb_z, df_s_0, df_s_1, dh_s_3, dh_s_4, df_0, df_1, \
                         dg_2, dg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_3[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_2[k];

        t_4[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_4[k]
                 + f_0 * df_1[k]
                 + pb_y[k] * dg_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pa_x, pb_z, ph_1, ph_5, df_s_2, dh_s_5, dh_s_6, \
                         dh_s_7, df_2, dg_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_5[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_4[k];

        t_6[k] = pa_x[k] * ph_1[k]
                 + f_2 * dh_s_6[k];

        t_7[k] = pa_x[k] * ph_5[k]
                 + f_2 * dh_s_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, df_s_6, df_s_7, df_s_8, dh_s_8, dh_s_9, \
                         dh_s_10, df_6, df_7, df_8, dg_9, dg_10, \
                         dg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = -f_1 * df_s_6[k]
                 + f_2 * dh_s_8[k]
                 + f_3 * df_6[k]
                 + pb_x[k] * dg_9[k];

        t_9[k] = -f_6 * df_s_7[k]
                 + f_2 * dh_s_9[k]
                 + f_0 * df_7[k]
                 + pb_x[k] * dg_10[k];

        t_10[k] = -f_4 * df_s_8[k]
                  + f_2 * dh_s_10[k]
                  + f_5 * df_8[k]
                  + pb_x[k] * dg_11[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_y, pb_z, pg_1, df_s_8, dh_s_11, dh_s_12, \
                         dh_s_13, df_8, dg_12, dg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * pg_1[k]
                  - f_1 * df_s_8[k]
                  + f_2 * dh_s_11[k]
                  + f_3 * df_8[k]
                  + pb_y[k] * dg_12[k];

        t_12[k] = f_2 * dh_s_12[k]
                  + pb_z[k] * dg_12[k];

        t_13[k] = -f_4 * df_s_8[k]
                  + f_2 * dh_s_13[k]
                  + f_5 * df_8[k]
                  + pb_z[k] * dg_13[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_y, pb_z, pg_2, pg_3, ph_2, ph_3, df_s_9, \
                         dh_s_14, dh_s_15, dh_s_16, df_9, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = -f_6 * df_s_9[k]
                  + f_2 * dh_s_14[k]
                  + f_0 * df_9[k]
                  + pb_z[k] * dg_14[k];

        t_15[k] = f_7 * pg_2[k]
                  + pa_y[k] * ph_2[k]
                  + f_2 * dh_s_15[k];

        t_16[k] = f_0 * pg_3[k]
                  + pa_y[k] * ph_3[k]
                  + f_2 * dh_s_16[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pa_y, pb_x, ph_5, df_s_12, df_s_14, dh_s_17, \
                         dh_s_18, dh_s_19, df_12, df_14, dg_18, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_y[k] * ph_5[k]
                  + f_2 * dh_s_17[k];

        t_18[k] = -f_1 * df_s_12[k]
                  + f_2 * dh_s_18[k]
                  + f_3 * df_12[k]
                  + pb_x[k] * dg_18[k];

        t_19[k] = -f_6 * df_s_14[k]
                  + f_2 * dh_s_19[k]
                  + f_0 * df_14[k]
                  + pb_x[k] * dg_20[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, df_s_15, df_s_16, df_s_18, dh_s_20, dh_s_21, \
                         dh_s_22, df_15, df_16, df_18, dg_21, dg_22, \
                         dg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_6 * df_s_15[k]
                  + f_2 * dh_s_20[k]
                  + f_0 * df_15[k]
                  + pb_x[k] * dg_21[k];

        t_21[k] = -f_4 * df_s_16[k]
                  + f_2 * dh_s_21[k]
                  + f_5 * df_16[k]
                  + pb_x[k] * dg_22[k];

        t_22[k] = -f_4 * df_s_18[k]
                  + f_2 * dh_s_22[k]
                  + f_5 * df_18[k]
                  + pb_x[k] * dg_23[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, df_s_16, df_s_17, df_s_18, dh_s_23, dh_s_24, \
                         dh_s_25, df_16, df_17, df_18, dg_24, dg_25, \
                         dg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_1 * df_s_16[k]
                  + f_2 * dh_s_23[k]
                  + f_3 * df_16[k]
                  + pb_y[k] * dg_24[k];

        t_24[k] = -f_6 * df_s_17[k]
                  + f_2 * dh_s_24[k]
                  + f_0 * df_17[k]
                  + pb_y[k] * dg_25[k];

        t_25[k] = -f_4 * df_s_18[k]
                  + f_2 * dh_s_25[k]
                  + f_5 * df_18[k]
                  + pb_y[k] * dg_26[k];
    }

#pragma omp simd aligned(t_26, t_27, pb_y, pb_z, pg_4, df_s_18, dh_s_26, dh_s_27, df_18, \
                         dg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * dh_s_26[k]
                  + pb_y[k] * dg_27[k];

        t_27[k] = f_0 * pg_4[k]
                  - f_1 * df_s_18[k]
                  + f_2 * dh_s_27[k]
                  + f_3 * df_18[k]
                  + pb_z[k] * dg_27[k];
    }
}

auto
compute_prim_dh_kinetic_energy_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_6 = buffer.data(df_s + 6);
    const auto *df_s_7 = buffer.data(df_s + 7);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_13 = buffer.data(df_s + 13);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, pg_0, df_s_0, df_s_1, dh_s_0, dh_s_1, \
                         dh_s_2, df_0, df_1, dg_0, dg_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_1[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_1[k];

        t_2[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_2[k]
                 + f_0 * df_1[k]
                 + pb_z[k] * dg_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_y, pg_1, ph_0, ph_3, df_s_6, dh_s_3, dh_s_4, \
                         dh_s_5, df_6, dg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = pa_x[k] * ph_0[k]
                 + f_2 * dh_s_3[k];

        t_4[k] = pa_x[k] * ph_3[k]
                 + f_2 * dh_s_4[k];

        t_5[k] = f_0 * pg_1[k]
                 - f_1 * df_s_6[k]
                 + f_2 * dh_s_5[k]
                 + f_3 * df_6[k]
                 + pb_y[k] * dg_6[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_y, pb_z, pg_2, ph_1, df_s_6, df_s_7, dh_s_6, \
                         dh_s_7, dh_s_8, df_6, df_7, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * df_s_6[k]
                 + f_2 * dh_s_6[k]
                 + f_5 * df_6[k]
                 + pb_z[k] * dg_7[k];

        t_7[k] = -f_6 * df_s_7[k]
                 + f_2 * dh_s_7[k]
                 + f_0 * df_7[k]
                 + pb_z[k] * dg_8[k];

        t_8[k] = f_7 * pg_2[k]
                 + pa_y[k] * ph_1[k]
                 + f_2 * dh_s_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_y, pg_3, ph_2, ph_3, df_s_11, dh_s_9, \
                         dh_s_10, dh_s_11, df_11, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * pg_3[k]
                 + pa_y[k] * ph_2[k]
                 + f_2 * dh_s_9[k];

        t_10[k] = pa_y[k] * ph_3[k]
                  + f_2 * dh_s_10[k];

        t_11[k] = -f_1 * df_s_11[k]
                  + f_2 * dh_s_11[k]
                  + f_3 * df_11[k]
                  + pb_y[k] * dg_12[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, pg_4, df_s_12, df_s_13, dh_s_12, \
                         dh_s_13, dh_s_14, df_12, df_13, dg_13, dg_14, \
                         dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_6 * df_s_12[k]
                  + f_2 * dh_s_12[k]
                  + f_0 * df_12[k]
                  + pb_y[k] * dg_13[k];

        t_13[k] = -f_4 * df_s_13[k]
                  + f_2 * dh_s_13[k]
                  + f_5 * df_13[k]
                  + pb_y[k] * dg_14[k];

        t_14[k] = f_0 * pg_4[k]
                  - f_1 * df_s_13[k]
                  + f_2 * dh_s_14[k]
                  + f_3 * df_13[k]
                  + pb_z[k] * dg_15[k];
    }
}

auto
compute_prim_dh_kinetic_energy_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t pg, const size_t ph,
                                 const size_t df_s, const size_t dh_s, const size_t df,
                                 const size_t dg, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pg_0 = buffer.data(pg + 0);
    const auto *pg_1 = buffer.data(pg + 1);
    const auto *pg_2 = buffer.data(pg + 2);
    const auto *pg_3 = buffer.data(pg + 3);
    const auto *pg_4 = buffer.data(pg + 4);

    const auto *ph_1 = buffer.data(ph + 1);
    const auto *ph_2 = buffer.data(ph + 2);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_5 = buffer.data(ph + 5);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_5 = buffer.data(df_s + 5);
    const auto *df_s_6 = buffer.data(df_s + 6);
    const auto *df_s_7 = buffer.data(df_s + 7);
    const auto *df_s_10 = buffer.data(df_s + 10);
    const auto *df_s_11 = buffer.data(df_s + 11);
    const auto *df_s_12 = buffer.data(df_s + 12);
    const auto *df_s_13 = buffer.data(df_s + 13);

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_1 = buffer.data(dh_s + 1);
    const auto *dh_s_2 = buffer.data(dh_s + 2);
    const auto *dh_s_3 = buffer.data(dh_s + 3);
    const auto *dh_s_4 = buffer.data(dh_s + 4);
    const auto *dh_s_5 = buffer.data(dh_s + 5);
    const auto *dh_s_6 = buffer.data(dh_s + 6);
    const auto *dh_s_7 = buffer.data(dh_s + 7);
    const auto *dh_s_8 = buffer.data(dh_s + 8);
    const auto *dh_s_9 = buffer.data(dh_s + 9);
    const auto *dh_s_10 = buffer.data(dh_s + 10);
    const auto *dh_s_11 = buffer.data(dh_s + 11);
    const auto *dh_s_12 = buffer.data(dh_s + 12);
    const auto *dh_s_13 = buffer.data(dh_s + 13);
    const auto *dh_s_14 = buffer.data(dh_s + 14);
    const auto *dh_s_15 = buffer.data(dh_s + 15);
    const auto *dh_s_16 = buffer.data(dh_s + 16);
    const auto *dh_s_17 = buffer.data(dh_s + 17);
    const auto *dh_s_18 = buffer.data(dh_s + 18);
    const auto *dh_s_19 = buffer.data(dh_s + 19);
    const auto *dh_s_20 = buffer.data(dh_s + 20);
    const auto *dh_s_21 = buffer.data(dh_s + 21);

    const auto *df_0 = buffer.data(df + 0);
    const auto *df_1 = buffer.data(df + 1);
    const auto *df_5 = buffer.data(df + 5);
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_7 = buffer.data(df + 7);
    const auto *df_10 = buffer.data(df + 10);
    const auto *df_11 = buffer.data(df + 11);
    const auto *df_12 = buffer.data(df + 12);
    const auto *df_13 = buffer.data(df + 13);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_7 = buffer.data(dg + 7);
    const auto *dg_8 = buffer.data(dg + 8);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_17 = buffer.data(dg + 17);
    const auto *dg_18 = buffer.data(dg + 18);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_z, pg_0, df_s_0, dh_s_0, dh_s_1, dh_s_2, \
                         df_0, dg_0, dg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pg_0[k]
                 - f_1 * df_s_0[k]
                 + f_2 * dh_s_0[k]
                 + f_3 * df_0[k]
                 + pb_x[k] * dg_0[k];

        t_1[k] = f_2 * dh_s_1[k]
                 + pb_z[k] * dg_0[k];

        t_2[k] = -f_4 * df_s_0[k]
                 + f_2 * dh_s_2[k]
                 + f_5 * df_0[k]
                 + pb_z[k] * dg_1[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, pb_z, ph_1, ph_5, df_s_1, dh_s_3, dh_s_4, \
                         dh_s_5, df_1, dg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_6 * df_s_1[k]
                 + f_2 * dh_s_3[k]
                 + f_0 * df_1[k]
                 + pb_z[k] * dg_2[k];

        t_4[k] = pa_x[k] * ph_1[k]
                 + f_2 * dh_s_4[k];

        t_5[k] = pa_x[k] * ph_5[k]
                 + f_2 * dh_s_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pg_1, df_s_5, df_s_6, dh_s_6, dh_s_7, \
                         dh_s_8, df_5, df_6, dg_6, dg_7, dg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_6 * df_s_5[k]
                 + f_2 * dh_s_6[k]
                 + f_0 * df_5[k]
                 + pb_x[k] * dg_6[k];

        t_7[k] = -f_4 * df_s_6[k]
                 + f_2 * dh_s_7[k]
                 + f_5 * df_6[k]
                 + pb_x[k] * dg_7[k];

        t_8[k] = f_0 * pg_1[k]
                 - f_1 * df_s_6[k]
                 + f_2 * dh_s_8[k]
                 + f_3 * df_6[k]
                 + pb_y[k] * dg_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_z, df_s_6, df_s_7, dh_s_9, dh_s_10, dh_s_11, \
                         df_6, df_7, dg_8, dg_9, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * dh_s_9[k]
                 + pb_z[k] * dg_8[k];

        t_10[k] = -f_4 * df_s_6[k]
                  + f_2 * dh_s_10[k]
                  + f_5 * df_6[k]
                  + pb_z[k] * dg_9[k];

        t_11[k] = -f_6 * df_s_7[k]
                  + f_2 * dh_s_11[k]
                  + f_0 * df_7[k]
                  + pb_z[k] * dg_10[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pg_2, pg_3, ph_2, ph_3, ph_5, dh_s_12, \
                         dh_s_13, dh_s_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * pg_2[k]
                  + pa_y[k] * ph_2[k]
                  + f_2 * dh_s_12[k];

        t_13[k] = f_0 * pg_3[k]
                  + pa_y[k] * ph_3[k]
                  + f_2 * dh_s_13[k];

        t_14[k] = pa_y[k] * ph_5[k]
                  + f_2 * dh_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, df_s_10, df_s_13, dh_s_15, dh_s_16, df_10, df_13, \
                         dg_13, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -f_6 * df_s_10[k]
                  + f_2 * dh_s_15[k]
                  + f_0 * df_10[k]
                  + pb_x[k] * dg_13[k];

        t_16[k] = -f_4 * df_s_13[k]
                  + f_2 * dh_s_16[k]
                  + f_5 * df_13[k]
                  + pb_x[k] * dg_14[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, df_s_11, df_s_12, df_s_13, dh_s_17, dh_s_18, \
                         dh_s_19, df_11, df_12, df_13, dg_15, dg_16, \
                         dg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_1 * df_s_11[k]
                  + f_2 * dh_s_17[k]
                  + f_3 * df_11[k]
                  + pb_y[k] * dg_15[k];

        t_18[k] = -f_6 * df_s_12[k]
                  + f_2 * dh_s_18[k]
                  + f_0 * df_12[k]
                  + pb_y[k] * dg_16[k];

        t_19[k] = -f_4 * df_s_13[k]
                  + f_2 * dh_s_19[k]
                  + f_5 * df_13[k]
                  + pb_y[k] * dg_17[k];
    }

#pragma omp simd aligned(t_20, t_21, pb_y, pb_z, pg_4, df_s_13, dh_s_20, dh_s_21, df_13, \
                         dg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * dh_s_20[k]
                  + pb_y[k] * dg_18[k];

        t_21[k] = f_0 * pg_4[k]
                  - f_1 * df_s_13[k]
                  + f_2 * dh_s_21[k]
                  + f_3 * df_13[k]
                  + pb_z[k] * dg_18[k];
    }
}

}  // namespace simdkin
