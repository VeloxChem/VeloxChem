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
    const auto *pg_5 = buffer.data(pg + 5);
    const auto *pg_10 = buffer.data(pg + 10);
    const auto *pg_12 = buffer.data(pg + 12);
    const auto *pg_14 = buffer.data(pg + 14);
    const auto *pg_18 = buffer.data(pg + 18);
    const auto *pg_21 = buffer.data(pg + 21);
    const auto *pg_25 = buffer.data(pg + 25);
    const auto *pg_27 = buffer.data(pg + 27);
    const auto *pg_28 = buffer.data(pg + 28);
    const auto *pg_29 = buffer.data(pg + 29);
    const auto *pg_32 = buffer.data(pg + 32);
    const auto *pg_34 = buffer.data(pg + 34);
    const auto *pg_35 = buffer.data(pg + 35);
    const auto *pg_37 = buffer.data(pg + 37);
    const auto *pg_39 = buffer.data(pg + 39);
    const auto *pg_41 = buffer.data(pg + 41);
    const auto *pg_42 = buffer.data(pg + 42);
    const auto *pg_43 = buffer.data(pg + 43);
    const auto *pg_44 = buffer.data(pg + 44);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_10 = buffer.data(ph + 10);
    const auto *ph_14 = buffer.data(ph + 14);
    const auto *ph_22 = buffer.data(ph + 22);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_41 = buffer.data(ph + 41);
    const auto *ph_42 = buffer.data(ph + 42);
    const auto *ph_44 = buffer.data(ph + 44);
    const auto *ph_46 = buffer.data(ph + 46);
    const auto *ph_47 = buffer.data(ph + 47);
    const auto *ph_49 = buffer.data(ph + 49);
    const auto *ph_50 = buffer.data(ph + 50);
    const auto *ph_51 = buffer.data(ph + 51);
    const auto *ph_57 = buffer.data(ph + 57);
    const auto *ph_58 = buffer.data(ph + 58);
    const auto *ph_59 = buffer.data(ph + 59);
    const auto *ph_60 = buffer.data(ph + 60);
    const auto *ph_62 = buffer.data(ph + 62);

    const auto *df_s_0 = buffer.data(df_s + 0);
    const auto *df_s_1 = buffer.data(df_s + 1);
    const auto *df_s_2 = buffer.data(df_s + 2);
    const auto *df_s_6 = buffer.data(df_s + 6);
    const auto *df_s_8 = buffer.data(df_s + 8);
    const auto *df_s_9 = buffer.data(df_s + 9);
    const auto *df_s_30 = buffer.data(df_s + 30);
    const auto *df_s_31 = buffer.data(df_s + 31);
    const auto *df_s_33 = buffer.data(df_s + 33);
    const auto *df_s_35 = buffer.data(df_s + 35);
    const auto *df_s_36 = buffer.data(df_s + 36);
    const auto *df_s_37 = buffer.data(df_s + 37);
    const auto *df_s_38 = buffer.data(df_s + 38);
    const auto *df_s_39 = buffer.data(df_s + 39);
    const auto *df_s_50 = buffer.data(df_s + 50);
    const auto *df_s_52 = buffer.data(df_s + 52);
    const auto *df_s_53 = buffer.data(df_s + 53);
    const auto *df_s_55 = buffer.data(df_s + 55);
    const auto *df_s_56 = buffer.data(df_s + 56);
    const auto *df_s_57 = buffer.data(df_s + 57);
    const auto *df_s_58 = buffer.data(df_s + 58);
    const auto *df_s_59 = buffer.data(df_s + 59);

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
    const auto *df_6 = buffer.data(df + 6);
    const auto *df_8 = buffer.data(df + 8);
    const auto *df_9 = buffer.data(df + 9);
    const auto *df_30 = buffer.data(df + 30);
    const auto *df_31 = buffer.data(df + 31);
    const auto *df_33 = buffer.data(df + 33);
    const auto *df_35 = buffer.data(df + 35);
    const auto *df_36 = buffer.data(df + 36);
    const auto *df_37 = buffer.data(df + 37);
    const auto *df_38 = buffer.data(df + 38);
    const auto *df_39 = buffer.data(df + 39);
    const auto *df_50 = buffer.data(df + 50);
    const auto *df_52 = buffer.data(df + 52);
    const auto *df_53 = buffer.data(df + 53);
    const auto *df_55 = buffer.data(df + 55);
    const auto *df_56 = buffer.data(df + 56);
    const auto *df_57 = buffer.data(df + 57);
    const auto *df_58 = buffer.data(df + 58);
    const auto *df_59 = buffer.data(df + 59);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_1 = buffer.data(dg + 1);
    const auto *dg_2 = buffer.data(dg + 2);
    const auto *dg_3 = buffer.data(dg + 3);
    const auto *dg_5 = buffer.data(dg + 5);
    const auto *dg_6 = buffer.data(dg + 6);
    const auto *dg_9 = buffer.data(dg + 9);
    const auto *dg_10 = buffer.data(dg + 10);
    const auto *dg_12 = buffer.data(dg + 12);
    const auto *dg_13 = buffer.data(dg + 13);
    const auto *dg_14 = buffer.data(dg + 14);
    const auto *dg_15 = buffer.data(dg + 15);
    const auto *dg_16 = buffer.data(dg + 16);
    const auto *dg_18 = buffer.data(dg + 18);
    const auto *dg_20 = buffer.data(dg + 20);
    const auto *dg_21 = buffer.data(dg + 21);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_27 = buffer.data(dg + 27);
    const auto *dg_28 = buffer.data(dg + 28);
    const auto *dg_30 = buffer.data(dg + 30);
    const auto *dg_32 = buffer.data(dg + 32);
    const auto *dg_35 = buffer.data(dg + 35);
    const auto *dg_39 = buffer.data(dg + 39);
    const auto *dg_41 = buffer.data(dg + 41);
    const auto *dg_42 = buffer.data(dg + 42);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_45 = buffer.data(dg + 45);
    const auto *dg_46 = buffer.data(dg + 46);
    const auto *dg_48 = buffer.data(dg + 48);
    const auto *dg_50 = buffer.data(dg + 50);
    const auto *dg_51 = buffer.data(dg + 51);
    const auto *dg_53 = buffer.data(dg + 53);
    const auto *dg_54 = buffer.data(dg + 54);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_56 = buffer.data(dg + 56);
    const auto *dg_57 = buffer.data(dg + 57);
    const auto *dg_58 = buffer.data(dg + 58);
    const auto *dg_59 = buffer.data(dg + 59);
    const auto *dg_70 = buffer.data(dg + 70);
    const auto *dg_71 = buffer.data(dg + 71);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_73 = buffer.data(dg + 73);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_75 = buffer.data(dg + 75);
    const auto *dg_77 = buffer.data(dg + 77);
    const auto *dg_78 = buffer.data(dg + 78);
    const auto *dg_80 = buffer.data(dg + 80);
    const auto *dg_81 = buffer.data(dg + 81);
    const auto *dg_82 = buffer.data(dg + 82);
    const auto *dg_84 = buffer.data(dg + 84);
    const auto *dg_85 = buffer.data(dg + 85);
    const auto *dg_86 = buffer.data(dg + 86);
    const auto *dg_87 = buffer.data(dg + 87);
    const auto *dg_88 = buffer.data(dg + 88);
    const auto *dg_89 = buffer.data(dg + 89);

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

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, pg_10, df_s_2, dh_s_8, dh_s_9, \
                         dh_s_10, df_2, dg_5, dg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * dh_s_8[k]
                 + pb_y[k] * dg_5[k];

        t_9[k] = -f_6 * df_s_2[k]
                 + f_2 * dh_s_9[k]
                 + f_0 * df_2[k]
                 + pb_z[k] * dg_5[k];

        t_10[k] = f_0 * pg_10[k]
                  + f_2 * dh_s_10[k]
                  + pb_x[k] * dg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, pg_12, dh_s_11, dh_s_12, dh_s_13, \
                         dg_6, dg_9, dg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * dh_s_11[k]
                  + pb_z[k] * dg_6[k];

        t_12[k] = f_0 * pg_12[k]
                  + f_2 * dh_s_12[k]
                  + pb_x[k] * dg_12[k];

        t_13[k] = f_2 * dh_s_13[k]
                  + pb_y[k] * dg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, pg_14, df_s_6, dh_s_14, dh_s_15, \
                         dh_s_16, df_6, dg_10, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * pg_14[k]
                  + f_2 * dh_s_14[k]
                  + pb_x[k] * dg_14[k];

        t_15[k] = -f_1 * df_s_6[k]
                  + f_2 * dh_s_15[k]
                  + f_3 * df_6[k]
                  + pb_y[k] * dg_10[k];

        t_16[k] = f_2 * dh_s_16[k]
                  + pb_z[k] * dg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, df_s_8, df_s_9, dh_s_17, dh_s_18, dh_s_19, \
                         df_8, df_9, dg_12, dg_13, dg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_6 * df_s_8[k]
                  + f_2 * dh_s_17[k]
                  + f_0 * df_8[k]
                  + pb_y[k] * dg_12[k];

        t_18[k] = -f_4 * df_s_9[k]
                  + f_2 * dh_s_18[k]
                  + f_5 * df_9[k]
                  + pb_y[k] * dg_13[k];

        t_19[k] = f_2 * dh_s_19[k]
                  + pb_y[k] * dg_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, pg_0, ph_0, df_s_9, dh_s_20, \
                         dh_s_21, dh_s_22, df_9, dg_14, dg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * df_s_9[k]
                  + f_2 * dh_s_20[k]
                  + f_3 * df_9[k]
                  + pb_z[k] * dg_14[k];

        t_21[k] = pa_y[k] * ph_0[k]
                  + f_2 * dh_s_21[k];

        t_22[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_22[k]
                  + pb_y[k] * dg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, pb_z, pg_18, ph_5, ph_24, \
                         dh_s_23, dh_s_24, dh_s_25, dh_s_26, dg_15, \
                         dg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * dh_s_23[k]
                  + pb_z[k] * dg_15[k];

        t_24[k] = f_7 * pg_18[k]
                  + pa_x[k] * ph_24[k]
                  + f_2 * dh_s_24[k];

        t_25[k] = f_2 * dh_s_25[k]
                  + pb_z[k] * dg_16[k];

        t_26[k] = pa_y[k] * ph_5[k]
                  + f_2 * dh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_y, pb_z, pg_5, pg_21, ph_27, dh_s_27, \
                         dh_s_28, dh_s_29, dg_18, dg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pg_21[k]
                  + pa_x[k] * ph_27[k]
                  + f_2 * dh_s_27[k];

        t_28[k] = f_2 * dh_s_28[k]
                  + pb_z[k] * dg_18[k];

        t_29[k] = f_5 * pg_5[k]
                  + f_2 * dh_s_29[k]
                  + pb_y[k] * dg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, pg_25, ph_9, dh_s_30, dh_s_31, \
                         dh_s_32, dg_21, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * ph_9[k]
                  + f_2 * dh_s_30[k];

        t_31[k] = f_5 * pg_25[k]
                  + f_2 * dh_s_31[k]
                  + pb_x[k] * dg_25[k];

        t_32[k] = f_2 * dh_s_32[k]
                  + pb_z[k] * dg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, pg_27, pg_28, ph_14, dh_s_33, dh_s_34, \
                         dh_s_35, dg_27, dg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_5 * pg_27[k]
                  + f_2 * dh_s_33[k]
                  + pb_x[k] * dg_27[k];

        t_34[k] = f_5 * pg_28[k]
                  + f_2 * dh_s_34[k]
                  + pb_x[k] * dg_28[k];

        t_35[k] = pa_y[k] * ph_14[k]
                  + f_2 * dh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_z, ph_36, ph_38, ph_39, dh_s_36, \
                         dh_s_37, dh_s_38, dh_s_39, dg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * ph_36[k]
                  + f_2 * dh_s_36[k];

        t_37[k] = f_2 * dh_s_37[k]
                  + pb_z[k] * dg_25[k];

        t_38[k] = pa_x[k] * ph_38[k]
                  + f_2 * dh_s_38[k];

        t_39[k] = pa_x[k] * ph_39[k]
                  + f_2 * dh_s_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_z, pb_y, ph_0, ph_40, ph_41, \
                         dh_s_40, dh_s_41, dh_s_42, dh_s_43, dg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_x[k] * ph_40[k]
                  + f_2 * dh_s_40[k];

        t_41[k] = pa_x[k] * ph_41[k]
                  + f_2 * dh_s_41[k];

        t_42[k] = pa_z[k] * ph_0[k]
                  + f_2 * dh_s_42[k];

        t_43[k] = f_2 * dh_s_43[k]
                  + pb_y[k] * dg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pa_z, pb_y, pb_z, pg_0, ph_3, dh_s_44, dh_s_45, \
                         dh_s_46, dg_30, dg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * pg_0[k]
                  + f_2 * dh_s_44[k]
                  + pb_z[k] * dg_30[k];

        t_45[k] = pa_z[k] * ph_3[k]
                  + f_2 * dh_s_45[k];

        t_46[k] = f_2 * dh_s_46[k]
                  + pb_y[k] * dg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_z, pg_35, pg_37, ph_6, ph_47, ph_49, \
                         dh_s_47, dh_s_48, dh_s_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_7 * pg_35[k]
                  + pa_x[k] * ph_47[k]
                  + f_2 * dh_s_47[k];

        t_48[k] = pa_z[k] * ph_6[k]
                  + f_2 * dh_s_48[k];

        t_49[k] = f_0 * pg_37[k]
                  + pa_x[k] * ph_49[k]
                  + f_2 * dh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pa_z, pb_y, pg_39, ph_10, ph_51, dh_s_50, \
                         dh_s_51, dh_s_52, dg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * dh_s_50[k]
                  + pb_y[k] * dg_35[k];

        t_51[k] = f_0 * pg_39[k]
                  + pa_x[k] * ph_51[k]
                  + f_2 * dh_s_51[k];

        t_52[k] = pa_z[k] * ph_10[k]
                  + f_2 * dh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, pg_41, pg_42, dh_s_53, dh_s_54, \
                         dh_s_55, dg_39, dg_41, dg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_5 * pg_41[k]
                  + f_2 * dh_s_53[k]
                  + pb_x[k] * dg_41[k];

        t_54[k] = f_5 * pg_42[k]
                  + f_2 * dh_s_54[k]
                  + pb_x[k] * dg_42[k];

        t_55[k] = f_2 * dh_s_55[k]
                  + pb_y[k] * dg_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_x, pg_44, ph_57, ph_58, ph_59, \
                         dh_s_56, dh_s_57, dh_s_58, dh_s_59, dg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_5 * pg_44[k]
                  + f_2 * dh_s_56[k]
                  + pb_x[k] * dg_44[k];

        t_57[k] = pa_x[k] * ph_57[k]
                  + f_2 * dh_s_57[k];

        t_58[k] = pa_x[k] * ph_58[k]
                  + f_2 * dh_s_58[k];

        t_59[k] = pa_x[k] * ph_59[k]
                  + f_2 * dh_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_x, pb_y, ph_60, ph_62, dh_s_60, dh_s_61, \
                         dh_s_62, dg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * ph_60[k]
                  + f_2 * dh_s_60[k];

        t_61[k] = f_2 * dh_s_61[k]
                  + pb_y[k] * dg_44[k];

        t_62[k] = pa_x[k] * ph_62[k]
                  + f_2 * dh_s_62[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_x, pb_z, df_s_30, df_s_31, dh_s_63, dh_s_64, \
                         dh_s_65, df_30, df_31, dg_45, dg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -f_1 * df_s_30[k]
                  + f_2 * dh_s_63[k]
                  + f_3 * df_30[k]
                  + pb_x[k] * dg_45[k];

        t_64[k] = -f_8 * df_s_31[k]
                  + f_2 * dh_s_64[k]
                  + f_7 * df_31[k]
                  + pb_x[k] * dg_46[k];

        t_65[k] = f_2 * dh_s_65[k]
                  + pb_z[k] * dg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, df_s_33, df_s_35, dh_s_66, dh_s_67, \
                         dh_s_68, df_33, df_35, dg_46, dg_48, dg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = -f_6 * df_s_33[k]
                  + f_2 * dh_s_66[k]
                  + f_0 * df_33[k]
                  + pb_x[k] * dg_48[k];

        t_67[k] = f_2 * dh_s_67[k]
                  + pb_z[k] * dg_46[k];

        t_68[k] = -f_6 * df_s_35[k]
                  + f_2 * dh_s_68[k]
                  + f_0 * df_35[k]
                  + pb_x[k] * dg_50[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_z, df_s_36, df_s_38, dh_s_69, dh_s_70, \
                         dh_s_71, df_36, df_38, dg_48, dg_51, dg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = -f_4 * df_s_36[k]
                  + f_2 * dh_s_69[k]
                  + f_5 * df_36[k]
                  + pb_x[k] * dg_51[k];

        t_70[k] = f_2 * dh_s_70[k]
                  + pb_z[k] * dg_48[k];

        t_71[k] = -f_4 * df_s_38[k]
                  + f_2 * dh_s_71[k]
                  + f_5 * df_38[k]
                  + pb_x[k] * dg_53[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pb_x, df_s_39, dh_s_72, dh_s_73, dh_s_74, \
                         dh_s_75, df_39, dg_54, dg_55, dg_56, dg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = -f_4 * df_s_39[k]
                  + f_2 * dh_s_72[k]
                  + f_5 * df_39[k]
                  + pb_x[k] * dg_54[k];

        t_73[k] = f_2 * dh_s_73[k]
                  + pb_x[k] * dg_55[k];

        t_74[k] = f_2 * dh_s_74[k]
                  + pb_x[k] * dg_56[k];

        t_75[k] = f_2 * dh_s_75[k]
                  + pb_x[k] * dg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_x, pb_y, pg_25, df_s_36, dh_s_76, dh_s_77, \
                         dh_s_78, df_36, dg_55, dg_58, dg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_2 * dh_s_76[k]
                  + pb_x[k] * dg_58[k];

        t_77[k] = f_2 * dh_s_77[k]
                  + pb_x[k] * dg_59[k];

        t_78[k] = f_0 * pg_25[k]
                  - f_1 * df_s_36[k]
                  + f_2 * dh_s_78[k]
                  + f_3 * df_36[k]
                  + pb_y[k] * dg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_z, df_s_36, df_s_37, dh_s_79, dh_s_80, dh_s_81, \
                         df_36, df_37, dg_55, dg_56, dg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_2 * dh_s_79[k]
                  + pb_z[k] * dg_55[k];

        t_80[k] = -f_4 * df_s_36[k]
                  + f_2 * dh_s_80[k]
                  + f_5 * df_36[k]
                  + pb_z[k] * dg_56[k];

        t_81[k] = -f_6 * df_s_37[k]
                  + f_2 * dh_s_81[k]
                  + f_0 * df_37[k]
                  + pb_z[k] * dg_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pb_y, pb_z, pg_29, ph_42, df_s_39, dh_s_82, \
                         dh_s_83, dh_s_84, df_39, dg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * pg_29[k]
                  + f_2 * dh_s_82[k]
                  + pb_y[k] * dg_59[k];

        t_83[k] = -f_1 * df_s_39[k]
                  + f_2 * dh_s_83[k]
                  + f_3 * df_39[k]
                  + pb_z[k] * dg_59[k];

        t_84[k] = pa_y[k] * ph_42[k]
                  + f_2 * dh_s_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pa_z, pg_32, ph_22, ph_24, ph_44, \
                         ph_46, dh_s_85, dh_s_86, dh_s_87, dh_s_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_z[k] * ph_22[k]
                  + f_2 * dh_s_85[k];

        t_86[k] = pa_y[k] * ph_44[k]
                  + f_2 * dh_s_86[k];

        t_87[k] = pa_z[k] * ph_24[k]
                  + f_2 * dh_s_87[k];

        t_88[k] = f_5 * pg_32[k]
                  + pa_y[k] * ph_46[k]
                  + f_2 * dh_s_88[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pa_z, pg_34, pg_35, ph_27, ph_47, \
                         ph_49, ph_50, dh_s_89, dh_s_90, dh_s_91, \
                         dh_s_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * ph_47[k]
                  + f_2 * dh_s_89[k];

        t_90[k] = pa_z[k] * ph_27[k]
                  + f_2 * dh_s_90[k];

        t_91[k] = f_0 * pg_34[k]
                  + pa_y[k] * ph_49[k]
                  + f_2 * dh_s_91[k];

        t_92[k] = f_5 * pg_35[k]
                  + pa_y[k] * ph_50[k]
                  + f_2 * dh_s_92[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_x, ph_51, dh_s_93, dh_s_94, dh_s_95, \
                         dh_s_96, dg_70, dg_71, dg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * ph_51[k]
                  + f_2 * dh_s_93[k];

        t_94[k] = f_2 * dh_s_94[k]
                  + pb_x[k] * dg_70[k];

        t_95[k] = f_2 * dh_s_95[k]
                  + pb_x[k] * dg_71[k];

        t_96[k] = f_2 * dh_s_96[k]
                  + pb_x[k] * dg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_z, pb_x, pb_z, pg_25, ph_36, dh_s_97, \
                         dh_s_98, dh_s_99, dh_s_100, dg_70, dg_73, \
                         dg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * dh_s_97[k]
                  + pb_x[k] * dg_73[k];

        t_98[k] = f_2 * dh_s_98[k]
                  + pb_x[k] * dg_74[k];

        t_99[k] = pa_z[k] * ph_36[k]
                  + f_2 * dh_s_99[k];

        t_100[k] = f_5 * pg_25[k]
                   + f_2 * dh_s_100[k]
                   + pb_z[k] * dg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_y, pb_y, pg_42, pg_43, pg_44, ph_59, ph_60, \
                         dh_s_101, dh_s_102, dh_s_103, dg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_7 * pg_42[k]
                   + pa_y[k] * ph_59[k]
                   + f_2 * dh_s_101[k];

        t_102[k] = f_0 * pg_43[k]
                   + pa_y[k] * ph_60[k]
                   + f_2 * dh_s_102[k];

        t_103[k] = f_5 * pg_44[k]
                   + f_2 * dh_s_103[k]
                   + pb_y[k] * dg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pb_x, pb_y, ph_62, df_s_50, dh_s_104, \
                         dh_s_105, dh_s_106, df_50, dg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * ph_62[k]
                   + f_2 * dh_s_104[k];

        t_105[k] = -f_1 * df_s_50[k]
                   + f_2 * dh_s_105[k]
                   + f_3 * df_50[k]
                   + pb_x[k] * dg_75[k];

        t_106[k] = f_2 * dh_s_106[k]
                   + pb_y[k] * dg_75[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_x, pb_y, df_s_52, df_s_53, dh_s_107, \
                         dh_s_108, dh_s_109, df_52, df_53, dg_77, \
                         dg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_8 * df_s_52[k]
                   + f_2 * dh_s_107[k]
                   + f_7 * df_52[k]
                   + pb_x[k] * dg_77[k];

        t_108[k] = -f_6 * df_s_53[k]
                   + f_2 * dh_s_108[k]
                   + f_0 * df_53[k]
                   + pb_x[k] * dg_78[k];

        t_109[k] = f_2 * dh_s_109[k]
                   + pb_y[k] * dg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pb_x, df_s_55, df_s_56, df_s_57, dh_s_110, \
                         dh_s_111, dh_s_112, df_55, df_56, df_57, dg_80, dg_81, \
                         dg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = -f_6 * df_s_55[k]
                   + f_2 * dh_s_110[k]
                   + f_0 * df_55[k]
                   + pb_x[k] * dg_80[k];

        t_111[k] = -f_4 * df_s_56[k]
                   + f_2 * dh_s_111[k]
                   + f_5 * df_56[k]
                   + pb_x[k] * dg_81[k];

        t_112[k] = -f_4 * df_s_57[k]
                   + f_2 * dh_s_112[k]
                   + f_5 * df_57[k]
                   + pb_x[k] * dg_82[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pb_y, df_s_59, dh_s_113, dh_s_114, \
                         dh_s_115, dh_s_116, df_59, dg_80, dg_84, dg_85, \
                         dg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_2 * dh_s_113[k]
                   + pb_y[k] * dg_80[k];

        t_114[k] = -f_4 * df_s_59[k]
                   + f_2 * dh_s_114[k]
                   + f_5 * df_59[k]
                   + pb_x[k] * dg_84[k];

        t_115[k] = f_2 * dh_s_115[k]
                   + pb_x[k] * dg_85[k];

        t_116[k] = f_2 * dh_s_116[k]
                   + pb_x[k] * dg_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, df_s_56, dh_s_117, dh_s_118, \
                         dh_s_119, dh_s_120, df_56, dg_85, dg_87, dg_88, \
                         dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_2 * dh_s_117[k]
                   + pb_x[k] * dg_87[k];

        t_118[k] = f_2 * dh_s_118[k]
                   + pb_x[k] * dg_88[k];

        t_119[k] = f_2 * dh_s_119[k]
                   + pb_x[k] * dg_89[k];

        t_120[k] = -f_1 * df_s_56[k]
                   + f_2 * dh_s_120[k]
                   + f_3 * df_56[k]
                   + pb_y[k] * dg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_y, df_s_57, df_s_58, df_s_59, dh_s_121, \
                         dh_s_122, dh_s_123, df_57, df_58, df_59, dg_86, dg_87, \
                         dg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_8 * df_s_57[k]
                   + f_2 * dh_s_121[k]
                   + f_7 * df_57[k]
                   + pb_y[k] * dg_86[k];

        t_122[k] = -f_6 * df_s_58[k]
                   + f_2 * dh_s_122[k]
                   + f_0 * df_58[k]
                   + pb_y[k] * dg_87[k];

        t_123[k] = -f_4 * df_s_59[k]
                   + f_2 * dh_s_123[k]
                   + f_5 * df_59[k]
                   + pb_y[k] * dg_88[k];
    }

#pragma omp simd aligned(t_124, t_125, pb_y, pb_z, pg_44, df_s_59, dh_s_124, dh_s_125, df_59, \
                         dg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * dh_s_124[k]
                   + pb_y[k] * dg_89[k];

        t_125[k] = f_0 * pg_44[k]
                   - f_1 * df_s_59[k]
                   + f_2 * dh_s_125[k]
                   + f_3 * df_59[k]
                   + pb_z[k] * dg_89[k];
    }
}

}  // namespace simdkin
