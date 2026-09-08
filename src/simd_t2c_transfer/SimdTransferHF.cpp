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


#include "SimdTransferHF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hf_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t hd, const size_t id,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.40625 * std::sqrt(35.0);
    const auto f_1 = 0.46875 * std::sqrt(35.0);
    const auto f_2 = 2.8125 * std::sqrt(35.0);
    const auto f_3 = 0.9375 * std::sqrt(35.0);
    const auto f_4 = 0.28125 * std::sqrt(35.0);
    const auto f_5 = 0.09375 * std::sqrt(35.0);
    const auto f_6 = 0.9375 * std::sqrt(210.0);
    const auto f_7 = 1.875 * std::sqrt(210.0);
    const auto f_8 = 0.1875 * std::sqrt(210.0);
    const auto f_9 = 0.46875 * std::sqrt(21.0);
    const auto f_10 = 1.875 * std::sqrt(21.0);
    const auto f_11 = 0.9375 * std::sqrt(21.0);
    const auto f_12 = 3.75 * std::sqrt(21.0);
    const auto f_13 = 0.09375 * std::sqrt(21.0);
    const auto f_14 = 0.375 * std::sqrt(21.0);
    const auto f_15 = 1.40625 * std::sqrt(14.0);
    const auto f_16 = 0.9375 * std::sqrt(14.0);
    const auto f_17 = 2.8125 * std::sqrt(14.0);
    const auto f_18 = 1.875 * std::sqrt(14.0);
    const auto f_19 = 0.28125 * std::sqrt(14.0);
    const auto f_20 = 0.1875 * std::sqrt(14.0);
    const auto f_21 = 0.46875 * std::sqrt(210.0);
    const auto f_22 = 0.09375 * std::sqrt(210.0);
    const auto f_23 = 5.625 * std::sqrt(14.0);
    const auto f_24 = 7.5 * std::sqrt(21.0);
    const auto f_25 = 0.375 * std::sqrt(210.0);
    const auto f_26 = 1.5 * std::sqrt(210.0);
    const auto f_27 = 2.25 * std::sqrt(35.0);
    const auto f_28 = 1.5 * std::sqrt(35.0);
    const auto f_29 = 1.40625 * std::sqrt(7.0);
    const auto f_30 = 0.46875 * std::sqrt(7.0);
    const auto f_31 = 0.9375 * std::sqrt(7.0);
    const auto f_32 = 0.3125 * std::sqrt(7.0);
    const auto f_33 = 11.25 * std::sqrt(7.0);
    const auto f_34 = 3.75 * std::sqrt(7.0);
    const auto f_35 = 0.15625 * std::sqrt(7.0);
    const auto f_36 = 1.25 * std::sqrt(7.0);
    const auto f_37 = 0.9375 * std::sqrt(42.0);
    const auto f_38 = 0.625 * std::sqrt(42.0);
    const auto f_39 = 7.5 * std::sqrt(42.0);
    const auto f_40 = 0.3125 * std::sqrt(42.0);
    const auto f_41 = 2.5 * std::sqrt(42.0);
    const auto f_42 = 0.09375 * std::sqrt(105.0);
    const auto f_43 = 0.375 * std::sqrt(105.0);
    const auto f_44 = 0.0625 * std::sqrt(105.0);
    const auto f_45 = 0.25 * std::sqrt(105.0);
    const auto f_46 = 0.75 * std::sqrt(105.0);
    const auto f_47 = 3.0 * std::sqrt(105.0);
    const auto f_48 = 0.03125 * std::sqrt(105.0);
    const auto f_49 = 0.125 * std::sqrt(105.0);
    const auto f_50 = std::sqrt(105.0);
    const auto f_51 = 0.28125 * std::sqrt(70.0);
    const auto f_52 = 0.1875 * std::sqrt(70.0);
    const auto f_53 = 0.125 * std::sqrt(70.0);
    const auto f_54 = 2.25 * std::sqrt(70.0);
    const auto f_55 = 1.5 * std::sqrt(70.0);
    const auto f_56 = 0.09375 * std::sqrt(70.0);
    const auto f_57 = 0.0625 * std::sqrt(70.0);
    const auto f_58 = 0.75 * std::sqrt(70.0);
    const auto f_59 = 0.5 * std::sqrt(70.0);
    const auto f_60 = 0.46875 * std::sqrt(42.0);
    const auto f_61 = 3.75 * std::sqrt(42.0);
    const auto f_62 = 0.15625 * std::sqrt(42.0);
    const auto f_63 = 1.25 * std::sqrt(42.0);
    const auto f_64 = 1.875 * std::sqrt(42.0);
    const auto f_65 = 7.5 * std::sqrt(7.0);
    const auto f_66 = 15.0 * std::sqrt(7.0);
    const auto f_67 = 0.375 * std::sqrt(70.0);
    const auto f_68 = 3.0 * std::sqrt(70.0);
    const auto f_69 = 0.5 * std::sqrt(105.0);
    const auto f_70 = 1.5 * std::sqrt(105.0);
    const auto f_71 = 0.46875 * std::sqrt(6.0);
    const auto f_72 = 0.15625 * std::sqrt(6.0);
    const auto f_73 = 0.9375 * std::sqrt(6.0);
    const auto f_74 = 0.3125 * std::sqrt(6.0);
    const auto f_75 = 5.625 * std::sqrt(6.0);
    const auto f_76 = 1.875 * std::sqrt(6.0);
    const auto f_77 = 3.75 * std::sqrt(6.0);
    const auto f_78 = 1.25 * std::sqrt(6.0);
    const auto f_79 = 0.09375 * std::sqrt(10.0);
    const auto f_80 = 0.375 * std::sqrt(10.0);
    const auto f_81 = 0.1875 * std::sqrt(10.0);
    const auto f_82 = 0.75 * std::sqrt(10.0);
    const auto f_83 = 1.125 * std::sqrt(10.0);
    const auto f_84 = 4.5 * std::sqrt(10.0);
    const auto f_85 = 3.0 * std::sqrt(10.0);
    const auto f_86 = 0.1875 * std::sqrt(15.0);
    const auto f_87 = 0.125 * std::sqrt(15.0);
    const auto f_88 = 0.375 * std::sqrt(15.0);
    const auto f_89 = 0.25 * std::sqrt(15.0);
    const auto f_90 = 2.25 * std::sqrt(15.0);
    const auto f_91 = 1.5 * std::sqrt(15.0);
    const auto f_92 = std::sqrt(15.0);
    const auto f_93 = 1.40625 * std::sqrt(10.0);
    const auto f_94 = 0.46875 * std::sqrt(10.0);
    const auto f_95 = 2.8125 * std::sqrt(10.0);
    const auto f_96 = 0.9375 * std::sqrt(10.0);
    const auto f_97 = 3.75 * std::sqrt(10.0);
    const auto f_98 = 1.25 * std::sqrt(10.0);
    const auto f_99 = 0.25 * std::sqrt(10.0);
    const auto f_100 = 1.875 * std::sqrt(15.0);
    const auto f_101 = 3.75 * std::sqrt(15.0);
    const auto f_102 = 5.0 * std::sqrt(15.0);
    const auto f_103 = 5.0 * std::sqrt(6.0);
    const auto f_104 = 0.25 * std::sqrt(6.0);
    const auto f_105 = std::sqrt(6.0);
    const auto f_106 = 0.9375 * std::sqrt(15.0);
    const auto f_107 = 2.5 * std::sqrt(15.0);
    const auto f_108 = 0.5 * std::sqrt(15.0);
    const auto f_109 = 1.875 * std::sqrt(7.0);
    const auto f_110 = 0.46875 * std::sqrt(14.0);
    const auto f_111 = 8.4375 * std::sqrt(14.0);
    const auto f_112 = 11.25 * std::sqrt(21.0);
    const auto f_113 = 0.5625 * std::sqrt(210.0);
    const auto f_114 = 2.25 * std::sqrt(210.0);
    const auto f_115 = 0.5625 * std::sqrt(35.0);
    const auto f_116 = 0.375 * std::sqrt(35.0);
    const auto f_117 = 3.375 * std::sqrt(35.0);
    const auto f_118 = 5.625 * std::sqrt(21.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;
    auto *g_45 = values + 45 * nvalues;
    auto *g_46 = values + 46 * nvalues;
    auto *g_47 = values + 47 * nvalues;
    auto *g_48 = values + 48 * nvalues;
    auto *g_49 = values + 49 * nvalues;
    auto *g_50 = values + 50 * nvalues;
    auto *g_51 = values + 51 * nvalues;
    auto *g_52 = values + 52 * nvalues;
    auto *g_53 = values + 53 * nvalues;
    auto *g_54 = values + 54 * nvalues;
    auto *g_55 = values + 55 * nvalues;
    auto *g_56 = values + 56 * nvalues;
    auto *g_57 = values + 57 * nvalues;
    auto *g_58 = values + 58 * nvalues;
    auto *g_59 = values + 59 * nvalues;
    auto *g_60 = values + 60 * nvalues;
    auto *g_61 = values + 61 * nvalues;
    auto *g_62 = values + 62 * nvalues;
    auto *g_63 = values + 63 * nvalues;
    auto *g_64 = values + 64 * nvalues;
    auto *g_65 = values + 65 * nvalues;
    auto *g_66 = values + 66 * nvalues;
    auto *g_67 = values + 67 * nvalues;
    auto *g_68 = values + 68 * nvalues;
    auto *g_69 = values + 69 * nvalues;
    auto *g_70 = values + 70 * nvalues;
    auto *g_71 = values + 71 * nvalues;
    auto *g_72 = values + 72 * nvalues;
    auto *g_73 = values + 73 * nvalues;
    auto *g_74 = values + 74 * nvalues;
    auto *g_75 = values + 75 * nvalues;
    auto *g_76 = values + 76 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

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
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_167 = buffer.data(id + 167);

#pragma omp simd aligned(ab_x, ab_y, hd_7, hd_9, hd_37, hd_39, hd_91, hd_93, id_7, id_21, \
                         id_37, id_63, id_91, id_129 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * hd_7[k]
                 - f_1 * ab_y[k] * hd_9[k]
                 - f_2 * ab_x[k] * hd_37[k]
                 + f_3 * ab_y[k] * hd_39[k]
                 + f_4 * ab_x[k] * hd_91[k]
                 - f_5 * ab_y[k] * hd_93[k]
                 + f_0 * id_7[k]
                 - f_1 * id_21[k]
                 - f_2 * id_37[k]
                 + f_3 * id_63[k]
                 + f_4 * id_91[k]
                 - f_5 * id_129[k];
    }

#pragma omp simd aligned(ab_x, hd_10, hd_40, hd_94, id_10, id_40, \
                         id_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_6 * ab_x[k] * hd_10[k]
                 - f_7 * ab_x[k] * hd_40[k]
                 + f_8 * ab_x[k] * hd_94[k]
                 + f_6 * id_10[k]
                 - f_7 * id_40[k]
                 + f_8 * id_94[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_7, hd_9, hd_11, hd_37, hd_39, hd_41, hd_91, hd_93, \
                         hd_95, id_7, id_21, id_23, id_37, id_63, id_65, id_91, id_129, \
                         id_131 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_9 * ab_x[k] * hd_7[k]
                 - f_9 * ab_y[k] * hd_9[k]
                 + f_10 * ab_y[k] * hd_11[k]
                 + f_11 * ab_x[k] * hd_37[k]
                 + f_11 * ab_y[k] * hd_39[k]
                 - f_12 * ab_y[k] * hd_41[k]
                 - f_13 * ab_x[k] * hd_91[k]
                 - f_13 * ab_y[k] * hd_93[k]
                 + f_14 * ab_y[k] * hd_95[k]
                 - f_9 * id_7[k]
                 - f_9 * id_21[k]
                 + f_10 * id_23[k]
                 + f_11 * id_37[k]
                 + f_11 * id_63[k]
                 - f_12 * id_65[k]
                 - f_13 * id_91[k]
                 - f_13 * id_129[k]
                 + f_14 * id_131[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_8, hd_10, hd_11, hd_38, hd_40, hd_41, hd_92, \
                         hd_94, hd_95, id_8, id_22, id_29, id_38, id_64, id_71, id_92, id_130, \
                         id_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_15 * ab_x[k] * hd_8[k]
                 - f_15 * ab_y[k] * hd_10[k]
                 + f_16 * ab_z[k] * hd_11[k]
                 + f_17 * ab_x[k] * hd_38[k]
                 + f_17 * ab_y[k] * hd_40[k]
                 - f_18 * ab_z[k] * hd_41[k]
                 - f_19 * ab_x[k] * hd_92[k]
                 - f_19 * ab_y[k] * hd_94[k]
                 + f_20 * ab_z[k] * hd_95[k]
                 - f_15 * id_8[k]
                 - f_15 * id_22[k]
                 + f_16 * id_29[k]
                 + f_17 * id_38[k]
                 + f_17 * id_64[k]
                 - f_18 * id_71[k]
                 - f_19 * id_92[k]
                 - f_19 * id_130[k]
                 + f_20 * id_137[k];
    }

#pragma omp simd aligned(ab_x, hd_6, hd_9, hd_11, hd_36, hd_39, hd_41, hd_90, hd_93, hd_95, \
                         id_6, id_9, id_11, id_36, id_39, id_41, id_90, id_93, \
                         id_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_9 * ab_x[k] * hd_6[k]
                 - f_9 * ab_x[k] * hd_9[k]
                 + f_10 * ab_x[k] * hd_11[k]
                 + f_11 * ab_x[k] * hd_36[k]
                 + f_11 * ab_x[k] * hd_39[k]
                 - f_12 * ab_x[k] * hd_41[k]
                 - f_13 * ab_x[k] * hd_90[k]
                 - f_13 * ab_x[k] * hd_93[k]
                 + f_14 * ab_x[k] * hd_95[k]
                 - f_9 * id_6[k]
                 - f_9 * id_9[k]
                 + f_10 * id_11[k]
                 + f_11 * id_36[k]
                 + f_11 * id_39[k]
                 - f_12 * id_41[k]
                 - f_13 * id_90[k]
                 - f_13 * id_93[k]
                 + f_14 * id_95[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_8, hd_10, hd_38, hd_40, hd_92, hd_94, id_8, id_22, \
                         id_38, id_64, id_92, id_130 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_21 * ab_x[k] * hd_8[k]
                 - f_21 * ab_y[k] * hd_10[k]
                 - f_6 * ab_x[k] * hd_38[k]
                 + f_6 * ab_y[k] * hd_40[k]
                 + f_22 * ab_x[k] * hd_92[k]
                 - f_22 * ab_y[k] * hd_94[k]
                 + f_21 * id_8[k]
                 - f_21 * id_22[k]
                 - f_6 * id_38[k]
                 + f_6 * id_64[k]
                 + f_22 * id_92[k]
                 - f_22 * id_130[k];
    }

#pragma omp simd aligned(ab_x, hd_6, hd_9, hd_36, hd_39, hd_90, hd_93, id_6, id_9, id_36, \
                         id_39, id_90, id_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * ab_x[k] * hd_6[k]
                 - f_0 * ab_x[k] * hd_9[k]
                 - f_3 * ab_x[k] * hd_36[k]
                 + f_2 * ab_x[k] * hd_39[k]
                 + f_5 * ab_x[k] * hd_90[k]
                 - f_4 * ab_x[k] * hd_93[k]
                 + f_1 * id_6[k]
                 - f_0 * id_9[k]
                 - f_3 * id_36[k]
                 + f_2 * id_39[k]
                 + f_5 * id_90[k]
                 - f_4 * id_93[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_25, hd_27, hd_28, hd_67, hd_69, hd_70, id_25, id_28, \
                         id_45, id_67, id_70, id_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_23 * ab_x[k] * hd_25[k]
                 - f_18 * ab_y[k] * hd_27[k]
                 - f_23 * ab_x[k] * hd_67[k]
                 + f_18 * ab_y[k] * hd_69[k]
                 + f_23 * id_25[k]
                 - f_18 * id_45[k]
                 - f_23 * id_67[k]
                 + f_18 * id_99[k];

        g_8[k] = f_24 * ab_x[k] * hd_28[k]
                 - f_24 * ab_x[k] * hd_70[k]
                 + f_24 * id_28[k]
                 - f_24 * id_70[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_25, hd_27, hd_29, hd_67, hd_69, hd_71, id_25, id_45, \
                         id_47, id_67, id_99, id_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_25 * ab_x[k] * hd_25[k]
                 - f_25 * ab_y[k] * hd_27[k]
                 + f_26 * ab_y[k] * hd_29[k]
                 + f_25 * ab_x[k] * hd_67[k]
                 + f_25 * ab_y[k] * hd_69[k]
                 - f_26 * ab_y[k] * hd_71[k]
                 - f_25 * id_25[k]
                 - f_25 * id_45[k]
                 + f_26 * id_47[k]
                 + f_25 * id_67[k]
                 + f_25 * id_99[k]
                 - f_26 * id_101[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_26, hd_28, hd_29, hd_68, hd_70, hd_71, id_26, \
                         id_46, id_53, id_68, id_100, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_27 * ab_x[k] * hd_26[k]
                  - f_27 * ab_y[k] * hd_28[k]
                  + f_28 * ab_z[k] * hd_29[k]
                  + f_27 * ab_x[k] * hd_68[k]
                  + f_27 * ab_y[k] * hd_70[k]
                  - f_28 * ab_z[k] * hd_71[k]
                  - f_27 * id_26[k]
                  - f_27 * id_46[k]
                  + f_28 * id_53[k]
                  + f_27 * id_68[k]
                  + f_27 * id_100[k]
                  - f_28 * id_107[k];
    }

#pragma omp simd aligned(ab_x, hd_24, hd_27, hd_29, hd_66, hd_69, hd_71, id_24, id_27, id_29, \
                         id_66, id_69, id_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_25 * ab_x[k] * hd_24[k]
                  - f_25 * ab_x[k] * hd_27[k]
                  + f_26 * ab_x[k] * hd_29[k]
                  + f_25 * ab_x[k] * hd_66[k]
                  + f_25 * ab_x[k] * hd_69[k]
                  - f_26 * ab_x[k] * hd_71[k]
                  - f_25 * id_24[k]
                  - f_25 * id_27[k]
                  + f_26 * id_29[k]
                  + f_25 * id_66[k]
                  + f_25 * id_69[k]
                  - f_26 * id_71[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_26, hd_28, hd_68, hd_70, id_26, id_46, id_68, \
                         id_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_12 * ab_x[k] * hd_26[k]
                  - f_12 * ab_y[k] * hd_28[k]
                  - f_12 * ab_x[k] * hd_68[k]
                  + f_12 * ab_y[k] * hd_70[k]
                  + f_12 * id_26[k]
                  - f_12 * id_46[k]
                  - f_12 * id_68[k]
                  + f_12 * id_100[k];
    }

#pragma omp simd aligned(ab_x, hd_24, hd_27, hd_66, hd_69, id_24, id_27, id_66, \
                         id_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_18 * ab_x[k] * hd_24[k]
                  - f_23 * ab_x[k] * hd_27[k]
                  - f_18 * ab_x[k] * hd_66[k]
                  + f_23 * ab_x[k] * hd_69[k]
                  + f_18 * id_24[k]
                  - f_23 * id_27[k]
                  - f_18 * id_66[k]
                  + f_23 * id_69[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_7, hd_9, hd_37, hd_39, hd_49, hd_51, hd_91, hd_93, \
                         hd_103, hd_105, id_7, id_21, id_37, id_49, id_63, id_75, id_91, \
                         id_103, id_129, id_141 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_29 * ab_x[k] * hd_7[k]
                  + f_30 * ab_y[k] * hd_9[k]
                  - f_31 * ab_x[k] * hd_37[k]
                  + f_32 * ab_y[k] * hd_39[k]
                  + f_33 * ab_x[k] * hd_49[k]
                  - f_34 * ab_y[k] * hd_51[k]
                  + f_30 * ab_x[k] * hd_91[k]
                  - f_35 * ab_y[k] * hd_93[k]
                  - f_34 * ab_x[k] * hd_103[k]
                  + f_36 * ab_y[k] * hd_105[k]
                  - f_29 * id_7[k]
                  + f_30 * id_21[k]
                  - f_31 * id_37[k]
                  + f_33 * id_49[k]
                  + f_32 * id_63[k]
                  - f_34 * id_75[k]
                  + f_30 * id_91[k]
                  - f_34 * id_103[k]
                  - f_35 * id_129[k]
                  + f_36 * id_141[k];
    }

#pragma omp simd aligned(ab_x, hd_10, hd_40, hd_52, hd_94, hd_106, id_10, id_40, id_52, id_94, \
                         id_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_37 * ab_x[k] * hd_10[k]
                  - f_38 * ab_x[k] * hd_40[k]
                  + f_39 * ab_x[k] * hd_52[k]
                  + f_40 * ab_x[k] * hd_94[k]
                  - f_41 * ab_x[k] * hd_106[k]
                  - f_37 * id_10[k]
                  - f_38 * id_40[k]
                  + f_39 * id_52[k]
                  + f_40 * id_94[k]
                  - f_41 * id_106[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_7, hd_9, hd_11, hd_37, hd_39, hd_41, hd_49, hd_51, \
                         hd_53, hd_91, hd_93, hd_95, hd_103, hd_105, hd_107, id_7, id_21, \
                         id_23, id_37, id_49, id_63, id_65, id_75, id_77, id_91, id_103, \
                         id_129, id_131, id_141, id_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_42 * ab_x[k] * hd_7[k]
                  + f_42 * ab_y[k] * hd_9[k]
                  - f_43 * ab_y[k] * hd_11[k]
                  + f_44 * ab_x[k] * hd_37[k]
                  + f_44 * ab_y[k] * hd_39[k]
                  - f_45 * ab_y[k] * hd_41[k]
                  - f_46 * ab_x[k] * hd_49[k]
                  - f_46 * ab_y[k] * hd_51[k]
                  + f_47 * ab_y[k] * hd_53[k]
                  - f_48 * ab_x[k] * hd_91[k]
                  - f_48 * ab_y[k] * hd_93[k]
                  + f_49 * ab_y[k] * hd_95[k]
                  + f_45 * ab_x[k] * hd_103[k]
                  + f_45 * ab_y[k] * hd_105[k]
                  - f_50 * ab_y[k] * hd_107[k]
                  + f_42 * id_7[k]
                  + f_42 * id_21[k]
                  - f_43 * id_23[k]
                  + f_44 * id_37[k]
                  - f_46 * id_49[k]
                  + f_44 * id_63[k]
                  - f_45 * id_65[k]
                  - f_46 * id_75[k]
                  + f_47 * id_77[k]
                  - f_48 * id_91[k]
                  + f_45 * id_103[k]
                  - f_48 * id_129[k]
                  + f_49 * id_131[k]
                  + f_45 * id_141[k]
                  - f_50 * id_143[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_8, hd_10, hd_11, hd_38, hd_40, hd_41, hd_50, \
                         hd_52, hd_53, hd_92, hd_94, hd_95, hd_104, hd_106, hd_107, id_8, \
                         id_22, id_29, id_38, id_50, id_64, id_71, id_76, id_83, id_92, \
                         id_104, id_130, id_137, id_142, id_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_51 * ab_x[k] * hd_8[k]
                  + f_51 * ab_y[k] * hd_10[k]
                  - f_52 * ab_z[k] * hd_11[k]
                  + f_52 * ab_x[k] * hd_38[k]
                  + f_52 * ab_y[k] * hd_40[k]
                  - f_53 * ab_z[k] * hd_41[k]
                  - f_54 * ab_x[k] * hd_50[k]
                  - f_54 * ab_y[k] * hd_52[k]
                  + f_55 * ab_z[k] * hd_53[k]
                  - f_56 * ab_x[k] * hd_92[k]
                  - f_56 * ab_y[k] * hd_94[k]
                  + f_57 * ab_z[k] * hd_95[k]
                  + f_58 * ab_x[k] * hd_104[k]
                  + f_58 * ab_y[k] * hd_106[k]
                  - f_59 * ab_z[k] * hd_107[k]
                  + f_51 * id_8[k]
                  + f_51 * id_22[k]
                  - f_52 * id_29[k]
                  + f_52 * id_38[k]
                  - f_54 * id_50[k]
                  + f_52 * id_64[k]
                  - f_53 * id_71[k]
                  - f_54 * id_76[k]
                  + f_55 * id_83[k]
                  - f_56 * id_92[k]
                  + f_58 * id_104[k]
                  - f_56 * id_130[k]
                  + f_57 * id_137[k]
                  + f_58 * id_142[k]
                  - f_59 * id_149[k];
    }

#pragma omp simd aligned(ab_x, hd_6, hd_9, hd_11, hd_36, hd_39, hd_41, hd_48, hd_51, hd_53, \
                         hd_90, hd_93, hd_95, hd_102, hd_105, hd_107, id_6, id_9, id_11, \
                         id_36, id_39, id_41, id_48, id_51, id_53, id_90, id_93, id_95, \
                         id_102, id_105, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_42 * ab_x[k] * hd_6[k]
                  + f_42 * ab_x[k] * hd_9[k]
                  - f_43 * ab_x[k] * hd_11[k]
                  + f_44 * ab_x[k] * hd_36[k]
                  + f_44 * ab_x[k] * hd_39[k]
                  - f_45 * ab_x[k] * hd_41[k]
                  - f_46 * ab_x[k] * hd_48[k]
                  - f_46 * ab_x[k] * hd_51[k]
                  + f_47 * ab_x[k] * hd_53[k]
                  - f_48 * ab_x[k] * hd_90[k]
                  - f_48 * ab_x[k] * hd_93[k]
                  + f_49 * ab_x[k] * hd_95[k]
                  + f_45 * ab_x[k] * hd_102[k]
                  + f_45 * ab_x[k] * hd_105[k]
                  - f_50 * ab_x[k] * hd_107[k]
                  + f_42 * id_6[k]
                  + f_42 * id_9[k]
                  - f_43 * id_11[k]
                  + f_44 * id_36[k]
                  + f_44 * id_39[k]
                  - f_45 * id_41[k]
                  - f_46 * id_48[k]
                  - f_46 * id_51[k]
                  + f_47 * id_53[k]
                  - f_48 * id_90[k]
                  - f_48 * id_93[k]
                  + f_49 * id_95[k]
                  + f_45 * id_102[k]
                  + f_45 * id_105[k]
                  - f_50 * id_107[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_8, hd_10, hd_38, hd_40, hd_50, hd_52, hd_92, hd_94, \
                         hd_104, hd_106, id_8, id_22, id_38, id_50, id_64, id_76, id_92, \
                         id_104, id_130, id_142 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_60 * ab_x[k] * hd_8[k]
                  + f_60 * ab_y[k] * hd_10[k]
                  - f_40 * ab_x[k] * hd_38[k]
                  + f_40 * ab_y[k] * hd_40[k]
                  + f_61 * ab_x[k] * hd_50[k]
                  - f_61 * ab_y[k] * hd_52[k]
                  + f_62 * ab_x[k] * hd_92[k]
                  - f_62 * ab_y[k] * hd_94[k]
                  - f_63 * ab_x[k] * hd_104[k]
                  + f_63 * ab_y[k] * hd_106[k]
                  - f_60 * id_8[k]
                  + f_60 * id_22[k]
                  - f_40 * id_38[k]
                  + f_61 * id_50[k]
                  + f_40 * id_64[k]
                  - f_61 * id_76[k]
                  + f_62 * id_92[k]
                  - f_63 * id_104[k]
                  - f_62 * id_130[k]
                  + f_63 * id_142[k];
    }

#pragma omp simd aligned(ab_x, hd_6, hd_9, hd_36, hd_39, hd_48, hd_51, hd_90, hd_93, hd_102, \
                         hd_105, id_6, id_9, id_36, id_39, id_48, id_51, id_90, id_93, id_102, \
                         id_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_30 * ab_x[k] * hd_6[k]
                  + f_29 * ab_x[k] * hd_9[k]
                  - f_32 * ab_x[k] * hd_36[k]
                  + f_31 * ab_x[k] * hd_39[k]
                  + f_34 * ab_x[k] * hd_48[k]
                  - f_33 * ab_x[k] * hd_51[k]
                  + f_35 * ab_x[k] * hd_90[k]
                  - f_30 * ab_x[k] * hd_93[k]
                  - f_36 * ab_x[k] * hd_102[k]
                  + f_34 * ab_x[k] * hd_105[k]
                  - f_30 * id_6[k]
                  + f_29 * id_9[k]
                  - f_32 * id_36[k]
                  + f_31 * id_39[k]
                  + f_34 * id_48[k]
                  - f_33 * id_51[k]
                  + f_35 * id_90[k]
                  - f_30 * id_93[k]
                  - f_36 * id_102[k]
                  + f_34 * id_105[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_25, hd_27, hd_67, hd_69, hd_79, hd_81, id_25, id_45, \
                         id_67, id_79, id_99, id_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_64 * ab_x[k] * hd_25[k]
                  + f_38 * ab_y[k] * hd_27[k]
                  - f_64 * ab_x[k] * hd_67[k]
                  + f_38 * ab_y[k] * hd_69[k]
                  + f_61 * ab_x[k] * hd_79[k]
                  - f_63 * ab_y[k] * hd_81[k]
                  - f_64 * id_25[k]
                  + f_38 * id_45[k]
                  - f_64 * id_67[k]
                  + f_61 * id_79[k]
                  + f_38 * id_99[k]
                  - f_63 * id_111[k];
    }

#pragma omp simd aligned(ab_x, hd_28, hd_70, hd_82, id_28, id_70, \
                         id_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_65 * ab_x[k] * hd_28[k]
                  - f_65 * ab_x[k] * hd_70[k]
                  + f_66 * ab_x[k] * hd_82[k]
                  - f_65 * id_28[k]
                  - f_65 * id_70[k]
                  + f_66 * id_82[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_25, hd_27, hd_29, hd_67, hd_69, hd_71, hd_79, hd_81, \
                         hd_83, id_25, id_45, id_47, id_67, id_79, id_99, id_101, id_111, \
                         id_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_67 * ab_x[k] * hd_25[k]
                  + f_67 * ab_y[k] * hd_27[k]
                  - f_55 * ab_y[k] * hd_29[k]
                  + f_67 * ab_x[k] * hd_67[k]
                  + f_67 * ab_y[k] * hd_69[k]
                  - f_55 * ab_y[k] * hd_71[k]
                  - f_58 * ab_x[k] * hd_79[k]
                  - f_58 * ab_y[k] * hd_81[k]
                  + f_68 * ab_y[k] * hd_83[k]
                  + f_67 * id_25[k]
                  + f_67 * id_45[k]
                  - f_55 * id_47[k]
                  + f_67 * id_67[k]
                  - f_58 * id_79[k]
                  + f_67 * id_99[k]
                  - f_55 * id_101[k]
                  - f_58 * id_111[k]
                  + f_68 * id_113[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_26, hd_28, hd_29, hd_68, hd_70, hd_71, hd_80, \
                         hd_82, hd_83, id_26, id_46, id_53, id_68, id_80, id_100, id_107, \
                         id_112, id_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_46 * ab_x[k] * hd_26[k]
                  + f_46 * ab_y[k] * hd_28[k]
                  - f_69 * ab_z[k] * hd_29[k]
                  + f_46 * ab_x[k] * hd_68[k]
                  + f_46 * ab_y[k] * hd_70[k]
                  - f_69 * ab_z[k] * hd_71[k]
                  - f_70 * ab_x[k] * hd_80[k]
                  - f_70 * ab_y[k] * hd_82[k]
                  + f_50 * ab_z[k] * hd_83[k]
                  + f_46 * id_26[k]
                  + f_46 * id_46[k]
                  - f_69 * id_53[k]
                  + f_46 * id_68[k]
                  - f_70 * id_80[k]
                  + f_46 * id_100[k]
                  - f_69 * id_107[k]
                  - f_70 * id_112[k]
                  + f_50 * id_119[k];
    }

#pragma omp simd aligned(ab_x, hd_24, hd_27, hd_29, hd_66, hd_69, hd_71, hd_78, hd_81, hd_83, \
                         id_24, id_27, id_29, id_66, id_69, id_71, id_78, id_81, \
                         id_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_67 * ab_x[k] * hd_24[k]
                  + f_67 * ab_x[k] * hd_27[k]
                  - f_55 * ab_x[k] * hd_29[k]
                  + f_67 * ab_x[k] * hd_66[k]
                  + f_67 * ab_x[k] * hd_69[k]
                  - f_55 * ab_x[k] * hd_71[k]
                  - f_58 * ab_x[k] * hd_78[k]
                  - f_58 * ab_x[k] * hd_81[k]
                  + f_68 * ab_x[k] * hd_83[k]
                  + f_67 * id_24[k]
                  + f_67 * id_27[k]
                  - f_55 * id_29[k]
                  + f_67 * id_66[k]
                  + f_67 * id_69[k]
                  - f_55 * id_71[k]
                  - f_58 * id_78[k]
                  - f_58 * id_81[k]
                  + f_68 * id_83[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_26, hd_28, hd_68, hd_70, hd_80, hd_82, id_26, id_46, \
                         id_68, id_80, id_100, id_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_34 * ab_x[k] * hd_26[k]
                  + f_34 * ab_y[k] * hd_28[k]
                  - f_34 * ab_x[k] * hd_68[k]
                  + f_34 * ab_y[k] * hd_70[k]
                  + f_65 * ab_x[k] * hd_80[k]
                  - f_65 * ab_y[k] * hd_82[k]
                  - f_34 * id_26[k]
                  + f_34 * id_46[k]
                  - f_34 * id_68[k]
                  + f_65 * id_80[k]
                  + f_34 * id_100[k]
                  - f_65 * id_112[k];
    }

#pragma omp simd aligned(ab_x, hd_24, hd_27, hd_66, hd_69, hd_78, hd_81, id_24, id_27, id_66, \
                         id_69, id_78, id_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_38 * ab_x[k] * hd_24[k]
                  + f_64 * ab_x[k] * hd_27[k]
                  - f_38 * ab_x[k] * hd_66[k]
                  + f_64 * ab_x[k] * hd_69[k]
                  + f_63 * ab_x[k] * hd_78[k]
                  - f_61 * ab_x[k] * hd_81[k]
                  - f_38 * id_24[k]
                  + f_64 * id_27[k]
                  - f_38 * id_66[k]
                  + f_64 * id_69[k]
                  + f_63 * id_78[k]
                  - f_61 * id_81[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_7, hd_9, hd_37, hd_39, hd_49, hd_51, hd_91, hd_93, \
                         hd_103, hd_105, hd_115, hd_117, id_7, id_21, id_37, id_49, id_63, \
                         id_75, id_91, id_103, id_115, id_129, id_141, \
                         id_153 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_71 * ab_x[k] * hd_7[k]
                  - f_72 * ab_y[k] * hd_9[k]
                  + f_73 * ab_x[k] * hd_37[k]
                  - f_74 * ab_y[k] * hd_39[k]
                  - f_75 * ab_x[k] * hd_49[k]
                  + f_76 * ab_y[k] * hd_51[k]
                  + f_71 * ab_x[k] * hd_91[k]
                  - f_72 * ab_y[k] * hd_93[k]
                  - f_75 * ab_x[k] * hd_103[k]
                  + f_76 * ab_y[k] * hd_105[k]
                  + f_77 * ab_x[k] * hd_115[k]
                  - f_78 * ab_y[k] * hd_117[k]
                  + f_71 * id_7[k]
                  - f_72 * id_21[k]
                  + f_73 * id_37[k]
                  - f_75 * id_49[k]
                  - f_74 * id_63[k]
                  + f_76 * id_75[k]
                  + f_71 * id_91[k]
                  - f_75 * id_103[k]
                  + f_77 * id_115[k]
                  - f_72 * id_129[k]
                  + f_76 * id_141[k]
                  - f_78 * id_153[k];
    }

#pragma omp simd aligned(ab_x, hd_10, hd_40, hd_52, hd_94, hd_106, hd_118, id_10, id_40, \
                         id_52, id_94, id_106, id_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = 1.875 * ab_x[k] * hd_10[k]
                  + 3.75 * ab_x[k] * hd_40[k]
                  - 22.5 * ab_x[k] * hd_52[k]
                  + 1.875 * ab_x[k] * hd_94[k]
                  - 22.5 * ab_x[k] * hd_106[k]
                  + 15.0 * ab_x[k] * hd_118[k]
                  + 1.875 * id_10[k]
                  + 3.75 * id_40[k]
                  - 22.5 * id_52[k]
                  + 1.875 * id_94[k]
                  - 22.5 * id_106[k]
                  + 15.0 * id_118[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_7, hd_9, hd_11, hd_37, hd_39, hd_41, hd_49, hd_51, \
                         hd_53, hd_91, hd_93, hd_95, hd_103, hd_105, hd_107, hd_115, hd_117, \
                         hd_119, id_7, id_21, id_23, id_37, id_49, id_63, id_65, id_75, id_77, \
                         id_91, id_103, id_115, id_129, id_131, id_141, id_143, id_153, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_79 * ab_x[k] * hd_7[k]
                  - f_79 * ab_y[k] * hd_9[k]
                  + f_80 * ab_y[k] * hd_11[k]
                  - f_81 * ab_x[k] * hd_37[k]
                  - f_81 * ab_y[k] * hd_39[k]
                  + f_82 * ab_y[k] * hd_41[k]
                  + f_83 * ab_x[k] * hd_49[k]
                  + f_83 * ab_y[k] * hd_51[k]
                  - f_84 * ab_y[k] * hd_53[k]
                  - f_79 * ab_x[k] * hd_91[k]
                  - f_79 * ab_y[k] * hd_93[k]
                  + f_80 * ab_y[k] * hd_95[k]
                  + f_83 * ab_x[k] * hd_103[k]
                  + f_83 * ab_y[k] * hd_105[k]
                  - f_84 * ab_y[k] * hd_107[k]
                  - f_82 * ab_x[k] * hd_115[k]
                  - f_82 * ab_y[k] * hd_117[k]
                  + f_85 * ab_y[k] * hd_119[k]
                  - f_79 * id_7[k]
                  - f_79 * id_21[k]
                  + f_80 * id_23[k]
                  - f_81 * id_37[k]
                  + f_83 * id_49[k]
                  - f_81 * id_63[k]
                  + f_82 * id_65[k]
                  + f_83 * id_75[k]
                  - f_84 * id_77[k]
                  - f_79 * id_91[k]
                  + f_83 * id_103[k]
                  - f_82 * id_115[k]
                  - f_79 * id_129[k]
                  + f_80 * id_131[k]
                  + f_83 * id_141[k]
                  - f_84 * id_143[k]
                  - f_82 * id_153[k]
                  + f_85 * id_155[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_8, hd_10, hd_11, hd_38, hd_40, hd_41, hd_50, \
                         hd_52, hd_53, hd_92, hd_94, hd_95, hd_104, hd_106, hd_107, hd_116, \
                         hd_118, hd_119, id_8, id_22, id_29, id_38, id_50, id_64, id_71, \
                         id_76, id_83, id_92, id_104, id_116, id_130, id_137, id_142, id_149, \
                         id_154, id_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_86 * ab_x[k] * hd_8[k]
                  - f_86 * ab_y[k] * hd_10[k]
                  + f_87 * ab_z[k] * hd_11[k]
                  - f_88 * ab_x[k] * hd_38[k]
                  - f_88 * ab_y[k] * hd_40[k]
                  + f_89 * ab_z[k] * hd_41[k]
                  + f_90 * ab_x[k] * hd_50[k]
                  + f_90 * ab_y[k] * hd_52[k]
                  - f_91 * ab_z[k] * hd_53[k]
                  - f_86 * ab_x[k] * hd_92[k]
                  - f_86 * ab_y[k] * hd_94[k]
                  + f_87 * ab_z[k] * hd_95[k]
                  + f_90 * ab_x[k] * hd_104[k]
                  + f_90 * ab_y[k] * hd_106[k]
                  - f_91 * ab_z[k] * hd_107[k]
                  - f_91 * ab_x[k] * hd_116[k]
                  - f_91 * ab_y[k] * hd_118[k]
                  + f_92 * ab_z[k] * hd_119[k]
                  - f_86 * id_8[k]
                  - f_86 * id_22[k]
                  + f_87 * id_29[k]
                  - f_88 * id_38[k]
                  + f_90 * id_50[k]
                  - f_88 * id_64[k]
                  + f_89 * id_71[k]
                  + f_90 * id_76[k]
                  - f_91 * id_83[k]
                  - f_86 * id_92[k]
                  + f_90 * id_104[k]
                  - f_91 * id_116[k]
                  - f_86 * id_130[k]
                  + f_87 * id_137[k]
                  + f_90 * id_142[k]
                  - f_91 * id_149[k]
                  - f_91 * id_154[k]
                  + f_92 * id_161[k];
    }

#pragma omp simd aligned(ab_x, hd_6, hd_9, hd_11, hd_36, hd_39, hd_41, hd_48, hd_51, hd_53, \
                         hd_90, hd_93, hd_95, hd_102, hd_105, hd_107, hd_114, hd_117, hd_119, \
                         id_6, id_9, id_11, id_36, id_39, id_41, id_48, id_51, id_53, id_90, \
                         id_93, id_95, id_102, id_105, id_107, id_114, id_117, \
                         id_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_79 * ab_x[k] * hd_6[k]
                  - f_79 * ab_x[k] * hd_9[k]
                  + f_80 * ab_x[k] * hd_11[k]
                  - f_81 * ab_x[k] * hd_36[k]
                  - f_81 * ab_x[k] * hd_39[k]
                  + f_82 * ab_x[k] * hd_41[k]
                  + f_83 * ab_x[k] * hd_48[k]
                  + f_83 * ab_x[k] * hd_51[k]
                  - f_84 * ab_x[k] * hd_53[k]
                  - f_79 * ab_x[k] * hd_90[k]
                  - f_79 * ab_x[k] * hd_93[k]
                  + f_80 * ab_x[k] * hd_95[k]
                  + f_83 * ab_x[k] * hd_102[k]
                  + f_83 * ab_x[k] * hd_105[k]
                  - f_84 * ab_x[k] * hd_107[k]
                  - f_82 * ab_x[k] * hd_114[k]
                  - f_82 * ab_x[k] * hd_117[k]
                  + f_85 * ab_x[k] * hd_119[k]
                  - f_79 * id_6[k]
                  - f_79 * id_9[k]
                  + f_80 * id_11[k]
                  - f_81 * id_36[k]
                  - f_81 * id_39[k]
                  + f_82 * id_41[k]
                  + f_83 * id_48[k]
                  + f_83 * id_51[k]
                  - f_84 * id_53[k]
                  - f_79 * id_90[k]
                  - f_79 * id_93[k]
                  + f_80 * id_95[k]
                  + f_83 * id_102[k]
                  + f_83 * id_105[k]
                  - f_84 * id_107[k]
                  - f_82 * id_114[k]
                  - f_82 * id_117[k]
                  + f_85 * id_119[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_8, hd_10, hd_38, hd_40, hd_50, hd_52, hd_92, hd_94, \
                         hd_104, hd_106, hd_116, hd_118, id_8, id_22, id_38, id_50, id_64, \
                         id_76, id_92, id_104, id_116, id_130, id_142, \
                         id_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = 0.9375 * ab_x[k] * hd_8[k]
                  - 0.9375 * ab_y[k] * hd_10[k]
                  + 1.875 * ab_x[k] * hd_38[k]
                  - 1.875 * ab_y[k] * hd_40[k]
                  - 11.25 * ab_x[k] * hd_50[k]
                  + 11.25 * ab_y[k] * hd_52[k]
                  + 0.9375 * ab_x[k] * hd_92[k]
                  - 0.9375 * ab_y[k] * hd_94[k]
                  - 11.25 * ab_x[k] * hd_104[k]
                  + 11.25 * ab_y[k] * hd_106[k]
                  + 7.5 * ab_x[k] * hd_116[k]
                  - 7.5 * ab_y[k] * hd_118[k]
                  + 0.9375 * id_8[k]
                  - 0.9375 * id_22[k]
                  + 1.875 * id_38[k]
                  - 11.25 * id_50[k]
                  - 1.875 * id_64[k]
                  + 11.25 * id_76[k]
                  + 0.9375 * id_92[k]
                  - 11.25 * id_104[k]
                  + 7.5 * id_116[k]
                  - 0.9375 * id_130[k]
                  + 11.25 * id_142[k]
                  - 7.5 * id_154[k];
    }

#pragma omp simd aligned(ab_x, hd_6, hd_9, hd_36, hd_39, hd_48, hd_51, hd_90, hd_93, hd_102, \
                         hd_105, hd_114, hd_117, id_6, id_9, id_36, id_39, id_48, id_51, \
                         id_90, id_93, id_102, id_105, id_114, id_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_72 * ab_x[k] * hd_6[k]
                  - f_71 * ab_x[k] * hd_9[k]
                  + f_74 * ab_x[k] * hd_36[k]
                  - f_73 * ab_x[k] * hd_39[k]
                  - f_76 * ab_x[k] * hd_48[k]
                  + f_75 * ab_x[k] * hd_51[k]
                  + f_72 * ab_x[k] * hd_90[k]
                  - f_71 * ab_x[k] * hd_93[k]
                  - f_76 * ab_x[k] * hd_102[k]
                  + f_75 * ab_x[k] * hd_105[k]
                  + f_78 * ab_x[k] * hd_114[k]
                  - f_77 * ab_x[k] * hd_117[k]
                  + f_72 * id_6[k]
                  - f_71 * id_9[k]
                  + f_74 * id_36[k]
                  - f_73 * id_39[k]
                  - f_76 * id_48[k]
                  + f_75 * id_51[k]
                  + f_72 * id_90[k]
                  - f_71 * id_93[k]
                  - f_76 * id_102[k]
                  + f_75 * id_105[k]
                  + f_78 * id_114[k]
                  - f_77 * id_117[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_13, hd_15, hd_43, hd_45, hd_55, hd_57, hd_97, hd_99, \
                         hd_109, hd_111, hd_121, hd_123, id_13, id_27, id_43, id_55, id_69, \
                         id_81, id_97, id_109, id_121, id_135, id_147, \
                         id_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_93 * ab_x[k] * hd_13[k]
                  - f_94 * ab_y[k] * hd_15[k]
                  + f_95 * ab_x[k] * hd_43[k]
                  - f_96 * ab_y[k] * hd_45[k]
                  - f_97 * ab_x[k] * hd_55[k]
                  + f_98 * ab_y[k] * hd_57[k]
                  + f_93 * ab_x[k] * hd_97[k]
                  - f_94 * ab_y[k] * hd_99[k]
                  - f_97 * ab_x[k] * hd_109[k]
                  + f_98 * ab_y[k] * hd_111[k]
                  + f_82 * ab_x[k] * hd_121[k]
                  - f_99 * ab_y[k] * hd_123[k]
                  + f_93 * id_13[k]
                  - f_94 * id_27[k]
                  + f_95 * id_43[k]
                  - f_97 * id_55[k]
                  - f_96 * id_69[k]
                  + f_98 * id_81[k]
                  + f_93 * id_97[k]
                  - f_97 * id_109[k]
                  + f_82 * id_121[k]
                  - f_94 * id_135[k]
                  + f_98 * id_147[k]
                  - f_99 * id_159[k];
    }

#pragma omp simd aligned(ab_x, hd_16, hd_46, hd_58, hd_100, hd_112, hd_124, id_16, id_46, \
                         id_58, id_100, id_112, id_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_100 * ab_x[k] * hd_16[k]
                  + f_101 * ab_x[k] * hd_46[k]
                  - f_102 * ab_x[k] * hd_58[k]
                  + f_100 * ab_x[k] * hd_100[k]
                  - f_102 * ab_x[k] * hd_112[k]
                  + f_92 * ab_x[k] * hd_124[k]
                  + f_100 * id_16[k]
                  + f_101 * id_46[k]
                  - f_102 * id_58[k]
                  + f_100 * id_100[k]
                  - f_102 * id_112[k]
                  + f_92 * id_124[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_13, hd_15, hd_17, hd_43, hd_45, hd_47, hd_55, hd_57, \
                         hd_59, hd_97, hd_99, hd_101, hd_109, hd_111, hd_113, hd_121, hd_123, \
                         hd_125, id_13, id_27, id_29, id_43, id_55, id_69, id_71, id_81, \
                         id_83, id_97, id_109, id_121, id_135, id_137, id_147, id_149, id_159, \
                         id_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_71 * ab_x[k] * hd_13[k]
                  - f_71 * ab_y[k] * hd_15[k]
                  + f_76 * ab_y[k] * hd_17[k]
                  - f_73 * ab_x[k] * hd_43[k]
                  - f_73 * ab_y[k] * hd_45[k]
                  + f_77 * ab_y[k] * hd_47[k]
                  + f_78 * ab_x[k] * hd_55[k]
                  + f_78 * ab_y[k] * hd_57[k]
                  - f_103 * ab_y[k] * hd_59[k]
                  - f_71 * ab_x[k] * hd_97[k]
                  - f_71 * ab_y[k] * hd_99[k]
                  + f_76 * ab_y[k] * hd_101[k]
                  + f_78 * ab_x[k] * hd_109[k]
                  + f_78 * ab_y[k] * hd_111[k]
                  - f_103 * ab_y[k] * hd_113[k]
                  - f_104 * ab_x[k] * hd_121[k]
                  - f_104 * ab_y[k] * hd_123[k]
                  + f_105 * ab_y[k] * hd_125[k]
                  - f_71 * id_13[k]
                  - f_71 * id_27[k]
                  + f_76 * id_29[k]
                  - f_73 * id_43[k]
                  + f_78 * id_55[k]
                  - f_73 * id_69[k]
                  + f_77 * id_71[k]
                  + f_78 * id_81[k]
                  - f_103 * id_83[k]
                  - f_71 * id_97[k]
                  + f_78 * id_109[k]
                  - f_104 * id_121[k]
                  - f_71 * id_135[k]
                  + f_76 * id_137[k]
                  + f_78 * id_147[k]
                  - f_103 * id_149[k]
                  - f_104 * id_159[k]
                  + f_105 * id_161[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_14, hd_16, hd_17, hd_44, hd_46, hd_47, hd_56, \
                         hd_58, hd_59, hd_98, hd_100, hd_101, hd_110, hd_112, hd_113, hd_122, \
                         hd_124, hd_125, id_14, id_28, id_35, id_44, id_56, id_70, id_77, \
                         id_82, id_89, id_98, id_110, id_122, id_136, id_143, id_148, id_155, \
                         id_160, id_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -2.8125 * ab_x[k] * hd_14[k]
                  - 2.8125 * ab_y[k] * hd_16[k]
                  + 1.875 * ab_z[k] * hd_17[k]
                  - 5.625 * ab_x[k] * hd_44[k]
                  - 5.625 * ab_y[k] * hd_46[k]
                  + 3.75 * ab_z[k] * hd_47[k]
                  + 7.5 * ab_x[k] * hd_56[k]
                  + 7.5 * ab_y[k] * hd_58[k]
                  - 5.0 * ab_z[k] * hd_59[k]
                  - 2.8125 * ab_x[k] * hd_98[k]
                  - 2.8125 * ab_y[k] * hd_100[k]
                  + 1.875 * ab_z[k] * hd_101[k]
                  + 7.5 * ab_x[k] * hd_110[k]
                  + 7.5 * ab_y[k] * hd_112[k]
                  - 5.0 * ab_z[k] * hd_113[k]
                  - 1.5 * ab_x[k] * hd_122[k]
                  - 1.5 * ab_y[k] * hd_124[k]
                  + ab_z[k] * hd_125[k]
                  - 2.8125 * id_14[k]
                  - 2.8125 * id_28[k]
                  + 1.875 * id_35[k]
                  - 5.625 * id_44[k]
                  + 7.5 * id_56[k]
                  - 5.625 * id_70[k]
                  + 3.75 * id_77[k]
                  + 7.5 * id_82[k]
                  - 5.0 * id_89[k]
                  - 2.8125 * id_98[k]
                  + 7.5 * id_110[k]
                  - 1.5 * id_122[k]
                  - 2.8125 * id_136[k]
                  + 1.875 * id_143[k]
                  + 7.5 * id_148[k]
                  - 5.0 * id_155[k]
                  - 1.5 * id_160[k]
                  + id_167[k];
    }

#pragma omp simd aligned(ab_x, hd_12, hd_15, hd_17, hd_42, hd_45, hd_47, hd_54, hd_57, hd_59, \
                         hd_96, hd_99, hd_101, hd_108, hd_111, hd_113, hd_120, hd_123, hd_125, \
                         id_12, id_15, id_17, id_42, id_45, id_47, id_54, id_57, id_59, id_96, \
                         id_99, id_101, id_108, id_111, id_113, id_120, id_123, \
                         id_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_71 * ab_x[k] * hd_12[k]
                  - f_71 * ab_x[k] * hd_15[k]
                  + f_76 * ab_x[k] * hd_17[k]
                  - f_73 * ab_x[k] * hd_42[k]
                  - f_73 * ab_x[k] * hd_45[k]
                  + f_77 * ab_x[k] * hd_47[k]
                  + f_78 * ab_x[k] * hd_54[k]
                  + f_78 * ab_x[k] * hd_57[k]
                  - f_103 * ab_x[k] * hd_59[k]
                  - f_71 * ab_x[k] * hd_96[k]
                  - f_71 * ab_x[k] * hd_99[k]
                  + f_76 * ab_x[k] * hd_101[k]
                  + f_78 * ab_x[k] * hd_108[k]
                  + f_78 * ab_x[k] * hd_111[k]
                  - f_103 * ab_x[k] * hd_113[k]
                  - f_104 * ab_x[k] * hd_120[k]
                  - f_104 * ab_x[k] * hd_123[k]
                  + f_105 * ab_x[k] * hd_125[k]
                  - f_71 * id_12[k]
                  - f_71 * id_15[k]
                  + f_76 * id_17[k]
                  - f_73 * id_42[k]
                  - f_73 * id_45[k]
                  + f_77 * id_47[k]
                  + f_78 * id_54[k]
                  + f_78 * id_57[k]
                  - f_103 * id_59[k]
                  - f_71 * id_96[k]
                  - f_71 * id_99[k]
                  + f_76 * id_101[k]
                  + f_78 * id_108[k]
                  + f_78 * id_111[k]
                  - f_103 * id_113[k]
                  - f_104 * id_120[k]
                  - f_104 * id_123[k]
                  + f_105 * id_125[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_14, hd_16, hd_44, hd_46, hd_56, hd_58, hd_98, hd_100, \
                         hd_110, hd_112, hd_122, hd_124, id_14, id_28, id_44, id_56, id_70, \
                         id_82, id_98, id_110, id_122, id_136, id_148, \
                         id_160 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_106 * ab_x[k] * hd_14[k]
                  - f_106 * ab_y[k] * hd_16[k]
                  + f_100 * ab_x[k] * hd_44[k]
                  - f_100 * ab_y[k] * hd_46[k]
                  - f_107 * ab_x[k] * hd_56[k]
                  + f_107 * ab_y[k] * hd_58[k]
                  + f_106 * ab_x[k] * hd_98[k]
                  - f_106 * ab_y[k] * hd_100[k]
                  - f_107 * ab_x[k] * hd_110[k]
                  + f_107 * ab_y[k] * hd_112[k]
                  + f_108 * ab_x[k] * hd_122[k]
                  - f_108 * ab_y[k] * hd_124[k]
                  + f_106 * id_14[k]
                  - f_106 * id_28[k]
                  + f_100 * id_44[k]
                  - f_107 * id_56[k]
                  - f_100 * id_70[k]
                  + f_107 * id_82[k]
                  + f_106 * id_98[k]
                  - f_107 * id_110[k]
                  + f_108 * id_122[k]
                  - f_106 * id_136[k]
                  + f_107 * id_148[k]
                  - f_108 * id_160[k];
    }

#pragma omp simd aligned(ab_x, hd_12, hd_15, hd_42, hd_45, hd_54, hd_57, hd_96, hd_99, hd_108, \
                         hd_111, hd_120, hd_123, id_12, id_15, id_42, id_45, id_54, id_57, \
                         id_96, id_99, id_108, id_111, id_120, id_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_94 * ab_x[k] * hd_12[k]
                  - f_93 * ab_x[k] * hd_15[k]
                  + f_96 * ab_x[k] * hd_42[k]
                  - f_95 * ab_x[k] * hd_45[k]
                  - f_98 * ab_x[k] * hd_54[k]
                  + f_97 * ab_x[k] * hd_57[k]
                  + f_94 * ab_x[k] * hd_96[k]
                  - f_93 * ab_x[k] * hd_99[k]
                  - f_98 * ab_x[k] * hd_108[k]
                  + f_97 * ab_x[k] * hd_111[k]
                  + f_99 * ab_x[k] * hd_120[k]
                  - f_82 * ab_x[k] * hd_123[k]
                  + f_94 * id_12[k]
                  - f_93 * id_15[k]
                  + f_96 * id_42[k]
                  - f_95 * id_45[k]
                  - f_98 * id_54[k]
                  + f_97 * id_57[k]
                  + f_94 * id_96[k]
                  - f_93 * id_99[k]
                  - f_98 * id_108[k]
                  + f_97 * id_111[k]
                  + f_99 * id_120[k]
                  - f_82 * id_123[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_1, hd_3, hd_19, hd_21, hd_31, hd_33, hd_61, hd_63, \
                         hd_73, hd_75, hd_85, hd_87, id_1, id_9, id_19, id_31, id_39, id_51, \
                         id_61, id_73, id_85, id_93, id_105, id_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_71 * ab_x[k] * hd_1[k]
                  - f_72 * ab_y[k] * hd_3[k]
                  + f_73 * ab_x[k] * hd_19[k]
                  - f_74 * ab_y[k] * hd_21[k]
                  - f_75 * ab_x[k] * hd_31[k]
                  + f_76 * ab_y[k] * hd_33[k]
                  + f_71 * ab_x[k] * hd_61[k]
                  - f_72 * ab_y[k] * hd_63[k]
                  - f_75 * ab_x[k] * hd_73[k]
                  + f_76 * ab_y[k] * hd_75[k]
                  + f_77 * ab_x[k] * hd_85[k]
                  - f_78 * ab_y[k] * hd_87[k]
                  + f_71 * id_1[k]
                  - f_72 * id_9[k]
                  + f_73 * id_19[k]
                  - f_75 * id_31[k]
                  - f_74 * id_39[k]
                  + f_76 * id_51[k]
                  + f_71 * id_61[k]
                  - f_75 * id_73[k]
                  + f_77 * id_85[k]
                  - f_72 * id_93[k]
                  + f_76 * id_105[k]
                  - f_78 * id_117[k];
    }

#pragma omp simd aligned(ab_x, hd_4, hd_22, hd_34, hd_64, hd_76, hd_88, id_4, id_22, id_34, \
                         id_64, id_76, id_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = 1.875 * ab_x[k] * hd_4[k]
                  + 3.75 * ab_x[k] * hd_22[k]
                  - 22.5 * ab_x[k] * hd_34[k]
                  + 1.875 * ab_x[k] * hd_64[k]
                  - 22.5 * ab_x[k] * hd_76[k]
                  + 15.0 * ab_x[k] * hd_88[k]
                  + 1.875 * id_4[k]
                  + 3.75 * id_22[k]
                  - 22.5 * id_34[k]
                  + 1.875 * id_64[k]
                  - 22.5 * id_76[k]
                  + 15.0 * id_88[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_1, hd_3, hd_5, hd_19, hd_21, hd_23, hd_31, hd_33, \
                         hd_35, hd_61, hd_63, hd_65, hd_73, hd_75, hd_77, hd_85, hd_87, hd_89, \
                         id_1, id_9, id_11, id_19, id_31, id_39, id_41, id_51, id_53, id_61, \
                         id_73, id_85, id_93, id_95, id_105, id_107, id_117, \
                         id_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_79 * ab_x[k] * hd_1[k]
                  - f_79 * ab_y[k] * hd_3[k]
                  + f_80 * ab_y[k] * hd_5[k]
                  - f_81 * ab_x[k] * hd_19[k]
                  - f_81 * ab_y[k] * hd_21[k]
                  + f_82 * ab_y[k] * hd_23[k]
                  + f_83 * ab_x[k] * hd_31[k]
                  + f_83 * ab_y[k] * hd_33[k]
                  - f_84 * ab_y[k] * hd_35[k]
                  - f_79 * ab_x[k] * hd_61[k]
                  - f_79 * ab_y[k] * hd_63[k]
                  + f_80 * ab_y[k] * hd_65[k]
                  + f_83 * ab_x[k] * hd_73[k]
                  + f_83 * ab_y[k] * hd_75[k]
                  - f_84 * ab_y[k] * hd_77[k]
                  - f_82 * ab_x[k] * hd_85[k]
                  - f_82 * ab_y[k] * hd_87[k]
                  + f_85 * ab_y[k] * hd_89[k]
                  - f_79 * id_1[k]
                  - f_79 * id_9[k]
                  + f_80 * id_11[k]
                  - f_81 * id_19[k]
                  + f_83 * id_31[k]
                  - f_81 * id_39[k]
                  + f_82 * id_41[k]
                  + f_83 * id_51[k]
                  - f_84 * id_53[k]
                  - f_79 * id_61[k]
                  + f_83 * id_73[k]
                  - f_82 * id_85[k]
                  - f_79 * id_93[k]
                  + f_80 * id_95[k]
                  + f_83 * id_105[k]
                  - f_84 * id_107[k]
                  - f_82 * id_117[k]
                  + f_85 * id_119[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_2, hd_4, hd_5, hd_20, hd_22, hd_23, hd_32, \
                         hd_34, hd_35, hd_62, hd_64, hd_65, hd_74, hd_76, hd_77, hd_86, hd_88, \
                         hd_89, id_2, id_10, id_17, id_20, id_32, id_40, id_47, id_52, id_59, \
                         id_62, id_74, id_86, id_94, id_101, id_106, id_113, id_118, \
                         id_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_86 * ab_x[k] * hd_2[k]
                  - f_86 * ab_y[k] * hd_4[k]
                  + f_87 * ab_z[k] * hd_5[k]
                  - f_88 * ab_x[k] * hd_20[k]
                  - f_88 * ab_y[k] * hd_22[k]
                  + f_89 * ab_z[k] * hd_23[k]
                  + f_90 * ab_x[k] * hd_32[k]
                  + f_90 * ab_y[k] * hd_34[k]
                  - f_91 * ab_z[k] * hd_35[k]
                  - f_86 * ab_x[k] * hd_62[k]
                  - f_86 * ab_y[k] * hd_64[k]
                  + f_87 * ab_z[k] * hd_65[k]
                  + f_90 * ab_x[k] * hd_74[k]
                  + f_90 * ab_y[k] * hd_76[k]
                  - f_91 * ab_z[k] * hd_77[k]
                  - f_91 * ab_x[k] * hd_86[k]
                  - f_91 * ab_y[k] * hd_88[k]
                  + f_92 * ab_z[k] * hd_89[k]
                  - f_86 * id_2[k]
                  - f_86 * id_10[k]
                  + f_87 * id_17[k]
                  - f_88 * id_20[k]
                  + f_90 * id_32[k]
                  - f_88 * id_40[k]
                  + f_89 * id_47[k]
                  + f_90 * id_52[k]
                  - f_91 * id_59[k]
                  - f_86 * id_62[k]
                  + f_90 * id_74[k]
                  - f_91 * id_86[k]
                  - f_86 * id_94[k]
                  + f_87 * id_101[k]
                  + f_90 * id_106[k]
                  - f_91 * id_113[k]
                  - f_91 * id_118[k]
                  + f_92 * id_125[k];
    }

#pragma omp simd aligned(ab_x, hd_0, hd_3, hd_5, hd_18, hd_21, hd_23, hd_30, hd_33, hd_35, \
                         hd_60, hd_63, hd_65, hd_72, hd_75, hd_77, hd_84, hd_87, hd_89, id_0, \
                         id_3, id_5, id_18, id_21, id_23, id_30, id_33, id_35, id_60, id_63, \
                         id_65, id_72, id_75, id_77, id_84, id_87, \
                         id_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_79 * ab_x[k] * hd_0[k]
                  - f_79 * ab_x[k] * hd_3[k]
                  + f_80 * ab_x[k] * hd_5[k]
                  - f_81 * ab_x[k] * hd_18[k]
                  - f_81 * ab_x[k] * hd_21[k]
                  + f_82 * ab_x[k] * hd_23[k]
                  + f_83 * ab_x[k] * hd_30[k]
                  + f_83 * ab_x[k] * hd_33[k]
                  - f_84 * ab_x[k] * hd_35[k]
                  - f_79 * ab_x[k] * hd_60[k]
                  - f_79 * ab_x[k] * hd_63[k]
                  + f_80 * ab_x[k] * hd_65[k]
                  + f_83 * ab_x[k] * hd_72[k]
                  + f_83 * ab_x[k] * hd_75[k]
                  - f_84 * ab_x[k] * hd_77[k]
                  - f_82 * ab_x[k] * hd_84[k]
                  - f_82 * ab_x[k] * hd_87[k]
                  + f_85 * ab_x[k] * hd_89[k]
                  - f_79 * id_0[k]
                  - f_79 * id_3[k]
                  + f_80 * id_5[k]
                  - f_81 * id_18[k]
                  - f_81 * id_21[k]
                  + f_82 * id_23[k]
                  + f_83 * id_30[k]
                  + f_83 * id_33[k]
                  - f_84 * id_35[k]
                  - f_79 * id_60[k]
                  - f_79 * id_63[k]
                  + f_80 * id_65[k]
                  + f_83 * id_72[k]
                  + f_83 * id_75[k]
                  - f_84 * id_77[k]
                  - f_82 * id_84[k]
                  - f_82 * id_87[k]
                  + f_85 * id_89[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_2, hd_4, hd_20, hd_22, hd_32, hd_34, hd_62, hd_64, \
                         hd_74, hd_76, hd_86, hd_88, id_2, id_10, id_20, id_32, id_40, id_52, \
                         id_62, id_74, id_86, id_94, id_106, id_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = 0.9375 * ab_x[k] * hd_2[k]
                  - 0.9375 * ab_y[k] * hd_4[k]
                  + 1.875 * ab_x[k] * hd_20[k]
                  - 1.875 * ab_y[k] * hd_22[k]
                  - 11.25 * ab_x[k] * hd_32[k]
                  + 11.25 * ab_y[k] * hd_34[k]
                  + 0.9375 * ab_x[k] * hd_62[k]
                  - 0.9375 * ab_y[k] * hd_64[k]
                  - 11.25 * ab_x[k] * hd_74[k]
                  + 11.25 * ab_y[k] * hd_76[k]
                  + 7.5 * ab_x[k] * hd_86[k]
                  - 7.5 * ab_y[k] * hd_88[k]
                  + 0.9375 * id_2[k]
                  - 0.9375 * id_10[k]
                  + 1.875 * id_20[k]
                  - 11.25 * id_32[k]
                  - 1.875 * id_40[k]
                  + 11.25 * id_52[k]
                  + 0.9375 * id_62[k]
                  - 11.25 * id_74[k]
                  + 7.5 * id_86[k]
                  - 0.9375 * id_94[k]
                  + 11.25 * id_106[k]
                  - 7.5 * id_118[k];
    }

#pragma omp simd aligned(ab_x, hd_0, hd_3, hd_18, hd_21, hd_30, hd_33, hd_60, hd_63, hd_72, \
                         hd_75, hd_84, hd_87, id_0, id_3, id_18, id_21, id_30, id_33, id_60, \
                         id_63, id_72, id_75, id_84, id_87 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_72 * ab_x[k] * hd_0[k]
                  - f_71 * ab_x[k] * hd_3[k]
                  + f_74 * ab_x[k] * hd_18[k]
                  - f_73 * ab_x[k] * hd_21[k]
                  - f_76 * ab_x[k] * hd_30[k]
                  + f_75 * ab_x[k] * hd_33[k]
                  + f_72 * ab_x[k] * hd_60[k]
                  - f_71 * ab_x[k] * hd_63[k]
                  - f_76 * ab_x[k] * hd_72[k]
                  + f_75 * ab_x[k] * hd_75[k]
                  + f_78 * ab_x[k] * hd_84[k]
                  - f_77 * ab_x[k] * hd_87[k]
                  + f_72 * id_0[k]
                  - f_71 * id_3[k]
                  + f_74 * id_18[k]
                  - f_73 * id_21[k]
                  - f_76 * id_30[k]
                  + f_75 * id_33[k]
                  + f_72 * id_60[k]
                  - f_71 * id_63[k]
                  - f_76 * id_72[k]
                  + f_75 * id_75[k]
                  + f_78 * id_84[k]
                  - f_77 * id_87[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_13, hd_15, hd_55, hd_57, hd_97, hd_99, hd_109, hd_111, \
                         id_13, id_27, id_55, id_81, id_97, id_109, id_135, \
                         id_147 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_37 * ab_x[k] * hd_13[k]
                  + f_40 * ab_y[k] * hd_15[k]
                  + f_64 * ab_x[k] * hd_55[k]
                  - f_38 * ab_y[k] * hd_57[k]
                  + f_37 * ab_x[k] * hd_97[k]
                  - f_40 * ab_y[k] * hd_99[k]
                  - f_64 * ab_x[k] * hd_109[k]
                  + f_38 * ab_y[k] * hd_111[k]
                  - f_37 * id_13[k]
                  + f_40 * id_27[k]
                  + f_64 * id_55[k]
                  - f_38 * id_81[k]
                  + f_37 * id_97[k]
                  - f_64 * id_109[k]
                  - f_40 * id_135[k]
                  + f_38 * id_147[k];
    }

#pragma omp simd aligned(ab_x, hd_16, hd_58, hd_100, hd_112, id_16, id_58, id_100, \
                         id_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_34 * ab_x[k] * hd_16[k]
                  + f_65 * ab_x[k] * hd_58[k]
                  + f_34 * ab_x[k] * hd_100[k]
                  - f_65 * ab_x[k] * hd_112[k]
                  - f_34 * id_16[k]
                  + f_65 * id_58[k]
                  + f_34 * id_100[k]
                  - f_65 * id_112[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_13, hd_15, hd_17, hd_55, hd_57, hd_59, hd_97, hd_99, \
                         hd_101, hd_109, hd_111, hd_113, id_13, id_27, id_29, id_55, id_81, \
                         id_83, id_97, id_109, id_135, id_137, id_147, \
                         id_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_52 * ab_x[k] * hd_13[k]
                  + f_52 * ab_y[k] * hd_15[k]
                  - f_58 * ab_y[k] * hd_17[k]
                  - f_67 * ab_x[k] * hd_55[k]
                  - f_67 * ab_y[k] * hd_57[k]
                  + f_55 * ab_y[k] * hd_59[k]
                  - f_52 * ab_x[k] * hd_97[k]
                  - f_52 * ab_y[k] * hd_99[k]
                  + f_58 * ab_y[k] * hd_101[k]
                  + f_67 * ab_x[k] * hd_109[k]
                  + f_67 * ab_y[k] * hd_111[k]
                  - f_55 * ab_y[k] * hd_113[k]
                  + f_52 * id_13[k]
                  + f_52 * id_27[k]
                  - f_58 * id_29[k]
                  - f_67 * id_55[k]
                  - f_67 * id_81[k]
                  + f_55 * id_83[k]
                  - f_52 * id_97[k]
                  + f_67 * id_109[k]
                  - f_52 * id_135[k]
                  + f_58 * id_137[k]
                  + f_67 * id_147[k]
                  - f_55 * id_149[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_14, hd_16, hd_17, hd_56, hd_58, hd_59, hd_98, \
                         hd_100, hd_101, hd_110, hd_112, hd_113, id_14, id_28, id_35, id_56, \
                         id_82, id_89, id_98, id_110, id_136, id_143, id_148, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_43 * ab_x[k] * hd_14[k]
                  + f_43 * ab_y[k] * hd_16[k]
                  - f_45 * ab_z[k] * hd_17[k]
                  - f_46 * ab_x[k] * hd_56[k]
                  - f_46 * ab_y[k] * hd_58[k]
                  + f_69 * ab_z[k] * hd_59[k]
                  - f_43 * ab_x[k] * hd_98[k]
                  - f_43 * ab_y[k] * hd_100[k]
                  + f_45 * ab_z[k] * hd_101[k]
                  + f_46 * ab_x[k] * hd_110[k]
                  + f_46 * ab_y[k] * hd_112[k]
                  - f_69 * ab_z[k] * hd_113[k]
                  + f_43 * id_14[k]
                  + f_43 * id_28[k]
                  - f_45 * id_35[k]
                  - f_46 * id_56[k]
                  - f_46 * id_82[k]
                  + f_69 * id_89[k]
                  - f_43 * id_98[k]
                  + f_46 * id_110[k]
                  - f_43 * id_136[k]
                  + f_45 * id_143[k]
                  + f_46 * id_148[k]
                  - f_69 * id_155[k];
    }

#pragma omp simd aligned(ab_x, hd_12, hd_15, hd_17, hd_54, hd_57, hd_59, hd_96, hd_99, hd_101, \
                         hd_108, hd_111, hd_113, id_12, id_15, id_17, id_54, id_57, id_59, \
                         id_96, id_99, id_101, id_108, id_111, id_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_52 * ab_x[k] * hd_12[k]
                  + f_52 * ab_x[k] * hd_15[k]
                  - f_58 * ab_x[k] * hd_17[k]
                  - f_67 * ab_x[k] * hd_54[k]
                  - f_67 * ab_x[k] * hd_57[k]
                  + f_55 * ab_x[k] * hd_59[k]
                  - f_52 * ab_x[k] * hd_96[k]
                  - f_52 * ab_x[k] * hd_99[k]
                  + f_58 * ab_x[k] * hd_101[k]
                  + f_67 * ab_x[k] * hd_108[k]
                  + f_67 * ab_x[k] * hd_111[k]
                  - f_55 * ab_x[k] * hd_113[k]
                  + f_52 * id_12[k]
                  + f_52 * id_15[k]
                  - f_58 * id_17[k]
                  - f_67 * id_54[k]
                  - f_67 * id_57[k]
                  + f_55 * id_59[k]
                  - f_52 * id_96[k]
                  - f_52 * id_99[k]
                  + f_58 * id_101[k]
                  + f_67 * id_108[k]
                  + f_67 * id_111[k]
                  - f_55 * id_113[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_14, hd_16, hd_56, hd_58, hd_98, hd_100, hd_110, \
                         hd_112, id_14, id_28, id_56, id_82, id_98, id_110, id_136, \
                         id_148 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_109 * ab_x[k] * hd_14[k]
                  + f_109 * ab_y[k] * hd_16[k]
                  + f_34 * ab_x[k] * hd_56[k]
                  - f_34 * ab_y[k] * hd_58[k]
                  + f_109 * ab_x[k] * hd_98[k]
                  - f_109 * ab_y[k] * hd_100[k]
                  - f_34 * ab_x[k] * hd_110[k]
                  + f_34 * ab_y[k] * hd_112[k]
                  - f_109 * id_14[k]
                  + f_109 * id_28[k]
                  + f_34 * id_56[k]
                  - f_34 * id_82[k]
                  + f_109 * id_98[k]
                  - f_34 * id_110[k]
                  - f_109 * id_136[k]
                  + f_34 * id_148[k];
    }

#pragma omp simd aligned(ab_x, hd_12, hd_15, hd_54, hd_57, hd_96, hd_99, hd_108, hd_111, \
                         id_12, id_15, id_54, id_57, id_96, id_99, id_108, \
                         id_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_40 * ab_x[k] * hd_12[k]
                  + f_37 * ab_x[k] * hd_15[k]
                  + f_38 * ab_x[k] * hd_54[k]
                  - f_64 * ab_x[k] * hd_57[k]
                  + f_40 * ab_x[k] * hd_96[k]
                  - f_37 * ab_x[k] * hd_99[k]
                  - f_38 * ab_x[k] * hd_108[k]
                  + f_64 * ab_x[k] * hd_111[k]
                  - f_40 * id_12[k]
                  + f_37 * id_15[k]
                  + f_38 * id_54[k]
                  - f_64 * id_57[k]
                  + f_40 * id_96[k]
                  - f_37 * id_99[k]
                  - f_38 * id_108[k]
                  + f_64 * id_111[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_1, hd_3, hd_19, hd_21, hd_31, hd_33, hd_61, hd_63, \
                         hd_73, hd_75, id_1, id_9, id_19, id_31, id_39, id_51, id_61, id_73, \
                         id_93, id_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_30 * ab_x[k] * hd_1[k]
                  + f_35 * ab_y[k] * hd_3[k]
                  + f_31 * ab_x[k] * hd_19[k]
                  - f_32 * ab_y[k] * hd_21[k]
                  + f_34 * ab_x[k] * hd_31[k]
                  - f_36 * ab_y[k] * hd_33[k]
                  + f_29 * ab_x[k] * hd_61[k]
                  - f_30 * ab_y[k] * hd_63[k]
                  - f_33 * ab_x[k] * hd_73[k]
                  + f_34 * ab_y[k] * hd_75[k]
                  - f_30 * id_1[k]
                  + f_35 * id_9[k]
                  + f_31 * id_19[k]
                  + f_34 * id_31[k]
                  - f_32 * id_39[k]
                  - f_36 * id_51[k]
                  + f_29 * id_61[k]
                  - f_33 * id_73[k]
                  - f_30 * id_93[k]
                  + f_34 * id_105[k];
    }

#pragma omp simd aligned(ab_x, hd_4, hd_22, hd_34, hd_64, hd_76, id_4, id_22, id_34, id_64, \
                         id_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_40 * ab_x[k] * hd_4[k]
                  + f_38 * ab_x[k] * hd_22[k]
                  + f_41 * ab_x[k] * hd_34[k]
                  + f_37 * ab_x[k] * hd_64[k]
                  - f_39 * ab_x[k] * hd_76[k]
                  - f_40 * id_4[k]
                  + f_38 * id_22[k]
                  + f_41 * id_34[k]
                  + f_37 * id_64[k]
                  - f_39 * id_76[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_1, hd_3, hd_5, hd_19, hd_21, hd_23, hd_31, hd_33, \
                         hd_35, hd_61, hd_63, hd_65, hd_73, hd_75, hd_77, id_1, id_9, id_11, \
                         id_19, id_31, id_39, id_41, id_51, id_53, id_61, id_73, id_93, id_95, \
                         id_105, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_48 * ab_x[k] * hd_1[k]
                  + f_48 * ab_y[k] * hd_3[k]
                  - f_49 * ab_y[k] * hd_5[k]
                  - f_44 * ab_x[k] * hd_19[k]
                  - f_44 * ab_y[k] * hd_21[k]
                  + f_45 * ab_y[k] * hd_23[k]
                  - f_45 * ab_x[k] * hd_31[k]
                  - f_45 * ab_y[k] * hd_33[k]
                  + f_50 * ab_y[k] * hd_35[k]
                  - f_42 * ab_x[k] * hd_61[k]
                  - f_42 * ab_y[k] * hd_63[k]
                  + f_43 * ab_y[k] * hd_65[k]
                  + f_46 * ab_x[k] * hd_73[k]
                  + f_46 * ab_y[k] * hd_75[k]
                  - f_47 * ab_y[k] * hd_77[k]
                  + f_48 * id_1[k]
                  + f_48 * id_9[k]
                  - f_49 * id_11[k]
                  - f_44 * id_19[k]
                  - f_45 * id_31[k]
                  - f_44 * id_39[k]
                  + f_45 * id_41[k]
                  - f_45 * id_51[k]
                  + f_50 * id_53[k]
                  - f_42 * id_61[k]
                  + f_46 * id_73[k]
                  - f_42 * id_93[k]
                  + f_43 * id_95[k]
                  + f_46 * id_105[k]
                  - f_47 * id_107[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_2, hd_4, hd_5, hd_20, hd_22, hd_23, hd_32, \
                         hd_34, hd_35, hd_62, hd_64, hd_65, hd_74, hd_76, hd_77, id_2, id_10, \
                         id_17, id_20, id_32, id_40, id_47, id_52, id_59, id_62, id_74, id_94, \
                         id_101, id_106, id_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_56 * ab_x[k] * hd_2[k]
                  + f_56 * ab_y[k] * hd_4[k]
                  - f_57 * ab_z[k] * hd_5[k]
                  - f_52 * ab_x[k] * hd_20[k]
                  - f_52 * ab_y[k] * hd_22[k]
                  + f_53 * ab_z[k] * hd_23[k]
                  - f_58 * ab_x[k] * hd_32[k]
                  - f_58 * ab_y[k] * hd_34[k]
                  + f_59 * ab_z[k] * hd_35[k]
                  - f_51 * ab_x[k] * hd_62[k]
                  - f_51 * ab_y[k] * hd_64[k]
                  + f_52 * ab_z[k] * hd_65[k]
                  + f_54 * ab_x[k] * hd_74[k]
                  + f_54 * ab_y[k] * hd_76[k]
                  - f_55 * ab_z[k] * hd_77[k]
                  + f_56 * id_2[k]
                  + f_56 * id_10[k]
                  - f_57 * id_17[k]
                  - f_52 * id_20[k]
                  - f_58 * id_32[k]
                  - f_52 * id_40[k]
                  + f_53 * id_47[k]
                  - f_58 * id_52[k]
                  + f_59 * id_59[k]
                  - f_51 * id_62[k]
                  + f_54 * id_74[k]
                  - f_51 * id_94[k]
                  + f_52 * id_101[k]
                  + f_54 * id_106[k]
                  - f_55 * id_113[k];
    }

#pragma omp simd aligned(ab_x, hd_0, hd_3, hd_5, hd_18, hd_21, hd_23, hd_30, hd_33, hd_35, \
                         hd_60, hd_63, hd_65, hd_72, hd_75, hd_77, id_0, id_3, id_5, id_18, \
                         id_21, id_23, id_30, id_33, id_35, id_60, id_63, id_65, id_72, id_75, \
                         id_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_48 * ab_x[k] * hd_0[k]
                  + f_48 * ab_x[k] * hd_3[k]
                  - f_49 * ab_x[k] * hd_5[k]
                  - f_44 * ab_x[k] * hd_18[k]
                  - f_44 * ab_x[k] * hd_21[k]
                  + f_45 * ab_x[k] * hd_23[k]
                  - f_45 * ab_x[k] * hd_30[k]
                  - f_45 * ab_x[k] * hd_33[k]
                  + f_50 * ab_x[k] * hd_35[k]
                  - f_42 * ab_x[k] * hd_60[k]
                  - f_42 * ab_x[k] * hd_63[k]
                  + f_43 * ab_x[k] * hd_65[k]
                  + f_46 * ab_x[k] * hd_72[k]
                  + f_46 * ab_x[k] * hd_75[k]
                  - f_47 * ab_x[k] * hd_77[k]
                  + f_48 * id_0[k]
                  + f_48 * id_3[k]
                  - f_49 * id_5[k]
                  - f_44 * id_18[k]
                  - f_44 * id_21[k]
                  + f_45 * id_23[k]
                  - f_45 * id_30[k]
                  - f_45 * id_33[k]
                  + f_50 * id_35[k]
                  - f_42 * id_60[k]
                  - f_42 * id_63[k]
                  + f_43 * id_65[k]
                  + f_46 * id_72[k]
                  + f_46 * id_75[k]
                  - f_47 * id_77[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_2, hd_4, hd_20, hd_22, hd_32, hd_34, hd_62, hd_64, \
                         hd_74, hd_76, id_2, id_10, id_20, id_32, id_40, id_52, id_62, id_74, \
                         id_94, id_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_62 * ab_x[k] * hd_2[k]
                  + f_62 * ab_y[k] * hd_4[k]
                  + f_40 * ab_x[k] * hd_20[k]
                  - f_40 * ab_y[k] * hd_22[k]
                  + f_63 * ab_x[k] * hd_32[k]
                  - f_63 * ab_y[k] * hd_34[k]
                  + f_60 * ab_x[k] * hd_62[k]
                  - f_60 * ab_y[k] * hd_64[k]
                  - f_61 * ab_x[k] * hd_74[k]
                  + f_61 * ab_y[k] * hd_76[k]
                  - f_62 * id_2[k]
                  + f_62 * id_10[k]
                  + f_40 * id_20[k]
                  + f_63 * id_32[k]
                  - f_40 * id_40[k]
                  - f_63 * id_52[k]
                  + f_60 * id_62[k]
                  - f_61 * id_74[k]
                  - f_60 * id_94[k]
                  + f_61 * id_106[k];
    }

#pragma omp simd aligned(ab_x, hd_0, hd_3, hd_18, hd_21, hd_30, hd_33, hd_60, hd_63, hd_72, \
                         hd_75, id_0, id_3, id_18, id_21, id_30, id_33, id_60, id_63, id_72, \
                         id_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_35 * ab_x[k] * hd_0[k]
                  + f_30 * ab_x[k] * hd_3[k]
                  + f_32 * ab_x[k] * hd_18[k]
                  - f_31 * ab_x[k] * hd_21[k]
                  + f_36 * ab_x[k] * hd_30[k]
                  - f_34 * ab_x[k] * hd_33[k]
                  + f_30 * ab_x[k] * hd_60[k]
                  - f_29 * ab_x[k] * hd_63[k]
                  - f_34 * ab_x[k] * hd_72[k]
                  + f_33 * ab_x[k] * hd_75[k]
                  - f_35 * id_0[k]
                  + f_30 * id_3[k]
                  + f_32 * id_18[k]
                  - f_31 * id_21[k]
                  + f_36 * id_30[k]
                  - f_34 * id_33[k]
                  + f_30 * id_60[k]
                  - f_29 * id_63[k]
                  - f_34 * id_72[k]
                  + f_33 * id_75[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_13, hd_15, hd_43, hd_45, hd_97, hd_99, id_13, id_27, \
                         id_43, id_69, id_97, id_135 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_15 * ab_x[k] * hd_13[k]
                  - f_110 * ab_y[k] * hd_15[k]
                  - f_111 * ab_x[k] * hd_43[k]
                  + f_17 * ab_y[k] * hd_45[k]
                  + f_15 * ab_x[k] * hd_97[k]
                  - f_110 * ab_y[k] * hd_99[k]
                  + f_15 * id_13[k]
                  - f_110 * id_27[k]
                  - f_111 * id_43[k]
                  + f_17 * id_69[k]
                  + f_15 * id_97[k]
                  - f_110 * id_135[k];
    }

#pragma omp simd aligned(ab_x, hd_16, hd_46, hd_100, id_16, id_46, \
                         id_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_10 * ab_x[k] * hd_16[k]
                  - f_112 * ab_x[k] * hd_46[k]
                  + f_10 * ab_x[k] * hd_100[k]
                  + f_10 * id_16[k]
                  - f_112 * id_46[k]
                  + f_10 * id_100[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_13, hd_15, hd_17, hd_43, hd_45, hd_47, hd_97, hd_99, \
                         hd_101, id_13, id_27, id_29, id_43, id_69, id_71, id_97, id_135, \
                         id_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_22 * ab_x[k] * hd_13[k]
                  - f_22 * ab_y[k] * hd_15[k]
                  + f_25 * ab_y[k] * hd_17[k]
                  + f_113 * ab_x[k] * hd_43[k]
                  + f_113 * ab_y[k] * hd_45[k]
                  - f_114 * ab_y[k] * hd_47[k]
                  - f_22 * ab_x[k] * hd_97[k]
                  - f_22 * ab_y[k] * hd_99[k]
                  + f_25 * ab_y[k] * hd_101[k]
                  - f_22 * id_13[k]
                  - f_22 * id_27[k]
                  + f_25 * id_29[k]
                  + f_113 * id_43[k]
                  + f_113 * id_69[k]
                  - f_114 * id_71[k]
                  - f_22 * id_97[k]
                  - f_22 * id_135[k]
                  + f_25 * id_137[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_14, hd_16, hd_17, hd_44, hd_46, hd_47, hd_98, \
                         hd_100, hd_101, id_14, id_28, id_35, id_44, id_70, id_77, id_98, \
                         id_136, id_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_115 * ab_x[k] * hd_14[k]
                  - f_115 * ab_y[k] * hd_16[k]
                  + f_116 * ab_z[k] * hd_17[k]
                  + f_117 * ab_x[k] * hd_44[k]
                  + f_117 * ab_y[k] * hd_46[k]
                  - f_27 * ab_z[k] * hd_47[k]
                  - f_115 * ab_x[k] * hd_98[k]
                  - f_115 * ab_y[k] * hd_100[k]
                  + f_116 * ab_z[k] * hd_101[k]
                  - f_115 * id_14[k]
                  - f_115 * id_28[k]
                  + f_116 * id_35[k]
                  + f_117 * id_44[k]
                  + f_117 * id_70[k]
                  - f_27 * id_77[k]
                  - f_115 * id_98[k]
                  - f_115 * id_136[k]
                  + f_116 * id_143[k];
    }

#pragma omp simd aligned(ab_x, hd_12, hd_15, hd_17, hd_42, hd_45, hd_47, hd_96, hd_99, hd_101, \
                         id_12, id_15, id_17, id_42, id_45, id_47, id_96, id_99, \
                         id_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_22 * ab_x[k] * hd_12[k]
                  - f_22 * ab_x[k] * hd_15[k]
                  + f_25 * ab_x[k] * hd_17[k]
                  + f_113 * ab_x[k] * hd_42[k]
                  + f_113 * ab_x[k] * hd_45[k]
                  - f_114 * ab_x[k] * hd_47[k]
                  - f_22 * ab_x[k] * hd_96[k]
                  - f_22 * ab_x[k] * hd_99[k]
                  + f_25 * ab_x[k] * hd_101[k]
                  - f_22 * id_12[k]
                  - f_22 * id_15[k]
                  + f_25 * id_17[k]
                  + f_113 * id_42[k]
                  + f_113 * id_45[k]
                  - f_114 * id_47[k]
                  - f_22 * id_96[k]
                  - f_22 * id_99[k]
                  + f_25 * id_101[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_14, hd_16, hd_44, hd_46, hd_98, hd_100, id_14, id_28, \
                         id_44, id_70, id_98, id_136 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_11 * ab_x[k] * hd_14[k]
                  - f_11 * ab_y[k] * hd_16[k]
                  - f_118 * ab_x[k] * hd_44[k]
                  + f_118 * ab_y[k] * hd_46[k]
                  + f_11 * ab_x[k] * hd_98[k]
                  - f_11 * ab_y[k] * hd_100[k]
                  + f_11 * id_14[k]
                  - f_11 * id_28[k]
                  - f_118 * id_44[k]
                  + f_118 * id_70[k]
                  + f_11 * id_98[k]
                  - f_11 * id_136[k];
    }

#pragma omp simd aligned(ab_x, hd_12, hd_15, hd_42, hd_45, hd_96, hd_99, id_12, id_15, id_42, \
                         id_45, id_96, id_99 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_110 * ab_x[k] * hd_12[k]
                  - f_15 * ab_x[k] * hd_15[k]
                  - f_17 * ab_x[k] * hd_42[k]
                  + f_111 * ab_x[k] * hd_45[k]
                  + f_110 * ab_x[k] * hd_96[k]
                  - f_15 * ab_x[k] * hd_99[k]
                  + f_110 * id_12[k]
                  - f_15 * id_15[k]
                  - f_17 * id_42[k]
                  + f_111 * id_45[k]
                  + f_110 * id_96[k]
                  - f_15 * id_99[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_1, hd_3, hd_19, hd_21, hd_61, hd_63, id_1, id_9, \
                         id_19, id_39, id_61, id_93 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_4 * ab_x[k] * hd_1[k]
                  - f_5 * ab_y[k] * hd_3[k]
                  - f_2 * ab_x[k] * hd_19[k]
                  + f_3 * ab_y[k] * hd_21[k]
                  + f_0 * ab_x[k] * hd_61[k]
                  - f_1 * ab_y[k] * hd_63[k]
                  + f_4 * id_1[k]
                  - f_5 * id_9[k]
                  - f_2 * id_19[k]
                  + f_3 * id_39[k]
                  + f_0 * id_61[k]
                  - f_1 * id_93[k];
    }

#pragma omp simd aligned(ab_x, hd_4, hd_22, hd_64, id_4, id_22, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_8 * ab_x[k] * hd_4[k]
                  - f_7 * ab_x[k] * hd_22[k]
                  + f_6 * ab_x[k] * hd_64[k]
                  + f_8 * id_4[k]
                  - f_7 * id_22[k]
                  + f_6 * id_64[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_1, hd_3, hd_5, hd_19, hd_21, hd_23, hd_61, hd_63, \
                         hd_65, id_1, id_9, id_11, id_19, id_39, id_41, id_61, id_93, \
                         id_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_13 * ab_x[k] * hd_1[k]
                  - f_13 * ab_y[k] * hd_3[k]
                  + f_14 * ab_y[k] * hd_5[k]
                  + f_11 * ab_x[k] * hd_19[k]
                  + f_11 * ab_y[k] * hd_21[k]
                  - f_12 * ab_y[k] * hd_23[k]
                  - f_9 * ab_x[k] * hd_61[k]
                  - f_9 * ab_y[k] * hd_63[k]
                  + f_10 * ab_y[k] * hd_65[k]
                  - f_13 * id_1[k]
                  - f_13 * id_9[k]
                  + f_14 * id_11[k]
                  + f_11 * id_19[k]
                  + f_11 * id_39[k]
                  - f_12 * id_41[k]
                  - f_9 * id_61[k]
                  - f_9 * id_93[k]
                  + f_10 * id_95[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, hd_2, hd_4, hd_5, hd_20, hd_22, hd_23, hd_62, \
                         hd_64, hd_65, id_2, id_10, id_17, id_20, id_40, id_47, id_62, id_94, \
                         id_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_19 * ab_x[k] * hd_2[k]
                  - f_19 * ab_y[k] * hd_4[k]
                  + f_20 * ab_z[k] * hd_5[k]
                  + f_17 * ab_x[k] * hd_20[k]
                  + f_17 * ab_y[k] * hd_22[k]
                  - f_18 * ab_z[k] * hd_23[k]
                  - f_15 * ab_x[k] * hd_62[k]
                  - f_15 * ab_y[k] * hd_64[k]
                  + f_16 * ab_z[k] * hd_65[k]
                  - f_19 * id_2[k]
                  - f_19 * id_10[k]
                  + f_20 * id_17[k]
                  + f_17 * id_20[k]
                  + f_17 * id_40[k]
                  - f_18 * id_47[k]
                  - f_15 * id_62[k]
                  - f_15 * id_94[k]
                  + f_16 * id_101[k];
    }

#pragma omp simd aligned(ab_x, hd_0, hd_3, hd_5, hd_18, hd_21, hd_23, hd_60, hd_63, hd_65, \
                         id_0, id_3, id_5, id_18, id_21, id_23, id_60, id_63, \
                         id_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_13 * ab_x[k] * hd_0[k]
                  - f_13 * ab_x[k] * hd_3[k]
                  + f_14 * ab_x[k] * hd_5[k]
                  + f_11 * ab_x[k] * hd_18[k]
                  + f_11 * ab_x[k] * hd_21[k]
                  - f_12 * ab_x[k] * hd_23[k]
                  - f_9 * ab_x[k] * hd_60[k]
                  - f_9 * ab_x[k] * hd_63[k]
                  + f_10 * ab_x[k] * hd_65[k]
                  - f_13 * id_0[k]
                  - f_13 * id_3[k]
                  + f_14 * id_5[k]
                  + f_11 * id_18[k]
                  + f_11 * id_21[k]
                  - f_12 * id_23[k]
                  - f_9 * id_60[k]
                  - f_9 * id_63[k]
                  + f_10 * id_65[k];
    }

#pragma omp simd aligned(ab_x, ab_y, hd_2, hd_4, hd_20, hd_22, hd_62, hd_64, id_2, id_10, \
                         id_20, id_40, id_62, id_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_22 * ab_x[k] * hd_2[k]
                  - f_22 * ab_y[k] * hd_4[k]
                  - f_6 * ab_x[k] * hd_20[k]
                  + f_6 * ab_y[k] * hd_22[k]
                  + f_21 * ab_x[k] * hd_62[k]
                  - f_21 * ab_y[k] * hd_64[k]
                  + f_22 * id_2[k]
                  - f_22 * id_10[k]
                  - f_6 * id_20[k]
                  + f_6 * id_40[k]
                  + f_21 * id_62[k]
                  - f_21 * id_94[k];
    }

#pragma omp simd aligned(ab_x, hd_0, hd_3, hd_18, hd_21, hd_60, hd_63, id_0, id_3, id_18, \
                         id_21, id_60, id_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_5 * ab_x[k] * hd_0[k]
                  - f_4 * ab_x[k] * hd_3[k]
                  - f_3 * ab_x[k] * hd_18[k]
                  + f_2 * ab_x[k] * hd_21[k]
                  + f_1 * ab_x[k] * hd_60[k]
                  - f_0 * ab_x[k] * hd_63[k]
                  + f_5 * id_0[k]
                  - f_4 * id_3[k]
                  - f_3 * id_18[k]
                  + f_2 * id_21[k]
                  + f_1 * id_60[k]
                  - f_0 * id_63[k];
    }
}

auto
compute_hrr_hf(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t hd, const size_t id, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_1 = buffer.data(hd + 1);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_4 = buffer.data(hd + 4);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_8 = buffer.data(hd + 8);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_10 = buffer.data(hd + 10);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_13 = buffer.data(hd + 13);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_15 = buffer.data(hd + 15);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_20 = buffer.data(hd + 20);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_22 = buffer.data(hd + 22);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_24 = buffer.data(hd + 24);
    const auto *hd_25 = buffer.data(hd + 25);
    const auto *hd_26 = buffer.data(hd + 26);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_31 = buffer.data(hd + 31);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_38 = buffer.data(hd + 38);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_40 = buffer.data(hd + 40);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_43 = buffer.data(hd + 43);
    const auto *hd_44 = buffer.data(hd + 44);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_49 = buffer.data(hd + 49);
    const auto *hd_50 = buffer.data(hd + 50);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_55 = buffer.data(hd + 55);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_62 = buffer.data(hd + 62);
    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_64 = buffer.data(hd + 64);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_67 = buffer.data(hd + 67);
    const auto *hd_68 = buffer.data(hd + 68);
    const auto *hd_69 = buffer.data(hd + 69);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_73 = buffer.data(hd + 73);
    const auto *hd_74 = buffer.data(hd + 74);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_79 = buffer.data(hd + 79);
    const auto *hd_80 = buffer.data(hd + 80);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_83 = buffer.data(hd + 83);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_85 = buffer.data(hd + 85);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_88 = buffer.data(hd + 88);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_92 = buffer.data(hd + 92);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_96 = buffer.data(hd + 96);
    const auto *hd_97 = buffer.data(hd + 97);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_114 = buffer.data(hd + 114);
    const auto *hd_115 = buffer.data(hd + 115);
    const auto *hd_116 = buffer.data(hd + 116);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_121 = buffer.data(hd + 121);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

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
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_167 = buffer.data(id + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, hd_0, hd_1, hd_2, hd_3, hd_4, id_0, \
                         id_1, id_2, id_3, id_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * hd_0[k]
                 + id_0[k];

        t_1[k] = ab_x[k] * hd_1[k]
                 + id_1[k];

        t_2[k] = ab_x[k] * hd_2[k]
                 + id_2[k];

        t_3[k] = ab_x[k] * hd_3[k]
                 + id_3[k];

        t_4[k] = ab_x[k] * hd_4[k]
                 + id_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, ab_y, ab_z, hd_3, hd_4, hd_5, id_5, \
                         id_9, id_10, id_11, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_x[k] * hd_5[k]
                 + id_5[k];

        t_6[k] = ab_y[k] * hd_3[k]
                 + id_9[k];

        t_7[k] = ab_y[k] * hd_4[k]
                 + id_10[k];

        t_8[k] = ab_y[k] * hd_5[k]
                 + id_11[k];

        t_9[k] = ab_z[k] * hd_5[k]
                 + id_17[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, hd_6, hd_7, hd_8, hd_9, hd_10, \
                         id_6, id_7, id_8, id_9, id_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = ab_x[k] * hd_6[k]
                  + id_6[k];

        t_11[k] = ab_x[k] * hd_7[k]
                  + id_7[k];

        t_12[k] = ab_x[k] * hd_8[k]
                  + id_8[k];

        t_13[k] = ab_x[k] * hd_9[k]
                  + id_9[k];

        t_14[k] = ab_x[k] * hd_10[k]
                  + id_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, ab_y, ab_z, hd_9, hd_10, hd_11, \
                         id_11, id_21, id_22, id_23, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = ab_x[k] * hd_11[k]
                  + id_11[k];

        t_16[k] = ab_y[k] * hd_9[k]
                  + id_21[k];

        t_17[k] = ab_y[k] * hd_10[k]
                  + id_22[k];

        t_18[k] = ab_y[k] * hd_11[k]
                  + id_23[k];

        t_19[k] = ab_z[k] * hd_11[k]
                  + id_29[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, hd_12, hd_13, hd_14, hd_15, \
                         hd_16, id_12, id_13, id_14, id_15, id_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = ab_x[k] * hd_12[k]
                  + id_12[k];

        t_21[k] = ab_x[k] * hd_13[k]
                  + id_13[k];

        t_22[k] = ab_x[k] * hd_14[k]
                  + id_14[k];

        t_23[k] = ab_x[k] * hd_15[k]
                  + id_15[k];

        t_24[k] = ab_x[k] * hd_16[k]
                  + id_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, ab_y, ab_z, hd_15, hd_16, hd_17, \
                         id_17, id_27, id_28, id_29, id_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = ab_x[k] * hd_17[k]
                  + id_17[k];

        t_26[k] = ab_y[k] * hd_15[k]
                  + id_27[k];

        t_27[k] = ab_y[k] * hd_16[k]
                  + id_28[k];

        t_28[k] = ab_y[k] * hd_17[k]
                  + id_29[k];

        t_29[k] = ab_z[k] * hd_17[k]
                  + id_35[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, hd_18, hd_19, hd_20, hd_21, \
                         hd_22, id_18, id_19, id_20, id_21, id_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = ab_x[k] * hd_18[k]
                  + id_18[k];

        t_31[k] = ab_x[k] * hd_19[k]
                  + id_19[k];

        t_32[k] = ab_x[k] * hd_20[k]
                  + id_20[k];

        t_33[k] = ab_x[k] * hd_21[k]
                  + id_21[k];

        t_34[k] = ab_x[k] * hd_22[k]
                  + id_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, ab_y, ab_z, hd_21, hd_22, hd_23, \
                         id_23, id_39, id_40, id_41, id_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = ab_x[k] * hd_23[k]
                  + id_23[k];

        t_36[k] = ab_y[k] * hd_21[k]
                  + id_39[k];

        t_37[k] = ab_y[k] * hd_22[k]
                  + id_40[k];

        t_38[k] = ab_y[k] * hd_23[k]
                  + id_41[k];

        t_39[k] = ab_z[k] * hd_23[k]
                  + id_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, hd_24, hd_25, hd_26, hd_27, \
                         hd_28, id_24, id_25, id_26, id_27, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = ab_x[k] * hd_24[k]
                  + id_24[k];

        t_41[k] = ab_x[k] * hd_25[k]
                  + id_25[k];

        t_42[k] = ab_x[k] * hd_26[k]
                  + id_26[k];

        t_43[k] = ab_x[k] * hd_27[k]
                  + id_27[k];

        t_44[k] = ab_x[k] * hd_28[k]
                  + id_28[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, ab_y, ab_z, hd_27, hd_28, hd_29, \
                         id_29, id_45, id_46, id_47, id_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_x[k] * hd_29[k]
                  + id_29[k];

        t_46[k] = ab_y[k] * hd_27[k]
                  + id_45[k];

        t_47[k] = ab_y[k] * hd_28[k]
                  + id_46[k];

        t_48[k] = ab_y[k] * hd_29[k]
                  + id_47[k];

        t_49[k] = ab_z[k] * hd_29[k]
                  + id_53[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, hd_30, hd_31, hd_32, hd_33, \
                         hd_34, id_30, id_31, id_32, id_33, id_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = ab_x[k] * hd_30[k]
                  + id_30[k];

        t_51[k] = ab_x[k] * hd_31[k]
                  + id_31[k];

        t_52[k] = ab_x[k] * hd_32[k]
                  + id_32[k];

        t_53[k] = ab_x[k] * hd_33[k]
                  + id_33[k];

        t_54[k] = ab_x[k] * hd_34[k]
                  + id_34[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, ab_y, ab_z, hd_33, hd_34, hd_35, \
                         id_35, id_51, id_52, id_53, id_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = ab_x[k] * hd_35[k]
                  + id_35[k];

        t_56[k] = ab_y[k] * hd_33[k]
                  + id_51[k];

        t_57[k] = ab_y[k] * hd_34[k]
                  + id_52[k];

        t_58[k] = ab_y[k] * hd_35[k]
                  + id_53[k];

        t_59[k] = ab_z[k] * hd_35[k]
                  + id_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, hd_36, hd_37, hd_38, hd_39, \
                         hd_40, id_36, id_37, id_38, id_39, id_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = ab_x[k] * hd_36[k]
                  + id_36[k];

        t_61[k] = ab_x[k] * hd_37[k]
                  + id_37[k];

        t_62[k] = ab_x[k] * hd_38[k]
                  + id_38[k];

        t_63[k] = ab_x[k] * hd_39[k]
                  + id_39[k];

        t_64[k] = ab_x[k] * hd_40[k]
                  + id_40[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, ab_y, ab_z, hd_39, hd_40, hd_41, \
                         id_41, id_63, id_64, id_65, id_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = ab_x[k] * hd_41[k]
                  + id_41[k];

        t_66[k] = ab_y[k] * hd_39[k]
                  + id_63[k];

        t_67[k] = ab_y[k] * hd_40[k]
                  + id_64[k];

        t_68[k] = ab_y[k] * hd_41[k]
                  + id_65[k];

        t_69[k] = ab_z[k] * hd_41[k]
                  + id_71[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, hd_42, hd_43, hd_44, hd_45, \
                         hd_46, id_42, id_43, id_44, id_45, id_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = ab_x[k] * hd_42[k]
                  + id_42[k];

        t_71[k] = ab_x[k] * hd_43[k]
                  + id_43[k];

        t_72[k] = ab_x[k] * hd_44[k]
                  + id_44[k];

        t_73[k] = ab_x[k] * hd_45[k]
                  + id_45[k];

        t_74[k] = ab_x[k] * hd_46[k]
                  + id_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, ab_y, ab_z, hd_45, hd_46, hd_47, \
                         id_47, id_69, id_70, id_71, id_77 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = ab_x[k] * hd_47[k]
                  + id_47[k];

        t_76[k] = ab_y[k] * hd_45[k]
                  + id_69[k];

        t_77[k] = ab_y[k] * hd_46[k]
                  + id_70[k];

        t_78[k] = ab_y[k] * hd_47[k]
                  + id_71[k];

        t_79[k] = ab_z[k] * hd_47[k]
                  + id_77[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, hd_48, hd_49, hd_50, hd_51, \
                         hd_52, id_48, id_49, id_50, id_51, id_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = ab_x[k] * hd_48[k]
                  + id_48[k];

        t_81[k] = ab_x[k] * hd_49[k]
                  + id_49[k];

        t_82[k] = ab_x[k] * hd_50[k]
                  + id_50[k];

        t_83[k] = ab_x[k] * hd_51[k]
                  + id_51[k];

        t_84[k] = ab_x[k] * hd_52[k]
                  + id_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, hd_51, hd_52, hd_53, \
                         id_53, id_75, id_76, id_77, id_83 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * hd_53[k]
                  + id_53[k];

        t_86[k] = ab_y[k] * hd_51[k]
                  + id_75[k];

        t_87[k] = ab_y[k] * hd_52[k]
                  + id_76[k];

        t_88[k] = ab_y[k] * hd_53[k]
                  + id_77[k];

        t_89[k] = ab_z[k] * hd_53[k]
                  + id_83[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, hd_54, hd_55, hd_56, hd_57, \
                         hd_58, id_54, id_55, id_56, id_57, id_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = ab_x[k] * hd_54[k]
                  + id_54[k];

        t_91[k] = ab_x[k] * hd_55[k]
                  + id_55[k];

        t_92[k] = ab_x[k] * hd_56[k]
                  + id_56[k];

        t_93[k] = ab_x[k] * hd_57[k]
                  + id_57[k];

        t_94[k] = ab_x[k] * hd_58[k]
                  + id_58[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, ab_y, ab_z, hd_57, hd_58, hd_59, \
                         id_59, id_81, id_82, id_83, id_89 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = ab_x[k] * hd_59[k]
                  + id_59[k];

        t_96[k] = ab_y[k] * hd_57[k]
                  + id_81[k];

        t_97[k] = ab_y[k] * hd_58[k]
                  + id_82[k];

        t_98[k] = ab_y[k] * hd_59[k]
                  + id_83[k];

        t_99[k] = ab_z[k] * hd_59[k]
                  + id_89[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, hd_60, hd_61, hd_62, hd_63, \
                         hd_64, id_60, id_61, id_62, id_63, id_64 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = ab_x[k] * hd_60[k]
                   + id_60[k];

        t_101[k] = ab_x[k] * hd_61[k]
                   + id_61[k];

        t_102[k] = ab_x[k] * hd_62[k]
                   + id_62[k];

        t_103[k] = ab_x[k] * hd_63[k]
                   + id_63[k];

        t_104[k] = ab_x[k] * hd_64[k]
                   + id_64[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, ab_y, ab_z, hd_63, hd_64, \
                         hd_65, id_65, id_93, id_94, id_95, id_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = ab_x[k] * hd_65[k]
                   + id_65[k];

        t_106[k] = ab_y[k] * hd_63[k]
                   + id_93[k];

        t_107[k] = ab_y[k] * hd_64[k]
                   + id_94[k];

        t_108[k] = ab_y[k] * hd_65[k]
                   + id_95[k];

        t_109[k] = ab_z[k] * hd_65[k]
                   + id_101[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, hd_66, hd_67, hd_68, hd_69, \
                         hd_70, id_66, id_67, id_68, id_69, id_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = ab_x[k] * hd_66[k]
                   + id_66[k];

        t_111[k] = ab_x[k] * hd_67[k]
                   + id_67[k];

        t_112[k] = ab_x[k] * hd_68[k]
                   + id_68[k];

        t_113[k] = ab_x[k] * hd_69[k]
                   + id_69[k];

        t_114[k] = ab_x[k] * hd_70[k]
                   + id_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, ab_y, ab_z, hd_69, hd_70, \
                         hd_71, id_71, id_99, id_100, id_101, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = ab_x[k] * hd_71[k]
                   + id_71[k];

        t_116[k] = ab_y[k] * hd_69[k]
                   + id_99[k];

        t_117[k] = ab_y[k] * hd_70[k]
                   + id_100[k];

        t_118[k] = ab_y[k] * hd_71[k]
                   + id_101[k];

        t_119[k] = ab_z[k] * hd_71[k]
                   + id_107[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, hd_72, hd_73, hd_74, hd_75, \
                         hd_76, id_72, id_73, id_74, id_75, id_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = ab_x[k] * hd_72[k]
                   + id_72[k];

        t_121[k] = ab_x[k] * hd_73[k]
                   + id_73[k];

        t_122[k] = ab_x[k] * hd_74[k]
                   + id_74[k];

        t_123[k] = ab_x[k] * hd_75[k]
                   + id_75[k];

        t_124[k] = ab_x[k] * hd_76[k]
                   + id_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, ab_y, ab_z, hd_75, hd_76, \
                         hd_77, id_77, id_105, id_106, id_107, id_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = ab_x[k] * hd_77[k]
                   + id_77[k];

        t_126[k] = ab_y[k] * hd_75[k]
                   + id_105[k];

        t_127[k] = ab_y[k] * hd_76[k]
                   + id_106[k];

        t_128[k] = ab_y[k] * hd_77[k]
                   + id_107[k];

        t_129[k] = ab_z[k] * hd_77[k]
                   + id_113[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, hd_78, hd_79, hd_80, hd_81, \
                         hd_82, id_78, id_79, id_80, id_81, id_82 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = ab_x[k] * hd_78[k]
                   + id_78[k];

        t_131[k] = ab_x[k] * hd_79[k]
                   + id_79[k];

        t_132[k] = ab_x[k] * hd_80[k]
                   + id_80[k];

        t_133[k] = ab_x[k] * hd_81[k]
                   + id_81[k];

        t_134[k] = ab_x[k] * hd_82[k]
                   + id_82[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, ab_y, ab_z, hd_81, hd_82, \
                         hd_83, id_83, id_111, id_112, id_113, id_119 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = ab_x[k] * hd_83[k]
                   + id_83[k];

        t_136[k] = ab_y[k] * hd_81[k]
                   + id_111[k];

        t_137[k] = ab_y[k] * hd_82[k]
                   + id_112[k];

        t_138[k] = ab_y[k] * hd_83[k]
                   + id_113[k];

        t_139[k] = ab_z[k] * hd_83[k]
                   + id_119[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, hd_84, hd_85, hd_86, hd_87, \
                         hd_88, id_84, id_85, id_86, id_87, id_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = ab_x[k] * hd_84[k]
                   + id_84[k];

        t_141[k] = ab_x[k] * hd_85[k]
                   + id_85[k];

        t_142[k] = ab_x[k] * hd_86[k]
                   + id_86[k];

        t_143[k] = ab_x[k] * hd_87[k]
                   + id_87[k];

        t_144[k] = ab_x[k] * hd_88[k]
                   + id_88[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, ab_y, ab_z, hd_87, hd_88, \
                         hd_89, id_89, id_117, id_118, id_119, id_125 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = ab_x[k] * hd_89[k]
                   + id_89[k];

        t_146[k] = ab_y[k] * hd_87[k]
                   + id_117[k];

        t_147[k] = ab_y[k] * hd_88[k]
                   + id_118[k];

        t_148[k] = ab_y[k] * hd_89[k]
                   + id_119[k];

        t_149[k] = ab_z[k] * hd_89[k]
                   + id_125[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, hd_90, hd_91, hd_92, hd_93, \
                         hd_94, id_90, id_91, id_92, id_93, id_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = ab_x[k] * hd_90[k]
                   + id_90[k];

        t_151[k] = ab_x[k] * hd_91[k]
                   + id_91[k];

        t_152[k] = ab_x[k] * hd_92[k]
                   + id_92[k];

        t_153[k] = ab_x[k] * hd_93[k]
                   + id_93[k];

        t_154[k] = ab_x[k] * hd_94[k]
                   + id_94[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, ab_y, ab_z, hd_93, hd_94, \
                         hd_95, id_95, id_129, id_130, id_131, id_137 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = ab_x[k] * hd_95[k]
                   + id_95[k];

        t_156[k] = ab_y[k] * hd_93[k]
                   + id_129[k];

        t_157[k] = ab_y[k] * hd_94[k]
                   + id_130[k];

        t_158[k] = ab_y[k] * hd_95[k]
                   + id_131[k];

        t_159[k] = ab_z[k] * hd_95[k]
                   + id_137[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, hd_96, hd_97, hd_98, hd_99, \
                         hd_100, id_96, id_97, id_98, id_99, id_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = ab_x[k] * hd_96[k]
                   + id_96[k];

        t_161[k] = ab_x[k] * hd_97[k]
                   + id_97[k];

        t_162[k] = ab_x[k] * hd_98[k]
                   + id_98[k];

        t_163[k] = ab_x[k] * hd_99[k]
                   + id_99[k];

        t_164[k] = ab_x[k] * hd_100[k]
                   + id_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, ab_y, ab_z, hd_99, hd_100, \
                         hd_101, id_101, id_135, id_136, id_137, \
                         id_143 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = ab_x[k] * hd_101[k]
                   + id_101[k];

        t_166[k] = ab_y[k] * hd_99[k]
                   + id_135[k];

        t_167[k] = ab_y[k] * hd_100[k]
                   + id_136[k];

        t_168[k] = ab_y[k] * hd_101[k]
                   + id_137[k];

        t_169[k] = ab_z[k] * hd_101[k]
                   + id_143[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, hd_102, hd_103, hd_104, \
                         hd_105, hd_106, id_102, id_103, id_104, id_105, \
                         id_106 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = ab_x[k] * hd_102[k]
                   + id_102[k];

        t_171[k] = ab_x[k] * hd_103[k]
                   + id_103[k];

        t_172[k] = ab_x[k] * hd_104[k]
                   + id_104[k];

        t_173[k] = ab_x[k] * hd_105[k]
                   + id_105[k];

        t_174[k] = ab_x[k] * hd_106[k]
                   + id_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, ab_y, ab_z, hd_105, hd_106, \
                         hd_107, id_107, id_141, id_142, id_143, \
                         id_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = ab_x[k] * hd_107[k]
                   + id_107[k];

        t_176[k] = ab_y[k] * hd_105[k]
                   + id_141[k];

        t_177[k] = ab_y[k] * hd_106[k]
                   + id_142[k];

        t_178[k] = ab_y[k] * hd_107[k]
                   + id_143[k];

        t_179[k] = ab_z[k] * hd_107[k]
                   + id_149[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, hd_108, hd_109, hd_110, \
                         hd_111, hd_112, id_108, id_109, id_110, id_111, \
                         id_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = ab_x[k] * hd_108[k]
                   + id_108[k];

        t_181[k] = ab_x[k] * hd_109[k]
                   + id_109[k];

        t_182[k] = ab_x[k] * hd_110[k]
                   + id_110[k];

        t_183[k] = ab_x[k] * hd_111[k]
                   + id_111[k];

        t_184[k] = ab_x[k] * hd_112[k]
                   + id_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, ab_y, ab_z, hd_111, hd_112, \
                         hd_113, id_113, id_147, id_148, id_149, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = ab_x[k] * hd_113[k]
                   + id_113[k];

        t_186[k] = ab_y[k] * hd_111[k]
                   + id_147[k];

        t_187[k] = ab_y[k] * hd_112[k]
                   + id_148[k];

        t_188[k] = ab_y[k] * hd_113[k]
                   + id_149[k];

        t_189[k] = ab_z[k] * hd_113[k]
                   + id_155[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, hd_114, hd_115, hd_116, \
                         hd_117, hd_118, id_114, id_115, id_116, id_117, \
                         id_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = ab_x[k] * hd_114[k]
                   + id_114[k];

        t_191[k] = ab_x[k] * hd_115[k]
                   + id_115[k];

        t_192[k] = ab_x[k] * hd_116[k]
                   + id_116[k];

        t_193[k] = ab_x[k] * hd_117[k]
                   + id_117[k];

        t_194[k] = ab_x[k] * hd_118[k]
                   + id_118[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, ab_y, ab_z, hd_117, hd_118, \
                         hd_119, id_119, id_153, id_154, id_155, \
                         id_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = ab_x[k] * hd_119[k]
                   + id_119[k];

        t_196[k] = ab_y[k] * hd_117[k]
                   + id_153[k];

        t_197[k] = ab_y[k] * hd_118[k]
                   + id_154[k];

        t_198[k] = ab_y[k] * hd_119[k]
                   + id_155[k];

        t_199[k] = ab_z[k] * hd_119[k]
                   + id_161[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, hd_120, hd_121, hd_122, \
                         hd_123, hd_124, id_120, id_121, id_122, id_123, \
                         id_124 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = ab_x[k] * hd_120[k]
                   + id_120[k];

        t_201[k] = ab_x[k] * hd_121[k]
                   + id_121[k];

        t_202[k] = ab_x[k] * hd_122[k]
                   + id_122[k];

        t_203[k] = ab_x[k] * hd_123[k]
                   + id_123[k];

        t_204[k] = ab_x[k] * hd_124[k]
                   + id_124[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, ab_y, ab_z, hd_123, hd_124, \
                         hd_125, id_125, id_159, id_160, id_161, \
                         id_167 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = ab_x[k] * hd_125[k]
                   + id_125[k];

        t_206[k] = ab_y[k] * hd_123[k]
                   + id_159[k];

        t_207[k] = ab_y[k] * hd_124[k]
                   + id_160[k];

        t_208[k] = ab_y[k] * hd_125[k]
                   + id_161[k];

        t_209[k] = ab_z[k] * hd_125[k]
                   + id_167[k];
    }
}

}  // namespace simdtrf
