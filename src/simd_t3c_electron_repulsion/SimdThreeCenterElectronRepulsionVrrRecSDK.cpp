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


#include "SimdThreeCenterElectronRepulsionVrrRecSDK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sdk_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t spk0,
                                                          const size_t spi, const size_t spk1,
                                                          const size_t sdh0, const size_t sdh1,
                                                          const size_t sdi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spk0_0 = buffer.data(spk0 + 0);
    const auto *spk0_3 = buffer.data(spk0 + 3);
    const auto *spk0_5 = buffer.data(spk0 + 5);
    const auto *spk0_6 = buffer.data(spk0 + 6);
    const auto *spk0_9 = buffer.data(spk0 + 9);
    const auto *spk0_10 = buffer.data(spk0 + 10);
    const auto *spk0_14 = buffer.data(spk0 + 14);
    const auto *spk0_15 = buffer.data(spk0 + 15);
    const auto *spk0_20 = buffer.data(spk0 + 20);
    const auto *spk0_39 = buffer.data(spk0 + 39);
    const auto *spk0_42 = buffer.data(spk0 + 42);
    const auto *spk0_46 = buffer.data(spk0 + 46);
    const auto *spk0_48 = buffer.data(spk0 + 48);
    const auto *spk0_51 = buffer.data(spk0 + 51);
    const auto *spk0_53 = buffer.data(spk0 + 53);
    const auto *spk0_54 = buffer.data(spk0 + 54);
    const auto *spk0_64 = buffer.data(spk0 + 64);
    const auto *spk0_66 = buffer.data(spk0 + 66);
    const auto *spk0_67 = buffer.data(spk0 + 67);
    const auto *spk0_68 = buffer.data(spk0 + 68);
    const auto *spk0_69 = buffer.data(spk0 + 69);
    const auto *spk0_71 = buffer.data(spk0 + 71);
    const auto *spk0_77 = buffer.data(spk0 + 77);
    const auto *spk0_81 = buffer.data(spk0 + 81);
    const auto *spk0_84 = buffer.data(spk0 + 84);
    const auto *spk0_86 = buffer.data(spk0 + 86);
    const auto *spk0_89 = buffer.data(spk0 + 89);
    const auto *spk0_90 = buffer.data(spk0 + 90);
    const auto *spk0_92 = buffer.data(spk0 + 92);
    const auto *spk0_100 = buffer.data(spk0 + 100);
    const auto *spk0_102 = buffer.data(spk0 + 102);
    const auto *spk0_103 = buffer.data(spk0 + 103);
    const auto *spk0_104 = buffer.data(spk0 + 104);
    const auto *spk0_105 = buffer.data(spk0 + 105);
    const auto *spk0_107 = buffer.data(spk0 + 107);

    const auto *spi_0 = buffer.data(spi + 0);
    const auto *spi_2 = buffer.data(spi + 2);
    const auto *spi_3 = buffer.data(spi + 3);
    const auto *spi_5 = buffer.data(spi + 5);
    const auto *spi_6 = buffer.data(spi + 6);
    const auto *spi_9 = buffer.data(spi + 9);
    const auto *spi_10 = buffer.data(spi + 10);
    const auto *spi_12 = buffer.data(spi + 12);
    const auto *spi_14 = buffer.data(spi + 14);
    const auto *spi_15 = buffer.data(spi + 15);
    const auto *spi_17 = buffer.data(spi + 17);
    const auto *spi_18 = buffer.data(spi + 18);
    const auto *spi_20 = buffer.data(spi + 20);
    const auto *spi_21 = buffer.data(spi + 21);
    const auto *spi_22 = buffer.data(spi + 22);
    const auto *spi_23 = buffer.data(spi + 23);
    const auto *spi_24 = buffer.data(spi + 24);
    const auto *spi_25 = buffer.data(spi + 25);
    const auto *spi_26 = buffer.data(spi + 26);
    const auto *spi_27 = buffer.data(spi + 27);
    const auto *spi_28 = buffer.data(spi + 28);
    const auto *spi_30 = buffer.data(spi + 30);
    const auto *spi_31 = buffer.data(spi + 31);
    const auto *spi_33 = buffer.data(spi + 33);
    const auto *spi_34 = buffer.data(spi + 34);
    const auto *spi_37 = buffer.data(spi + 37);
    const auto *spi_38 = buffer.data(spi + 38);
    const auto *spi_40 = buffer.data(spi + 40);
    const auto *spi_43 = buffer.data(spi + 43);
    const auto *spi_45 = buffer.data(spi + 45);
    const auto *spi_46 = buffer.data(spi + 46);
    const auto *spi_49 = buffer.data(spi + 49);
    const auto *spi_50 = buffer.data(spi + 50);
    const auto *spi_51 = buffer.data(spi + 51);
    const auto *spi_52 = buffer.data(spi + 52);
    const auto *spi_53 = buffer.data(spi + 53);
    const auto *spi_54 = buffer.data(spi + 54);
    const auto *spi_55 = buffer.data(spi + 55);
    const auto *spi_61 = buffer.data(spi + 61);
    const auto *spi_65 = buffer.data(spi + 65);
    const auto *spi_68 = buffer.data(spi + 68);
    const auto *spi_70 = buffer.data(spi + 70);
    const auto *spi_73 = buffer.data(spi + 73);
    const auto *spi_74 = buffer.data(spi + 74);
    const auto *spi_76 = buffer.data(spi + 76);
    const auto *spi_77 = buffer.data(spi + 77);
    const auto *spi_78 = buffer.data(spi + 78);
    const auto *spi_79 = buffer.data(spi + 79);
    const auto *spi_80 = buffer.data(spi + 80);
    const auto *spi_81 = buffer.data(spi + 81);
    const auto *spi_82 = buffer.data(spi + 82);
    const auto *spi_83 = buffer.data(spi + 83);

    const auto *spk1_0 = buffer.data(spk1 + 0);
    const auto *spk1_3 = buffer.data(spk1 + 3);
    const auto *spk1_5 = buffer.data(spk1 + 5);
    const auto *spk1_6 = buffer.data(spk1 + 6);
    const auto *spk1_9 = buffer.data(spk1 + 9);
    const auto *spk1_10 = buffer.data(spk1 + 10);
    const auto *spk1_14 = buffer.data(spk1 + 14);
    const auto *spk1_15 = buffer.data(spk1 + 15);
    const auto *spk1_20 = buffer.data(spk1 + 20);
    const auto *spk1_39 = buffer.data(spk1 + 39);
    const auto *spk1_42 = buffer.data(spk1 + 42);
    const auto *spk1_46 = buffer.data(spk1 + 46);
    const auto *spk1_48 = buffer.data(spk1 + 48);
    const auto *spk1_51 = buffer.data(spk1 + 51);
    const auto *spk1_53 = buffer.data(spk1 + 53);
    const auto *spk1_54 = buffer.data(spk1 + 54);
    const auto *spk1_64 = buffer.data(spk1 + 64);
    const auto *spk1_66 = buffer.data(spk1 + 66);
    const auto *spk1_67 = buffer.data(spk1 + 67);
    const auto *spk1_68 = buffer.data(spk1 + 68);
    const auto *spk1_69 = buffer.data(spk1 + 69);
    const auto *spk1_71 = buffer.data(spk1 + 71);
    const auto *spk1_77 = buffer.data(spk1 + 77);
    const auto *spk1_81 = buffer.data(spk1 + 81);
    const auto *spk1_84 = buffer.data(spk1 + 84);
    const auto *spk1_86 = buffer.data(spk1 + 86);
    const auto *spk1_89 = buffer.data(spk1 + 89);
    const auto *spk1_90 = buffer.data(spk1 + 90);
    const auto *spk1_92 = buffer.data(spk1 + 92);
    const auto *spk1_100 = buffer.data(spk1 + 100);
    const auto *spk1_102 = buffer.data(spk1 + 102);
    const auto *spk1_103 = buffer.data(spk1 + 103);
    const auto *spk1_104 = buffer.data(spk1 + 104);
    const auto *spk1_105 = buffer.data(spk1 + 105);
    const auto *spk1_107 = buffer.data(spk1 + 107);

    const auto *sdh0_0 = buffer.data(sdh0 + 0);
    const auto *sdh0_3 = buffer.data(sdh0 + 3);
    const auto *sdh0_5 = buffer.data(sdh0 + 5);
    const auto *sdh0_6 = buffer.data(sdh0 + 6);
    const auto *sdh0_9 = buffer.data(sdh0 + 9);
    const auto *sdh0_10 = buffer.data(sdh0 + 10);
    const auto *sdh0_12 = buffer.data(sdh0 + 12);
    const auto *sdh0_14 = buffer.data(sdh0 + 14);
    const auto *sdh0_15 = buffer.data(sdh0 + 15);
    const auto *sdh0_17 = buffer.data(sdh0 + 17);
    const auto *sdh0_18 = buffer.data(sdh0 + 18);
    const auto *sdh0_19 = buffer.data(sdh0 + 19);
    const auto *sdh0_20 = buffer.data(sdh0 + 20);
    const auto *sdh0_63 = buffer.data(sdh0 + 63);
    const auto *sdh0_66 = buffer.data(sdh0 + 66);
    const auto *sdh0_68 = buffer.data(sdh0 + 68);
    const auto *sdh0_69 = buffer.data(sdh0 + 69);
    const auto *sdh0_72 = buffer.data(sdh0 + 72);
    const auto *sdh0_73 = buffer.data(sdh0 + 73);
    const auto *sdh0_75 = buffer.data(sdh0 + 75);
    const auto *sdh0_77 = buffer.data(sdh0 + 77);

    const auto *sdh1_0 = buffer.data(sdh1 + 0);
    const auto *sdh1_3 = buffer.data(sdh1 + 3);
    const auto *sdh1_5 = buffer.data(sdh1 + 5);
    const auto *sdh1_6 = buffer.data(sdh1 + 6);
    const auto *sdh1_9 = buffer.data(sdh1 + 9);
    const auto *sdh1_10 = buffer.data(sdh1 + 10);
    const auto *sdh1_12 = buffer.data(sdh1 + 12);
    const auto *sdh1_14 = buffer.data(sdh1 + 14);
    const auto *sdh1_15 = buffer.data(sdh1 + 15);
    const auto *sdh1_17 = buffer.data(sdh1 + 17);
    const auto *sdh1_18 = buffer.data(sdh1 + 18);
    const auto *sdh1_19 = buffer.data(sdh1 + 19);
    const auto *sdh1_20 = buffer.data(sdh1 + 20);
    const auto *sdh1_63 = buffer.data(sdh1 + 63);
    const auto *sdh1_66 = buffer.data(sdh1 + 66);
    const auto *sdh1_68 = buffer.data(sdh1 + 68);
    const auto *sdh1_69 = buffer.data(sdh1 + 69);
    const auto *sdh1_72 = buffer.data(sdh1 + 72);
    const auto *sdh1_73 = buffer.data(sdh1 + 73);
    const auto *sdh1_75 = buffer.data(sdh1 + 75);
    const auto *sdh1_77 = buffer.data(sdh1 + 77);

    const auto *sdi_0 = buffer.data(sdi + 0);
    const auto *sdi_2 = buffer.data(sdi + 2);
    const auto *sdi_3 = buffer.data(sdi + 3);
    const auto *sdi_5 = buffer.data(sdi + 5);
    const auto *sdi_6 = buffer.data(sdi + 6);
    const auto *sdi_9 = buffer.data(sdi + 9);
    const auto *sdi_10 = buffer.data(sdi + 10);
    const auto *sdi_12 = buffer.data(sdi + 12);
    const auto *sdi_14 = buffer.data(sdi + 14);
    const auto *sdi_15 = buffer.data(sdi + 15);
    const auto *sdi_17 = buffer.data(sdi + 17);
    const auto *sdi_18 = buffer.data(sdi + 18);
    const auto *sdi_20 = buffer.data(sdi + 20);
    const auto *sdi_21 = buffer.data(sdi + 21);
    const auto *sdi_22 = buffer.data(sdi + 22);
    const auto *sdi_23 = buffer.data(sdi + 23);
    const auto *sdi_24 = buffer.data(sdi + 24);
    const auto *sdi_25 = buffer.data(sdi + 25);
    const auto *sdi_26 = buffer.data(sdi + 26);
    const auto *sdi_27 = buffer.data(sdi + 27);
    const auto *sdi_28 = buffer.data(sdi + 28);
    const auto *sdi_30 = buffer.data(sdi + 30);
    const auto *sdi_31 = buffer.data(sdi + 31);
    const auto *sdi_33 = buffer.data(sdi + 33);
    const auto *sdi_34 = buffer.data(sdi + 34);
    const auto *sdi_37 = buffer.data(sdi + 37);
    const auto *sdi_38 = buffer.data(sdi + 38);
    const auto *sdi_42 = buffer.data(sdi + 42);
    const auto *sdi_49 = buffer.data(sdi + 49);
    const auto *sdi_50 = buffer.data(sdi + 50);
    const auto *sdi_51 = buffer.data(sdi + 51);
    const auto *sdi_52 = buffer.data(sdi + 52);
    const auto *sdi_53 = buffer.data(sdi + 53);
    const auto *sdi_54 = buffer.data(sdi + 54);
    const auto *sdi_55 = buffer.data(sdi + 55);
    const auto *sdi_56 = buffer.data(sdi + 56);
    const auto *sdi_58 = buffer.data(sdi + 58);
    const auto *sdi_59 = buffer.data(sdi + 59);
    const auto *sdi_61 = buffer.data(sdi + 61);
    const auto *sdi_62 = buffer.data(sdi + 62);
    const auto *sdi_65 = buffer.data(sdi + 65);
    const auto *sdi_66 = buffer.data(sdi + 66);
    const auto *sdi_70 = buffer.data(sdi + 70);
    const auto *sdi_77 = buffer.data(sdi + 77);
    const auto *sdi_78 = buffer.data(sdi + 78);
    const auto *sdi_79 = buffer.data(sdi + 79);
    const auto *sdi_80 = buffer.data(sdi + 80);
    const auto *sdi_81 = buffer.data(sdi + 81);
    const auto *sdi_82 = buffer.data(sdi + 82);
    const auto *sdi_83 = buffer.data(sdi + 83);
    const auto *sdi_84 = buffer.data(sdi + 84);
    const auto *sdi_86 = buffer.data(sdi + 86);
    const auto *sdi_87 = buffer.data(sdi + 87);
    const auto *sdi_89 = buffer.data(sdi + 89);
    const auto *sdi_90 = buffer.data(sdi + 90);
    const auto *sdi_93 = buffer.data(sdi + 93);
    const auto *sdi_94 = buffer.data(sdi + 94);
    const auto *sdi_96 = buffer.data(sdi + 96);
    const auto *sdi_98 = buffer.data(sdi + 98);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, spi_0, spi_3, sdh0_0, sdh0_3, \
                         sdh1_0, sdh1_3, sdi_0, sdi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spi_0[k]
                 + f_1 * sdh0_0[k]
                 - f_2 * sdh1_0[k]
                 + f_3 * pc_x[k] * sdi_0[k];

        t_1[k] = f_3 * pc_y[k] * sdi_0[k];

        t_2[k] = f_3 * pc_z[k] * sdi_0[k];

        t_3[k] = f_0 * spi_3[k]
                 + f_4 * sdh0_3[k]
                 - f_5 * sdh1_3[k]
                 + f_3 * pc_x[k] * sdi_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, spi_5, spi_6, sdh0_5, sdh0_6, sdh1_5, \
                         sdh1_6, sdi_2, sdi_5, sdi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sdi_2[k];

        t_5[k] = f_0 * spi_5[k]
                 + f_4 * sdh0_5[k]
                 - f_5 * sdh1_5[k]
                 + f_3 * pc_x[k] * sdi_5[k];

        t_6[k] = f_0 * spi_6[k]
                 + f_6 * sdh0_6[k]
                 - f_7 * sdh1_6[k]
                 + f_3 * pc_x[k] * sdi_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, spi_9, sdh0_9, sdh1_9, sdi_3, sdi_5, \
                         sdi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sdi_3[k];

        t_8[k] = f_3 * pc_y[k] * sdi_5[k];

        t_9[k] = f_0 * spi_9[k]
                 + f_6 * sdh0_9[k]
                 - f_7 * sdh1_9[k]
                 + f_3 * pc_x[k] * sdi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, spi_10, spi_12, sdh0_10, sdh0_12, \
                         sdh1_10, sdh1_12, sdi_6, sdi_10, sdi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * spi_10[k]
                  + f_8 * sdh0_10[k]
                  - f_9 * sdh1_10[k]
                  + f_3 * pc_x[k] * sdi_10[k];

        t_11[k] = f_3 * pc_z[k] * sdi_6[k];

        t_12[k] = f_0 * spi_12[k]
                  + f_8 * sdh0_12[k]
                  - f_9 * sdh1_12[k]
                  + f_3 * pc_x[k] * sdi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, spi_14, spi_15, sdh0_14, sdh0_15, \
                         sdh1_14, sdh1_15, sdi_9, sdi_14, sdi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sdi_9[k];

        t_14[k] = f_0 * spi_14[k]
                  + f_8 * sdh0_14[k]
                  - f_9 * sdh1_14[k]
                  + f_3 * pc_x[k] * sdi_14[k];

        t_15[k] = f_0 * spi_15[k]
                  + f_10 * sdh0_15[k]
                  - f_11 * sdh1_15[k]
                  + f_3 * pc_x[k] * sdi_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, spi_17, spi_18, sdh0_17, sdh0_18, \
                         sdh1_17, sdh1_18, sdi_10, sdi_17, sdi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sdi_10[k];

        t_17[k] = f_0 * spi_17[k]
                  + f_10 * sdh0_17[k]
                  - f_11 * sdh1_17[k]
                  + f_3 * pc_x[k] * sdi_17[k];

        t_18[k] = f_0 * spi_18[k]
                  + f_10 * sdh0_18[k]
                  - f_11 * sdh1_18[k]
                  + f_3 * pc_x[k] * sdi_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pc_x, pc_y, spi_20, spi_21, spi_22, sdh0_20, \
                         sdh1_20, sdi_14, sdi_20, sdi_21, sdi_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sdi_14[k];

        t_20[k] = f_0 * spi_20[k]
                  + f_10 * sdh0_20[k]
                  - f_11 * sdh1_20[k]
                  + f_3 * pc_x[k] * sdi_20[k];

        t_21[k] = f_0 * spi_21[k]
                  + f_3 * pc_x[k] * sdi_21[k];

        t_22[k] = f_0 * spi_22[k]
                  + f_3 * pc_x[k] * sdi_22[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, pc_x, spi_23, spi_24, spi_25, spi_26, \
                         spi_27, sdi_23, sdi_24, sdi_25, sdi_26, \
                         sdi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * spi_23[k]
                  + f_3 * pc_x[k] * sdi_23[k];

        t_24[k] = f_0 * spi_24[k]
                  + f_3 * pc_x[k] * sdi_24[k];

        t_25[k] = f_0 * spi_25[k]
                  + f_3 * pc_x[k] * sdi_25[k];

        t_26[k] = f_0 * spi_26[k]
                  + f_3 * pc_x[k] * sdi_26[k];

        t_27[k] = f_0 * spi_27[k]
                  + f_3 * pc_x[k] * sdi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, sdh0_15, sdh0_17, sdh0_18, \
                         sdh1_15, sdh1_17, sdh1_18, sdi_21, sdi_23, \
                         sdi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * sdh0_15[k]
                  - f_2 * sdh1_15[k]
                  + f_3 * pc_y[k] * sdi_21[k];

        t_29[k] = f_3 * pc_z[k] * sdi_21[k];

        t_30[k] = f_4 * sdh0_17[k]
                  - f_5 * sdh1_17[k]
                  + f_3 * pc_y[k] * sdi_23[k];

        t_31[k] = f_6 * sdh0_18[k]
                  - f_7 * sdh1_18[k]
                  + f_3 * pc_y[k] * sdi_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, sdh0_19, sdh0_20, sdh1_19, \
                         sdh1_20, sdi_25, sdi_26, sdi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_8 * sdh0_19[k]
                  - f_9 * sdh1_19[k]
                  + f_3 * pc_y[k] * sdi_25[k];

        t_33[k] = f_10 * sdh0_20[k]
                  - f_11 * sdh1_20[k]
                  + f_3 * pc_y[k] * sdi_26[k];

        t_34[k] = f_3 * pc_y[k] * sdi_27[k];

        t_35[k] = f_1 * sdh0_20[k]
                  - f_2 * sdh1_20[k]
                  + f_3 * pc_z[k] * sdi_27[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_x, pb_y, pc_x, pc_y, pc_z, spk0_0, \
                         spk0_39, spi_0, spi_31, spk1_0, spk1_39, \
                         sdi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * spk0_0[k]
                  - f_12 * pc_y[k] * spk1_0[k];

        t_37[k] = f_13 * spi_0[k]
                  + f_3 * pc_y[k] * sdi_28[k];

        t_38[k] = f_3 * pc_z[k] * sdi_28[k];

        t_39[k] = pb_x[k] * spk0_39[k]
                  + f_14 * spi_31[k]
                  - f_12 * pc_x[k] * spk1_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pb_x, pb_y, pc_x, pc_y, spk0_5, spk0_42, spi_2, \
                         spi_34, spk1_5, spk1_42, sdi_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_13 * spi_2[k]
                  + f_3 * pc_y[k] * sdi_30[k];

        t_41[k] = pb_y[k] * spk0_5[k]
                  - f_12 * pc_y[k] * spk1_5[k];

        t_42[k] = pb_x[k] * spk0_42[k]
                  + f_15 * spi_34[k]
                  - f_12 * pc_x[k] * spk1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_y, pc_y, pc_z, spk0_9, spi_5, spk1_9, sdi_31, \
                         sdi_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * pc_z[k] * sdi_31[k];

        t_44[k] = f_13 * spi_5[k]
                  + f_3 * pc_y[k] * sdi_33[k];

        t_45[k] = pb_y[k] * spk0_9[k]
                  - f_12 * pc_y[k] * spk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pc_x, pc_z, spk0_46, spk0_48, spi_38, spi_40, \
                         spk1_46, spk1_48, sdi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_x[k] * spk0_46[k]
                  + f_16 * spi_38[k]
                  - f_12 * pc_x[k] * spk1_46[k];

        t_47[k] = f_3 * pc_z[k] * sdi_34[k];

        t_48[k] = pb_x[k] * spk0_48[k]
                  + f_16 * spi_40[k]
                  - f_12 * pc_x[k] * spk1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_y, pc_x, pc_y, spk0_14, spk0_51, spi_9, \
                         spi_43, spk1_14, spk1_51, sdi_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_13 * spi_9[k]
                  + f_3 * pc_y[k] * sdi_37[k];

        t_50[k] = pb_y[k] * spk0_14[k]
                  - f_12 * pc_y[k] * spk1_14[k];

        t_51[k] = pb_x[k] * spk0_51[k]
                  + f_0 * spi_43[k]
                  - f_12 * pc_x[k] * spk1_51[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pc_x, pc_z, spk0_53, spk0_54, spi_45, spi_46, \
                         spk1_53, spk1_54, sdi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * pc_z[k] * sdi_38[k];

        t_53[k] = pb_x[k] * spk0_53[k]
                  + f_0 * spi_45[k]
                  - f_12 * pc_x[k] * spk1_53[k];

        t_54[k] = pb_x[k] * spk0_54[k]
                  + f_0 * spi_46[k]
                  - f_12 * pc_x[k] * spk1_54[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_y, pc_x, pc_y, spk0_20, spi_14, spi_49, \
                         spi_50, spk1_20, sdi_42, sdi_49, sdi_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_13 * spi_14[k]
                  + f_3 * pc_y[k] * sdi_42[k];

        t_56[k] = pb_y[k] * spk0_20[k]
                  - f_12 * pc_y[k] * spk1_20[k];

        t_57[k] = f_13 * spi_49[k]
                  + f_3 * pc_x[k] * sdi_49[k];

        t_58[k] = f_13 * spi_50[k]
                  + f_3 * pc_x[k] * sdi_50[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pc_x, spi_51, spi_52, spi_53, spi_54, \
                         spi_55, sdi_51, sdi_52, sdi_53, sdi_54, \
                         sdi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_13 * spi_51[k]
                  + f_3 * pc_x[k] * sdi_51[k];

        t_60[k] = f_13 * spi_52[k]
                  + f_3 * pc_x[k] * sdi_52[k];

        t_61[k] = f_13 * spi_53[k]
                  + f_3 * pc_x[k] * sdi_53[k];

        t_62[k] = f_13 * spi_54[k]
                  + f_3 * pc_x[k] * sdi_54[k];

        t_63[k] = f_13 * spi_55[k]
                  + f_3 * pc_x[k] * sdi_55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_x, pc_x, pc_z, spk0_64, spk0_66, spk0_67, \
                         spk1_64, spk1_66, spk1_67, sdi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_x[k] * spk0_64[k]
                  - f_12 * pc_x[k] * spk1_64[k];

        t_65[k] = f_3 * pc_z[k] * sdi_49[k];

        t_66[k] = pb_x[k] * spk0_66[k]
                  - f_12 * pc_x[k] * spk1_66[k];

        t_67[k] = pb_x[k] * spk0_67[k]
                  - f_12 * pc_x[k] * spk1_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_x, pc_x, pc_y, spk0_68, spk0_69, spk0_71, \
                         spi_27, spk1_68, spk1_69, spk1_71, sdi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * spk0_68[k]
                  - f_12 * pc_x[k] * spk1_68[k];

        t_69[k] = pb_x[k] * spk0_69[k]
                  - f_12 * pc_x[k] * spk1_69[k];

        t_70[k] = f_13 * spi_27[k]
                  + f_3 * pc_y[k] * sdi_55[k];

        t_71[k] = pb_x[k] * spk0_71[k]
                  - f_12 * pc_x[k] * spk1_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_z, pc_y, pc_z, spk0_0, spk0_3, \
                         spi_0, spk1_0, spk1_3, sdi_56, sdi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pb_z[k] * spk0_0[k]
                  - f_12 * pc_z[k] * spk1_0[k];

        t_73[k] = f_3 * pc_y[k] * sdi_56[k];

        t_74[k] = f_13 * spi_0[k]
                  + f_3 * pc_z[k] * sdi_56[k];

        t_75[k] = pb_z[k] * spk0_3[k]
                  - f_12 * pc_z[k] * spk1_3[k];

        t_76[k] = f_3 * pc_y[k] * sdi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_x, pb_z, pc_x, pc_z, spk0_6, spk0_77, spi_3, \
                         spi_61, spk1_6, spk1_77, sdi_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pb_x[k] * spk0_77[k]
                  + f_14 * spi_61[k]
                  - f_12 * pc_x[k] * spk1_77[k];

        t_78[k] = pb_z[k] * spk0_6[k]
                  - f_12 * pc_z[k] * spk1_6[k];

        t_79[k] = f_13 * spi_3[k]
                  + f_3 * pc_z[k] * sdi_59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, pb_z, pc_x, pc_y, pc_z, spk0_10, spk0_81, \
                         spi_65, spk1_10, spk1_81, sdi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_3 * pc_y[k] * sdi_61[k];

        t_81[k] = pb_x[k] * spk0_81[k]
                  + f_15 * spi_65[k]
                  - f_12 * pc_x[k] * spk1_81[k];

        t_82[k] = pb_z[k] * spk0_10[k]
                  - f_12 * pc_z[k] * spk1_10[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, pc_x, pc_y, pc_z, spk0_84, spi_6, spi_68, \
                         spk1_84, sdi_62, sdi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_13 * spi_6[k]
                  + f_3 * pc_z[k] * sdi_62[k];

        t_84[k] = pb_x[k] * spk0_84[k]
                  + f_16 * spi_68[k]
                  - f_12 * pc_x[k] * spk1_84[k];

        t_85[k] = f_3 * pc_y[k] * sdi_65[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_x, pb_z, pc_x, pc_z, spk0_15, spk0_86, spi_10, \
                         spi_70, spk1_15, spk1_86, sdi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * spk0_86[k]
                  + f_16 * spi_70[k]
                  - f_12 * pc_x[k] * spk1_86[k];

        t_87[k] = pb_z[k] * spk0_15[k]
                  - f_12 * pc_z[k] * spk1_15[k];

        t_88[k] = f_13 * spi_10[k]
                  + f_3 * pc_z[k] * sdi_66[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pb_x, pc_x, pc_y, spk0_89, spk0_90, spi_73, spi_74, \
                         spk1_89, spk1_90, sdi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_x[k] * spk0_89[k]
                  + f_0 * spi_73[k]
                  - f_12 * pc_x[k] * spk1_89[k];

        t_90[k] = pb_x[k] * spk0_90[k]
                  + f_0 * spi_74[k]
                  - f_12 * pc_x[k] * spk1_90[k];

        t_91[k] = f_3 * pc_y[k] * sdi_70[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pc_x, spk0_92, spi_76, spi_77, spi_78, \
                         spi_79, spk1_92, sdi_77, sdi_78, sdi_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_x[k] * spk0_92[k]
                  + f_0 * spi_76[k]
                  - f_12 * pc_x[k] * spk1_92[k];

        t_93[k] = f_13 * spi_77[k]
                  + f_3 * pc_x[k] * sdi_77[k];

        t_94[k] = f_13 * spi_78[k]
                  + f_3 * pc_x[k] * sdi_78[k];

        t_95[k] = f_13 * spi_79[k]
                  + f_3 * pc_x[k] * sdi_79[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, spi_80, spi_81, spi_82, spi_83, sdi_80, \
                         sdi_81, sdi_82, sdi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_13 * spi_80[k]
                  + f_3 * pc_x[k] * sdi_80[k];

        t_97[k] = f_13 * spi_81[k]
                  + f_3 * pc_x[k] * sdi_81[k];

        t_98[k] = f_13 * spi_82[k]
                  + f_3 * pc_x[k] * sdi_82[k];

        t_99[k] = f_13 * spi_83[k]
                  + f_3 * pc_x[k] * sdi_83[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_x, pc_x, pc_z, spk0_100, spk0_102, \
                         spk0_103, spi_21, spk1_100, spk1_102, spk1_103, \
                         sdi_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_x[k] * spk0_100[k]
                   - f_12 * pc_x[k] * spk1_100[k];

        t_101[k] = f_13 * spi_21[k]
                   + f_3 * pc_z[k] * sdi_77[k];

        t_102[k] = pb_x[k] * spk0_102[k]
                   - f_12 * pc_x[k] * spk1_102[k];

        t_103[k] = pb_x[k] * spk0_103[k]
                   - f_12 * pc_x[k] * spk1_103[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_x, pc_x, pc_y, spk0_104, spk0_105, \
                         spk0_107, spk1_104, spk1_105, spk1_107, \
                         sdi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * spk0_104[k]
                   - f_12 * pc_x[k] * spk1_104[k];

        t_105[k] = pb_x[k] * spk0_105[k]
                   - f_12 * pc_x[k] * spk1_105[k];

        t_106[k] = f_3 * pc_y[k] * sdi_83[k];

        t_107[k] = pb_x[k] * spk0_107[k]
                   - f_12 * pc_x[k] * spk1_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pc_x, pc_y, pc_z, spi_28, sdh0_63, \
                         sdh0_66, sdh1_63, sdh1_66, sdi_84, sdi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * sdh0_63[k]
                   - f_2 * sdh1_63[k]
                   + f_3 * pc_x[k] * sdi_84[k];

        t_109[k] = f_0 * spi_28[k]
                   + f_3 * pc_y[k] * sdi_84[k];

        t_110[k] = f_3 * pc_z[k] * sdi_84[k];

        t_111[k] = f_4 * sdh0_66[k]
                   - f_5 * sdh1_66[k]
                   + f_3 * pc_x[k] * sdi_87[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pc_x, pc_y, pc_z, spi_30, sdh0_68, \
                         sdh0_69, sdh1_68, sdh1_69, sdi_86, sdi_87, sdi_89, \
                         sdi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_0 * spi_30[k]
                   + f_3 * pc_y[k] * sdi_86[k];

        t_113[k] = f_4 * sdh0_68[k]
                   - f_5 * sdh1_68[k]
                   + f_3 * pc_x[k] * sdi_89[k];

        t_114[k] = f_6 * sdh0_69[k]
                   - f_7 * sdh1_69[k]
                   + f_3 * pc_x[k] * sdi_90[k];

        t_115[k] = f_3 * pc_z[k] * sdi_87[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pc_x, pc_y, pc_z, spi_33, sdh0_72, \
                         sdh0_73, sdh1_72, sdh1_73, sdi_89, sdi_90, sdi_93, \
                         sdi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_0 * spi_33[k]
                   + f_3 * pc_y[k] * sdi_89[k];

        t_117[k] = f_6 * sdh0_72[k]
                   - f_7 * sdh1_72[k]
                   + f_3 * pc_x[k] * sdi_93[k];

        t_118[k] = f_8 * sdh0_73[k]
                   - f_9 * sdh1_73[k]
                   + f_3 * pc_x[k] * sdi_94[k];

        t_119[k] = f_3 * pc_z[k] * sdi_90[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pc_x, pc_y, spi_37, sdh0_75, sdh0_77, sdh1_75, \
                         sdh1_77, sdi_93, sdi_96, sdi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * sdh0_75[k]
                   - f_9 * sdh1_75[k]
                   + f_3 * pc_x[k] * sdi_96[k];

        t_121[k] = f_0 * spi_37[k]
                   + f_3 * pc_y[k] * sdi_93[k];

        t_122[k] = f_8 * sdh0_77[k]
                   - f_9 * sdh1_77[k]
                   + f_3 * pc_x[k] * sdi_98[k];
    }
}

static auto
compute_prim_sdk_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t spk0,
                                                          const size_t spi, const size_t spk1,
                                                          const size_t sdh0, const size_t sdh1,
                                                          const size_t sdi, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 3.0 / gamma;
    const auto f_2 = 3.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.0 / gamma;
    const auto f_5 = 2.0 * p / (gamma * q);
    const auto f_6 = 1.5 / gamma;
    const auto f_7 = 1.5 * p / (gamma * q);
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);
    const auto f_10 = 0.5 / gamma;
    const auto f_11 = 0.5 * p / (gamma * q);
    const auto f_12 = gamma / q;
    const auto f_13 = 0.5 / q;
    const auto f_14 = 2.5 / q;
    const auto f_15 = 2.0 / q;
    const auto f_16 = 1.5 / q;

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
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spk0_39 = buffer.data(spk0 + 39);
    const auto *spk0_42 = buffer.data(spk0 + 42);
    const auto *spk0_46 = buffer.data(spk0 + 46);
    const auto *spk0_51 = buffer.data(spk0 + 51);
    const auto *spk0_64 = buffer.data(spk0 + 64);
    const auto *spk0_72 = buffer.data(spk0 + 72);
    const auto *spk0_77 = buffer.data(spk0 + 77);
    const auto *spk0_81 = buffer.data(spk0 + 81);
    const auto *spk0_86 = buffer.data(spk0 + 86);
    const auto *spk0_92 = buffer.data(spk0 + 92);
    const auto *spk0_102 = buffer.data(spk0 + 102);
    const auto *spk0_103 = buffer.data(spk0 + 103);
    const auto *spk0_104 = buffer.data(spk0 + 104);
    const auto *spk0_105 = buffer.data(spk0 + 105);
    const auto *spk0_107 = buffer.data(spk0 + 107);

    const auto *spi_28 = buffer.data(spi + 28);
    const auto *spi_31 = buffer.data(spi + 31);
    const auto *spi_34 = buffer.data(spi + 34);
    const auto *spi_38 = buffer.data(spi + 38);
    const auto *spi_42 = buffer.data(spi + 42);
    const auto *spi_49 = buffer.data(spi + 49);
    const auto *spi_51 = buffer.data(spi + 51);
    const auto *spi_52 = buffer.data(spi + 52);
    const auto *spi_53 = buffer.data(spi + 53);
    const auto *spi_54 = buffer.data(spi + 54);
    const auto *spi_55 = buffer.data(spi + 55);
    const auto *spi_56 = buffer.data(spi + 56);
    const auto *spi_58 = buffer.data(spi + 58);
    const auto *spi_59 = buffer.data(spi + 59);
    const auto *spi_61 = buffer.data(spi + 61);
    const auto *spi_62 = buffer.data(spi + 62);
    const auto *spi_65 = buffer.data(spi + 65);
    const auto *spi_66 = buffer.data(spi + 66);
    const auto *spi_70 = buffer.data(spi + 70);
    const auto *spi_77 = buffer.data(spi + 77);
    const auto *spi_79 = buffer.data(spi + 79);
    const auto *spi_80 = buffer.data(spi + 80);
    const auto *spi_81 = buffer.data(spi + 81);
    const auto *spi_82 = buffer.data(spi + 82);
    const auto *spi_83 = buffer.data(spi + 83);

    const auto *spk1_39 = buffer.data(spk1 + 39);
    const auto *spk1_42 = buffer.data(spk1 + 42);
    const auto *spk1_46 = buffer.data(spk1 + 46);
    const auto *spk1_51 = buffer.data(spk1 + 51);
    const auto *spk1_64 = buffer.data(spk1 + 64);
    const auto *spk1_72 = buffer.data(spk1 + 72);
    const auto *spk1_77 = buffer.data(spk1 + 77);
    const auto *spk1_81 = buffer.data(spk1 + 81);
    const auto *spk1_86 = buffer.data(spk1 + 86);
    const auto *spk1_92 = buffer.data(spk1 + 92);
    const auto *spk1_102 = buffer.data(spk1 + 102);
    const auto *spk1_103 = buffer.data(spk1 + 103);
    const auto *spk1_104 = buffer.data(spk1 + 104);
    const auto *spk1_105 = buffer.data(spk1 + 105);
    const auto *spk1_107 = buffer.data(spk1 + 107);

    const auto *sdh0_78 = buffer.data(sdh0 + 78);
    const auto *sdh0_80 = buffer.data(sdh0 + 80);
    const auto *sdh0_81 = buffer.data(sdh0 + 81);
    const auto *sdh0_82 = buffer.data(sdh0 + 82);
    const auto *sdh0_83 = buffer.data(sdh0 + 83);
    const auto *sdh0_96 = buffer.data(sdh0 + 96);
    const auto *sdh0_101 = buffer.data(sdh0 + 101);
    const auto *sdh0_102 = buffer.data(sdh0 + 102);
    const auto *sdh0_105 = buffer.data(sdh0 + 105);
    const auto *sdh0_108 = buffer.data(sdh0 + 108);
    const auto *sdh0_110 = buffer.data(sdh0 + 110);
    const auto *sdh0_111 = buffer.data(sdh0 + 111);
    const auto *sdh0_114 = buffer.data(sdh0 + 114);
    const auto *sdh0_115 = buffer.data(sdh0 + 115);
    const auto *sdh0_117 = buffer.data(sdh0 + 117);
    const auto *sdh0_119 = buffer.data(sdh0 + 119);
    const auto *sdh0_120 = buffer.data(sdh0 + 120);
    const auto *sdh0_122 = buffer.data(sdh0 + 122);
    const auto *sdh0_123 = buffer.data(sdh0 + 123);
    const auto *sdh0_124 = buffer.data(sdh0 + 124);
    const auto *sdh0_125 = buffer.data(sdh0 + 125);

    const auto *sdh1_78 = buffer.data(sdh1 + 78);
    const auto *sdh1_80 = buffer.data(sdh1 + 80);
    const auto *sdh1_81 = buffer.data(sdh1 + 81);
    const auto *sdh1_82 = buffer.data(sdh1 + 82);
    const auto *sdh1_83 = buffer.data(sdh1 + 83);
    const auto *sdh1_96 = buffer.data(sdh1 + 96);
    const auto *sdh1_101 = buffer.data(sdh1 + 101);
    const auto *sdh1_102 = buffer.data(sdh1 + 102);
    const auto *sdh1_105 = buffer.data(sdh1 + 105);
    const auto *sdh1_108 = buffer.data(sdh1 + 108);
    const auto *sdh1_110 = buffer.data(sdh1 + 110);
    const auto *sdh1_111 = buffer.data(sdh1 + 111);
    const auto *sdh1_114 = buffer.data(sdh1 + 114);
    const auto *sdh1_115 = buffer.data(sdh1 + 115);
    const auto *sdh1_117 = buffer.data(sdh1 + 117);
    const auto *sdh1_119 = buffer.data(sdh1 + 119);
    const auto *sdh1_120 = buffer.data(sdh1 + 120);
    const auto *sdh1_122 = buffer.data(sdh1 + 122);
    const auto *sdh1_123 = buffer.data(sdh1 + 123);
    const auto *sdh1_124 = buffer.data(sdh1 + 124);
    const auto *sdh1_125 = buffer.data(sdh1 + 125);

    const auto *sdi_94 = buffer.data(sdi + 94);
    const auto *sdi_98 = buffer.data(sdi + 98);
    const auto *sdi_99 = buffer.data(sdi + 99);
    const auto *sdi_101 = buffer.data(sdi + 101);
    const auto *sdi_102 = buffer.data(sdi + 102);
    const auto *sdi_104 = buffer.data(sdi + 104);
    const auto *sdi_105 = buffer.data(sdi + 105);
    const auto *sdi_106 = buffer.data(sdi + 106);
    const auto *sdi_107 = buffer.data(sdi + 107);
    const auto *sdi_108 = buffer.data(sdi + 108);
    const auto *sdi_109 = buffer.data(sdi + 109);
    const auto *sdi_110 = buffer.data(sdi + 110);
    const auto *sdi_111 = buffer.data(sdi + 111);
    const auto *sdi_112 = buffer.data(sdi + 112);
    const auto *sdi_114 = buffer.data(sdi + 114);
    const auto *sdi_115 = buffer.data(sdi + 115);
    const auto *sdi_117 = buffer.data(sdi + 117);
    const auto *sdi_118 = buffer.data(sdi + 118);
    const auto *sdi_121 = buffer.data(sdi + 121);
    const auto *sdi_122 = buffer.data(sdi + 122);
    const auto *sdi_124 = buffer.data(sdi + 124);
    const auto *sdi_126 = buffer.data(sdi + 126);
    const auto *sdi_129 = buffer.data(sdi + 129);
    const auto *sdi_130 = buffer.data(sdi + 130);
    const auto *sdi_133 = buffer.data(sdi + 133);
    const auto *sdi_134 = buffer.data(sdi + 134);
    const auto *sdi_135 = buffer.data(sdi + 135);
    const auto *sdi_136 = buffer.data(sdi + 136);
    const auto *sdi_137 = buffer.data(sdi + 137);
    const auto *sdi_138 = buffer.data(sdi + 138);
    const auto *sdi_139 = buffer.data(sdi + 139);
    const auto *sdi_140 = buffer.data(sdi + 140);
    const auto *sdi_142 = buffer.data(sdi + 142);
    const auto *sdi_143 = buffer.data(sdi + 143);
    const auto *sdi_145 = buffer.data(sdi + 145);
    const auto *sdi_146 = buffer.data(sdi + 146);
    const auto *sdi_149 = buffer.data(sdi + 149);
    const auto *sdi_150 = buffer.data(sdi + 150);
    const auto *sdi_152 = buffer.data(sdi + 152);
    const auto *sdi_154 = buffer.data(sdi + 154);
    const auto *sdi_155 = buffer.data(sdi + 155);
    const auto *sdi_157 = buffer.data(sdi + 157);
    const auto *sdi_158 = buffer.data(sdi + 158);
    const auto *sdi_160 = buffer.data(sdi + 160);
    const auto *sdi_161 = buffer.data(sdi + 161);
    const auto *sdi_162 = buffer.data(sdi + 162);
    const auto *sdi_163 = buffer.data(sdi + 163);
    const auto *sdi_164 = buffer.data(sdi + 164);
    const auto *sdi_165 = buffer.data(sdi + 165);
    const auto *sdi_166 = buffer.data(sdi + 166);
    const auto *sdi_167 = buffer.data(sdi + 167);

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pc_x, pc_z, sdh0_78, sdh0_80, sdh0_81, \
                         sdh1_78, sdh1_80, sdh1_81, sdi_94, sdi_99, sdi_101, \
                         sdi_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_10 * sdh0_78[k]
                   - f_11 * sdh1_78[k]
                   + f_3 * pc_x[k] * sdi_99[k];

        t_124[k] = f_3 * pc_z[k] * sdi_94[k];

        t_125[k] = f_10 * sdh0_80[k]
                   - f_11 * sdh1_80[k]
                   + f_3 * pc_x[k] * sdi_101[k];

        t_126[k] = f_10 * sdh0_81[k]
                   - f_11 * sdh1_81[k]
                   + f_3 * pc_x[k] * sdi_102[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, spi_42, sdh0_83, \
                         sdh1_83, sdi_98, sdi_104, sdi_105, sdi_106, \
                         sdi_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_0 * spi_42[k]
                   + f_3 * pc_y[k] * sdi_98[k];

        t_128[k] = f_10 * sdh0_83[k]
                   - f_11 * sdh1_83[k]
                   + f_3 * pc_x[k] * sdi_104[k];

        t_129[k] = f_3 * pc_x[k] * sdi_105[k];

        t_130[k] = f_3 * pc_x[k] * sdi_106[k];

        t_131[k] = f_3 * pc_x[k] * sdi_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pc_x, pc_y, spi_49, sdh0_78, \
                         sdh1_78, sdi_105, sdi_108, sdi_109, sdi_110, \
                         sdi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_3 * pc_x[k] * sdi_108[k];

        t_133[k] = f_3 * pc_x[k] * sdi_109[k];

        t_134[k] = f_3 * pc_x[k] * sdi_110[k];

        t_135[k] = f_3 * pc_x[k] * sdi_111[k];

        t_136[k] = f_0 * spi_49[k]
                   + f_1 * sdh0_78[k]
                   - f_2 * sdh1_78[k]
                   + f_3 * pc_y[k] * sdi_105[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, spi_51, spi_52, sdh0_80, sdh0_81, \
                         sdh1_80, sdh1_81, sdi_105, sdi_107, sdi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_3 * pc_z[k] * sdi_105[k];

        t_138[k] = f_0 * spi_51[k]
                   + f_4 * sdh0_80[k]
                   - f_5 * sdh1_80[k]
                   + f_3 * pc_y[k] * sdi_107[k];

        t_139[k] = f_0 * spi_52[k]
                   + f_6 * sdh0_81[k]
                   - f_7 * sdh1_81[k]
                   + f_3 * pc_y[k] * sdi_108[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_y, pc_z, spi_53, spi_54, spi_55, \
                         sdh0_82, sdh0_83, sdh1_82, sdh1_83, sdi_109, sdi_110, \
                         sdi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * spi_53[k]
                   + f_8 * sdh0_82[k]
                   - f_9 * sdh1_82[k]
                   + f_3 * pc_y[k] * sdi_109[k];

        t_141[k] = f_0 * spi_54[k]
                   + f_10 * sdh0_83[k]
                   - f_11 * sdh1_83[k]
                   + f_3 * pc_y[k] * sdi_110[k];

        t_142[k] = f_0 * spi_55[k]
                   + f_3 * pc_y[k] * sdi_111[k];

        t_143[k] = f_1 * sdh0_83[k]
                   - f_2 * sdh1_83[k]
                   + f_3 * pc_z[k] * sdi_111[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_y, pb_z, pc_y, pc_z, spk0_39, spk0_72, \
                         spi_28, spi_56, spk1_39, spk1_72, sdi_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_y[k] * spk0_72[k]
                   - f_12 * pc_y[k] * spk1_72[k];

        t_145[k] = f_13 * spi_56[k]
                   + f_3 * pc_y[k] * sdi_112[k];

        t_146[k] = f_13 * spi_28[k]
                   + f_3 * pc_z[k] * sdi_112[k];

        t_147[k] = pb_z[k] * spk0_39[k]
                   - f_12 * pc_z[k] * spk1_39[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_y, pb_z, pc_y, pc_z, spk0_42, spk0_77, \
                         spi_31, spi_58, spk1_42, spk1_77, sdi_114, \
                         sdi_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_13 * spi_58[k]
                   + f_3 * pc_y[k] * sdi_114[k];

        t_149[k] = pb_y[k] * spk0_77[k]
                   - f_12 * pc_y[k] * spk1_77[k];

        t_150[k] = pb_z[k] * spk0_42[k]
                   - f_12 * pc_z[k] * spk1_42[k];

        t_151[k] = f_13 * spi_31[k]
                   + f_3 * pc_z[k] * sdi_115[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, pc_y, pc_z, spk0_46, spk0_81, \
                         spi_34, spi_61, spk1_46, spk1_81, sdi_117, \
                         sdi_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_13 * spi_61[k]
                   + f_3 * pc_y[k] * sdi_117[k];

        t_153[k] = pb_y[k] * spk0_81[k]
                   - f_12 * pc_y[k] * spk1_81[k];

        t_154[k] = pb_z[k] * spk0_46[k]
                   - f_12 * pc_z[k] * spk1_46[k];

        t_155[k] = f_13 * spi_34[k]
                   + f_3 * pc_z[k] * sdi_118[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_y, pc_x, pc_y, spk0_86, spi_65, spk1_86, \
                         sdh0_96, sdh1_96, sdi_121, sdi_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_8 * sdh0_96[k]
                   - f_9 * sdh1_96[k]
                   + f_3 * pc_x[k] * sdi_124[k];

        t_157[k] = f_13 * spi_65[k]
                   + f_3 * pc_y[k] * sdi_121[k];

        t_158[k] = pb_y[k] * spk0_86[k]
                   - f_12 * pc_y[k] * spk1_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_z, pc_x, pc_z, spk0_51, spi_38, spk1_51, \
                         sdh0_101, sdh1_101, sdi_122, sdi_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pb_z[k] * spk0_51[k]
                   - f_12 * pc_z[k] * spk1_51[k];

        t_160[k] = f_13 * spi_38[k]
                   + f_3 * pc_z[k] * sdi_122[k];

        t_161[k] = f_10 * sdh0_101[k]
                   - f_11 * sdh1_101[k]
                   + f_3 * pc_x[k] * sdi_129[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_y, pc_x, pc_y, spk0_92, spi_70, \
                         spk1_92, sdh0_102, sdh1_102, sdi_126, sdi_130, \
                         sdi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_10 * sdh0_102[k]
                   - f_11 * sdh1_102[k]
                   + f_3 * pc_x[k] * sdi_130[k];

        t_163[k] = f_13 * spi_70[k]
                   + f_3 * pc_y[k] * sdi_126[k];

        t_164[k] = pb_y[k] * spk0_92[k]
                   - f_12 * pc_y[k] * spk1_92[k];

        t_165[k] = f_3 * pc_x[k] * sdi_133[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, pc_x, sdi_134, sdi_135, \
                         sdi_136, sdi_137, sdi_138, sdi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_3 * pc_x[k] * sdi_134[k];

        t_167[k] = f_3 * pc_x[k] * sdi_135[k];

        t_168[k] = f_3 * pc_x[k] * sdi_136[k];

        t_169[k] = f_3 * pc_x[k] * sdi_137[k];

        t_170[k] = f_3 * pc_x[k] * sdi_138[k];

        t_171[k] = f_3 * pc_x[k] * sdi_139[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_y, pb_z, pc_y, pc_z, spk0_64, spk0_102, \
                         spi_49, spi_79, spk1_64, spk1_102, sdi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * spk0_64[k]
                   - f_12 * pc_z[k] * spk1_64[k];

        t_173[k] = f_13 * spi_49[k]
                   + f_3 * pc_z[k] * sdi_133[k];

        t_174[k] = pb_y[k] * spk0_102[k]
                   + f_14 * spi_79[k]
                   - f_12 * pc_y[k] * spk1_102[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_y, pc_y, spk0_103, spk0_104, spk0_105, \
                         spi_80, spi_81, spi_82, spk1_103, spk1_104, \
                         spk1_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pb_y[k] * spk0_103[k]
                   + f_15 * spi_80[k]
                   - f_12 * pc_y[k] * spk1_103[k];

        t_176[k] = pb_y[k] * spk0_104[k]
                   + f_16 * spi_81[k]
                   - f_12 * pc_y[k] * spk1_104[k];

        t_177[k] = pb_y[k] * spk0_105[k]
                   + f_0 * spi_82[k]
                   - f_12 * pc_y[k] * spk1_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pb_y, pc_x, pc_y, spk0_107, spi_83, \
                         spk1_107, sdh0_105, sdh1_105, sdi_139, \
                         sdi_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_13 * spi_83[k]
                   + f_3 * pc_y[k] * sdi_139[k];

        t_179[k] = pb_y[k] * spk0_107[k]
                   - f_12 * pc_y[k] * spk1_107[k];

        t_180[k] = f_1 * sdh0_105[k]
                   - f_2 * sdh1_105[k]
                   + f_3 * pc_x[k] * sdi_140[k];

        t_181[k] = f_3 * pc_y[k] * sdi_140[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, spi_56, sdh0_108, \
                         sdh0_110, sdh1_108, sdh1_110, sdi_140, sdi_142, sdi_143, \
                         sdi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_0 * spi_56[k]
                   + f_3 * pc_z[k] * sdi_140[k];

        t_183[k] = f_4 * sdh0_108[k]
                   - f_5 * sdh1_108[k]
                   + f_3 * pc_x[k] * sdi_143[k];

        t_184[k] = f_3 * pc_y[k] * sdi_142[k];

        t_185[k] = f_4 * sdh0_110[k]
                   - f_5 * sdh1_110[k]
                   + f_3 * pc_x[k] * sdi_145[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, pc_y, pc_z, spi_59, sdh0_111, \
                         sdh0_114, sdh1_111, sdh1_114, sdi_143, sdi_145, sdi_146, \
                         sdi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_6 * sdh0_111[k]
                   - f_7 * sdh1_111[k]
                   + f_3 * pc_x[k] * sdi_146[k];

        t_187[k] = f_0 * spi_59[k]
                   + f_3 * pc_z[k] * sdi_143[k];

        t_188[k] = f_3 * pc_y[k] * sdi_145[k];

        t_189[k] = f_6 * sdh0_114[k]
                   - f_7 * sdh1_114[k]
                   + f_3 * pc_x[k] * sdi_149[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, spi_62, sdh0_115, \
                         sdh0_117, sdh1_115, sdh1_117, sdi_146, sdi_149, sdi_150, \
                         sdi_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * sdh0_115[k]
                   - f_9 * sdh1_115[k]
                   + f_3 * pc_x[k] * sdi_150[k];

        t_191[k] = f_0 * spi_62[k]
                   + f_3 * pc_z[k] * sdi_146[k];

        t_192[k] = f_8 * sdh0_117[k]
                   - f_9 * sdh1_117[k]
                   + f_3 * pc_x[k] * sdi_152[k];

        t_193[k] = f_3 * pc_y[k] * sdi_149[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pc_x, pc_z, spi_66, sdh0_119, sdh0_120, \
                         sdh1_119, sdh1_120, sdi_150, sdi_154, \
                         sdi_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_8 * sdh0_119[k]
                   - f_9 * sdh1_119[k]
                   + f_3 * pc_x[k] * sdi_154[k];

        t_195[k] = f_10 * sdh0_120[k]
                   - f_11 * sdh1_120[k]
                   + f_3 * pc_x[k] * sdi_155[k];

        t_196[k] = f_0 * spi_66[k]
                   + f_3 * pc_z[k] * sdi_150[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pc_x, pc_y, sdh0_122, sdh0_123, sdh0_125, \
                         sdh1_122, sdh1_123, sdh1_125, sdi_154, sdi_157, sdi_158, \
                         sdi_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_10 * sdh0_122[k]
                   - f_11 * sdh1_122[k]
                   + f_3 * pc_x[k] * sdi_157[k];

        t_198[k] = f_10 * sdh0_123[k]
                   - f_11 * sdh1_123[k]
                   + f_3 * pc_x[k] * sdi_158[k];

        t_199[k] = f_3 * pc_y[k] * sdi_154[k];

        t_200[k] = f_10 * sdh0_125[k]
                   - f_11 * sdh1_125[k]
                   + f_3 * pc_x[k] * sdi_160[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, t_207, pc_x, sdi_161, \
                         sdi_162, sdi_163, sdi_164, sdi_165, sdi_166, \
                         sdi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_3 * pc_x[k] * sdi_161[k];

        t_202[k] = f_3 * pc_x[k] * sdi_162[k];

        t_203[k] = f_3 * pc_x[k] * sdi_163[k];

        t_204[k] = f_3 * pc_x[k] * sdi_164[k];

        t_205[k] = f_3 * pc_x[k] * sdi_165[k];

        t_206[k] = f_3 * pc_x[k] * sdi_166[k];

        t_207[k] = f_3 * pc_x[k] * sdi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pc_y, pc_z, spi_77, sdh0_120, sdh0_122, \
                         sdh0_123, sdh1_120, sdh1_122, sdh1_123, sdi_161, sdi_163, \
                         sdi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * sdh0_120[k]
                   - f_2 * sdh1_120[k]
                   + f_3 * pc_y[k] * sdi_161[k];

        t_209[k] = f_0 * spi_77[k]
                   + f_3 * pc_z[k] * sdi_161[k];

        t_210[k] = f_4 * sdh0_122[k]
                   - f_5 * sdh1_122[k]
                   + f_3 * pc_y[k] * sdi_163[k];

        t_211[k] = f_6 * sdh0_123[k]
                   - f_7 * sdh1_123[k]
                   + f_3 * pc_y[k] * sdi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_y, pc_z, spi_83, sdh0_124, sdh0_125, \
                         sdh1_124, sdh1_125, sdi_165, sdi_166, \
                         sdi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * sdh0_124[k]
                   - f_9 * sdh1_124[k]
                   + f_3 * pc_y[k] * sdi_165[k];

        t_213[k] = f_10 * sdh0_125[k]
                   - f_11 * sdh1_125[k]
                   + f_3 * pc_y[k] * sdi_166[k];

        t_214[k] = f_3 * pc_y[k] * sdi_167[k];

        t_215[k] = f_0 * spi_83[k]
                   + f_1 * sdh0_125[k]
                   - f_2 * sdh1_125[k]
                   + f_3 * pc_z[k] * sdi_167[k];
    }
}

auto
compute_prim_sdk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t spk0, const size_t spi,
                                                   const size_t spk1, const size_t sdh0,
                                                   const size_t sdh1, const size_t sdi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sdk_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, spk0, spi,
                                                              spk1, sdh0, sdh1, sdi, ncols,
                                                              gamma, p, q);

    compute_prim_sdk_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, spk0, spi,
                                                              spk1, sdh0, sdh1, sdi, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
