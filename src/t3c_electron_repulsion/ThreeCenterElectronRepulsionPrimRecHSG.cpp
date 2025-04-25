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

#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsg,
                                 size_t idx_eri_0_fsg,
                                 size_t idx_eri_1_fsg,
                                 size_t idx_eri_1_gsf,
                                 size_t idx_eri_1_gsg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : FSG

    auto g_xxx_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg);

    auto g_xxx_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 1);

    auto g_xxx_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 2);

    auto g_xxx_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 3);

    auto g_xxx_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 4);

    auto g_xxx_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 5);

    auto g_xxx_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 6);

    auto g_xxx_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 7);

    auto g_xxx_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 8);

    auto g_xxx_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 9);

    auto g_xxx_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 10);

    auto g_xxx_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 11);

    auto g_xxx_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 12);

    auto g_xxx_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 13);

    auto g_xxx_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 14);

    auto g_xxy_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 15);

    auto g_xxy_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 17);

    auto g_xxy_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 20);

    auto g_xxy_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 24);

    auto g_xxz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 30);

    auto g_xxz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 31);

    auto g_xxz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 33);

    auto g_xxz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 36);

    auto g_xyy_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 46);

    auto g_xyy_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 48);

    auto g_xyy_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 49);

    auto g_xyy_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 51);

    auto g_xyy_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 52);

    auto g_xyy_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 53);

    auto g_xyy_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 55);

    auto g_xyy_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 56);

    auto g_xyy_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 57);

    auto g_xyy_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 58);

    auto g_xyy_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 59);

    auto g_xzz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 77);

    auto g_xzz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 79);

    auto g_xzz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 80);

    auto g_xzz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 82);

    auto g_xzz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 83);

    auto g_xzz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 84);

    auto g_xzz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 85);

    auto g_xzz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 86);

    auto g_xzz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 87);

    auto g_xzz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 88);

    auto g_xzz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 89);

    auto g_yyy_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 90);

    auto g_yyy_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 91);

    auto g_yyy_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 92);

    auto g_yyy_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 93);

    auto g_yyy_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 94);

    auto g_yyy_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 95);

    auto g_yyy_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 96);

    auto g_yyy_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 97);

    auto g_yyy_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 98);

    auto g_yyy_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 99);

    auto g_yyy_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 100);

    auto g_yyy_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 101);

    auto g_yyy_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 102);

    auto g_yyy_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 103);

    auto g_yyy_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 104);

    auto g_yyz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 106);

    auto g_yyz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 108);

    auto g_yyz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 111);

    auto g_yyz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 115);

    auto g_yzz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 120);

    auto g_yzz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 122);

    auto g_yzz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 124);

    auto g_yzz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 125);

    auto g_yzz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 127);

    auto g_yzz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 128);

    auto g_yzz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 129);

    auto g_yzz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 131);

    auto g_yzz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 132);

    auto g_yzz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 133);

    auto g_yzz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 134);

    auto g_zzz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 135);

    auto g_zzz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 136);

    auto g_zzz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 137);

    auto g_zzz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 138);

    auto g_zzz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 139);

    auto g_zzz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 140);

    auto g_zzz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 141);

    auto g_zzz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 142);

    auto g_zzz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 143);

    auto g_zzz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 144);

    auto g_zzz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 145);

    auto g_zzz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 146);

    auto g_zzz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 147);

    auto g_zzz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 148);

    auto g_zzz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 149);

    /// Set up components of auxilary buffer : FSG

    auto g_xxx_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg);

    auto g_xxx_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 1);

    auto g_xxx_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 2);

    auto g_xxx_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 3);

    auto g_xxx_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 4);

    auto g_xxx_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 5);

    auto g_xxx_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 6);

    auto g_xxx_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 7);

    auto g_xxx_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 8);

    auto g_xxx_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 9);

    auto g_xxx_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 10);

    auto g_xxx_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 11);

    auto g_xxx_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 12);

    auto g_xxx_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 13);

    auto g_xxx_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 14);

    auto g_xxy_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 15);

    auto g_xxy_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 17);

    auto g_xxy_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 20);

    auto g_xxy_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 24);

    auto g_xxz_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 30);

    auto g_xxz_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 31);

    auto g_xxz_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 33);

    auto g_xxz_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 36);

    auto g_xyy_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 46);

    auto g_xyy_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 48);

    auto g_xyy_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 49);

    auto g_xyy_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 51);

    auto g_xyy_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 52);

    auto g_xyy_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 53);

    auto g_xyy_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 55);

    auto g_xyy_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 56);

    auto g_xyy_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 57);

    auto g_xyy_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 58);

    auto g_xyy_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 59);

    auto g_xzz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 77);

    auto g_xzz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 79);

    auto g_xzz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 80);

    auto g_xzz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 82);

    auto g_xzz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 83);

    auto g_xzz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 84);

    auto g_xzz_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 85);

    auto g_xzz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 86);

    auto g_xzz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 87);

    auto g_xzz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 88);

    auto g_xzz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 89);

    auto g_yyy_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 90);

    auto g_yyy_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 91);

    auto g_yyy_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 92);

    auto g_yyy_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 93);

    auto g_yyy_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 94);

    auto g_yyy_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 95);

    auto g_yyy_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 96);

    auto g_yyy_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 97);

    auto g_yyy_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 98);

    auto g_yyy_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 99);

    auto g_yyy_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 100);

    auto g_yyy_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 101);

    auto g_yyy_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 102);

    auto g_yyy_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 103);

    auto g_yyy_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 104);

    auto g_yyz_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 106);

    auto g_yyz_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 108);

    auto g_yyz_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 111);

    auto g_yyz_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 115);

    auto g_yzz_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 120);

    auto g_yzz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 122);

    auto g_yzz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 124);

    auto g_yzz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 125);

    auto g_yzz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 127);

    auto g_yzz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 128);

    auto g_yzz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 129);

    auto g_yzz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 131);

    auto g_yzz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 132);

    auto g_yzz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 133);

    auto g_yzz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 134);

    auto g_zzz_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 135);

    auto g_zzz_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 136);

    auto g_zzz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 137);

    auto g_zzz_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 138);

    auto g_zzz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 139);

    auto g_zzz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 140);

    auto g_zzz_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 141);

    auto g_zzz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 142);

    auto g_zzz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 143);

    auto g_zzz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 144);

    auto g_zzz_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 145);

    auto g_zzz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 146);

    auto g_zzz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 147);

    auto g_zzz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 148);

    auto g_zzz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 149);

    /// Set up components of auxilary buffer : GSF

    auto g_xxxx_0_xxx_1 = pbuffer.data(idx_eri_1_gsf);

    auto g_xxxx_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 1);

    auto g_xxxx_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 2);

    auto g_xxxx_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 3);

    auto g_xxxx_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 4);

    auto g_xxxx_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 5);

    auto g_xxxx_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 6);

    auto g_xxxx_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 7);

    auto g_xxxx_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 8);

    auto g_xxxx_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 9);

    auto g_xxxz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 22);

    auto g_xxxz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 24);

    auto g_xxxz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 25);

    auto g_xxxz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 27);

    auto g_xxxz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 28);

    auto g_xxxz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 29);

    auto g_xxyy_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 30);

    auto g_xxyy_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 31);

    auto g_xxyy_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 32);

    auto g_xxyy_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 33);

    auto g_xxyy_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 34);

    auto g_xxyy_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 35);

    auto g_xxyy_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 36);

    auto g_xxyy_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 37);

    auto g_xxyy_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 38);

    auto g_xxyy_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 39);

    auto g_xxzz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 50);

    auto g_xxzz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 51);

    auto g_xxzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 52);

    auto g_xxzz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 53);

    auto g_xxzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 54);

    auto g_xxzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 55);

    auto g_xxzz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 56);

    auto g_xxzz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 57);

    auto g_xxzz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 58);

    auto g_xxzz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 59);

    auto g_xyyy_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 61);

    auto g_xyyy_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 63);

    auto g_xyyy_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 64);

    auto g_xyyy_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 66);

    auto g_xyyy_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 67);

    auto g_xyyy_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 68);

    auto g_xzzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 92);

    auto g_xzzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 94);

    auto g_xzzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 95);

    auto g_xzzz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 97);

    auto g_xzzz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 98);

    auto g_xzzz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 99);

    auto g_yyyy_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 100);

    auto g_yyyy_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 101);

    auto g_yyyy_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 102);

    auto g_yyyy_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 103);

    auto g_yyyy_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 104);

    auto g_yyyy_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 105);

    auto g_yyyy_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 106);

    auto g_yyyy_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 107);

    auto g_yyyy_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 108);

    auto g_yyyy_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 109);

    auto g_yyyz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 112);

    auto g_yyyz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 114);

    auto g_yyyz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 115);

    auto g_yyyz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 117);

    auto g_yyyz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 118);

    auto g_yyyz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 119);

    auto g_yyzz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 120);

    auto g_yyzz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 121);

    auto g_yyzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 122);

    auto g_yyzz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 123);

    auto g_yyzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 124);

    auto g_yyzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 125);

    auto g_yyzz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 126);

    auto g_yyzz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 127);

    auto g_yyzz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 128);

    auto g_yyzz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 129);

    auto g_yzzz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 131);

    auto g_yzzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 132);

    auto g_yzzz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 133);

    auto g_yzzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 134);

    auto g_yzzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 135);

    auto g_yzzz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 136);

    auto g_yzzz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 137);

    auto g_yzzz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 138);

    auto g_yzzz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 139);

    auto g_zzzz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 140);

    auto g_zzzz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 141);

    auto g_zzzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 142);

    auto g_zzzz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 143);

    auto g_zzzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 144);

    auto g_zzzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 145);

    auto g_zzzz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 146);

    auto g_zzzz_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 147);

    auto g_zzzz_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 148);

    auto g_zzzz_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 149);

    /// Set up components of auxilary buffer : GSG

    auto g_xxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg);

    auto g_xxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 1);

    auto g_xxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 2);

    auto g_xxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 3);

    auto g_xxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 4);

    auto g_xxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 5);

    auto g_xxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 6);

    auto g_xxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 7);

    auto g_xxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 8);

    auto g_xxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 9);

    auto g_xxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 10);

    auto g_xxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 11);

    auto g_xxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 12);

    auto g_xxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 13);

    auto g_xxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 14);

    auto g_xxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 15);

    auto g_xxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 16);

    auto g_xxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 17);

    auto g_xxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 18);

    auto g_xxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 20);

    auto g_xxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 21);

    auto g_xxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 24);

    auto g_xxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 25);

    auto g_xxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 30);

    auto g_xxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 31);

    auto g_xxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 32);

    auto g_xxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 33);

    auto g_xxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 34);

    auto g_xxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 35);

    auto g_xxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 36);

    auto g_xxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 37);

    auto g_xxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 38);

    auto g_xxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 39);

    auto g_xxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 41);

    auto g_xxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 42);

    auto g_xxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 43);

    auto g_xxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 44);

    auto g_xxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 45);

    auto g_xxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 46);

    auto g_xxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 47);

    auto g_xxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 48);

    auto g_xxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 49);

    auto g_xxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 50);

    auto g_xxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 51);

    auto g_xxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 52);

    auto g_xxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 53);

    auto g_xxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 54);

    auto g_xxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 55);

    auto g_xxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 56);

    auto g_xxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 57);

    auto g_xxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 58);

    auto g_xxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 59);

    auto g_xxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 75);

    auto g_xxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 76);

    auto g_xxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 77);

    auto g_xxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 78);

    auto g_xxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 79);

    auto g_xxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 80);

    auto g_xxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 81);

    auto g_xxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 82);

    auto g_xxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 83);

    auto g_xxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 84);

    auto g_xxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 85);

    auto g_xxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 86);

    auto g_xxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 87);

    auto g_xxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 88);

    auto g_xxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 89);

    auto g_xyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 90);

    auto g_xyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 91);

    auto g_xyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 93);

    auto g_xyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 94);

    auto g_xyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 96);

    auto g_xyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 97);

    auto g_xyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 98);

    auto g_xyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 100);

    auto g_xyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 101);

    auto g_xyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 102);

    auto g_xyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 103);

    auto g_xyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 104);

    auto g_xzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 135);

    auto g_xzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 137);

    auto g_xzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 139);

    auto g_xzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 140);

    auto g_xzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 142);

    auto g_xzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 143);

    auto g_xzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 144);

    auto g_xzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 145);

    auto g_xzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 146);

    auto g_xzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 147);

    auto g_xzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 148);

    auto g_xzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 149);

    auto g_yyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 150);

    auto g_yyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 151);

    auto g_yyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 152);

    auto g_yyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 153);

    auto g_yyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 154);

    auto g_yyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 155);

    auto g_yyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 156);

    auto g_yyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 157);

    auto g_yyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 158);

    auto g_yyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 159);

    auto g_yyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 160);

    auto g_yyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 161);

    auto g_yyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 162);

    auto g_yyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 163);

    auto g_yyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 164);

    auto g_yyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 166);

    auto g_yyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 167);

    auto g_yyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 168);

    auto g_yyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 169);

    auto g_yyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 170);

    auto g_yyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 171);

    auto g_yyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 172);

    auto g_yyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 173);

    auto g_yyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 174);

    auto g_yyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 175);

    auto g_yyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 176);

    auto g_yyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 177);

    auto g_yyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 178);

    auto g_yyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 179);

    auto g_yyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 180);

    auto g_yyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 181);

    auto g_yyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 182);

    auto g_yyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 183);

    auto g_yyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 184);

    auto g_yyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 185);

    auto g_yyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 186);

    auto g_yyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 187);

    auto g_yyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 188);

    auto g_yyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 189);

    auto g_yyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 190);

    auto g_yyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 191);

    auto g_yyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 192);

    auto g_yyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 193);

    auto g_yyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 194);

    auto g_yzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 195);

    auto g_yzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 196);

    auto g_yzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 197);

    auto g_yzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 198);

    auto g_yzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 199);

    auto g_yzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 200);

    auto g_yzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 201);

    auto g_yzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 202);

    auto g_yzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 203);

    auto g_yzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 204);

    auto g_yzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 205);

    auto g_yzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 206);

    auto g_yzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 207);

    auto g_yzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 208);

    auto g_yzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 209);

    auto g_zzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_gsg + 210);

    auto g_zzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_gsg + 211);

    auto g_zzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_gsg + 212);

    auto g_zzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_gsg + 213);

    auto g_zzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_gsg + 214);

    auto g_zzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_gsg + 215);

    auto g_zzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_gsg + 216);

    auto g_zzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_gsg + 217);

    auto g_zzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_gsg + 218);

    auto g_zzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_gsg + 219);

    auto g_zzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_gsg + 220);

    auto g_zzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_gsg + 221);

    auto g_zzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_gsg + 222);

    auto g_zzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_gsg + 223);

    auto g_zzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_gsg + 224);

    /// Set up 0-15 components of targeted buffer : HSG

    auto g_xxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg);

    auto g_xxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 1);

    auto g_xxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 2);

    auto g_xxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 3);

    auto g_xxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 4);

    auto g_xxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 5);

    auto g_xxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 6);

    auto g_xxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 7);

    auto g_xxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 8);

    auto g_xxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 9);

    auto g_xxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 10);

    auto g_xxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 11);

    auto g_xxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 12);

    auto g_xxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 13);

    auto g_xxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 14);

    #pragma omp simd aligned(g_xxx_0_xxxx_0, g_xxx_0_xxxx_1, g_xxx_0_xxxy_0, g_xxx_0_xxxy_1, g_xxx_0_xxxz_0, g_xxx_0_xxxz_1, g_xxx_0_xxyy_0, g_xxx_0_xxyy_1, g_xxx_0_xxyz_0, g_xxx_0_xxyz_1, g_xxx_0_xxzz_0, g_xxx_0_xxzz_1, g_xxx_0_xyyy_0, g_xxx_0_xyyy_1, g_xxx_0_xyyz_0, g_xxx_0_xyyz_1, g_xxx_0_xyzz_0, g_xxx_0_xyzz_1, g_xxx_0_xzzz_0, g_xxx_0_xzzz_1, g_xxx_0_yyyy_0, g_xxx_0_yyyy_1, g_xxx_0_yyyz_0, g_xxx_0_yyyz_1, g_xxx_0_yyzz_0, g_xxx_0_yyzz_1, g_xxx_0_yzzz_0, g_xxx_0_yzzz_1, g_xxx_0_zzzz_0, g_xxx_0_zzzz_1, g_xxxx_0_xxx_1, g_xxxx_0_xxxx_1, g_xxxx_0_xxxy_1, g_xxxx_0_xxxz_1, g_xxxx_0_xxy_1, g_xxxx_0_xxyy_1, g_xxxx_0_xxyz_1, g_xxxx_0_xxz_1, g_xxxx_0_xxzz_1, g_xxxx_0_xyy_1, g_xxxx_0_xyyy_1, g_xxxx_0_xyyz_1, g_xxxx_0_xyz_1, g_xxxx_0_xyzz_1, g_xxxx_0_xzz_1, g_xxxx_0_xzzz_1, g_xxxx_0_yyy_1, g_xxxx_0_yyyy_1, g_xxxx_0_yyyz_1, g_xxxx_0_yyz_1, g_xxxx_0_yyzz_1, g_xxxx_0_yzz_1, g_xxxx_0_yzzz_1, g_xxxx_0_zzz_1, g_xxxx_0_zzzz_1, g_xxxxx_0_xxxx_0, g_xxxxx_0_xxxy_0, g_xxxxx_0_xxxz_0, g_xxxxx_0_xxyy_0, g_xxxxx_0_xxyz_0, g_xxxxx_0_xxzz_0, g_xxxxx_0_xyyy_0, g_xxxxx_0_xyyz_0, g_xxxxx_0_xyzz_0, g_xxxxx_0_xzzz_0, g_xxxxx_0_yyyy_0, g_xxxxx_0_yyyz_0, g_xxxxx_0_yyzz_0, g_xxxxx_0_yzzz_0, g_xxxxx_0_zzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xxxx_0[i] = 4.0 * g_xxx_0_xxxx_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxx_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxx_1[i] * fi_acd_0 + g_xxxx_0_xxxx_1[i] * wa_x[i];

        g_xxxxx_0_xxxy_0[i] = 4.0 * g_xxx_0_xxxy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxy_1[i] * fi_acd_0 + g_xxxx_0_xxxy_1[i] * wa_x[i];

        g_xxxxx_0_xxxz_0[i] = 4.0 * g_xxx_0_xxxz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxz_1[i] * fi_acd_0 + g_xxxx_0_xxxz_1[i] * wa_x[i];

        g_xxxxx_0_xxyy_0[i] = 4.0 * g_xxx_0_xxyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyy_1[i] * fi_acd_0 + g_xxxx_0_xxyy_1[i] * wa_x[i];

        g_xxxxx_0_xxyz_0[i] = 4.0 * g_xxx_0_xxyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyz_1[i] * fi_acd_0 + g_xxxx_0_xxyz_1[i] * wa_x[i];

        g_xxxxx_0_xxzz_0[i] = 4.0 * g_xxx_0_xxzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xzz_1[i] * fi_acd_0 + g_xxxx_0_xxzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyy_0[i] = 4.0 * g_xxx_0_xyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyy_1[i] * fz_be_0 + g_xxxx_0_yyy_1[i] * fi_acd_0 + g_xxxx_0_xyyy_1[i] * wa_x[i];

        g_xxxxx_0_xyyz_0[i] = 4.0 * g_xxx_0_xyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyz_1[i] * fz_be_0 + g_xxxx_0_yyz_1[i] * fi_acd_0 + g_xxxx_0_xyyz_1[i] * wa_x[i];

        g_xxxxx_0_xyzz_0[i] = 4.0 * g_xxx_0_xyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyzz_1[i] * fz_be_0 + g_xxxx_0_yzz_1[i] * fi_acd_0 + g_xxxx_0_xyzz_1[i] * wa_x[i];

        g_xxxxx_0_xzzz_0[i] = 4.0 * g_xxx_0_xzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xzzz_1[i] * fz_be_0 + g_xxxx_0_zzz_1[i] * fi_acd_0 + g_xxxx_0_xzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyy_0[i] = 4.0 * g_xxx_0_yyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyy_1[i] * fz_be_0 + g_xxxx_0_yyyy_1[i] * wa_x[i];

        g_xxxxx_0_yyyz_0[i] = 4.0 * g_xxx_0_yyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyz_1[i] * fz_be_0 + g_xxxx_0_yyyz_1[i] * wa_x[i];

        g_xxxxx_0_yyzz_0[i] = 4.0 * g_xxx_0_yyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyzz_1[i] * fz_be_0 + g_xxxx_0_yyzz_1[i] * wa_x[i];

        g_xxxxx_0_yzzz_0[i] = 4.0 * g_xxx_0_yzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yzzz_1[i] * fz_be_0 + g_xxxx_0_yzzz_1[i] * wa_x[i];

        g_xxxxx_0_zzzz_0[i] = 4.0 * g_xxx_0_zzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_zzzz_1[i] * fz_be_0 + g_xxxx_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : HSG

    auto g_xxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 15);

    auto g_xxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 16);

    auto g_xxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 17);

    auto g_xxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 18);

    auto g_xxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 19);

    auto g_xxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 20);

    auto g_xxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 21);

    auto g_xxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 22);

    auto g_xxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 23);

    auto g_xxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 24);

    auto g_xxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 25);

    auto g_xxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 26);

    auto g_xxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 27);

    auto g_xxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 28);

    auto g_xxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 29);

    #pragma omp simd aligned(g_xxxx_0_xxx_1, g_xxxx_0_xxxx_1, g_xxxx_0_xxxy_1, g_xxxx_0_xxxz_1, g_xxxx_0_xxy_1, g_xxxx_0_xxyy_1, g_xxxx_0_xxyz_1, g_xxxx_0_xxz_1, g_xxxx_0_xxzz_1, g_xxxx_0_xyy_1, g_xxxx_0_xyyy_1, g_xxxx_0_xyyz_1, g_xxxx_0_xyz_1, g_xxxx_0_xyzz_1, g_xxxx_0_xzz_1, g_xxxx_0_xzzz_1, g_xxxx_0_yyy_1, g_xxxx_0_yyyy_1, g_xxxx_0_yyyz_1, g_xxxx_0_yyz_1, g_xxxx_0_yyzz_1, g_xxxx_0_yzz_1, g_xxxx_0_yzzz_1, g_xxxx_0_zzz_1, g_xxxx_0_zzzz_1, g_xxxxy_0_xxxx_0, g_xxxxy_0_xxxy_0, g_xxxxy_0_xxxz_0, g_xxxxy_0_xxyy_0, g_xxxxy_0_xxyz_0, g_xxxxy_0_xxzz_0, g_xxxxy_0_xyyy_0, g_xxxxy_0_xyyz_0, g_xxxxy_0_xyzz_0, g_xxxxy_0_xzzz_0, g_xxxxy_0_yyyy_0, g_xxxxy_0_yyyz_0, g_xxxxy_0_yyzz_0, g_xxxxy_0_yzzz_0, g_xxxxy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xxxx_0[i] = g_xxxx_0_xxxx_1[i] * wa_y[i];

        g_xxxxy_0_xxxy_0[i] = g_xxxx_0_xxx_1[i] * fi_acd_0 + g_xxxx_0_xxxy_1[i] * wa_y[i];

        g_xxxxy_0_xxxz_0[i] = g_xxxx_0_xxxz_1[i] * wa_y[i];

        g_xxxxy_0_xxyy_0[i] = 2.0 * g_xxxx_0_xxy_1[i] * fi_acd_0 + g_xxxx_0_xxyy_1[i] * wa_y[i];

        g_xxxxy_0_xxyz_0[i] = g_xxxx_0_xxz_1[i] * fi_acd_0 + g_xxxx_0_xxyz_1[i] * wa_y[i];

        g_xxxxy_0_xxzz_0[i] = g_xxxx_0_xxzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyy_0[i] = 3.0 * g_xxxx_0_xyy_1[i] * fi_acd_0 + g_xxxx_0_xyyy_1[i] * wa_y[i];

        g_xxxxy_0_xyyz_0[i] = 2.0 * g_xxxx_0_xyz_1[i] * fi_acd_0 + g_xxxx_0_xyyz_1[i] * wa_y[i];

        g_xxxxy_0_xyzz_0[i] = g_xxxx_0_xzz_1[i] * fi_acd_0 + g_xxxx_0_xyzz_1[i] * wa_y[i];

        g_xxxxy_0_xzzz_0[i] = g_xxxx_0_xzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyy_0[i] = 4.0 * g_xxxx_0_yyy_1[i] * fi_acd_0 + g_xxxx_0_yyyy_1[i] * wa_y[i];

        g_xxxxy_0_yyyz_0[i] = 3.0 * g_xxxx_0_yyz_1[i] * fi_acd_0 + g_xxxx_0_yyyz_1[i] * wa_y[i];

        g_xxxxy_0_yyzz_0[i] = 2.0 * g_xxxx_0_yzz_1[i] * fi_acd_0 + g_xxxx_0_yyzz_1[i] * wa_y[i];

        g_xxxxy_0_yzzz_0[i] = g_xxxx_0_zzz_1[i] * fi_acd_0 + g_xxxx_0_yzzz_1[i] * wa_y[i];

        g_xxxxy_0_zzzz_0[i] = g_xxxx_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 30-45 components of targeted buffer : HSG

    auto g_xxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 30);

    auto g_xxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 31);

    auto g_xxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 32);

    auto g_xxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 33);

    auto g_xxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 34);

    auto g_xxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 35);

    auto g_xxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 36);

    auto g_xxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 37);

    auto g_xxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 38);

    auto g_xxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 39);

    auto g_xxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 40);

    auto g_xxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 41);

    auto g_xxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 42);

    auto g_xxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 43);

    auto g_xxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 44);

    #pragma omp simd aligned(g_xxxx_0_xxx_1, g_xxxx_0_xxxx_1, g_xxxx_0_xxxy_1, g_xxxx_0_xxxz_1, g_xxxx_0_xxy_1, g_xxxx_0_xxyy_1, g_xxxx_0_xxyz_1, g_xxxx_0_xxz_1, g_xxxx_0_xxzz_1, g_xxxx_0_xyy_1, g_xxxx_0_xyyy_1, g_xxxx_0_xyyz_1, g_xxxx_0_xyz_1, g_xxxx_0_xyzz_1, g_xxxx_0_xzz_1, g_xxxx_0_xzzz_1, g_xxxx_0_yyy_1, g_xxxx_0_yyyy_1, g_xxxx_0_yyyz_1, g_xxxx_0_yyz_1, g_xxxx_0_yyzz_1, g_xxxx_0_yzz_1, g_xxxx_0_yzzz_1, g_xxxx_0_zzz_1, g_xxxx_0_zzzz_1, g_xxxxz_0_xxxx_0, g_xxxxz_0_xxxy_0, g_xxxxz_0_xxxz_0, g_xxxxz_0_xxyy_0, g_xxxxz_0_xxyz_0, g_xxxxz_0_xxzz_0, g_xxxxz_0_xyyy_0, g_xxxxz_0_xyyz_0, g_xxxxz_0_xyzz_0, g_xxxxz_0_xzzz_0, g_xxxxz_0_yyyy_0, g_xxxxz_0_yyyz_0, g_xxxxz_0_yyzz_0, g_xxxxz_0_yzzz_0, g_xxxxz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xxxx_0[i] = g_xxxx_0_xxxx_1[i] * wa_z[i];

        g_xxxxz_0_xxxy_0[i] = g_xxxx_0_xxxy_1[i] * wa_z[i];

        g_xxxxz_0_xxxz_0[i] = g_xxxx_0_xxx_1[i] * fi_acd_0 + g_xxxx_0_xxxz_1[i] * wa_z[i];

        g_xxxxz_0_xxyy_0[i] = g_xxxx_0_xxyy_1[i] * wa_z[i];

        g_xxxxz_0_xxyz_0[i] = g_xxxx_0_xxy_1[i] * fi_acd_0 + g_xxxx_0_xxyz_1[i] * wa_z[i];

        g_xxxxz_0_xxzz_0[i] = 2.0 * g_xxxx_0_xxz_1[i] * fi_acd_0 + g_xxxx_0_xxzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyy_0[i] = g_xxxx_0_xyyy_1[i] * wa_z[i];

        g_xxxxz_0_xyyz_0[i] = g_xxxx_0_xyy_1[i] * fi_acd_0 + g_xxxx_0_xyyz_1[i] * wa_z[i];

        g_xxxxz_0_xyzz_0[i] = 2.0 * g_xxxx_0_xyz_1[i] * fi_acd_0 + g_xxxx_0_xyzz_1[i] * wa_z[i];

        g_xxxxz_0_xzzz_0[i] = 3.0 * g_xxxx_0_xzz_1[i] * fi_acd_0 + g_xxxx_0_xzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyy_0[i] = g_xxxx_0_yyyy_1[i] * wa_z[i];

        g_xxxxz_0_yyyz_0[i] = g_xxxx_0_yyy_1[i] * fi_acd_0 + g_xxxx_0_yyyz_1[i] * wa_z[i];

        g_xxxxz_0_yyzz_0[i] = 2.0 * g_xxxx_0_yyz_1[i] * fi_acd_0 + g_xxxx_0_yyzz_1[i] * wa_z[i];

        g_xxxxz_0_yzzz_0[i] = 3.0 * g_xxxx_0_yzz_1[i] * fi_acd_0 + g_xxxx_0_yzzz_1[i] * wa_z[i];

        g_xxxxz_0_zzzz_0[i] = 4.0 * g_xxxx_0_zzz_1[i] * fi_acd_0 + g_xxxx_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 45-60 components of targeted buffer : HSG

    auto g_xxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 45);

    auto g_xxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 46);

    auto g_xxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 47);

    auto g_xxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 48);

    auto g_xxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 49);

    auto g_xxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 50);

    auto g_xxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 51);

    auto g_xxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 52);

    auto g_xxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 53);

    auto g_xxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 54);

    auto g_xxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 55);

    auto g_xxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 56);

    auto g_xxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 57);

    auto g_xxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 58);

    auto g_xxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 59);

    #pragma omp simd aligned(g_xxx_0_xxxx_0, g_xxx_0_xxxx_1, g_xxx_0_xxxz_0, g_xxx_0_xxxz_1, g_xxx_0_xxzz_0, g_xxx_0_xxzz_1, g_xxx_0_xzzz_0, g_xxx_0_xzzz_1, g_xxxy_0_xxxx_1, g_xxxy_0_xxxz_1, g_xxxy_0_xxzz_1, g_xxxy_0_xzzz_1, g_xxxyy_0_xxxx_0, g_xxxyy_0_xxxy_0, g_xxxyy_0_xxxz_0, g_xxxyy_0_xxyy_0, g_xxxyy_0_xxyz_0, g_xxxyy_0_xxzz_0, g_xxxyy_0_xyyy_0, g_xxxyy_0_xyyz_0, g_xxxyy_0_xyzz_0, g_xxxyy_0_xzzz_0, g_xxxyy_0_yyyy_0, g_xxxyy_0_yyyz_0, g_xxxyy_0_yyzz_0, g_xxxyy_0_yzzz_0, g_xxxyy_0_zzzz_0, g_xxyy_0_xxxy_1, g_xxyy_0_xxy_1, g_xxyy_0_xxyy_1, g_xxyy_0_xxyz_1, g_xxyy_0_xyy_1, g_xxyy_0_xyyy_1, g_xxyy_0_xyyz_1, g_xxyy_0_xyz_1, g_xxyy_0_xyzz_1, g_xxyy_0_yyy_1, g_xxyy_0_yyyy_1, g_xxyy_0_yyyz_1, g_xxyy_0_yyz_1, g_xxyy_0_yyzz_1, g_xxyy_0_yzz_1, g_xxyy_0_yzzz_1, g_xxyy_0_zzzz_1, g_xyy_0_xxxy_0, g_xyy_0_xxxy_1, g_xyy_0_xxyy_0, g_xyy_0_xxyy_1, g_xyy_0_xxyz_0, g_xyy_0_xxyz_1, g_xyy_0_xyyy_0, g_xyy_0_xyyy_1, g_xyy_0_xyyz_0, g_xyy_0_xyyz_1, g_xyy_0_xyzz_0, g_xyy_0_xyzz_1, g_xyy_0_yyyy_0, g_xyy_0_yyyy_1, g_xyy_0_yyyz_0, g_xyy_0_yyyz_1, g_xyy_0_yyzz_0, g_xyy_0_yyzz_1, g_xyy_0_yzzz_0, g_xyy_0_yzzz_1, g_xyy_0_zzzz_0, g_xyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xxxx_0[i] = g_xxx_0_xxxx_0[i] * fbe_0 - g_xxx_0_xxxx_1[i] * fz_be_0 + g_xxxy_0_xxxx_1[i] * wa_y[i];

        g_xxxyy_0_xxxy_0[i] = 2.0 * g_xyy_0_xxxy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxy_1[i] * fi_acd_0 + g_xxyy_0_xxxy_1[i] * wa_x[i];

        g_xxxyy_0_xxxz_0[i] = g_xxx_0_xxxz_0[i] * fbe_0 - g_xxx_0_xxxz_1[i] * fz_be_0 + g_xxxy_0_xxxz_1[i] * wa_y[i];

        g_xxxyy_0_xxyy_0[i] = 2.0 * g_xyy_0_xxyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyy_1[i] * fi_acd_0 + g_xxyy_0_xxyy_1[i] * wa_x[i];

        g_xxxyy_0_xxyz_0[i] = 2.0 * g_xyy_0_xxyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyz_1[i] * fi_acd_0 + g_xxyy_0_xxyz_1[i] * wa_x[i];

        g_xxxyy_0_xxzz_0[i] = g_xxx_0_xxzz_0[i] * fbe_0 - g_xxx_0_xxzz_1[i] * fz_be_0 + g_xxxy_0_xxzz_1[i] * wa_y[i];

        g_xxxyy_0_xyyy_0[i] = 2.0 * g_xyy_0_xyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyy_1[i] * fz_be_0 + g_xxyy_0_yyy_1[i] * fi_acd_0 + g_xxyy_0_xyyy_1[i] * wa_x[i];

        g_xxxyy_0_xyyz_0[i] = 2.0 * g_xyy_0_xyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyz_1[i] * fz_be_0 + g_xxyy_0_yyz_1[i] * fi_acd_0 + g_xxyy_0_xyyz_1[i] * wa_x[i];

        g_xxxyy_0_xyzz_0[i] = 2.0 * g_xyy_0_xyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyzz_1[i] * fz_be_0 + g_xxyy_0_yzz_1[i] * fi_acd_0 + g_xxyy_0_xyzz_1[i] * wa_x[i];

        g_xxxyy_0_xzzz_0[i] = g_xxx_0_xzzz_0[i] * fbe_0 - g_xxx_0_xzzz_1[i] * fz_be_0 + g_xxxy_0_xzzz_1[i] * wa_y[i];

        g_xxxyy_0_yyyy_0[i] = 2.0 * g_xyy_0_yyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyy_1[i] * fz_be_0 + g_xxyy_0_yyyy_1[i] * wa_x[i];

        g_xxxyy_0_yyyz_0[i] = 2.0 * g_xyy_0_yyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyz_1[i] * fz_be_0 + g_xxyy_0_yyyz_1[i] * wa_x[i];

        g_xxxyy_0_yyzz_0[i] = 2.0 * g_xyy_0_yyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyzz_1[i] * fz_be_0 + g_xxyy_0_yyzz_1[i] * wa_x[i];

        g_xxxyy_0_yzzz_0[i] = 2.0 * g_xyy_0_yzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yzzz_1[i] * fz_be_0 + g_xxyy_0_yzzz_1[i] * wa_x[i];

        g_xxxyy_0_zzzz_0[i] = 2.0 * g_xyy_0_zzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_zzzz_1[i] * fz_be_0 + g_xxyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 60-75 components of targeted buffer : HSG

    auto g_xxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 60);

    auto g_xxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 61);

    auto g_xxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 62);

    auto g_xxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 63);

    auto g_xxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 64);

    auto g_xxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 65);

    auto g_xxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 66);

    auto g_xxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 67);

    auto g_xxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 68);

    auto g_xxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 69);

    auto g_xxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 70);

    auto g_xxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 71);

    auto g_xxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 72);

    auto g_xxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 73);

    auto g_xxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 74);

    #pragma omp simd aligned(g_xxxy_0_xxxy_1, g_xxxy_0_xxyy_1, g_xxxy_0_xyyy_1, g_xxxy_0_yyyy_1, g_xxxyz_0_xxxx_0, g_xxxyz_0_xxxy_0, g_xxxyz_0_xxxz_0, g_xxxyz_0_xxyy_0, g_xxxyz_0_xxyz_0, g_xxxyz_0_xxzz_0, g_xxxyz_0_xyyy_0, g_xxxyz_0_xyyz_0, g_xxxyz_0_xyzz_0, g_xxxyz_0_xzzz_0, g_xxxyz_0_yyyy_0, g_xxxyz_0_yyyz_0, g_xxxyz_0_yyzz_0, g_xxxyz_0_yzzz_0, g_xxxyz_0_zzzz_0, g_xxxz_0_xxxx_1, g_xxxz_0_xxxz_1, g_xxxz_0_xxyz_1, g_xxxz_0_xxz_1, g_xxxz_0_xxzz_1, g_xxxz_0_xyyz_1, g_xxxz_0_xyz_1, g_xxxz_0_xyzz_1, g_xxxz_0_xzz_1, g_xxxz_0_xzzz_1, g_xxxz_0_yyyz_1, g_xxxz_0_yyz_1, g_xxxz_0_yyzz_1, g_xxxz_0_yzz_1, g_xxxz_0_yzzz_1, g_xxxz_0_zzz_1, g_xxxz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xxxx_0[i] = g_xxxz_0_xxxx_1[i] * wa_y[i];

        g_xxxyz_0_xxxy_0[i] = g_xxxy_0_xxxy_1[i] * wa_z[i];

        g_xxxyz_0_xxxz_0[i] = g_xxxz_0_xxxz_1[i] * wa_y[i];

        g_xxxyz_0_xxyy_0[i] = g_xxxy_0_xxyy_1[i] * wa_z[i];

        g_xxxyz_0_xxyz_0[i] = g_xxxz_0_xxz_1[i] * fi_acd_0 + g_xxxz_0_xxyz_1[i] * wa_y[i];

        g_xxxyz_0_xxzz_0[i] = g_xxxz_0_xxzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyy_0[i] = g_xxxy_0_xyyy_1[i] * wa_z[i];

        g_xxxyz_0_xyyz_0[i] = 2.0 * g_xxxz_0_xyz_1[i] * fi_acd_0 + g_xxxz_0_xyyz_1[i] * wa_y[i];

        g_xxxyz_0_xyzz_0[i] = g_xxxz_0_xzz_1[i] * fi_acd_0 + g_xxxz_0_xyzz_1[i] * wa_y[i];

        g_xxxyz_0_xzzz_0[i] = g_xxxz_0_xzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyy_0[i] = g_xxxy_0_yyyy_1[i] * wa_z[i];

        g_xxxyz_0_yyyz_0[i] = 3.0 * g_xxxz_0_yyz_1[i] * fi_acd_0 + g_xxxz_0_yyyz_1[i] * wa_y[i];

        g_xxxyz_0_yyzz_0[i] = 2.0 * g_xxxz_0_yzz_1[i] * fi_acd_0 + g_xxxz_0_yyzz_1[i] * wa_y[i];

        g_xxxyz_0_yzzz_0[i] = g_xxxz_0_zzz_1[i] * fi_acd_0 + g_xxxz_0_yzzz_1[i] * wa_y[i];

        g_xxxyz_0_zzzz_0[i] = g_xxxz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 75-90 components of targeted buffer : HSG

    auto g_xxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 75);

    auto g_xxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 76);

    auto g_xxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 77);

    auto g_xxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 78);

    auto g_xxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 79);

    auto g_xxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 80);

    auto g_xxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 81);

    auto g_xxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 82);

    auto g_xxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 83);

    auto g_xxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 84);

    auto g_xxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 85);

    auto g_xxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 86);

    auto g_xxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 87);

    auto g_xxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 88);

    auto g_xxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 89);

    #pragma omp simd aligned(g_xxx_0_xxxx_0, g_xxx_0_xxxx_1, g_xxx_0_xxxy_0, g_xxx_0_xxxy_1, g_xxx_0_xxyy_0, g_xxx_0_xxyy_1, g_xxx_0_xyyy_0, g_xxx_0_xyyy_1, g_xxxz_0_xxxx_1, g_xxxz_0_xxxy_1, g_xxxz_0_xxyy_1, g_xxxz_0_xyyy_1, g_xxxzz_0_xxxx_0, g_xxxzz_0_xxxy_0, g_xxxzz_0_xxxz_0, g_xxxzz_0_xxyy_0, g_xxxzz_0_xxyz_0, g_xxxzz_0_xxzz_0, g_xxxzz_0_xyyy_0, g_xxxzz_0_xyyz_0, g_xxxzz_0_xyzz_0, g_xxxzz_0_xzzz_0, g_xxxzz_0_yyyy_0, g_xxxzz_0_yyyz_0, g_xxxzz_0_yyzz_0, g_xxxzz_0_yzzz_0, g_xxxzz_0_zzzz_0, g_xxzz_0_xxxz_1, g_xxzz_0_xxyz_1, g_xxzz_0_xxz_1, g_xxzz_0_xxzz_1, g_xxzz_0_xyyz_1, g_xxzz_0_xyz_1, g_xxzz_0_xyzz_1, g_xxzz_0_xzz_1, g_xxzz_0_xzzz_1, g_xxzz_0_yyyy_1, g_xxzz_0_yyyz_1, g_xxzz_0_yyz_1, g_xxzz_0_yyzz_1, g_xxzz_0_yzz_1, g_xxzz_0_yzzz_1, g_xxzz_0_zzz_1, g_xxzz_0_zzzz_1, g_xzz_0_xxxz_0, g_xzz_0_xxxz_1, g_xzz_0_xxyz_0, g_xzz_0_xxyz_1, g_xzz_0_xxzz_0, g_xzz_0_xxzz_1, g_xzz_0_xyyz_0, g_xzz_0_xyyz_1, g_xzz_0_xyzz_0, g_xzz_0_xyzz_1, g_xzz_0_xzzz_0, g_xzz_0_xzzz_1, g_xzz_0_yyyy_0, g_xzz_0_yyyy_1, g_xzz_0_yyyz_0, g_xzz_0_yyyz_1, g_xzz_0_yyzz_0, g_xzz_0_yyzz_1, g_xzz_0_yzzz_0, g_xzz_0_yzzz_1, g_xzz_0_zzzz_0, g_xzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xxxx_0[i] = g_xxx_0_xxxx_0[i] * fbe_0 - g_xxx_0_xxxx_1[i] * fz_be_0 + g_xxxz_0_xxxx_1[i] * wa_z[i];

        g_xxxzz_0_xxxy_0[i] = g_xxx_0_xxxy_0[i] * fbe_0 - g_xxx_0_xxxy_1[i] * fz_be_0 + g_xxxz_0_xxxy_1[i] * wa_z[i];

        g_xxxzz_0_xxxz_0[i] = 2.0 * g_xzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxz_1[i] * fi_acd_0 + g_xxzz_0_xxxz_1[i] * wa_x[i];

        g_xxxzz_0_xxyy_0[i] = g_xxx_0_xxyy_0[i] * fbe_0 - g_xxx_0_xxyy_1[i] * fz_be_0 + g_xxxz_0_xxyy_1[i] * wa_z[i];

        g_xxxzz_0_xxyz_0[i] = 2.0 * g_xzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyz_1[i] * fi_acd_0 + g_xxzz_0_xxyz_1[i] * wa_x[i];

        g_xxxzz_0_xxzz_0[i] = 2.0 * g_xzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xzz_1[i] * fi_acd_0 + g_xxzz_0_xxzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyy_0[i] = g_xxx_0_xyyy_0[i] * fbe_0 - g_xxx_0_xyyy_1[i] * fz_be_0 + g_xxxz_0_xyyy_1[i] * wa_z[i];

        g_xxxzz_0_xyyz_0[i] = 2.0 * g_xzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyz_1[i] * fz_be_0 + g_xxzz_0_yyz_1[i] * fi_acd_0 + g_xxzz_0_xyyz_1[i] * wa_x[i];

        g_xxxzz_0_xyzz_0[i] = 2.0 * g_xzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyzz_1[i] * fz_be_0 + g_xxzz_0_yzz_1[i] * fi_acd_0 + g_xxzz_0_xyzz_1[i] * wa_x[i];

        g_xxxzz_0_xzzz_0[i] = 2.0 * g_xzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xzzz_1[i] * fz_be_0 + g_xxzz_0_zzz_1[i] * fi_acd_0 + g_xxzz_0_xzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyy_0[i] = 2.0 * g_xzz_0_yyyy_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyy_1[i] * fz_be_0 + g_xxzz_0_yyyy_1[i] * wa_x[i];

        g_xxxzz_0_yyyz_0[i] = 2.0 * g_xzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyz_1[i] * fz_be_0 + g_xxzz_0_yyyz_1[i] * wa_x[i];

        g_xxxzz_0_yyzz_0[i] = 2.0 * g_xzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyzz_1[i] * fz_be_0 + g_xxzz_0_yyzz_1[i] * wa_x[i];

        g_xxxzz_0_yzzz_0[i] = 2.0 * g_xzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yzzz_1[i] * fz_be_0 + g_xxzz_0_yzzz_1[i] * wa_x[i];

        g_xxxzz_0_zzzz_0[i] = 2.0 * g_xzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_zzzz_1[i] * fz_be_0 + g_xxzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 90-105 components of targeted buffer : HSG

    auto g_xxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 90);

    auto g_xxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 91);

    auto g_xxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 92);

    auto g_xxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 93);

    auto g_xxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 94);

    auto g_xxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 95);

    auto g_xxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 96);

    auto g_xxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 97);

    auto g_xxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 98);

    auto g_xxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 99);

    auto g_xxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 100);

    auto g_xxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 101);

    auto g_xxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 102);

    auto g_xxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 103);

    auto g_xxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 104);

    #pragma omp simd aligned(g_xxy_0_xxxx_0, g_xxy_0_xxxx_1, g_xxy_0_xxxz_0, g_xxy_0_xxxz_1, g_xxy_0_xxzz_0, g_xxy_0_xxzz_1, g_xxy_0_xzzz_0, g_xxy_0_xzzz_1, g_xxyy_0_xxxx_1, g_xxyy_0_xxxz_1, g_xxyy_0_xxzz_1, g_xxyy_0_xzzz_1, g_xxyyy_0_xxxx_0, g_xxyyy_0_xxxy_0, g_xxyyy_0_xxxz_0, g_xxyyy_0_xxyy_0, g_xxyyy_0_xxyz_0, g_xxyyy_0_xxzz_0, g_xxyyy_0_xyyy_0, g_xxyyy_0_xyyz_0, g_xxyyy_0_xyzz_0, g_xxyyy_0_xzzz_0, g_xxyyy_0_yyyy_0, g_xxyyy_0_yyyz_0, g_xxyyy_0_yyzz_0, g_xxyyy_0_yzzz_0, g_xxyyy_0_zzzz_0, g_xyyy_0_xxxy_1, g_xyyy_0_xxy_1, g_xyyy_0_xxyy_1, g_xyyy_0_xxyz_1, g_xyyy_0_xyy_1, g_xyyy_0_xyyy_1, g_xyyy_0_xyyz_1, g_xyyy_0_xyz_1, g_xyyy_0_xyzz_1, g_xyyy_0_yyy_1, g_xyyy_0_yyyy_1, g_xyyy_0_yyyz_1, g_xyyy_0_yyz_1, g_xyyy_0_yyzz_1, g_xyyy_0_yzz_1, g_xyyy_0_yzzz_1, g_xyyy_0_zzzz_1, g_yyy_0_xxxy_0, g_yyy_0_xxxy_1, g_yyy_0_xxyy_0, g_yyy_0_xxyy_1, g_yyy_0_xxyz_0, g_yyy_0_xxyz_1, g_yyy_0_xyyy_0, g_yyy_0_xyyy_1, g_yyy_0_xyyz_0, g_yyy_0_xyyz_1, g_yyy_0_xyzz_0, g_yyy_0_xyzz_1, g_yyy_0_yyyy_0, g_yyy_0_yyyy_1, g_yyy_0_yyyz_0, g_yyy_0_yyyz_1, g_yyy_0_yyzz_0, g_yyy_0_yyzz_1, g_yyy_0_yzzz_0, g_yyy_0_yzzz_1, g_yyy_0_zzzz_0, g_yyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xxxx_0[i] = 2.0 * g_xxy_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxx_1[i] * fz_be_0 + g_xxyy_0_xxxx_1[i] * wa_y[i];

        g_xxyyy_0_xxxy_0[i] = g_yyy_0_xxxy_0[i] * fbe_0 - g_yyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxy_1[i] * fi_acd_0 + g_xyyy_0_xxxy_1[i] * wa_x[i];

        g_xxyyy_0_xxxz_0[i] = 2.0 * g_xxy_0_xxxz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxz_1[i] * fz_be_0 + g_xxyy_0_xxxz_1[i] * wa_y[i];

        g_xxyyy_0_xxyy_0[i] = g_yyy_0_xxyy_0[i] * fbe_0 - g_yyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyy_1[i] * fi_acd_0 + g_xyyy_0_xxyy_1[i] * wa_x[i];

        g_xxyyy_0_xxyz_0[i] = g_yyy_0_xxyz_0[i] * fbe_0 - g_yyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyz_1[i] * fi_acd_0 + g_xyyy_0_xxyz_1[i] * wa_x[i];

        g_xxyyy_0_xxzz_0[i] = 2.0 * g_xxy_0_xxzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxzz_1[i] * fz_be_0 + g_xxyy_0_xxzz_1[i] * wa_y[i];

        g_xxyyy_0_xyyy_0[i] = g_yyy_0_xyyy_0[i] * fbe_0 - g_yyy_0_xyyy_1[i] * fz_be_0 + g_xyyy_0_yyy_1[i] * fi_acd_0 + g_xyyy_0_xyyy_1[i] * wa_x[i];

        g_xxyyy_0_xyyz_0[i] = g_yyy_0_xyyz_0[i] * fbe_0 - g_yyy_0_xyyz_1[i] * fz_be_0 + g_xyyy_0_yyz_1[i] * fi_acd_0 + g_xyyy_0_xyyz_1[i] * wa_x[i];

        g_xxyyy_0_xyzz_0[i] = g_yyy_0_xyzz_0[i] * fbe_0 - g_yyy_0_xyzz_1[i] * fz_be_0 + g_xyyy_0_yzz_1[i] * fi_acd_0 + g_xyyy_0_xyzz_1[i] * wa_x[i];

        g_xxyyy_0_xzzz_0[i] = 2.0 * g_xxy_0_xzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xzzz_1[i] * fz_be_0 + g_xxyy_0_xzzz_1[i] * wa_y[i];

        g_xxyyy_0_yyyy_0[i] = g_yyy_0_yyyy_0[i] * fbe_0 - g_yyy_0_yyyy_1[i] * fz_be_0 + g_xyyy_0_yyyy_1[i] * wa_x[i];

        g_xxyyy_0_yyyz_0[i] = g_yyy_0_yyyz_0[i] * fbe_0 - g_yyy_0_yyyz_1[i] * fz_be_0 + g_xyyy_0_yyyz_1[i] * wa_x[i];

        g_xxyyy_0_yyzz_0[i] = g_yyy_0_yyzz_0[i] * fbe_0 - g_yyy_0_yyzz_1[i] * fz_be_0 + g_xyyy_0_yyzz_1[i] * wa_x[i];

        g_xxyyy_0_yzzz_0[i] = g_yyy_0_yzzz_0[i] * fbe_0 - g_yyy_0_yzzz_1[i] * fz_be_0 + g_xyyy_0_yzzz_1[i] * wa_x[i];

        g_xxyyy_0_zzzz_0[i] = g_yyy_0_zzzz_0[i] * fbe_0 - g_yyy_0_zzzz_1[i] * fz_be_0 + g_xyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 105-120 components of targeted buffer : HSG

    auto g_xxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 105);

    auto g_xxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 106);

    auto g_xxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 107);

    auto g_xxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 108);

    auto g_xxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 109);

    auto g_xxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 110);

    auto g_xxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 111);

    auto g_xxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 112);

    auto g_xxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 113);

    auto g_xxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 114);

    auto g_xxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 115);

    auto g_xxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 116);

    auto g_xxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 117);

    auto g_xxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 118);

    auto g_xxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 119);

    #pragma omp simd aligned(g_xxyy_0_xxx_1, g_xxyy_0_xxxx_1, g_xxyy_0_xxxy_1, g_xxyy_0_xxxz_1, g_xxyy_0_xxy_1, g_xxyy_0_xxyy_1, g_xxyy_0_xxyz_1, g_xxyy_0_xxz_1, g_xxyy_0_xxzz_1, g_xxyy_0_xyy_1, g_xxyy_0_xyyy_1, g_xxyy_0_xyyz_1, g_xxyy_0_xyz_1, g_xxyy_0_xyzz_1, g_xxyy_0_xzz_1, g_xxyy_0_xzzz_1, g_xxyy_0_yyy_1, g_xxyy_0_yyyy_1, g_xxyy_0_yyyz_1, g_xxyy_0_yyz_1, g_xxyy_0_yyzz_1, g_xxyy_0_yzz_1, g_xxyy_0_yzzz_1, g_xxyy_0_zzz_1, g_xxyy_0_zzzz_1, g_xxyyz_0_xxxx_0, g_xxyyz_0_xxxy_0, g_xxyyz_0_xxxz_0, g_xxyyz_0_xxyy_0, g_xxyyz_0_xxyz_0, g_xxyyz_0_xxzz_0, g_xxyyz_0_xyyy_0, g_xxyyz_0_xyyz_0, g_xxyyz_0_xyzz_0, g_xxyyz_0_xzzz_0, g_xxyyz_0_yyyy_0, g_xxyyz_0_yyyz_0, g_xxyyz_0_yyzz_0, g_xxyyz_0_yzzz_0, g_xxyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xxxx_0[i] = g_xxyy_0_xxxx_1[i] * wa_z[i];

        g_xxyyz_0_xxxy_0[i] = g_xxyy_0_xxxy_1[i] * wa_z[i];

        g_xxyyz_0_xxxz_0[i] = g_xxyy_0_xxx_1[i] * fi_acd_0 + g_xxyy_0_xxxz_1[i] * wa_z[i];

        g_xxyyz_0_xxyy_0[i] = g_xxyy_0_xxyy_1[i] * wa_z[i];

        g_xxyyz_0_xxyz_0[i] = g_xxyy_0_xxy_1[i] * fi_acd_0 + g_xxyy_0_xxyz_1[i] * wa_z[i];

        g_xxyyz_0_xxzz_0[i] = 2.0 * g_xxyy_0_xxz_1[i] * fi_acd_0 + g_xxyy_0_xxzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyy_0[i] = g_xxyy_0_xyyy_1[i] * wa_z[i];

        g_xxyyz_0_xyyz_0[i] = g_xxyy_0_xyy_1[i] * fi_acd_0 + g_xxyy_0_xyyz_1[i] * wa_z[i];

        g_xxyyz_0_xyzz_0[i] = 2.0 * g_xxyy_0_xyz_1[i] * fi_acd_0 + g_xxyy_0_xyzz_1[i] * wa_z[i];

        g_xxyyz_0_xzzz_0[i] = 3.0 * g_xxyy_0_xzz_1[i] * fi_acd_0 + g_xxyy_0_xzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyy_0[i] = g_xxyy_0_yyyy_1[i] * wa_z[i];

        g_xxyyz_0_yyyz_0[i] = g_xxyy_0_yyy_1[i] * fi_acd_0 + g_xxyy_0_yyyz_1[i] * wa_z[i];

        g_xxyyz_0_yyzz_0[i] = 2.0 * g_xxyy_0_yyz_1[i] * fi_acd_0 + g_xxyy_0_yyzz_1[i] * wa_z[i];

        g_xxyyz_0_yzzz_0[i] = 3.0 * g_xxyy_0_yzz_1[i] * fi_acd_0 + g_xxyy_0_yzzz_1[i] * wa_z[i];

        g_xxyyz_0_zzzz_0[i] = 4.0 * g_xxyy_0_zzz_1[i] * fi_acd_0 + g_xxyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 120-135 components of targeted buffer : HSG

    auto g_xxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 120);

    auto g_xxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 121);

    auto g_xxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 122);

    auto g_xxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 123);

    auto g_xxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 124);

    auto g_xxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 125);

    auto g_xxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 126);

    auto g_xxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 127);

    auto g_xxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 128);

    auto g_xxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 129);

    auto g_xxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 130);

    auto g_xxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 131);

    auto g_xxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 132);

    auto g_xxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 133);

    auto g_xxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 134);

    #pragma omp simd aligned(g_xxyzz_0_xxxx_0, g_xxyzz_0_xxxy_0, g_xxyzz_0_xxxz_0, g_xxyzz_0_xxyy_0, g_xxyzz_0_xxyz_0, g_xxyzz_0_xxzz_0, g_xxyzz_0_xyyy_0, g_xxyzz_0_xyyz_0, g_xxyzz_0_xyzz_0, g_xxyzz_0_xzzz_0, g_xxyzz_0_yyyy_0, g_xxyzz_0_yyyz_0, g_xxyzz_0_yyzz_0, g_xxyzz_0_yzzz_0, g_xxyzz_0_zzzz_0, g_xxzz_0_xxx_1, g_xxzz_0_xxxx_1, g_xxzz_0_xxxy_1, g_xxzz_0_xxxz_1, g_xxzz_0_xxy_1, g_xxzz_0_xxyy_1, g_xxzz_0_xxyz_1, g_xxzz_0_xxz_1, g_xxzz_0_xxzz_1, g_xxzz_0_xyy_1, g_xxzz_0_xyyy_1, g_xxzz_0_xyyz_1, g_xxzz_0_xyz_1, g_xxzz_0_xyzz_1, g_xxzz_0_xzz_1, g_xxzz_0_xzzz_1, g_xxzz_0_yyy_1, g_xxzz_0_yyyy_1, g_xxzz_0_yyyz_1, g_xxzz_0_yyz_1, g_xxzz_0_yyzz_1, g_xxzz_0_yzz_1, g_xxzz_0_yzzz_1, g_xxzz_0_zzz_1, g_xxzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xxxx_0[i] = g_xxzz_0_xxxx_1[i] * wa_y[i];

        g_xxyzz_0_xxxy_0[i] = g_xxzz_0_xxx_1[i] * fi_acd_0 + g_xxzz_0_xxxy_1[i] * wa_y[i];

        g_xxyzz_0_xxxz_0[i] = g_xxzz_0_xxxz_1[i] * wa_y[i];

        g_xxyzz_0_xxyy_0[i] = 2.0 * g_xxzz_0_xxy_1[i] * fi_acd_0 + g_xxzz_0_xxyy_1[i] * wa_y[i];

        g_xxyzz_0_xxyz_0[i] = g_xxzz_0_xxz_1[i] * fi_acd_0 + g_xxzz_0_xxyz_1[i] * wa_y[i];

        g_xxyzz_0_xxzz_0[i] = g_xxzz_0_xxzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyy_0[i] = 3.0 * g_xxzz_0_xyy_1[i] * fi_acd_0 + g_xxzz_0_xyyy_1[i] * wa_y[i];

        g_xxyzz_0_xyyz_0[i] = 2.0 * g_xxzz_0_xyz_1[i] * fi_acd_0 + g_xxzz_0_xyyz_1[i] * wa_y[i];

        g_xxyzz_0_xyzz_0[i] = g_xxzz_0_xzz_1[i] * fi_acd_0 + g_xxzz_0_xyzz_1[i] * wa_y[i];

        g_xxyzz_0_xzzz_0[i] = g_xxzz_0_xzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyy_0[i] = 4.0 * g_xxzz_0_yyy_1[i] * fi_acd_0 + g_xxzz_0_yyyy_1[i] * wa_y[i];

        g_xxyzz_0_yyyz_0[i] = 3.0 * g_xxzz_0_yyz_1[i] * fi_acd_0 + g_xxzz_0_yyyz_1[i] * wa_y[i];

        g_xxyzz_0_yyzz_0[i] = 2.0 * g_xxzz_0_yzz_1[i] * fi_acd_0 + g_xxzz_0_yyzz_1[i] * wa_y[i];

        g_xxyzz_0_yzzz_0[i] = g_xxzz_0_zzz_1[i] * fi_acd_0 + g_xxzz_0_yzzz_1[i] * wa_y[i];

        g_xxyzz_0_zzzz_0[i] = g_xxzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 135-150 components of targeted buffer : HSG

    auto g_xxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 135);

    auto g_xxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 136);

    auto g_xxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 137);

    auto g_xxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 138);

    auto g_xxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 139);

    auto g_xxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 140);

    auto g_xxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 141);

    auto g_xxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 142);

    auto g_xxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 143);

    auto g_xxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 144);

    auto g_xxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 145);

    auto g_xxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 146);

    auto g_xxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 147);

    auto g_xxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 148);

    auto g_xxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 149);

    #pragma omp simd aligned(g_xxz_0_xxxx_0, g_xxz_0_xxxx_1, g_xxz_0_xxxy_0, g_xxz_0_xxxy_1, g_xxz_0_xxyy_0, g_xxz_0_xxyy_1, g_xxz_0_xyyy_0, g_xxz_0_xyyy_1, g_xxzz_0_xxxx_1, g_xxzz_0_xxxy_1, g_xxzz_0_xxyy_1, g_xxzz_0_xyyy_1, g_xxzzz_0_xxxx_0, g_xxzzz_0_xxxy_0, g_xxzzz_0_xxxz_0, g_xxzzz_0_xxyy_0, g_xxzzz_0_xxyz_0, g_xxzzz_0_xxzz_0, g_xxzzz_0_xyyy_0, g_xxzzz_0_xyyz_0, g_xxzzz_0_xyzz_0, g_xxzzz_0_xzzz_0, g_xxzzz_0_yyyy_0, g_xxzzz_0_yyyz_0, g_xxzzz_0_yyzz_0, g_xxzzz_0_yzzz_0, g_xxzzz_0_zzzz_0, g_xzzz_0_xxxz_1, g_xzzz_0_xxyz_1, g_xzzz_0_xxz_1, g_xzzz_0_xxzz_1, g_xzzz_0_xyyz_1, g_xzzz_0_xyz_1, g_xzzz_0_xyzz_1, g_xzzz_0_xzz_1, g_xzzz_0_xzzz_1, g_xzzz_0_yyyy_1, g_xzzz_0_yyyz_1, g_xzzz_0_yyz_1, g_xzzz_0_yyzz_1, g_xzzz_0_yzz_1, g_xzzz_0_yzzz_1, g_xzzz_0_zzz_1, g_xzzz_0_zzzz_1, g_zzz_0_xxxz_0, g_zzz_0_xxxz_1, g_zzz_0_xxyz_0, g_zzz_0_xxyz_1, g_zzz_0_xxzz_0, g_zzz_0_xxzz_1, g_zzz_0_xyyz_0, g_zzz_0_xyyz_1, g_zzz_0_xyzz_0, g_zzz_0_xyzz_1, g_zzz_0_xzzz_0, g_zzz_0_xzzz_1, g_zzz_0_yyyy_0, g_zzz_0_yyyy_1, g_zzz_0_yyyz_0, g_zzz_0_yyyz_1, g_zzz_0_yyzz_0, g_zzz_0_yyzz_1, g_zzz_0_yzzz_0, g_zzz_0_yzzz_1, g_zzz_0_zzzz_0, g_zzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xxxx_0[i] = 2.0 * g_xxz_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxx_1[i] * fz_be_0 + g_xxzz_0_xxxx_1[i] * wa_z[i];

        g_xxzzz_0_xxxy_0[i] = 2.0 * g_xxz_0_xxxy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxy_1[i] * fz_be_0 + g_xxzz_0_xxxy_1[i] * wa_z[i];

        g_xxzzz_0_xxxz_0[i] = g_zzz_0_xxxz_0[i] * fbe_0 - g_zzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxz_1[i] * fi_acd_0 + g_xzzz_0_xxxz_1[i] * wa_x[i];

        g_xxzzz_0_xxyy_0[i] = 2.0 * g_xxz_0_xxyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxyy_1[i] * fz_be_0 + g_xxzz_0_xxyy_1[i] * wa_z[i];

        g_xxzzz_0_xxyz_0[i] = g_zzz_0_xxyz_0[i] * fbe_0 - g_zzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyz_1[i] * fi_acd_0 + g_xzzz_0_xxyz_1[i] * wa_x[i];

        g_xxzzz_0_xxzz_0[i] = g_zzz_0_xxzz_0[i] * fbe_0 - g_zzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xzz_1[i] * fi_acd_0 + g_xzzz_0_xxzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyy_0[i] = 2.0 * g_xxz_0_xyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xyyy_1[i] * fz_be_0 + g_xxzz_0_xyyy_1[i] * wa_z[i];

        g_xxzzz_0_xyyz_0[i] = g_zzz_0_xyyz_0[i] * fbe_0 - g_zzz_0_xyyz_1[i] * fz_be_0 + g_xzzz_0_yyz_1[i] * fi_acd_0 + g_xzzz_0_xyyz_1[i] * wa_x[i];

        g_xxzzz_0_xyzz_0[i] = g_zzz_0_xyzz_0[i] * fbe_0 - g_zzz_0_xyzz_1[i] * fz_be_0 + g_xzzz_0_yzz_1[i] * fi_acd_0 + g_xzzz_0_xyzz_1[i] * wa_x[i];

        g_xxzzz_0_xzzz_0[i] = g_zzz_0_xzzz_0[i] * fbe_0 - g_zzz_0_xzzz_1[i] * fz_be_0 + g_xzzz_0_zzz_1[i] * fi_acd_0 + g_xzzz_0_xzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyy_0[i] = g_zzz_0_yyyy_0[i] * fbe_0 - g_zzz_0_yyyy_1[i] * fz_be_0 + g_xzzz_0_yyyy_1[i] * wa_x[i];

        g_xxzzz_0_yyyz_0[i] = g_zzz_0_yyyz_0[i] * fbe_0 - g_zzz_0_yyyz_1[i] * fz_be_0 + g_xzzz_0_yyyz_1[i] * wa_x[i];

        g_xxzzz_0_yyzz_0[i] = g_zzz_0_yyzz_0[i] * fbe_0 - g_zzz_0_yyzz_1[i] * fz_be_0 + g_xzzz_0_yyzz_1[i] * wa_x[i];

        g_xxzzz_0_yzzz_0[i] = g_zzz_0_yzzz_0[i] * fbe_0 - g_zzz_0_yzzz_1[i] * fz_be_0 + g_xzzz_0_yzzz_1[i] * wa_x[i];

        g_xxzzz_0_zzzz_0[i] = g_zzz_0_zzzz_0[i] * fbe_0 - g_zzz_0_zzzz_1[i] * fz_be_0 + g_xzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 150-165 components of targeted buffer : HSG

    auto g_xyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 150);

    auto g_xyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 151);

    auto g_xyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 152);

    auto g_xyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 153);

    auto g_xyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 154);

    auto g_xyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 155);

    auto g_xyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 156);

    auto g_xyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 157);

    auto g_xyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 158);

    auto g_xyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 159);

    auto g_xyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 160);

    auto g_xyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 161);

    auto g_xyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 162);

    auto g_xyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 163);

    auto g_xyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 164);

    #pragma omp simd aligned(g_xyyyy_0_xxxx_0, g_xyyyy_0_xxxy_0, g_xyyyy_0_xxxz_0, g_xyyyy_0_xxyy_0, g_xyyyy_0_xxyz_0, g_xyyyy_0_xxzz_0, g_xyyyy_0_xyyy_0, g_xyyyy_0_xyyz_0, g_xyyyy_0_xyzz_0, g_xyyyy_0_xzzz_0, g_xyyyy_0_yyyy_0, g_xyyyy_0_yyyz_0, g_xyyyy_0_yyzz_0, g_xyyyy_0_yzzz_0, g_xyyyy_0_zzzz_0, g_yyyy_0_xxx_1, g_yyyy_0_xxxx_1, g_yyyy_0_xxxy_1, g_yyyy_0_xxxz_1, g_yyyy_0_xxy_1, g_yyyy_0_xxyy_1, g_yyyy_0_xxyz_1, g_yyyy_0_xxz_1, g_yyyy_0_xxzz_1, g_yyyy_0_xyy_1, g_yyyy_0_xyyy_1, g_yyyy_0_xyyz_1, g_yyyy_0_xyz_1, g_yyyy_0_xyzz_1, g_yyyy_0_xzz_1, g_yyyy_0_xzzz_1, g_yyyy_0_yyy_1, g_yyyy_0_yyyy_1, g_yyyy_0_yyyz_1, g_yyyy_0_yyz_1, g_yyyy_0_yyzz_1, g_yyyy_0_yzz_1, g_yyyy_0_yzzz_1, g_yyyy_0_zzz_1, g_yyyy_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xxxx_0[i] = 4.0 * g_yyyy_0_xxx_1[i] * fi_acd_0 + g_yyyy_0_xxxx_1[i] * wa_x[i];

        g_xyyyy_0_xxxy_0[i] = 3.0 * g_yyyy_0_xxy_1[i] * fi_acd_0 + g_yyyy_0_xxxy_1[i] * wa_x[i];

        g_xyyyy_0_xxxz_0[i] = 3.0 * g_yyyy_0_xxz_1[i] * fi_acd_0 + g_yyyy_0_xxxz_1[i] * wa_x[i];

        g_xyyyy_0_xxyy_0[i] = 2.0 * g_yyyy_0_xyy_1[i] * fi_acd_0 + g_yyyy_0_xxyy_1[i] * wa_x[i];

        g_xyyyy_0_xxyz_0[i] = 2.0 * g_yyyy_0_xyz_1[i] * fi_acd_0 + g_yyyy_0_xxyz_1[i] * wa_x[i];

        g_xyyyy_0_xxzz_0[i] = 2.0 * g_yyyy_0_xzz_1[i] * fi_acd_0 + g_yyyy_0_xxzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyy_0[i] = g_yyyy_0_yyy_1[i] * fi_acd_0 + g_yyyy_0_xyyy_1[i] * wa_x[i];

        g_xyyyy_0_xyyz_0[i] = g_yyyy_0_yyz_1[i] * fi_acd_0 + g_yyyy_0_xyyz_1[i] * wa_x[i];

        g_xyyyy_0_xyzz_0[i] = g_yyyy_0_yzz_1[i] * fi_acd_0 + g_yyyy_0_xyzz_1[i] * wa_x[i];

        g_xyyyy_0_xzzz_0[i] = g_yyyy_0_zzz_1[i] * fi_acd_0 + g_yyyy_0_xzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyy_0[i] = g_yyyy_0_yyyy_1[i] * wa_x[i];

        g_xyyyy_0_yyyz_0[i] = g_yyyy_0_yyyz_1[i] * wa_x[i];

        g_xyyyy_0_yyzz_0[i] = g_yyyy_0_yyzz_1[i] * wa_x[i];

        g_xyyyy_0_yzzz_0[i] = g_yyyy_0_yzzz_1[i] * wa_x[i];

        g_xyyyy_0_zzzz_0[i] = g_yyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 165-180 components of targeted buffer : HSG

    auto g_xyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 165);

    auto g_xyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 166);

    auto g_xyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 167);

    auto g_xyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 168);

    auto g_xyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 169);

    auto g_xyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 170);

    auto g_xyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 171);

    auto g_xyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 172);

    auto g_xyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 173);

    auto g_xyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 174);

    auto g_xyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 175);

    auto g_xyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 176);

    auto g_xyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 177);

    auto g_xyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 178);

    auto g_xyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 179);

    #pragma omp simd aligned(g_xyyy_0_xxxx_1, g_xyyy_0_xxxy_1, g_xyyy_0_xxyy_1, g_xyyy_0_xyyy_1, g_xyyyz_0_xxxx_0, g_xyyyz_0_xxxy_0, g_xyyyz_0_xxxz_0, g_xyyyz_0_xxyy_0, g_xyyyz_0_xxyz_0, g_xyyyz_0_xxzz_0, g_xyyyz_0_xyyy_0, g_xyyyz_0_xyyz_0, g_xyyyz_0_xyzz_0, g_xyyyz_0_xzzz_0, g_xyyyz_0_yyyy_0, g_xyyyz_0_yyyz_0, g_xyyyz_0_yyzz_0, g_xyyyz_0_yzzz_0, g_xyyyz_0_zzzz_0, g_yyyz_0_xxxz_1, g_yyyz_0_xxyz_1, g_yyyz_0_xxz_1, g_yyyz_0_xxzz_1, g_yyyz_0_xyyz_1, g_yyyz_0_xyz_1, g_yyyz_0_xyzz_1, g_yyyz_0_xzz_1, g_yyyz_0_xzzz_1, g_yyyz_0_yyyy_1, g_yyyz_0_yyyz_1, g_yyyz_0_yyz_1, g_yyyz_0_yyzz_1, g_yyyz_0_yzz_1, g_yyyz_0_yzzz_1, g_yyyz_0_zzz_1, g_yyyz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xxxx_0[i] = g_xyyy_0_xxxx_1[i] * wa_z[i];

        g_xyyyz_0_xxxy_0[i] = g_xyyy_0_xxxy_1[i] * wa_z[i];

        g_xyyyz_0_xxxz_0[i] = 3.0 * g_yyyz_0_xxz_1[i] * fi_acd_0 + g_yyyz_0_xxxz_1[i] * wa_x[i];

        g_xyyyz_0_xxyy_0[i] = g_xyyy_0_xxyy_1[i] * wa_z[i];

        g_xyyyz_0_xxyz_0[i] = 2.0 * g_yyyz_0_xyz_1[i] * fi_acd_0 + g_yyyz_0_xxyz_1[i] * wa_x[i];

        g_xyyyz_0_xxzz_0[i] = 2.0 * g_yyyz_0_xzz_1[i] * fi_acd_0 + g_yyyz_0_xxzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyy_0[i] = g_xyyy_0_xyyy_1[i] * wa_z[i];

        g_xyyyz_0_xyyz_0[i] = g_yyyz_0_yyz_1[i] * fi_acd_0 + g_yyyz_0_xyyz_1[i] * wa_x[i];

        g_xyyyz_0_xyzz_0[i] = g_yyyz_0_yzz_1[i] * fi_acd_0 + g_yyyz_0_xyzz_1[i] * wa_x[i];

        g_xyyyz_0_xzzz_0[i] = g_yyyz_0_zzz_1[i] * fi_acd_0 + g_yyyz_0_xzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyy_0[i] = g_yyyz_0_yyyy_1[i] * wa_x[i];

        g_xyyyz_0_yyyz_0[i] = g_yyyz_0_yyyz_1[i] * wa_x[i];

        g_xyyyz_0_yyzz_0[i] = g_yyyz_0_yyzz_1[i] * wa_x[i];

        g_xyyyz_0_yzzz_0[i] = g_yyyz_0_yzzz_1[i] * wa_x[i];

        g_xyyyz_0_zzzz_0[i] = g_yyyz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 180-195 components of targeted buffer : HSG

    auto g_xyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 180);

    auto g_xyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 181);

    auto g_xyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 182);

    auto g_xyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 183);

    auto g_xyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 184);

    auto g_xyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 185);

    auto g_xyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 186);

    auto g_xyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 187);

    auto g_xyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 188);

    auto g_xyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 189);

    auto g_xyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 190);

    auto g_xyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 191);

    auto g_xyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 192);

    auto g_xyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 193);

    auto g_xyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 194);

    #pragma omp simd aligned(g_xyyzz_0_xxxx_0, g_xyyzz_0_xxxy_0, g_xyyzz_0_xxxz_0, g_xyyzz_0_xxyy_0, g_xyyzz_0_xxyz_0, g_xyyzz_0_xxzz_0, g_xyyzz_0_xyyy_0, g_xyyzz_0_xyyz_0, g_xyyzz_0_xyzz_0, g_xyyzz_0_xzzz_0, g_xyyzz_0_yyyy_0, g_xyyzz_0_yyyz_0, g_xyyzz_0_yyzz_0, g_xyyzz_0_yzzz_0, g_xyyzz_0_zzzz_0, g_yyzz_0_xxx_1, g_yyzz_0_xxxx_1, g_yyzz_0_xxxy_1, g_yyzz_0_xxxz_1, g_yyzz_0_xxy_1, g_yyzz_0_xxyy_1, g_yyzz_0_xxyz_1, g_yyzz_0_xxz_1, g_yyzz_0_xxzz_1, g_yyzz_0_xyy_1, g_yyzz_0_xyyy_1, g_yyzz_0_xyyz_1, g_yyzz_0_xyz_1, g_yyzz_0_xyzz_1, g_yyzz_0_xzz_1, g_yyzz_0_xzzz_1, g_yyzz_0_yyy_1, g_yyzz_0_yyyy_1, g_yyzz_0_yyyz_1, g_yyzz_0_yyz_1, g_yyzz_0_yyzz_1, g_yyzz_0_yzz_1, g_yyzz_0_yzzz_1, g_yyzz_0_zzz_1, g_yyzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xxxx_0[i] = 4.0 * g_yyzz_0_xxx_1[i] * fi_acd_0 + g_yyzz_0_xxxx_1[i] * wa_x[i];

        g_xyyzz_0_xxxy_0[i] = 3.0 * g_yyzz_0_xxy_1[i] * fi_acd_0 + g_yyzz_0_xxxy_1[i] * wa_x[i];

        g_xyyzz_0_xxxz_0[i] = 3.0 * g_yyzz_0_xxz_1[i] * fi_acd_0 + g_yyzz_0_xxxz_1[i] * wa_x[i];

        g_xyyzz_0_xxyy_0[i] = 2.0 * g_yyzz_0_xyy_1[i] * fi_acd_0 + g_yyzz_0_xxyy_1[i] * wa_x[i];

        g_xyyzz_0_xxyz_0[i] = 2.0 * g_yyzz_0_xyz_1[i] * fi_acd_0 + g_yyzz_0_xxyz_1[i] * wa_x[i];

        g_xyyzz_0_xxzz_0[i] = 2.0 * g_yyzz_0_xzz_1[i] * fi_acd_0 + g_yyzz_0_xxzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyy_0[i] = g_yyzz_0_yyy_1[i] * fi_acd_0 + g_yyzz_0_xyyy_1[i] * wa_x[i];

        g_xyyzz_0_xyyz_0[i] = g_yyzz_0_yyz_1[i] * fi_acd_0 + g_yyzz_0_xyyz_1[i] * wa_x[i];

        g_xyyzz_0_xyzz_0[i] = g_yyzz_0_yzz_1[i] * fi_acd_0 + g_yyzz_0_xyzz_1[i] * wa_x[i];

        g_xyyzz_0_xzzz_0[i] = g_yyzz_0_zzz_1[i] * fi_acd_0 + g_yyzz_0_xzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyy_0[i] = g_yyzz_0_yyyy_1[i] * wa_x[i];

        g_xyyzz_0_yyyz_0[i] = g_yyzz_0_yyyz_1[i] * wa_x[i];

        g_xyyzz_0_yyzz_0[i] = g_yyzz_0_yyzz_1[i] * wa_x[i];

        g_xyyzz_0_yzzz_0[i] = g_yyzz_0_yzzz_1[i] * wa_x[i];

        g_xyyzz_0_zzzz_0[i] = g_yyzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 195-210 components of targeted buffer : HSG

    auto g_xyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 195);

    auto g_xyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 196);

    auto g_xyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 197);

    auto g_xyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 198);

    auto g_xyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 199);

    auto g_xyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 200);

    auto g_xyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 201);

    auto g_xyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 202);

    auto g_xyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 203);

    auto g_xyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 204);

    auto g_xyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 205);

    auto g_xyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 206);

    auto g_xyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 207);

    auto g_xyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 208);

    auto g_xyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 209);

    #pragma omp simd aligned(g_xyzzz_0_xxxx_0, g_xyzzz_0_xxxy_0, g_xyzzz_0_xxxz_0, g_xyzzz_0_xxyy_0, g_xyzzz_0_xxyz_0, g_xyzzz_0_xxzz_0, g_xyzzz_0_xyyy_0, g_xyzzz_0_xyyz_0, g_xyzzz_0_xyzz_0, g_xyzzz_0_xzzz_0, g_xyzzz_0_yyyy_0, g_xyzzz_0_yyyz_0, g_xyzzz_0_yyzz_0, g_xyzzz_0_yzzz_0, g_xyzzz_0_zzzz_0, g_xzzz_0_xxxx_1, g_xzzz_0_xxxz_1, g_xzzz_0_xxzz_1, g_xzzz_0_xzzz_1, g_yzzz_0_xxxy_1, g_yzzz_0_xxy_1, g_yzzz_0_xxyy_1, g_yzzz_0_xxyz_1, g_yzzz_0_xyy_1, g_yzzz_0_xyyy_1, g_yzzz_0_xyyz_1, g_yzzz_0_xyz_1, g_yzzz_0_xyzz_1, g_yzzz_0_yyy_1, g_yzzz_0_yyyy_1, g_yzzz_0_yyyz_1, g_yzzz_0_yyz_1, g_yzzz_0_yyzz_1, g_yzzz_0_yzz_1, g_yzzz_0_yzzz_1, g_yzzz_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xxxx_0[i] = g_xzzz_0_xxxx_1[i] * wa_y[i];

        g_xyzzz_0_xxxy_0[i] = 3.0 * g_yzzz_0_xxy_1[i] * fi_acd_0 + g_yzzz_0_xxxy_1[i] * wa_x[i];

        g_xyzzz_0_xxxz_0[i] = g_xzzz_0_xxxz_1[i] * wa_y[i];

        g_xyzzz_0_xxyy_0[i] = 2.0 * g_yzzz_0_xyy_1[i] * fi_acd_0 + g_yzzz_0_xxyy_1[i] * wa_x[i];

        g_xyzzz_0_xxyz_0[i] = 2.0 * g_yzzz_0_xyz_1[i] * fi_acd_0 + g_yzzz_0_xxyz_1[i] * wa_x[i];

        g_xyzzz_0_xxzz_0[i] = g_xzzz_0_xxzz_1[i] * wa_y[i];

        g_xyzzz_0_xyyy_0[i] = g_yzzz_0_yyy_1[i] * fi_acd_0 + g_yzzz_0_xyyy_1[i] * wa_x[i];

        g_xyzzz_0_xyyz_0[i] = g_yzzz_0_yyz_1[i] * fi_acd_0 + g_yzzz_0_xyyz_1[i] * wa_x[i];

        g_xyzzz_0_xyzz_0[i] = g_yzzz_0_yzz_1[i] * fi_acd_0 + g_yzzz_0_xyzz_1[i] * wa_x[i];

        g_xyzzz_0_xzzz_0[i] = g_xzzz_0_xzzz_1[i] * wa_y[i];

        g_xyzzz_0_yyyy_0[i] = g_yzzz_0_yyyy_1[i] * wa_x[i];

        g_xyzzz_0_yyyz_0[i] = g_yzzz_0_yyyz_1[i] * wa_x[i];

        g_xyzzz_0_yyzz_0[i] = g_yzzz_0_yyzz_1[i] * wa_x[i];

        g_xyzzz_0_yzzz_0[i] = g_yzzz_0_yzzz_1[i] * wa_x[i];

        g_xyzzz_0_zzzz_0[i] = g_yzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 210-225 components of targeted buffer : HSG

    auto g_xzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 210);

    auto g_xzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 211);

    auto g_xzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 212);

    auto g_xzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 213);

    auto g_xzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 214);

    auto g_xzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 215);

    auto g_xzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 216);

    auto g_xzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 217);

    auto g_xzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 218);

    auto g_xzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 219);

    auto g_xzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 220);

    auto g_xzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 221);

    auto g_xzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 222);

    auto g_xzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 223);

    auto g_xzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 224);

    #pragma omp simd aligned(g_xzzzz_0_xxxx_0, g_xzzzz_0_xxxy_0, g_xzzzz_0_xxxz_0, g_xzzzz_0_xxyy_0, g_xzzzz_0_xxyz_0, g_xzzzz_0_xxzz_0, g_xzzzz_0_xyyy_0, g_xzzzz_0_xyyz_0, g_xzzzz_0_xyzz_0, g_xzzzz_0_xzzz_0, g_xzzzz_0_yyyy_0, g_xzzzz_0_yyyz_0, g_xzzzz_0_yyzz_0, g_xzzzz_0_yzzz_0, g_xzzzz_0_zzzz_0, g_zzzz_0_xxx_1, g_zzzz_0_xxxx_1, g_zzzz_0_xxxy_1, g_zzzz_0_xxxz_1, g_zzzz_0_xxy_1, g_zzzz_0_xxyy_1, g_zzzz_0_xxyz_1, g_zzzz_0_xxz_1, g_zzzz_0_xxzz_1, g_zzzz_0_xyy_1, g_zzzz_0_xyyy_1, g_zzzz_0_xyyz_1, g_zzzz_0_xyz_1, g_zzzz_0_xyzz_1, g_zzzz_0_xzz_1, g_zzzz_0_xzzz_1, g_zzzz_0_yyy_1, g_zzzz_0_yyyy_1, g_zzzz_0_yyyz_1, g_zzzz_0_yyz_1, g_zzzz_0_yyzz_1, g_zzzz_0_yzz_1, g_zzzz_0_yzzz_1, g_zzzz_0_zzz_1, g_zzzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xxxx_0[i] = 4.0 * g_zzzz_0_xxx_1[i] * fi_acd_0 + g_zzzz_0_xxxx_1[i] * wa_x[i];

        g_xzzzz_0_xxxy_0[i] = 3.0 * g_zzzz_0_xxy_1[i] * fi_acd_0 + g_zzzz_0_xxxy_1[i] * wa_x[i];

        g_xzzzz_0_xxxz_0[i] = 3.0 * g_zzzz_0_xxz_1[i] * fi_acd_0 + g_zzzz_0_xxxz_1[i] * wa_x[i];

        g_xzzzz_0_xxyy_0[i] = 2.0 * g_zzzz_0_xyy_1[i] * fi_acd_0 + g_zzzz_0_xxyy_1[i] * wa_x[i];

        g_xzzzz_0_xxyz_0[i] = 2.0 * g_zzzz_0_xyz_1[i] * fi_acd_0 + g_zzzz_0_xxyz_1[i] * wa_x[i];

        g_xzzzz_0_xxzz_0[i] = 2.0 * g_zzzz_0_xzz_1[i] * fi_acd_0 + g_zzzz_0_xxzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyy_0[i] = g_zzzz_0_yyy_1[i] * fi_acd_0 + g_zzzz_0_xyyy_1[i] * wa_x[i];

        g_xzzzz_0_xyyz_0[i] = g_zzzz_0_yyz_1[i] * fi_acd_0 + g_zzzz_0_xyyz_1[i] * wa_x[i];

        g_xzzzz_0_xyzz_0[i] = g_zzzz_0_yzz_1[i] * fi_acd_0 + g_zzzz_0_xyzz_1[i] * wa_x[i];

        g_xzzzz_0_xzzz_0[i] = g_zzzz_0_zzz_1[i] * fi_acd_0 + g_zzzz_0_xzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyy_0[i] = g_zzzz_0_yyyy_1[i] * wa_x[i];

        g_xzzzz_0_yyyz_0[i] = g_zzzz_0_yyyz_1[i] * wa_x[i];

        g_xzzzz_0_yyzz_0[i] = g_zzzz_0_yyzz_1[i] * wa_x[i];

        g_xzzzz_0_yzzz_0[i] = g_zzzz_0_yzzz_1[i] * wa_x[i];

        g_xzzzz_0_zzzz_0[i] = g_zzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 225-240 components of targeted buffer : HSG

    auto g_yyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 225);

    auto g_yyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 226);

    auto g_yyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 227);

    auto g_yyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 228);

    auto g_yyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 229);

    auto g_yyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 230);

    auto g_yyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 231);

    auto g_yyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 232);

    auto g_yyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 233);

    auto g_yyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 234);

    auto g_yyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 235);

    auto g_yyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 236);

    auto g_yyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 237);

    auto g_yyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 238);

    auto g_yyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 239);

    #pragma omp simd aligned(g_yyy_0_xxxx_0, g_yyy_0_xxxx_1, g_yyy_0_xxxy_0, g_yyy_0_xxxy_1, g_yyy_0_xxxz_0, g_yyy_0_xxxz_1, g_yyy_0_xxyy_0, g_yyy_0_xxyy_1, g_yyy_0_xxyz_0, g_yyy_0_xxyz_1, g_yyy_0_xxzz_0, g_yyy_0_xxzz_1, g_yyy_0_xyyy_0, g_yyy_0_xyyy_1, g_yyy_0_xyyz_0, g_yyy_0_xyyz_1, g_yyy_0_xyzz_0, g_yyy_0_xyzz_1, g_yyy_0_xzzz_0, g_yyy_0_xzzz_1, g_yyy_0_yyyy_0, g_yyy_0_yyyy_1, g_yyy_0_yyyz_0, g_yyy_0_yyyz_1, g_yyy_0_yyzz_0, g_yyy_0_yyzz_1, g_yyy_0_yzzz_0, g_yyy_0_yzzz_1, g_yyy_0_zzzz_0, g_yyy_0_zzzz_1, g_yyyy_0_xxx_1, g_yyyy_0_xxxx_1, g_yyyy_0_xxxy_1, g_yyyy_0_xxxz_1, g_yyyy_0_xxy_1, g_yyyy_0_xxyy_1, g_yyyy_0_xxyz_1, g_yyyy_0_xxz_1, g_yyyy_0_xxzz_1, g_yyyy_0_xyy_1, g_yyyy_0_xyyy_1, g_yyyy_0_xyyz_1, g_yyyy_0_xyz_1, g_yyyy_0_xyzz_1, g_yyyy_0_xzz_1, g_yyyy_0_xzzz_1, g_yyyy_0_yyy_1, g_yyyy_0_yyyy_1, g_yyyy_0_yyyz_1, g_yyyy_0_yyz_1, g_yyyy_0_yyzz_1, g_yyyy_0_yzz_1, g_yyyy_0_yzzz_1, g_yyyy_0_zzz_1, g_yyyy_0_zzzz_1, g_yyyyy_0_xxxx_0, g_yyyyy_0_xxxy_0, g_yyyyy_0_xxxz_0, g_yyyyy_0_xxyy_0, g_yyyyy_0_xxyz_0, g_yyyyy_0_xxzz_0, g_yyyyy_0_xyyy_0, g_yyyyy_0_xyyz_0, g_yyyyy_0_xyzz_0, g_yyyyy_0_xzzz_0, g_yyyyy_0_yyyy_0, g_yyyyy_0_yyyz_0, g_yyyyy_0_yyzz_0, g_yyyyy_0_yzzz_0, g_yyyyy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xxxx_0[i] = 4.0 * g_yyy_0_xxxx_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxx_1[i] * fz_be_0 + g_yyyy_0_xxxx_1[i] * wa_y[i];

        g_yyyyy_0_xxxy_0[i] = 4.0 * g_yyy_0_xxxy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxy_1[i] * fz_be_0 + g_yyyy_0_xxx_1[i] * fi_acd_0 + g_yyyy_0_xxxy_1[i] * wa_y[i];

        g_yyyyy_0_xxxz_0[i] = 4.0 * g_yyy_0_xxxz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxz_1[i] * fz_be_0 + g_yyyy_0_xxxz_1[i] * wa_y[i];

        g_yyyyy_0_xxyy_0[i] = 4.0 * g_yyy_0_xxyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxy_1[i] * fi_acd_0 + g_yyyy_0_xxyy_1[i] * wa_y[i];

        g_yyyyy_0_xxyz_0[i] = 4.0 * g_yyy_0_xxyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyz_1[i] * fz_be_0 + g_yyyy_0_xxz_1[i] * fi_acd_0 + g_yyyy_0_xxyz_1[i] * wa_y[i];

        g_yyyyy_0_xxzz_0[i] = 4.0 * g_yyy_0_xxzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxzz_1[i] * fz_be_0 + g_yyyy_0_xxzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyy_0[i] = 4.0 * g_yyy_0_xyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xyy_1[i] * fi_acd_0 + g_yyyy_0_xyyy_1[i] * wa_y[i];

        g_yyyyy_0_xyyz_0[i] = 4.0 * g_yyy_0_xyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xyz_1[i] * fi_acd_0 + g_yyyy_0_xyyz_1[i] * wa_y[i];

        g_yyyyy_0_xyzz_0[i] = 4.0 * g_yyy_0_xyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyzz_1[i] * fz_be_0 + g_yyyy_0_xzz_1[i] * fi_acd_0 + g_yyyy_0_xyzz_1[i] * wa_y[i];

        g_yyyyy_0_xzzz_0[i] = 4.0 * g_yyy_0_xzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xzzz_1[i] * fz_be_0 + g_yyyy_0_xzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyy_0[i] = 4.0 * g_yyy_0_yyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_0_yyy_1[i] * fi_acd_0 + g_yyyy_0_yyyy_1[i] * wa_y[i];

        g_yyyyy_0_yyyz_0[i] = 4.0 * g_yyy_0_yyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_yyz_1[i] * fi_acd_0 + g_yyyy_0_yyyz_1[i] * wa_y[i];

        g_yyyyy_0_yyzz_0[i] = 4.0 * g_yyy_0_yyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_yzz_1[i] * fi_acd_0 + g_yyyy_0_yyzz_1[i] * wa_y[i];

        g_yyyyy_0_yzzz_0[i] = 4.0 * g_yyy_0_yzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yzzz_1[i] * fz_be_0 + g_yyyy_0_zzz_1[i] * fi_acd_0 + g_yyyy_0_yzzz_1[i] * wa_y[i];

        g_yyyyy_0_zzzz_0[i] = 4.0 * g_yyy_0_zzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_zzzz_1[i] * fz_be_0 + g_yyyy_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 240-255 components of targeted buffer : HSG

    auto g_yyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 240);

    auto g_yyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 241);

    auto g_yyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 242);

    auto g_yyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 243);

    auto g_yyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 244);

    auto g_yyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 245);

    auto g_yyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 246);

    auto g_yyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 247);

    auto g_yyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 248);

    auto g_yyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 249);

    auto g_yyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 250);

    auto g_yyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 251);

    auto g_yyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 252);

    auto g_yyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 253);

    auto g_yyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 254);

    #pragma omp simd aligned(g_yyyy_0_xxx_1, g_yyyy_0_xxxx_1, g_yyyy_0_xxxy_1, g_yyyy_0_xxxz_1, g_yyyy_0_xxy_1, g_yyyy_0_xxyy_1, g_yyyy_0_xxyz_1, g_yyyy_0_xxz_1, g_yyyy_0_xxzz_1, g_yyyy_0_xyy_1, g_yyyy_0_xyyy_1, g_yyyy_0_xyyz_1, g_yyyy_0_xyz_1, g_yyyy_0_xyzz_1, g_yyyy_0_xzz_1, g_yyyy_0_xzzz_1, g_yyyy_0_yyy_1, g_yyyy_0_yyyy_1, g_yyyy_0_yyyz_1, g_yyyy_0_yyz_1, g_yyyy_0_yyzz_1, g_yyyy_0_yzz_1, g_yyyy_0_yzzz_1, g_yyyy_0_zzz_1, g_yyyy_0_zzzz_1, g_yyyyz_0_xxxx_0, g_yyyyz_0_xxxy_0, g_yyyyz_0_xxxz_0, g_yyyyz_0_xxyy_0, g_yyyyz_0_xxyz_0, g_yyyyz_0_xxzz_0, g_yyyyz_0_xyyy_0, g_yyyyz_0_xyyz_0, g_yyyyz_0_xyzz_0, g_yyyyz_0_xzzz_0, g_yyyyz_0_yyyy_0, g_yyyyz_0_yyyz_0, g_yyyyz_0_yyzz_0, g_yyyyz_0_yzzz_0, g_yyyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xxxx_0[i] = g_yyyy_0_xxxx_1[i] * wa_z[i];

        g_yyyyz_0_xxxy_0[i] = g_yyyy_0_xxxy_1[i] * wa_z[i];

        g_yyyyz_0_xxxz_0[i] = g_yyyy_0_xxx_1[i] * fi_acd_0 + g_yyyy_0_xxxz_1[i] * wa_z[i];

        g_yyyyz_0_xxyy_0[i] = g_yyyy_0_xxyy_1[i] * wa_z[i];

        g_yyyyz_0_xxyz_0[i] = g_yyyy_0_xxy_1[i] * fi_acd_0 + g_yyyy_0_xxyz_1[i] * wa_z[i];

        g_yyyyz_0_xxzz_0[i] = 2.0 * g_yyyy_0_xxz_1[i] * fi_acd_0 + g_yyyy_0_xxzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyy_0[i] = g_yyyy_0_xyyy_1[i] * wa_z[i];

        g_yyyyz_0_xyyz_0[i] = g_yyyy_0_xyy_1[i] * fi_acd_0 + g_yyyy_0_xyyz_1[i] * wa_z[i];

        g_yyyyz_0_xyzz_0[i] = 2.0 * g_yyyy_0_xyz_1[i] * fi_acd_0 + g_yyyy_0_xyzz_1[i] * wa_z[i];

        g_yyyyz_0_xzzz_0[i] = 3.0 * g_yyyy_0_xzz_1[i] * fi_acd_0 + g_yyyy_0_xzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyy_0[i] = g_yyyy_0_yyyy_1[i] * wa_z[i];

        g_yyyyz_0_yyyz_0[i] = g_yyyy_0_yyy_1[i] * fi_acd_0 + g_yyyy_0_yyyz_1[i] * wa_z[i];

        g_yyyyz_0_yyzz_0[i] = 2.0 * g_yyyy_0_yyz_1[i] * fi_acd_0 + g_yyyy_0_yyzz_1[i] * wa_z[i];

        g_yyyyz_0_yzzz_0[i] = 3.0 * g_yyyy_0_yzz_1[i] * fi_acd_0 + g_yyyy_0_yzzz_1[i] * wa_z[i];

        g_yyyyz_0_zzzz_0[i] = 4.0 * g_yyyy_0_zzz_1[i] * fi_acd_0 + g_yyyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 255-270 components of targeted buffer : HSG

    auto g_yyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 255);

    auto g_yyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 256);

    auto g_yyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 257);

    auto g_yyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 258);

    auto g_yyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 259);

    auto g_yyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 260);

    auto g_yyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 261);

    auto g_yyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 262);

    auto g_yyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 263);

    auto g_yyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 264);

    auto g_yyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 265);

    auto g_yyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 266);

    auto g_yyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 267);

    auto g_yyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 268);

    auto g_yyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 269);

    #pragma omp simd aligned(g_yyy_0_xxxy_0, g_yyy_0_xxxy_1, g_yyy_0_xxyy_0, g_yyy_0_xxyy_1, g_yyy_0_xyyy_0, g_yyy_0_xyyy_1, g_yyy_0_yyyy_0, g_yyy_0_yyyy_1, g_yyyz_0_xxxy_1, g_yyyz_0_xxyy_1, g_yyyz_0_xyyy_1, g_yyyz_0_yyyy_1, g_yyyzz_0_xxxx_0, g_yyyzz_0_xxxy_0, g_yyyzz_0_xxxz_0, g_yyyzz_0_xxyy_0, g_yyyzz_0_xxyz_0, g_yyyzz_0_xxzz_0, g_yyyzz_0_xyyy_0, g_yyyzz_0_xyyz_0, g_yyyzz_0_xyzz_0, g_yyyzz_0_xzzz_0, g_yyyzz_0_yyyy_0, g_yyyzz_0_yyyz_0, g_yyyzz_0_yyzz_0, g_yyyzz_0_yzzz_0, g_yyyzz_0_zzzz_0, g_yyzz_0_xxxx_1, g_yyzz_0_xxxz_1, g_yyzz_0_xxyz_1, g_yyzz_0_xxz_1, g_yyzz_0_xxzz_1, g_yyzz_0_xyyz_1, g_yyzz_0_xyz_1, g_yyzz_0_xyzz_1, g_yyzz_0_xzz_1, g_yyzz_0_xzzz_1, g_yyzz_0_yyyz_1, g_yyzz_0_yyz_1, g_yyzz_0_yyzz_1, g_yyzz_0_yzz_1, g_yyzz_0_yzzz_1, g_yyzz_0_zzz_1, g_yyzz_0_zzzz_1, g_yzz_0_xxxx_0, g_yzz_0_xxxx_1, g_yzz_0_xxxz_0, g_yzz_0_xxxz_1, g_yzz_0_xxyz_0, g_yzz_0_xxyz_1, g_yzz_0_xxzz_0, g_yzz_0_xxzz_1, g_yzz_0_xyyz_0, g_yzz_0_xyyz_1, g_yzz_0_xyzz_0, g_yzz_0_xyzz_1, g_yzz_0_xzzz_0, g_yzz_0_xzzz_1, g_yzz_0_yyyz_0, g_yzz_0_yyyz_1, g_yzz_0_yyzz_0, g_yzz_0_yyzz_1, g_yzz_0_yzzz_0, g_yzz_0_yzzz_1, g_yzz_0_zzzz_0, g_yzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xxxx_0[i] = 2.0 * g_yzz_0_xxxx_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxx_1[i] * fz_be_0 + g_yyzz_0_xxxx_1[i] * wa_y[i];

        g_yyyzz_0_xxxy_0[i] = g_yyy_0_xxxy_0[i] * fbe_0 - g_yyy_0_xxxy_1[i] * fz_be_0 + g_yyyz_0_xxxy_1[i] * wa_z[i];

        g_yyyzz_0_xxxz_0[i] = 2.0 * g_yzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxz_1[i] * fz_be_0 + g_yyzz_0_xxxz_1[i] * wa_y[i];

        g_yyyzz_0_xxyy_0[i] = g_yyy_0_xxyy_0[i] * fbe_0 - g_yyy_0_xxyy_1[i] * fz_be_0 + g_yyyz_0_xxyy_1[i] * wa_z[i];

        g_yyyzz_0_xxyz_0[i] = 2.0 * g_yzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyz_1[i] * fz_be_0 + g_yyzz_0_xxz_1[i] * fi_acd_0 + g_yyzz_0_xxyz_1[i] * wa_y[i];

        g_yyyzz_0_xxzz_0[i] = 2.0 * g_yzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxzz_1[i] * fz_be_0 + g_yyzz_0_xxzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyy_0[i] = g_yyy_0_xyyy_0[i] * fbe_0 - g_yyy_0_xyyy_1[i] * fz_be_0 + g_yyyz_0_xyyy_1[i] * wa_z[i];

        g_yyyzz_0_xyyz_0[i] = 2.0 * g_yzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xyz_1[i] * fi_acd_0 + g_yyzz_0_xyyz_1[i] * wa_y[i];

        g_yyyzz_0_xyzz_0[i] = 2.0 * g_yzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyzz_1[i] * fz_be_0 + g_yyzz_0_xzz_1[i] * fi_acd_0 + g_yyzz_0_xyzz_1[i] * wa_y[i];

        g_yyyzz_0_xzzz_0[i] = 2.0 * g_yzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xzzz_1[i] * fz_be_0 + g_yyzz_0_xzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyy_0[i] = g_yyy_0_yyyy_0[i] * fbe_0 - g_yyy_0_yyyy_1[i] * fz_be_0 + g_yyyz_0_yyyy_1[i] * wa_z[i];

        g_yyyzz_0_yyyz_0[i] = 2.0 * g_yzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_yyz_1[i] * fi_acd_0 + g_yyzz_0_yyyz_1[i] * wa_y[i];

        g_yyyzz_0_yyzz_0[i] = 2.0 * g_yzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_yzz_1[i] * fi_acd_0 + g_yyzz_0_yyzz_1[i] * wa_y[i];

        g_yyyzz_0_yzzz_0[i] = 2.0 * g_yzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yzzz_1[i] * fz_be_0 + g_yyzz_0_zzz_1[i] * fi_acd_0 + g_yyzz_0_yzzz_1[i] * wa_y[i];

        g_yyyzz_0_zzzz_0[i] = 2.0 * g_yzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_zzzz_1[i] * fz_be_0 + g_yyzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 270-285 components of targeted buffer : HSG

    auto g_yyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 270);

    auto g_yyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 271);

    auto g_yyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 272);

    auto g_yyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 273);

    auto g_yyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 274);

    auto g_yyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 275);

    auto g_yyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 276);

    auto g_yyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 277);

    auto g_yyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 278);

    auto g_yyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 279);

    auto g_yyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 280);

    auto g_yyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 281);

    auto g_yyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 282);

    auto g_yyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 283);

    auto g_yyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 284);

    #pragma omp simd aligned(g_yyz_0_xxxy_0, g_yyz_0_xxxy_1, g_yyz_0_xxyy_0, g_yyz_0_xxyy_1, g_yyz_0_xyyy_0, g_yyz_0_xyyy_1, g_yyz_0_yyyy_0, g_yyz_0_yyyy_1, g_yyzz_0_xxxy_1, g_yyzz_0_xxyy_1, g_yyzz_0_xyyy_1, g_yyzz_0_yyyy_1, g_yyzzz_0_xxxx_0, g_yyzzz_0_xxxy_0, g_yyzzz_0_xxxz_0, g_yyzzz_0_xxyy_0, g_yyzzz_0_xxyz_0, g_yyzzz_0_xxzz_0, g_yyzzz_0_xyyy_0, g_yyzzz_0_xyyz_0, g_yyzzz_0_xyzz_0, g_yyzzz_0_xzzz_0, g_yyzzz_0_yyyy_0, g_yyzzz_0_yyyz_0, g_yyzzz_0_yyzz_0, g_yyzzz_0_yzzz_0, g_yyzzz_0_zzzz_0, g_yzzz_0_xxxx_1, g_yzzz_0_xxxz_1, g_yzzz_0_xxyz_1, g_yzzz_0_xxz_1, g_yzzz_0_xxzz_1, g_yzzz_0_xyyz_1, g_yzzz_0_xyz_1, g_yzzz_0_xyzz_1, g_yzzz_0_xzz_1, g_yzzz_0_xzzz_1, g_yzzz_0_yyyz_1, g_yzzz_0_yyz_1, g_yzzz_0_yyzz_1, g_yzzz_0_yzz_1, g_yzzz_0_yzzz_1, g_yzzz_0_zzz_1, g_yzzz_0_zzzz_1, g_zzz_0_xxxx_0, g_zzz_0_xxxx_1, g_zzz_0_xxxz_0, g_zzz_0_xxxz_1, g_zzz_0_xxyz_0, g_zzz_0_xxyz_1, g_zzz_0_xxzz_0, g_zzz_0_xxzz_1, g_zzz_0_xyyz_0, g_zzz_0_xyyz_1, g_zzz_0_xyzz_0, g_zzz_0_xyzz_1, g_zzz_0_xzzz_0, g_zzz_0_xzzz_1, g_zzz_0_yyyz_0, g_zzz_0_yyyz_1, g_zzz_0_yyzz_0, g_zzz_0_yyzz_1, g_zzz_0_yzzz_0, g_zzz_0_yzzz_1, g_zzz_0_zzzz_0, g_zzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xxxx_0[i] = g_zzz_0_xxxx_0[i] * fbe_0 - g_zzz_0_xxxx_1[i] * fz_be_0 + g_yzzz_0_xxxx_1[i] * wa_y[i];

        g_yyzzz_0_xxxy_0[i] = 2.0 * g_yyz_0_xxxy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxy_1[i] * fz_be_0 + g_yyzz_0_xxxy_1[i] * wa_z[i];

        g_yyzzz_0_xxxz_0[i] = g_zzz_0_xxxz_0[i] * fbe_0 - g_zzz_0_xxxz_1[i] * fz_be_0 + g_yzzz_0_xxxz_1[i] * wa_y[i];

        g_yyzzz_0_xxyy_0[i] = 2.0 * g_yyz_0_xxyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxyy_1[i] * fz_be_0 + g_yyzz_0_xxyy_1[i] * wa_z[i];

        g_yyzzz_0_xxyz_0[i] = g_zzz_0_xxyz_0[i] * fbe_0 - g_zzz_0_xxyz_1[i] * fz_be_0 + g_yzzz_0_xxz_1[i] * fi_acd_0 + g_yzzz_0_xxyz_1[i] * wa_y[i];

        g_yyzzz_0_xxzz_0[i] = g_zzz_0_xxzz_0[i] * fbe_0 - g_zzz_0_xxzz_1[i] * fz_be_0 + g_yzzz_0_xxzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyy_0[i] = 2.0 * g_yyz_0_xyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xyyy_1[i] * fz_be_0 + g_yyzz_0_xyyy_1[i] * wa_z[i];

        g_yyzzz_0_xyyz_0[i] = g_zzz_0_xyyz_0[i] * fbe_0 - g_zzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xyz_1[i] * fi_acd_0 + g_yzzz_0_xyyz_1[i] * wa_y[i];

        g_yyzzz_0_xyzz_0[i] = g_zzz_0_xyzz_0[i] * fbe_0 - g_zzz_0_xyzz_1[i] * fz_be_0 + g_yzzz_0_xzz_1[i] * fi_acd_0 + g_yzzz_0_xyzz_1[i] * wa_y[i];

        g_yyzzz_0_xzzz_0[i] = g_zzz_0_xzzz_0[i] * fbe_0 - g_zzz_0_xzzz_1[i] * fz_be_0 + g_yzzz_0_xzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyy_0[i] = 2.0 * g_yyz_0_yyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_yyyy_1[i] * fz_be_0 + g_yyzz_0_yyyy_1[i] * wa_z[i];

        g_yyzzz_0_yyyz_0[i] = g_zzz_0_yyyz_0[i] * fbe_0 - g_zzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_yyz_1[i] * fi_acd_0 + g_yzzz_0_yyyz_1[i] * wa_y[i];

        g_yyzzz_0_yyzz_0[i] = g_zzz_0_yyzz_0[i] * fbe_0 - g_zzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_yzz_1[i] * fi_acd_0 + g_yzzz_0_yyzz_1[i] * wa_y[i];

        g_yyzzz_0_yzzz_0[i] = g_zzz_0_yzzz_0[i] * fbe_0 - g_zzz_0_yzzz_1[i] * fz_be_0 + g_yzzz_0_zzz_1[i] * fi_acd_0 + g_yzzz_0_yzzz_1[i] * wa_y[i];

        g_yyzzz_0_zzzz_0[i] = g_zzz_0_zzzz_0[i] * fbe_0 - g_zzz_0_zzzz_1[i] * fz_be_0 + g_yzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 285-300 components of targeted buffer : HSG

    auto g_yzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 285);

    auto g_yzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 286);

    auto g_yzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 287);

    auto g_yzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 288);

    auto g_yzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 289);

    auto g_yzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 290);

    auto g_yzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 291);

    auto g_yzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 292);

    auto g_yzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 293);

    auto g_yzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 294);

    auto g_yzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 295);

    auto g_yzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 296);

    auto g_yzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 297);

    auto g_yzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 298);

    auto g_yzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 299);

    #pragma omp simd aligned(g_yzzzz_0_xxxx_0, g_yzzzz_0_xxxy_0, g_yzzzz_0_xxxz_0, g_yzzzz_0_xxyy_0, g_yzzzz_0_xxyz_0, g_yzzzz_0_xxzz_0, g_yzzzz_0_xyyy_0, g_yzzzz_0_xyyz_0, g_yzzzz_0_xyzz_0, g_yzzzz_0_xzzz_0, g_yzzzz_0_yyyy_0, g_yzzzz_0_yyyz_0, g_yzzzz_0_yyzz_0, g_yzzzz_0_yzzz_0, g_yzzzz_0_zzzz_0, g_zzzz_0_xxx_1, g_zzzz_0_xxxx_1, g_zzzz_0_xxxy_1, g_zzzz_0_xxxz_1, g_zzzz_0_xxy_1, g_zzzz_0_xxyy_1, g_zzzz_0_xxyz_1, g_zzzz_0_xxz_1, g_zzzz_0_xxzz_1, g_zzzz_0_xyy_1, g_zzzz_0_xyyy_1, g_zzzz_0_xyyz_1, g_zzzz_0_xyz_1, g_zzzz_0_xyzz_1, g_zzzz_0_xzz_1, g_zzzz_0_xzzz_1, g_zzzz_0_yyy_1, g_zzzz_0_yyyy_1, g_zzzz_0_yyyz_1, g_zzzz_0_yyz_1, g_zzzz_0_yyzz_1, g_zzzz_0_yzz_1, g_zzzz_0_yzzz_1, g_zzzz_0_zzz_1, g_zzzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xxxx_0[i] = g_zzzz_0_xxxx_1[i] * wa_y[i];

        g_yzzzz_0_xxxy_0[i] = g_zzzz_0_xxx_1[i] * fi_acd_0 + g_zzzz_0_xxxy_1[i] * wa_y[i];

        g_yzzzz_0_xxxz_0[i] = g_zzzz_0_xxxz_1[i] * wa_y[i];

        g_yzzzz_0_xxyy_0[i] = 2.0 * g_zzzz_0_xxy_1[i] * fi_acd_0 + g_zzzz_0_xxyy_1[i] * wa_y[i];

        g_yzzzz_0_xxyz_0[i] = g_zzzz_0_xxz_1[i] * fi_acd_0 + g_zzzz_0_xxyz_1[i] * wa_y[i];

        g_yzzzz_0_xxzz_0[i] = g_zzzz_0_xxzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyy_0[i] = 3.0 * g_zzzz_0_xyy_1[i] * fi_acd_0 + g_zzzz_0_xyyy_1[i] * wa_y[i];

        g_yzzzz_0_xyyz_0[i] = 2.0 * g_zzzz_0_xyz_1[i] * fi_acd_0 + g_zzzz_0_xyyz_1[i] * wa_y[i];

        g_yzzzz_0_xyzz_0[i] = g_zzzz_0_xzz_1[i] * fi_acd_0 + g_zzzz_0_xyzz_1[i] * wa_y[i];

        g_yzzzz_0_xzzz_0[i] = g_zzzz_0_xzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyy_0[i] = 4.0 * g_zzzz_0_yyy_1[i] * fi_acd_0 + g_zzzz_0_yyyy_1[i] * wa_y[i];

        g_yzzzz_0_yyyz_0[i] = 3.0 * g_zzzz_0_yyz_1[i] * fi_acd_0 + g_zzzz_0_yyyz_1[i] * wa_y[i];

        g_yzzzz_0_yyzz_0[i] = 2.0 * g_zzzz_0_yzz_1[i] * fi_acd_0 + g_zzzz_0_yyzz_1[i] * wa_y[i];

        g_yzzzz_0_yzzz_0[i] = g_zzzz_0_zzz_1[i] * fi_acd_0 + g_zzzz_0_yzzz_1[i] * wa_y[i];

        g_yzzzz_0_zzzz_0[i] = g_zzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 300-315 components of targeted buffer : HSG

    auto g_zzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 300);

    auto g_zzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 301);

    auto g_zzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 302);

    auto g_zzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 303);

    auto g_zzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 304);

    auto g_zzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 305);

    auto g_zzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 306);

    auto g_zzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 307);

    auto g_zzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 308);

    auto g_zzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 309);

    auto g_zzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 310);

    auto g_zzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 311);

    auto g_zzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 312);

    auto g_zzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 313);

    auto g_zzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 314);

    #pragma omp simd aligned(g_zzz_0_xxxx_0, g_zzz_0_xxxx_1, g_zzz_0_xxxy_0, g_zzz_0_xxxy_1, g_zzz_0_xxxz_0, g_zzz_0_xxxz_1, g_zzz_0_xxyy_0, g_zzz_0_xxyy_1, g_zzz_0_xxyz_0, g_zzz_0_xxyz_1, g_zzz_0_xxzz_0, g_zzz_0_xxzz_1, g_zzz_0_xyyy_0, g_zzz_0_xyyy_1, g_zzz_0_xyyz_0, g_zzz_0_xyyz_1, g_zzz_0_xyzz_0, g_zzz_0_xyzz_1, g_zzz_0_xzzz_0, g_zzz_0_xzzz_1, g_zzz_0_yyyy_0, g_zzz_0_yyyy_1, g_zzz_0_yyyz_0, g_zzz_0_yyyz_1, g_zzz_0_yyzz_0, g_zzz_0_yyzz_1, g_zzz_0_yzzz_0, g_zzz_0_yzzz_1, g_zzz_0_zzzz_0, g_zzz_0_zzzz_1, g_zzzz_0_xxx_1, g_zzzz_0_xxxx_1, g_zzzz_0_xxxy_1, g_zzzz_0_xxxz_1, g_zzzz_0_xxy_1, g_zzzz_0_xxyy_1, g_zzzz_0_xxyz_1, g_zzzz_0_xxz_1, g_zzzz_0_xxzz_1, g_zzzz_0_xyy_1, g_zzzz_0_xyyy_1, g_zzzz_0_xyyz_1, g_zzzz_0_xyz_1, g_zzzz_0_xyzz_1, g_zzzz_0_xzz_1, g_zzzz_0_xzzz_1, g_zzzz_0_yyy_1, g_zzzz_0_yyyy_1, g_zzzz_0_yyyz_1, g_zzzz_0_yyz_1, g_zzzz_0_yyzz_1, g_zzzz_0_yzz_1, g_zzzz_0_yzzz_1, g_zzzz_0_zzz_1, g_zzzz_0_zzzz_1, g_zzzzz_0_xxxx_0, g_zzzzz_0_xxxy_0, g_zzzzz_0_xxxz_0, g_zzzzz_0_xxyy_0, g_zzzzz_0_xxyz_0, g_zzzzz_0_xxzz_0, g_zzzzz_0_xyyy_0, g_zzzzz_0_xyyz_0, g_zzzzz_0_xyzz_0, g_zzzzz_0_xzzz_0, g_zzzzz_0_yyyy_0, g_zzzzz_0_yyyz_0, g_zzzzz_0_yyzz_0, g_zzzzz_0_yzzz_0, g_zzzzz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xxxx_0[i] = 4.0 * g_zzz_0_xxxx_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxx_1[i] * fz_be_0 + g_zzzz_0_xxxx_1[i] * wa_z[i];

        g_zzzzz_0_xxxy_0[i] = 4.0 * g_zzz_0_xxxy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxy_1[i] * fz_be_0 + g_zzzz_0_xxxy_1[i] * wa_z[i];

        g_zzzzz_0_xxxz_0[i] = 4.0 * g_zzz_0_xxxz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxz_1[i] * fz_be_0 + g_zzzz_0_xxx_1[i] * fi_acd_0 + g_zzzz_0_xxxz_1[i] * wa_z[i];

        g_zzzzz_0_xxyy_0[i] = 4.0 * g_zzz_0_xxyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyy_1[i] * fz_be_0 + g_zzzz_0_xxyy_1[i] * wa_z[i];

        g_zzzzz_0_xxyz_0[i] = 4.0 * g_zzz_0_xxyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyz_1[i] * fz_be_0 + g_zzzz_0_xxy_1[i] * fi_acd_0 + g_zzzz_0_xxyz_1[i] * wa_z[i];

        g_zzzzz_0_xxzz_0[i] = 4.0 * g_zzz_0_xxzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxz_1[i] * fi_acd_0 + g_zzzz_0_xxzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyy_0[i] = 4.0 * g_zzz_0_xyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyy_1[i] * fz_be_0 + g_zzzz_0_xyyy_1[i] * wa_z[i];

        g_zzzzz_0_xyyz_0[i] = 4.0 * g_zzz_0_xyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyz_1[i] * fz_be_0 + g_zzzz_0_xyy_1[i] * fi_acd_0 + g_zzzz_0_xyyz_1[i] * wa_z[i];

        g_zzzzz_0_xyzz_0[i] = 4.0 * g_zzz_0_xyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xyz_1[i] * fi_acd_0 + g_zzzz_0_xyzz_1[i] * wa_z[i];

        g_zzzzz_0_xzzz_0[i] = 4.0 * g_zzz_0_xzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xzz_1[i] * fi_acd_0 + g_zzzz_0_xzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyy_0[i] = 4.0 * g_zzz_0_yyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyy_1[i] * fz_be_0 + g_zzzz_0_yyyy_1[i] * wa_z[i];

        g_zzzzz_0_yyyz_0[i] = 4.0 * g_zzz_0_yyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyz_1[i] * fz_be_0 + g_zzzz_0_yyy_1[i] * fi_acd_0 + g_zzzz_0_yyyz_1[i] * wa_z[i];

        g_zzzzz_0_yyzz_0[i] = 4.0 * g_zzz_0_yyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_yyz_1[i] * fi_acd_0 + g_zzzz_0_yyzz_1[i] * wa_z[i];

        g_zzzzz_0_yzzz_0[i] = 4.0 * g_zzz_0_yzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_yzz_1[i] * fi_acd_0 + g_zzzz_0_yzzz_1[i] * wa_z[i];

        g_zzzzz_0_zzzz_0[i] = 4.0 * g_zzz_0_zzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_zzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_zzz_1[i] * fi_acd_0 + g_zzzz_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

