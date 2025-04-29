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

#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsh,
                                 size_t idx_eri_0_dsh,
                                 size_t idx_eri_1_dsh,
                                 size_t idx_eri_1_fsg,
                                 size_t idx_eri_1_fsh,
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

    /// Set up components of auxilary buffer : DSH

    auto g_xx_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh);

    auto g_xx_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 1);

    auto g_xx_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 2);

    auto g_xx_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 3);

    auto g_xx_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 4);

    auto g_xx_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 5);

    auto g_xx_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 6);

    auto g_xx_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 7);

    auto g_xx_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 8);

    auto g_xx_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 9);

    auto g_xx_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 10);

    auto g_xx_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 11);

    auto g_xx_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 12);

    auto g_xx_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 13);

    auto g_xx_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 14);

    auto g_xx_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 15);

    auto g_xx_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 16);

    auto g_xx_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 17);

    auto g_xx_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 18);

    auto g_xx_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 19);

    auto g_xx_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 20);

    auto g_yy_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 63);

    auto g_yy_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 64);

    auto g_yy_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 65);

    auto g_yy_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 66);

    auto g_yy_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 67);

    auto g_yy_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 68);

    auto g_yy_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 69);

    auto g_yy_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 70);

    auto g_yy_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 71);

    auto g_yy_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 72);

    auto g_yy_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 73);

    auto g_yy_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 74);

    auto g_yy_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 75);

    auto g_yy_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 76);

    auto g_yy_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 77);

    auto g_yy_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 78);

    auto g_yy_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 79);

    auto g_yy_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 80);

    auto g_yy_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 81);

    auto g_yy_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 82);

    auto g_yy_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 83);

    auto g_zz_0_xxxxx_0 = pbuffer.data(idx_eri_0_dsh + 105);

    auto g_zz_0_xxxxy_0 = pbuffer.data(idx_eri_0_dsh + 106);

    auto g_zz_0_xxxxz_0 = pbuffer.data(idx_eri_0_dsh + 107);

    auto g_zz_0_xxxyy_0 = pbuffer.data(idx_eri_0_dsh + 108);

    auto g_zz_0_xxxyz_0 = pbuffer.data(idx_eri_0_dsh + 109);

    auto g_zz_0_xxxzz_0 = pbuffer.data(idx_eri_0_dsh + 110);

    auto g_zz_0_xxyyy_0 = pbuffer.data(idx_eri_0_dsh + 111);

    auto g_zz_0_xxyyz_0 = pbuffer.data(idx_eri_0_dsh + 112);

    auto g_zz_0_xxyzz_0 = pbuffer.data(idx_eri_0_dsh + 113);

    auto g_zz_0_xxzzz_0 = pbuffer.data(idx_eri_0_dsh + 114);

    auto g_zz_0_xyyyy_0 = pbuffer.data(idx_eri_0_dsh + 115);

    auto g_zz_0_xyyyz_0 = pbuffer.data(idx_eri_0_dsh + 116);

    auto g_zz_0_xyyzz_0 = pbuffer.data(idx_eri_0_dsh + 117);

    auto g_zz_0_xyzzz_0 = pbuffer.data(idx_eri_0_dsh + 118);

    auto g_zz_0_xzzzz_0 = pbuffer.data(idx_eri_0_dsh + 119);

    auto g_zz_0_yyyyy_0 = pbuffer.data(idx_eri_0_dsh + 120);

    auto g_zz_0_yyyyz_0 = pbuffer.data(idx_eri_0_dsh + 121);

    auto g_zz_0_yyyzz_0 = pbuffer.data(idx_eri_0_dsh + 122);

    auto g_zz_0_yyzzz_0 = pbuffer.data(idx_eri_0_dsh + 123);

    auto g_zz_0_yzzzz_0 = pbuffer.data(idx_eri_0_dsh + 124);

    auto g_zz_0_zzzzz_0 = pbuffer.data(idx_eri_0_dsh + 125);

    /// Set up components of auxilary buffer : DSH

    auto g_xx_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh);

    auto g_xx_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 1);

    auto g_xx_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 2);

    auto g_xx_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 3);

    auto g_xx_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 4);

    auto g_xx_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 5);

    auto g_xx_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 6);

    auto g_xx_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 7);

    auto g_xx_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 8);

    auto g_xx_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 9);

    auto g_xx_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 10);

    auto g_xx_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 11);

    auto g_xx_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 12);

    auto g_xx_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 13);

    auto g_xx_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 14);

    auto g_xx_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 15);

    auto g_xx_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 16);

    auto g_xx_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 17);

    auto g_xx_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 18);

    auto g_xx_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 19);

    auto g_xx_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 20);

    auto g_yy_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 63);

    auto g_yy_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 64);

    auto g_yy_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 65);

    auto g_yy_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 66);

    auto g_yy_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 67);

    auto g_yy_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 68);

    auto g_yy_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 69);

    auto g_yy_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 70);

    auto g_yy_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 71);

    auto g_yy_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 72);

    auto g_yy_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 73);

    auto g_yy_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 74);

    auto g_yy_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 75);

    auto g_yy_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 76);

    auto g_yy_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 77);

    auto g_yy_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 78);

    auto g_yy_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 79);

    auto g_yy_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 80);

    auto g_yy_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 81);

    auto g_yy_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 82);

    auto g_yy_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 83);

    auto g_zz_0_xxxxx_1 = pbuffer.data(idx_eri_1_dsh + 105);

    auto g_zz_0_xxxxy_1 = pbuffer.data(idx_eri_1_dsh + 106);

    auto g_zz_0_xxxxz_1 = pbuffer.data(idx_eri_1_dsh + 107);

    auto g_zz_0_xxxyy_1 = pbuffer.data(idx_eri_1_dsh + 108);

    auto g_zz_0_xxxyz_1 = pbuffer.data(idx_eri_1_dsh + 109);

    auto g_zz_0_xxxzz_1 = pbuffer.data(idx_eri_1_dsh + 110);

    auto g_zz_0_xxyyy_1 = pbuffer.data(idx_eri_1_dsh + 111);

    auto g_zz_0_xxyyz_1 = pbuffer.data(idx_eri_1_dsh + 112);

    auto g_zz_0_xxyzz_1 = pbuffer.data(idx_eri_1_dsh + 113);

    auto g_zz_0_xxzzz_1 = pbuffer.data(idx_eri_1_dsh + 114);

    auto g_zz_0_xyyyy_1 = pbuffer.data(idx_eri_1_dsh + 115);

    auto g_zz_0_xyyyz_1 = pbuffer.data(idx_eri_1_dsh + 116);

    auto g_zz_0_xyyzz_1 = pbuffer.data(idx_eri_1_dsh + 117);

    auto g_zz_0_xyzzz_1 = pbuffer.data(idx_eri_1_dsh + 118);

    auto g_zz_0_xzzzz_1 = pbuffer.data(idx_eri_1_dsh + 119);

    auto g_zz_0_yyyyy_1 = pbuffer.data(idx_eri_1_dsh + 120);

    auto g_zz_0_yyyyz_1 = pbuffer.data(idx_eri_1_dsh + 121);

    auto g_zz_0_yyyzz_1 = pbuffer.data(idx_eri_1_dsh + 122);

    auto g_zz_0_yyzzz_1 = pbuffer.data(idx_eri_1_dsh + 123);

    auto g_zz_0_yzzzz_1 = pbuffer.data(idx_eri_1_dsh + 124);

    auto g_zz_0_zzzzz_1 = pbuffer.data(idx_eri_1_dsh + 125);

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

    auto g_xxz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 32);

    auto g_xxz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 34);

    auto g_xxz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 35);

    auto g_xxz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 37);

    auto g_xxz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 38);

    auto g_xxz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 39);

    auto g_xxz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 41);

    auto g_xxz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 42);

    auto g_xxz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 43);

    auto g_xxz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 44);

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

    auto g_xzz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 77);

    auto g_xzz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 79);

    auto g_xzz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 80);

    auto g_xzz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 82);

    auto g_xzz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 83);

    auto g_xzz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 84);

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

    auto g_yyz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 107);

    auto g_yyz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 109);

    auto g_yyz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 110);

    auto g_yyz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 112);

    auto g_yyz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 113);

    auto g_yyz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 114);

    auto g_yyz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 116);

    auto g_yyz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 117);

    auto g_yyz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 118);

    auto g_yyz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 119);

    auto g_yzz_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 121);

    auto g_yzz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 122);

    auto g_yzz_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 123);

    auto g_yzz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 124);

    auto g_yzz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 125);

    auto g_yzz_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 126);

    auto g_yzz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 127);

    auto g_yzz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 128);

    auto g_yzz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 129);

    auto g_yzz_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 130);

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

    /// Set up components of auxilary buffer : FSH

    auto g_xxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh);

    auto g_xxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 1);

    auto g_xxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 2);

    auto g_xxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 3);

    auto g_xxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 4);

    auto g_xxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 5);

    auto g_xxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 6);

    auto g_xxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 7);

    auto g_xxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 8);

    auto g_xxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 9);

    auto g_xxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 10);

    auto g_xxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 11);

    auto g_xxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 12);

    auto g_xxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 13);

    auto g_xxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 14);

    auto g_xxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 15);

    auto g_xxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 16);

    auto g_xxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 17);

    auto g_xxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 18);

    auto g_xxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 19);

    auto g_xxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 20);

    auto g_xxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 21);

    auto g_xxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 22);

    auto g_xxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 23);

    auto g_xxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 24);

    auto g_xxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 26);

    auto g_xxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 27);

    auto g_xxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 30);

    auto g_xxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 31);

    auto g_xxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 35);

    auto g_xxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 36);

    auto g_xxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 42);

    auto g_xxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 43);

    auto g_xxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 44);

    auto g_xxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 45);

    auto g_xxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 46);

    auto g_xxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 47);

    auto g_xxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 48);

    auto g_xxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 49);

    auto g_xxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 50);

    auto g_xxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 51);

    auto g_xxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 52);

    auto g_xxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 53);

    auto g_xxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 54);

    auto g_xxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 55);

    auto g_xxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 56);

    auto g_xxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 58);

    auto g_xxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 59);

    auto g_xxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 60);

    auto g_xxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 61);

    auto g_xxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 62);

    auto g_xyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 63);

    auto g_xyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 64);

    auto g_xyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 66);

    auto g_xyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 67);

    auto g_xyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 69);

    auto g_xyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 70);

    auto g_xyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 71);

    auto g_xyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 73);

    auto g_xyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 74);

    auto g_xyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 75);

    auto g_xyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 76);

    auto g_xyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 78);

    auto g_xyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 79);

    auto g_xyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 80);

    auto g_xyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 81);

    auto g_xyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 82);

    auto g_xyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 83);

    auto g_xzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 105);

    auto g_xzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 107);

    auto g_xzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 109);

    auto g_xzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 110);

    auto g_xzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 112);

    auto g_xzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 113);

    auto g_xzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 114);

    auto g_xzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 116);

    auto g_xzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 117);

    auto g_xzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 118);

    auto g_xzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 119);

    auto g_xzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 120);

    auto g_xzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 121);

    auto g_xzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 122);

    auto g_xzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 123);

    auto g_xzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 124);

    auto g_xzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 125);

    auto g_yyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 126);

    auto g_yyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 127);

    auto g_yyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 128);

    auto g_yyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 129);

    auto g_yyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 130);

    auto g_yyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 131);

    auto g_yyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 132);

    auto g_yyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 133);

    auto g_yyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 134);

    auto g_yyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 135);

    auto g_yyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 136);

    auto g_yyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 137);

    auto g_yyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 138);

    auto g_yyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 139);

    auto g_yyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 140);

    auto g_yyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 141);

    auto g_yyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 142);

    auto g_yyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 143);

    auto g_yyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 144);

    auto g_yyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 145);

    auto g_yyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 146);

    auto g_yyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 148);

    auto g_yyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 149);

    auto g_yyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 150);

    auto g_yyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 151);

    auto g_yyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 152);

    auto g_yyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 153);

    auto g_yyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 154);

    auto g_yyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 155);

    auto g_yyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 156);

    auto g_yyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 157);

    auto g_yyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 158);

    auto g_yyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 159);

    auto g_yyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 160);

    auto g_yyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 161);

    auto g_yyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 162);

    auto g_yyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 163);

    auto g_yyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 164);

    auto g_yyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 165);

    auto g_yyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 166);

    auto g_yyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 167);

    auto g_yzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 168);

    auto g_yzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 169);

    auto g_yzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 170);

    auto g_yzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 171);

    auto g_yzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 172);

    auto g_yzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 173);

    auto g_yzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 174);

    auto g_yzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 175);

    auto g_yzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 176);

    auto g_yzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 177);

    auto g_yzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 178);

    auto g_yzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 179);

    auto g_yzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 180);

    auto g_yzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 181);

    auto g_yzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 182);

    auto g_yzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 183);

    auto g_yzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 184);

    auto g_yzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 185);

    auto g_yzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 186);

    auto g_yzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 187);

    auto g_yzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 188);

    auto g_zzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 189);

    auto g_zzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 190);

    auto g_zzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 191);

    auto g_zzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 192);

    auto g_zzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 193);

    auto g_zzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 194);

    auto g_zzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 195);

    auto g_zzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 196);

    auto g_zzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 197);

    auto g_zzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 198);

    auto g_zzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 199);

    auto g_zzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 200);

    auto g_zzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 201);

    auto g_zzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 202);

    auto g_zzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 203);

    auto g_zzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 204);

    auto g_zzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 205);

    auto g_zzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 206);

    auto g_zzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 207);

    auto g_zzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 208);

    auto g_zzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 209);

    /// Set up 0-21 components of targeted buffer : GSH

    auto g_xxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh);

    auto g_xxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 1);

    auto g_xxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 2);

    auto g_xxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 3);

    auto g_xxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 4);

    auto g_xxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 5);

    auto g_xxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 6);

    auto g_xxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 7);

    auto g_xxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 8);

    auto g_xxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 9);

    auto g_xxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 10);

    auto g_xxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 11);

    auto g_xxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 12);

    auto g_xxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 13);

    auto g_xxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 14);

    auto g_xxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 15);

    auto g_xxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 16);

    auto g_xxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 17);

    auto g_xxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 18);

    auto g_xxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 19);

    auto g_xxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 20);

    #pragma omp simd aligned(g_xx_0_xxxxx_0, g_xx_0_xxxxx_1, g_xx_0_xxxxy_0, g_xx_0_xxxxy_1, g_xx_0_xxxxz_0, g_xx_0_xxxxz_1, g_xx_0_xxxyy_0, g_xx_0_xxxyy_1, g_xx_0_xxxyz_0, g_xx_0_xxxyz_1, g_xx_0_xxxzz_0, g_xx_0_xxxzz_1, g_xx_0_xxyyy_0, g_xx_0_xxyyy_1, g_xx_0_xxyyz_0, g_xx_0_xxyyz_1, g_xx_0_xxyzz_0, g_xx_0_xxyzz_1, g_xx_0_xxzzz_0, g_xx_0_xxzzz_1, g_xx_0_xyyyy_0, g_xx_0_xyyyy_1, g_xx_0_xyyyz_0, g_xx_0_xyyyz_1, g_xx_0_xyyzz_0, g_xx_0_xyyzz_1, g_xx_0_xyzzz_0, g_xx_0_xyzzz_1, g_xx_0_xzzzz_0, g_xx_0_xzzzz_1, g_xx_0_yyyyy_0, g_xx_0_yyyyy_1, g_xx_0_yyyyz_0, g_xx_0_yyyyz_1, g_xx_0_yyyzz_0, g_xx_0_yyyzz_1, g_xx_0_yyzzz_0, g_xx_0_yyzzz_1, g_xx_0_yzzzz_0, g_xx_0_yzzzz_1, g_xx_0_zzzzz_0, g_xx_0_zzzzz_1, g_xxx_0_xxxx_1, g_xxx_0_xxxxx_1, g_xxx_0_xxxxy_1, g_xxx_0_xxxxz_1, g_xxx_0_xxxy_1, g_xxx_0_xxxyy_1, g_xxx_0_xxxyz_1, g_xxx_0_xxxz_1, g_xxx_0_xxxzz_1, g_xxx_0_xxyy_1, g_xxx_0_xxyyy_1, g_xxx_0_xxyyz_1, g_xxx_0_xxyz_1, g_xxx_0_xxyzz_1, g_xxx_0_xxzz_1, g_xxx_0_xxzzz_1, g_xxx_0_xyyy_1, g_xxx_0_xyyyy_1, g_xxx_0_xyyyz_1, g_xxx_0_xyyz_1, g_xxx_0_xyyzz_1, g_xxx_0_xyzz_1, g_xxx_0_xyzzz_1, g_xxx_0_xzzz_1, g_xxx_0_xzzzz_1, g_xxx_0_yyyy_1, g_xxx_0_yyyyy_1, g_xxx_0_yyyyz_1, g_xxx_0_yyyz_1, g_xxx_0_yyyzz_1, g_xxx_0_yyzz_1, g_xxx_0_yyzzz_1, g_xxx_0_yzzz_1, g_xxx_0_yzzzz_1, g_xxx_0_zzzz_1, g_xxx_0_zzzzz_1, g_xxxx_0_xxxxx_0, g_xxxx_0_xxxxy_0, g_xxxx_0_xxxxz_0, g_xxxx_0_xxxyy_0, g_xxxx_0_xxxyz_0, g_xxxx_0_xxxzz_0, g_xxxx_0_xxyyy_0, g_xxxx_0_xxyyz_0, g_xxxx_0_xxyzz_0, g_xxxx_0_xxzzz_0, g_xxxx_0_xyyyy_0, g_xxxx_0_xyyyz_0, g_xxxx_0_xyyzz_0, g_xxxx_0_xyzzz_0, g_xxxx_0_xzzzz_0, g_xxxx_0_yyyyy_0, g_xxxx_0_yyyyz_0, g_xxxx_0_yyyzz_0, g_xxxx_0_yyzzz_0, g_xxxx_0_yzzzz_0, g_xxxx_0_zzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxxxx_0[i] = 3.0 * g_xx_0_xxxxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxx_1[i] * wa_x[i];

        g_xxxx_0_xxxxy_0[i] = 3.0 * g_xx_0_xxxxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxy_1[i] * wa_x[i];

        g_xxxx_0_xxxxz_0[i] = 3.0 * g_xx_0_xxxxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxz_1[i] * wa_x[i];

        g_xxxx_0_xxxyy_0[i] = 3.0 * g_xx_0_xxxyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyy_1[i] * fi_acd_0 + g_xxx_0_xxxyy_1[i] * wa_x[i];

        g_xxxx_0_xxxyz_0[i] = 3.0 * g_xx_0_xxxyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyz_1[i] * fi_acd_0 + g_xxx_0_xxxyz_1[i] * wa_x[i];

        g_xxxx_0_xxxzz_0[i] = 3.0 * g_xx_0_xxxzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxzz_1[i] * fi_acd_0 + g_xxx_0_xxxzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyy_0[i] = 3.0 * g_xx_0_xxyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyy_1[i] * wa_x[i];

        g_xxxx_0_xxyyz_0[i] = 3.0 * g_xx_0_xxyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyz_1[i] * wa_x[i];

        g_xxxx_0_xxyzz_0[i] = 3.0 * g_xx_0_xxyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyzz_1[i] * fi_acd_0 + g_xxx_0_xxyzz_1[i] * wa_x[i];

        g_xxxx_0_xxzzz_0[i] = 3.0 * g_xx_0_xxzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyy_0[i] = 3.0 * g_xx_0_xyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyy_1[i] * fz_be_0 + g_xxx_0_yyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyy_1[i] * wa_x[i];

        g_xxxx_0_xyyyz_0[i] = 3.0 * g_xx_0_xyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyz_1[i] * fz_be_0 + g_xxx_0_yyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyz_1[i] * wa_x[i];

        g_xxxx_0_xyyzz_0[i] = 3.0 * g_xx_0_xyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyzz_1[i] * fz_be_0 + g_xxx_0_yyzz_1[i] * fi_acd_0 + g_xxx_0_xyyzz_1[i] * wa_x[i];

        g_xxxx_0_xyzzz_0[i] = 3.0 * g_xx_0_xyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyzzz_1[i] * fz_be_0 + g_xxx_0_yzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzz_1[i] * wa_x[i];

        g_xxxx_0_xzzzz_0[i] = 3.0 * g_xx_0_xzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzzzz_1[i] * fz_be_0 + g_xxx_0_zzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyy_0[i] = 3.0 * g_xx_0_yyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyy_1[i] * wa_x[i];

        g_xxxx_0_yyyyz_0[i] = 3.0 * g_xx_0_yyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyz_1[i] * wa_x[i];

        g_xxxx_0_yyyzz_0[i] = 3.0 * g_xx_0_yyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyzz_1[i] * fz_be_0 + g_xxx_0_yyyzz_1[i] * wa_x[i];

        g_xxxx_0_yyzzz_0[i] = 3.0 * g_xx_0_yyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyzzz_1[i] * fz_be_0 + g_xxx_0_yyzzz_1[i] * wa_x[i];

        g_xxxx_0_yzzzz_0[i] = 3.0 * g_xx_0_yzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzz_1[i] * wa_x[i];

        g_xxxx_0_zzzzz_0[i] = 3.0 * g_xx_0_zzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : GSH

    auto g_xxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 21);

    auto g_xxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 22);

    auto g_xxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 23);

    auto g_xxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 24);

    auto g_xxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 25);

    auto g_xxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 26);

    auto g_xxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 27);

    auto g_xxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 28);

    auto g_xxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 29);

    auto g_xxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 30);

    auto g_xxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 31);

    auto g_xxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 32);

    auto g_xxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 33);

    auto g_xxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 34);

    auto g_xxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 35);

    auto g_xxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 36);

    auto g_xxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 37);

    auto g_xxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 38);

    auto g_xxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 39);

    auto g_xxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 40);

    auto g_xxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 41);

    #pragma omp simd aligned(g_xxx_0_xxxx_1, g_xxx_0_xxxxx_1, g_xxx_0_xxxxy_1, g_xxx_0_xxxxz_1, g_xxx_0_xxxy_1, g_xxx_0_xxxyy_1, g_xxx_0_xxxyz_1, g_xxx_0_xxxz_1, g_xxx_0_xxxzz_1, g_xxx_0_xxyy_1, g_xxx_0_xxyyy_1, g_xxx_0_xxyyz_1, g_xxx_0_xxyz_1, g_xxx_0_xxyzz_1, g_xxx_0_xxzz_1, g_xxx_0_xxzzz_1, g_xxx_0_xyyy_1, g_xxx_0_xyyyy_1, g_xxx_0_xyyyz_1, g_xxx_0_xyyz_1, g_xxx_0_xyyzz_1, g_xxx_0_xyzz_1, g_xxx_0_xyzzz_1, g_xxx_0_xzzz_1, g_xxx_0_xzzzz_1, g_xxx_0_yyyy_1, g_xxx_0_yyyyy_1, g_xxx_0_yyyyz_1, g_xxx_0_yyyz_1, g_xxx_0_yyyzz_1, g_xxx_0_yyzz_1, g_xxx_0_yyzzz_1, g_xxx_0_yzzz_1, g_xxx_0_yzzzz_1, g_xxx_0_zzzz_1, g_xxx_0_zzzzz_1, g_xxxy_0_xxxxx_0, g_xxxy_0_xxxxy_0, g_xxxy_0_xxxxz_0, g_xxxy_0_xxxyy_0, g_xxxy_0_xxxyz_0, g_xxxy_0_xxxzz_0, g_xxxy_0_xxyyy_0, g_xxxy_0_xxyyz_0, g_xxxy_0_xxyzz_0, g_xxxy_0_xxzzz_0, g_xxxy_0_xyyyy_0, g_xxxy_0_xyyyz_0, g_xxxy_0_xyyzz_0, g_xxxy_0_xyzzz_0, g_xxxy_0_xzzzz_0, g_xxxy_0_yyyyy_0, g_xxxy_0_yyyyz_0, g_xxxy_0_yyyzz_0, g_xxxy_0_yyzzz_0, g_xxxy_0_yzzzz_0, g_xxxy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxxxx_0[i] = g_xxx_0_xxxxx_1[i] * wa_y[i];

        g_xxxy_0_xxxxy_0[i] = g_xxx_0_xxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxy_1[i] * wa_y[i];

        g_xxxy_0_xxxxz_0[i] = g_xxx_0_xxxxz_1[i] * wa_y[i];

        g_xxxy_0_xxxyy_0[i] = 2.0 * g_xxx_0_xxxy_1[i] * fi_acd_0 + g_xxx_0_xxxyy_1[i] * wa_y[i];

        g_xxxy_0_xxxyz_0[i] = g_xxx_0_xxxz_1[i] * fi_acd_0 + g_xxx_0_xxxyz_1[i] * wa_y[i];

        g_xxxy_0_xxxzz_0[i] = g_xxx_0_xxxzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyy_0[i] = 3.0 * g_xxx_0_xxyy_1[i] * fi_acd_0 + g_xxx_0_xxyyy_1[i] * wa_y[i];

        g_xxxy_0_xxyyz_0[i] = 2.0 * g_xxx_0_xxyz_1[i] * fi_acd_0 + g_xxx_0_xxyyz_1[i] * wa_y[i];

        g_xxxy_0_xxyzz_0[i] = g_xxx_0_xxzz_1[i] * fi_acd_0 + g_xxx_0_xxyzz_1[i] * wa_y[i];

        g_xxxy_0_xxzzz_0[i] = g_xxx_0_xxzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyy_0[i] = 4.0 * g_xxx_0_xyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyy_1[i] * wa_y[i];

        g_xxxy_0_xyyyz_0[i] = 3.0 * g_xxx_0_xyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyz_1[i] * wa_y[i];

        g_xxxy_0_xyyzz_0[i] = 2.0 * g_xxx_0_xyzz_1[i] * fi_acd_0 + g_xxx_0_xyyzz_1[i] * wa_y[i];

        g_xxxy_0_xyzzz_0[i] = g_xxx_0_xzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzz_1[i] * wa_y[i];

        g_xxxy_0_xzzzz_0[i] = g_xxx_0_xzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyy_0[i] = 5.0 * g_xxx_0_yyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyy_1[i] * wa_y[i];

        g_xxxy_0_yyyyz_0[i] = 4.0 * g_xxx_0_yyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyz_1[i] * wa_y[i];

        g_xxxy_0_yyyzz_0[i] = 3.0 * g_xxx_0_yyzz_1[i] * fi_acd_0 + g_xxx_0_yyyzz_1[i] * wa_y[i];

        g_xxxy_0_yyzzz_0[i] = 2.0 * g_xxx_0_yzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzz_1[i] * wa_y[i];

        g_xxxy_0_yzzzz_0[i] = g_xxx_0_zzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzz_1[i] * wa_y[i];

        g_xxxy_0_zzzzz_0[i] = g_xxx_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 42-63 components of targeted buffer : GSH

    auto g_xxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 42);

    auto g_xxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 43);

    auto g_xxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 44);

    auto g_xxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 45);

    auto g_xxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 46);

    auto g_xxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 47);

    auto g_xxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 48);

    auto g_xxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 49);

    auto g_xxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 50);

    auto g_xxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 51);

    auto g_xxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 52);

    auto g_xxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 53);

    auto g_xxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 54);

    auto g_xxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 55);

    auto g_xxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 56);

    auto g_xxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 57);

    auto g_xxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 58);

    auto g_xxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 59);

    auto g_xxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 60);

    auto g_xxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 61);

    auto g_xxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 62);

    #pragma omp simd aligned(g_xxx_0_xxxx_1, g_xxx_0_xxxxx_1, g_xxx_0_xxxxy_1, g_xxx_0_xxxxz_1, g_xxx_0_xxxy_1, g_xxx_0_xxxyy_1, g_xxx_0_xxxyz_1, g_xxx_0_xxxz_1, g_xxx_0_xxxzz_1, g_xxx_0_xxyy_1, g_xxx_0_xxyyy_1, g_xxx_0_xxyyz_1, g_xxx_0_xxyz_1, g_xxx_0_xxyzz_1, g_xxx_0_xxzz_1, g_xxx_0_xxzzz_1, g_xxx_0_xyyy_1, g_xxx_0_xyyyy_1, g_xxx_0_xyyyz_1, g_xxx_0_xyyz_1, g_xxx_0_xyyzz_1, g_xxx_0_xyzz_1, g_xxx_0_xyzzz_1, g_xxx_0_xzzz_1, g_xxx_0_xzzzz_1, g_xxx_0_yyyy_1, g_xxx_0_yyyyy_1, g_xxx_0_yyyyz_1, g_xxx_0_yyyz_1, g_xxx_0_yyyzz_1, g_xxx_0_yyzz_1, g_xxx_0_yyzzz_1, g_xxx_0_yzzz_1, g_xxx_0_yzzzz_1, g_xxx_0_zzzz_1, g_xxx_0_zzzzz_1, g_xxxz_0_xxxxx_0, g_xxxz_0_xxxxy_0, g_xxxz_0_xxxxz_0, g_xxxz_0_xxxyy_0, g_xxxz_0_xxxyz_0, g_xxxz_0_xxxzz_0, g_xxxz_0_xxyyy_0, g_xxxz_0_xxyyz_0, g_xxxz_0_xxyzz_0, g_xxxz_0_xxzzz_0, g_xxxz_0_xyyyy_0, g_xxxz_0_xyyyz_0, g_xxxz_0_xyyzz_0, g_xxxz_0_xyzzz_0, g_xxxz_0_xzzzz_0, g_xxxz_0_yyyyy_0, g_xxxz_0_yyyyz_0, g_xxxz_0_yyyzz_0, g_xxxz_0_yyzzz_0, g_xxxz_0_yzzzz_0, g_xxxz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxxxx_0[i] = g_xxx_0_xxxxx_1[i] * wa_z[i];

        g_xxxz_0_xxxxy_0[i] = g_xxx_0_xxxxy_1[i] * wa_z[i];

        g_xxxz_0_xxxxz_0[i] = g_xxx_0_xxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxz_1[i] * wa_z[i];

        g_xxxz_0_xxxyy_0[i] = g_xxx_0_xxxyy_1[i] * wa_z[i];

        g_xxxz_0_xxxyz_0[i] = g_xxx_0_xxxy_1[i] * fi_acd_0 + g_xxx_0_xxxyz_1[i] * wa_z[i];

        g_xxxz_0_xxxzz_0[i] = 2.0 * g_xxx_0_xxxz_1[i] * fi_acd_0 + g_xxx_0_xxxzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyy_0[i] = g_xxx_0_xxyyy_1[i] * wa_z[i];

        g_xxxz_0_xxyyz_0[i] = g_xxx_0_xxyy_1[i] * fi_acd_0 + g_xxx_0_xxyyz_1[i] * wa_z[i];

        g_xxxz_0_xxyzz_0[i] = 2.0 * g_xxx_0_xxyz_1[i] * fi_acd_0 + g_xxx_0_xxyzz_1[i] * wa_z[i];

        g_xxxz_0_xxzzz_0[i] = 3.0 * g_xxx_0_xxzz_1[i] * fi_acd_0 + g_xxx_0_xxzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyy_0[i] = g_xxx_0_xyyyy_1[i] * wa_z[i];

        g_xxxz_0_xyyyz_0[i] = g_xxx_0_xyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyz_1[i] * wa_z[i];

        g_xxxz_0_xyyzz_0[i] = 2.0 * g_xxx_0_xyyz_1[i] * fi_acd_0 + g_xxx_0_xyyzz_1[i] * wa_z[i];

        g_xxxz_0_xyzzz_0[i] = 3.0 * g_xxx_0_xyzz_1[i] * fi_acd_0 + g_xxx_0_xyzzz_1[i] * wa_z[i];

        g_xxxz_0_xzzzz_0[i] = 4.0 * g_xxx_0_xzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyy_0[i] = g_xxx_0_yyyyy_1[i] * wa_z[i];

        g_xxxz_0_yyyyz_0[i] = g_xxx_0_yyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyz_1[i] * wa_z[i];

        g_xxxz_0_yyyzz_0[i] = 2.0 * g_xxx_0_yyyz_1[i] * fi_acd_0 + g_xxx_0_yyyzz_1[i] * wa_z[i];

        g_xxxz_0_yyzzz_0[i] = 3.0 * g_xxx_0_yyzz_1[i] * fi_acd_0 + g_xxx_0_yyzzz_1[i] * wa_z[i];

        g_xxxz_0_yzzzz_0[i] = 4.0 * g_xxx_0_yzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzz_1[i] * wa_z[i];

        g_xxxz_0_zzzzz_0[i] = 5.0 * g_xxx_0_zzzz_1[i] * fi_acd_0 + g_xxx_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 63-84 components of targeted buffer : GSH

    auto g_xxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 63);

    auto g_xxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 64);

    auto g_xxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 65);

    auto g_xxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 66);

    auto g_xxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 67);

    auto g_xxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 68);

    auto g_xxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 69);

    auto g_xxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 70);

    auto g_xxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 71);

    auto g_xxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 72);

    auto g_xxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 73);

    auto g_xxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 74);

    auto g_xxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 75);

    auto g_xxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 76);

    auto g_xxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 77);

    auto g_xxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 78);

    auto g_xxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 79);

    auto g_xxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 80);

    auto g_xxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 81);

    auto g_xxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 82);

    auto g_xxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 83);

    #pragma omp simd aligned(g_xx_0_xxxxx_0, g_xx_0_xxxxx_1, g_xx_0_xxxxz_0, g_xx_0_xxxxz_1, g_xx_0_xxxzz_0, g_xx_0_xxxzz_1, g_xx_0_xxzzz_0, g_xx_0_xxzzz_1, g_xx_0_xzzzz_0, g_xx_0_xzzzz_1, g_xxy_0_xxxxx_1, g_xxy_0_xxxxz_1, g_xxy_0_xxxzz_1, g_xxy_0_xxzzz_1, g_xxy_0_xzzzz_1, g_xxyy_0_xxxxx_0, g_xxyy_0_xxxxy_0, g_xxyy_0_xxxxz_0, g_xxyy_0_xxxyy_0, g_xxyy_0_xxxyz_0, g_xxyy_0_xxxzz_0, g_xxyy_0_xxyyy_0, g_xxyy_0_xxyyz_0, g_xxyy_0_xxyzz_0, g_xxyy_0_xxzzz_0, g_xxyy_0_xyyyy_0, g_xxyy_0_xyyyz_0, g_xxyy_0_xyyzz_0, g_xxyy_0_xyzzz_0, g_xxyy_0_xzzzz_0, g_xxyy_0_yyyyy_0, g_xxyy_0_yyyyz_0, g_xxyy_0_yyyzz_0, g_xxyy_0_yyzzz_0, g_xxyy_0_yzzzz_0, g_xxyy_0_zzzzz_0, g_xyy_0_xxxxy_1, g_xyy_0_xxxy_1, g_xyy_0_xxxyy_1, g_xyy_0_xxxyz_1, g_xyy_0_xxyy_1, g_xyy_0_xxyyy_1, g_xyy_0_xxyyz_1, g_xyy_0_xxyz_1, g_xyy_0_xxyzz_1, g_xyy_0_xyyy_1, g_xyy_0_xyyyy_1, g_xyy_0_xyyyz_1, g_xyy_0_xyyz_1, g_xyy_0_xyyzz_1, g_xyy_0_xyzz_1, g_xyy_0_xyzzz_1, g_xyy_0_yyyy_1, g_xyy_0_yyyyy_1, g_xyy_0_yyyyz_1, g_xyy_0_yyyz_1, g_xyy_0_yyyzz_1, g_xyy_0_yyzz_1, g_xyy_0_yyzzz_1, g_xyy_0_yzzz_1, g_xyy_0_yzzzz_1, g_xyy_0_zzzzz_1, g_yy_0_xxxxy_0, g_yy_0_xxxxy_1, g_yy_0_xxxyy_0, g_yy_0_xxxyy_1, g_yy_0_xxxyz_0, g_yy_0_xxxyz_1, g_yy_0_xxyyy_0, g_yy_0_xxyyy_1, g_yy_0_xxyyz_0, g_yy_0_xxyyz_1, g_yy_0_xxyzz_0, g_yy_0_xxyzz_1, g_yy_0_xyyyy_0, g_yy_0_xyyyy_1, g_yy_0_xyyyz_0, g_yy_0_xyyyz_1, g_yy_0_xyyzz_0, g_yy_0_xyyzz_1, g_yy_0_xyzzz_0, g_yy_0_xyzzz_1, g_yy_0_yyyyy_0, g_yy_0_yyyyy_1, g_yy_0_yyyyz_0, g_yy_0_yyyyz_1, g_yy_0_yyyzz_0, g_yy_0_yyyzz_1, g_yy_0_yyzzz_0, g_yy_0_yyzzz_1, g_yy_0_yzzzz_0, g_yy_0_yzzzz_1, g_yy_0_zzzzz_0, g_yy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxxxx_0[i] = g_xx_0_xxxxx_0[i] * fbe_0 - g_xx_0_xxxxx_1[i] * fz_be_0 + g_xxy_0_xxxxx_1[i] * wa_y[i];

        g_xxyy_0_xxxxy_0[i] = g_yy_0_xxxxy_0[i] * fbe_0 - g_yy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxy_1[i] * fi_acd_0 + g_xyy_0_xxxxy_1[i] * wa_x[i];

        g_xxyy_0_xxxxz_0[i] = g_xx_0_xxxxz_0[i] * fbe_0 - g_xx_0_xxxxz_1[i] * fz_be_0 + g_xxy_0_xxxxz_1[i] * wa_y[i];

        g_xxyy_0_xxxyy_0[i] = g_yy_0_xxxyy_0[i] * fbe_0 - g_yy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyy_1[i] * fi_acd_0 + g_xyy_0_xxxyy_1[i] * wa_x[i];

        g_xxyy_0_xxxyz_0[i] = g_yy_0_xxxyz_0[i] * fbe_0 - g_yy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyz_1[i] * fi_acd_0 + g_xyy_0_xxxyz_1[i] * wa_x[i];

        g_xxyy_0_xxxzz_0[i] = g_xx_0_xxxzz_0[i] * fbe_0 - g_xx_0_xxxzz_1[i] * fz_be_0 + g_xxy_0_xxxzz_1[i] * wa_y[i];

        g_xxyy_0_xxyyy_0[i] = g_yy_0_xxyyy_0[i] * fbe_0 - g_yy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyy_1[i] * fi_acd_0 + g_xyy_0_xxyyy_1[i] * wa_x[i];

        g_xxyy_0_xxyyz_0[i] = g_yy_0_xxyyz_0[i] * fbe_0 - g_yy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyz_1[i] * fi_acd_0 + g_xyy_0_xxyyz_1[i] * wa_x[i];

        g_xxyy_0_xxyzz_0[i] = g_yy_0_xxyzz_0[i] * fbe_0 - g_yy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyzz_1[i] * fi_acd_0 + g_xyy_0_xxyzz_1[i] * wa_x[i];

        g_xxyy_0_xxzzz_0[i] = g_xx_0_xxzzz_0[i] * fbe_0 - g_xx_0_xxzzz_1[i] * fz_be_0 + g_xxy_0_xxzzz_1[i] * wa_y[i];

        g_xxyy_0_xyyyy_0[i] = g_yy_0_xyyyy_0[i] * fbe_0 - g_yy_0_xyyyy_1[i] * fz_be_0 + g_xyy_0_yyyy_1[i] * fi_acd_0 + g_xyy_0_xyyyy_1[i] * wa_x[i];

        g_xxyy_0_xyyyz_0[i] = g_yy_0_xyyyz_0[i] * fbe_0 - g_yy_0_xyyyz_1[i] * fz_be_0 + g_xyy_0_yyyz_1[i] * fi_acd_0 + g_xyy_0_xyyyz_1[i] * wa_x[i];

        g_xxyy_0_xyyzz_0[i] = g_yy_0_xyyzz_0[i] * fbe_0 - g_yy_0_xyyzz_1[i] * fz_be_0 + g_xyy_0_yyzz_1[i] * fi_acd_0 + g_xyy_0_xyyzz_1[i] * wa_x[i];

        g_xxyy_0_xyzzz_0[i] = g_yy_0_xyzzz_0[i] * fbe_0 - g_yy_0_xyzzz_1[i] * fz_be_0 + g_xyy_0_yzzz_1[i] * fi_acd_0 + g_xyy_0_xyzzz_1[i] * wa_x[i];

        g_xxyy_0_xzzzz_0[i] = g_xx_0_xzzzz_0[i] * fbe_0 - g_xx_0_xzzzz_1[i] * fz_be_0 + g_xxy_0_xzzzz_1[i] * wa_y[i];

        g_xxyy_0_yyyyy_0[i] = g_yy_0_yyyyy_0[i] * fbe_0 - g_yy_0_yyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyy_1[i] * wa_x[i];

        g_xxyy_0_yyyyz_0[i] = g_yy_0_yyyyz_0[i] * fbe_0 - g_yy_0_yyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyz_1[i] * wa_x[i];

        g_xxyy_0_yyyzz_0[i] = g_yy_0_yyyzz_0[i] * fbe_0 - g_yy_0_yyyzz_1[i] * fz_be_0 + g_xyy_0_yyyzz_1[i] * wa_x[i];

        g_xxyy_0_yyzzz_0[i] = g_yy_0_yyzzz_0[i] * fbe_0 - g_yy_0_yyzzz_1[i] * fz_be_0 + g_xyy_0_yyzzz_1[i] * wa_x[i];

        g_xxyy_0_yzzzz_0[i] = g_yy_0_yzzzz_0[i] * fbe_0 - g_yy_0_yzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzz_1[i] * wa_x[i];

        g_xxyy_0_zzzzz_0[i] = g_yy_0_zzzzz_0[i] * fbe_0 - g_yy_0_zzzzz_1[i] * fz_be_0 + g_xyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 84-105 components of targeted buffer : GSH

    auto g_xxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 84);

    auto g_xxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 85);

    auto g_xxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 86);

    auto g_xxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 87);

    auto g_xxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 88);

    auto g_xxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 89);

    auto g_xxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 90);

    auto g_xxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 91);

    auto g_xxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 92);

    auto g_xxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 93);

    auto g_xxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 94);

    auto g_xxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 95);

    auto g_xxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 96);

    auto g_xxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 97);

    auto g_xxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 98);

    auto g_xxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 99);

    auto g_xxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 100);

    auto g_xxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 101);

    auto g_xxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 102);

    auto g_xxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 103);

    auto g_xxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 104);

    #pragma omp simd aligned(g_xxy_0_xxxxy_1, g_xxy_0_xxxyy_1, g_xxy_0_xxyyy_1, g_xxy_0_xyyyy_1, g_xxy_0_yyyyy_1, g_xxyz_0_xxxxx_0, g_xxyz_0_xxxxy_0, g_xxyz_0_xxxxz_0, g_xxyz_0_xxxyy_0, g_xxyz_0_xxxyz_0, g_xxyz_0_xxxzz_0, g_xxyz_0_xxyyy_0, g_xxyz_0_xxyyz_0, g_xxyz_0_xxyzz_0, g_xxyz_0_xxzzz_0, g_xxyz_0_xyyyy_0, g_xxyz_0_xyyyz_0, g_xxyz_0_xyyzz_0, g_xxyz_0_xyzzz_0, g_xxyz_0_xzzzz_0, g_xxyz_0_yyyyy_0, g_xxyz_0_yyyyz_0, g_xxyz_0_yyyzz_0, g_xxyz_0_yyzzz_0, g_xxyz_0_yzzzz_0, g_xxyz_0_zzzzz_0, g_xxz_0_xxxxx_1, g_xxz_0_xxxxz_1, g_xxz_0_xxxyz_1, g_xxz_0_xxxz_1, g_xxz_0_xxxzz_1, g_xxz_0_xxyyz_1, g_xxz_0_xxyz_1, g_xxz_0_xxyzz_1, g_xxz_0_xxzz_1, g_xxz_0_xxzzz_1, g_xxz_0_xyyyz_1, g_xxz_0_xyyz_1, g_xxz_0_xyyzz_1, g_xxz_0_xyzz_1, g_xxz_0_xyzzz_1, g_xxz_0_xzzz_1, g_xxz_0_xzzzz_1, g_xxz_0_yyyyz_1, g_xxz_0_yyyz_1, g_xxz_0_yyyzz_1, g_xxz_0_yyzz_1, g_xxz_0_yyzzz_1, g_xxz_0_yzzz_1, g_xxz_0_yzzzz_1, g_xxz_0_zzzz_1, g_xxz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxxxx_0[i] = g_xxz_0_xxxxx_1[i] * wa_y[i];

        g_xxyz_0_xxxxy_0[i] = g_xxy_0_xxxxy_1[i] * wa_z[i];

        g_xxyz_0_xxxxz_0[i] = g_xxz_0_xxxxz_1[i] * wa_y[i];

        g_xxyz_0_xxxyy_0[i] = g_xxy_0_xxxyy_1[i] * wa_z[i];

        g_xxyz_0_xxxyz_0[i] = g_xxz_0_xxxz_1[i] * fi_acd_0 + g_xxz_0_xxxyz_1[i] * wa_y[i];

        g_xxyz_0_xxxzz_0[i] = g_xxz_0_xxxzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyy_0[i] = g_xxy_0_xxyyy_1[i] * wa_z[i];

        g_xxyz_0_xxyyz_0[i] = 2.0 * g_xxz_0_xxyz_1[i] * fi_acd_0 + g_xxz_0_xxyyz_1[i] * wa_y[i];

        g_xxyz_0_xxyzz_0[i] = g_xxz_0_xxzz_1[i] * fi_acd_0 + g_xxz_0_xxyzz_1[i] * wa_y[i];

        g_xxyz_0_xxzzz_0[i] = g_xxz_0_xxzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyy_0[i] = g_xxy_0_xyyyy_1[i] * wa_z[i];

        g_xxyz_0_xyyyz_0[i] = 3.0 * g_xxz_0_xyyz_1[i] * fi_acd_0 + g_xxz_0_xyyyz_1[i] * wa_y[i];

        g_xxyz_0_xyyzz_0[i] = 2.0 * g_xxz_0_xyzz_1[i] * fi_acd_0 + g_xxz_0_xyyzz_1[i] * wa_y[i];

        g_xxyz_0_xyzzz_0[i] = g_xxz_0_xzzz_1[i] * fi_acd_0 + g_xxz_0_xyzzz_1[i] * wa_y[i];

        g_xxyz_0_xzzzz_0[i] = g_xxz_0_xzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyy_0[i] = g_xxy_0_yyyyy_1[i] * wa_z[i];

        g_xxyz_0_yyyyz_0[i] = 4.0 * g_xxz_0_yyyz_1[i] * fi_acd_0 + g_xxz_0_yyyyz_1[i] * wa_y[i];

        g_xxyz_0_yyyzz_0[i] = 3.0 * g_xxz_0_yyzz_1[i] * fi_acd_0 + g_xxz_0_yyyzz_1[i] * wa_y[i];

        g_xxyz_0_yyzzz_0[i] = 2.0 * g_xxz_0_yzzz_1[i] * fi_acd_0 + g_xxz_0_yyzzz_1[i] * wa_y[i];

        g_xxyz_0_yzzzz_0[i] = g_xxz_0_zzzz_1[i] * fi_acd_0 + g_xxz_0_yzzzz_1[i] * wa_y[i];

        g_xxyz_0_zzzzz_0[i] = g_xxz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 105-126 components of targeted buffer : GSH

    auto g_xxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 105);

    auto g_xxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 106);

    auto g_xxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 107);

    auto g_xxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 108);

    auto g_xxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 109);

    auto g_xxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 110);

    auto g_xxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 111);

    auto g_xxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 112);

    auto g_xxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 113);

    auto g_xxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 114);

    auto g_xxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 115);

    auto g_xxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 116);

    auto g_xxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 117);

    auto g_xxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 118);

    auto g_xxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 119);

    auto g_xxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 120);

    auto g_xxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 121);

    auto g_xxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 122);

    auto g_xxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 123);

    auto g_xxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 124);

    auto g_xxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 125);

    #pragma omp simd aligned(g_xx_0_xxxxx_0, g_xx_0_xxxxx_1, g_xx_0_xxxxy_0, g_xx_0_xxxxy_1, g_xx_0_xxxyy_0, g_xx_0_xxxyy_1, g_xx_0_xxyyy_0, g_xx_0_xxyyy_1, g_xx_0_xyyyy_0, g_xx_0_xyyyy_1, g_xxz_0_xxxxx_1, g_xxz_0_xxxxy_1, g_xxz_0_xxxyy_1, g_xxz_0_xxyyy_1, g_xxz_0_xyyyy_1, g_xxzz_0_xxxxx_0, g_xxzz_0_xxxxy_0, g_xxzz_0_xxxxz_0, g_xxzz_0_xxxyy_0, g_xxzz_0_xxxyz_0, g_xxzz_0_xxxzz_0, g_xxzz_0_xxyyy_0, g_xxzz_0_xxyyz_0, g_xxzz_0_xxyzz_0, g_xxzz_0_xxzzz_0, g_xxzz_0_xyyyy_0, g_xxzz_0_xyyyz_0, g_xxzz_0_xyyzz_0, g_xxzz_0_xyzzz_0, g_xxzz_0_xzzzz_0, g_xxzz_0_yyyyy_0, g_xxzz_0_yyyyz_0, g_xxzz_0_yyyzz_0, g_xxzz_0_yyzzz_0, g_xxzz_0_yzzzz_0, g_xxzz_0_zzzzz_0, g_xzz_0_xxxxz_1, g_xzz_0_xxxyz_1, g_xzz_0_xxxz_1, g_xzz_0_xxxzz_1, g_xzz_0_xxyyz_1, g_xzz_0_xxyz_1, g_xzz_0_xxyzz_1, g_xzz_0_xxzz_1, g_xzz_0_xxzzz_1, g_xzz_0_xyyyz_1, g_xzz_0_xyyz_1, g_xzz_0_xyyzz_1, g_xzz_0_xyzz_1, g_xzz_0_xyzzz_1, g_xzz_0_xzzz_1, g_xzz_0_xzzzz_1, g_xzz_0_yyyyy_1, g_xzz_0_yyyyz_1, g_xzz_0_yyyz_1, g_xzz_0_yyyzz_1, g_xzz_0_yyzz_1, g_xzz_0_yyzzz_1, g_xzz_0_yzzz_1, g_xzz_0_yzzzz_1, g_xzz_0_zzzz_1, g_xzz_0_zzzzz_1, g_zz_0_xxxxz_0, g_zz_0_xxxxz_1, g_zz_0_xxxyz_0, g_zz_0_xxxyz_1, g_zz_0_xxxzz_0, g_zz_0_xxxzz_1, g_zz_0_xxyyz_0, g_zz_0_xxyyz_1, g_zz_0_xxyzz_0, g_zz_0_xxyzz_1, g_zz_0_xxzzz_0, g_zz_0_xxzzz_1, g_zz_0_xyyyz_0, g_zz_0_xyyyz_1, g_zz_0_xyyzz_0, g_zz_0_xyyzz_1, g_zz_0_xyzzz_0, g_zz_0_xyzzz_1, g_zz_0_xzzzz_0, g_zz_0_xzzzz_1, g_zz_0_yyyyy_0, g_zz_0_yyyyy_1, g_zz_0_yyyyz_0, g_zz_0_yyyyz_1, g_zz_0_yyyzz_0, g_zz_0_yyyzz_1, g_zz_0_yyzzz_0, g_zz_0_yyzzz_1, g_zz_0_yzzzz_0, g_zz_0_yzzzz_1, g_zz_0_zzzzz_0, g_zz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxxxx_0[i] = g_xx_0_xxxxx_0[i] * fbe_0 - g_xx_0_xxxxx_1[i] * fz_be_0 + g_xxz_0_xxxxx_1[i] * wa_z[i];

        g_xxzz_0_xxxxy_0[i] = g_xx_0_xxxxy_0[i] * fbe_0 - g_xx_0_xxxxy_1[i] * fz_be_0 + g_xxz_0_xxxxy_1[i] * wa_z[i];

        g_xxzz_0_xxxxz_0[i] = g_zz_0_xxxxz_0[i] * fbe_0 - g_zz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxz_1[i] * fi_acd_0 + g_xzz_0_xxxxz_1[i] * wa_x[i];

        g_xxzz_0_xxxyy_0[i] = g_xx_0_xxxyy_0[i] * fbe_0 - g_xx_0_xxxyy_1[i] * fz_be_0 + g_xxz_0_xxxyy_1[i] * wa_z[i];

        g_xxzz_0_xxxyz_0[i] = g_zz_0_xxxyz_0[i] * fbe_0 - g_zz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyz_1[i] * fi_acd_0 + g_xzz_0_xxxyz_1[i] * wa_x[i];

        g_xxzz_0_xxxzz_0[i] = g_zz_0_xxxzz_0[i] * fbe_0 - g_zz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxzz_1[i] * fi_acd_0 + g_xzz_0_xxxzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyy_0[i] = g_xx_0_xxyyy_0[i] * fbe_0 - g_xx_0_xxyyy_1[i] * fz_be_0 + g_xxz_0_xxyyy_1[i] * wa_z[i];

        g_xxzz_0_xxyyz_0[i] = g_zz_0_xxyyz_0[i] * fbe_0 - g_zz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyz_1[i] * fi_acd_0 + g_xzz_0_xxyyz_1[i] * wa_x[i];

        g_xxzz_0_xxyzz_0[i] = g_zz_0_xxyzz_0[i] * fbe_0 - g_zz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyzz_1[i] * fi_acd_0 + g_xzz_0_xxyzz_1[i] * wa_x[i];

        g_xxzz_0_xxzzz_0[i] = g_zz_0_xxzzz_0[i] * fbe_0 - g_zz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xzzz_1[i] * fi_acd_0 + g_xzz_0_xxzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyy_0[i] = g_xx_0_xyyyy_0[i] * fbe_0 - g_xx_0_xyyyy_1[i] * fz_be_0 + g_xxz_0_xyyyy_1[i] * wa_z[i];

        g_xxzz_0_xyyyz_0[i] = g_zz_0_xyyyz_0[i] * fbe_0 - g_zz_0_xyyyz_1[i] * fz_be_0 + g_xzz_0_yyyz_1[i] * fi_acd_0 + g_xzz_0_xyyyz_1[i] * wa_x[i];

        g_xxzz_0_xyyzz_0[i] = g_zz_0_xyyzz_0[i] * fbe_0 - g_zz_0_xyyzz_1[i] * fz_be_0 + g_xzz_0_yyzz_1[i] * fi_acd_0 + g_xzz_0_xyyzz_1[i] * wa_x[i];

        g_xxzz_0_xyzzz_0[i] = g_zz_0_xyzzz_0[i] * fbe_0 - g_zz_0_xyzzz_1[i] * fz_be_0 + g_xzz_0_yzzz_1[i] * fi_acd_0 + g_xzz_0_xyzzz_1[i] * wa_x[i];

        g_xxzz_0_xzzzz_0[i] = g_zz_0_xzzzz_0[i] * fbe_0 - g_zz_0_xzzzz_1[i] * fz_be_0 + g_xzz_0_zzzz_1[i] * fi_acd_0 + g_xzz_0_xzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyy_0[i] = g_zz_0_yyyyy_0[i] * fbe_0 - g_zz_0_yyyyy_1[i] * fz_be_0 + g_xzz_0_yyyyy_1[i] * wa_x[i];

        g_xxzz_0_yyyyz_0[i] = g_zz_0_yyyyz_0[i] * fbe_0 - g_zz_0_yyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyz_1[i] * wa_x[i];

        g_xxzz_0_yyyzz_0[i] = g_zz_0_yyyzz_0[i] * fbe_0 - g_zz_0_yyyzz_1[i] * fz_be_0 + g_xzz_0_yyyzz_1[i] * wa_x[i];

        g_xxzz_0_yyzzz_0[i] = g_zz_0_yyzzz_0[i] * fbe_0 - g_zz_0_yyzzz_1[i] * fz_be_0 + g_xzz_0_yyzzz_1[i] * wa_x[i];

        g_xxzz_0_yzzzz_0[i] = g_zz_0_yzzzz_0[i] * fbe_0 - g_zz_0_yzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzz_1[i] * wa_x[i];

        g_xxzz_0_zzzzz_0[i] = g_zz_0_zzzzz_0[i] * fbe_0 - g_zz_0_zzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 126-147 components of targeted buffer : GSH

    auto g_xyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 126);

    auto g_xyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 127);

    auto g_xyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 128);

    auto g_xyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 129);

    auto g_xyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 130);

    auto g_xyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 131);

    auto g_xyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 132);

    auto g_xyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 133);

    auto g_xyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 134);

    auto g_xyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 135);

    auto g_xyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 136);

    auto g_xyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 137);

    auto g_xyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 138);

    auto g_xyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 139);

    auto g_xyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 140);

    auto g_xyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 141);

    auto g_xyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 142);

    auto g_xyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 143);

    auto g_xyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 144);

    auto g_xyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 145);

    auto g_xyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 146);

    #pragma omp simd aligned(g_xyyy_0_xxxxx_0, g_xyyy_0_xxxxy_0, g_xyyy_0_xxxxz_0, g_xyyy_0_xxxyy_0, g_xyyy_0_xxxyz_0, g_xyyy_0_xxxzz_0, g_xyyy_0_xxyyy_0, g_xyyy_0_xxyyz_0, g_xyyy_0_xxyzz_0, g_xyyy_0_xxzzz_0, g_xyyy_0_xyyyy_0, g_xyyy_0_xyyyz_0, g_xyyy_0_xyyzz_0, g_xyyy_0_xyzzz_0, g_xyyy_0_xzzzz_0, g_xyyy_0_yyyyy_0, g_xyyy_0_yyyyz_0, g_xyyy_0_yyyzz_0, g_xyyy_0_yyzzz_0, g_xyyy_0_yzzzz_0, g_xyyy_0_zzzzz_0, g_yyy_0_xxxx_1, g_yyy_0_xxxxx_1, g_yyy_0_xxxxy_1, g_yyy_0_xxxxz_1, g_yyy_0_xxxy_1, g_yyy_0_xxxyy_1, g_yyy_0_xxxyz_1, g_yyy_0_xxxz_1, g_yyy_0_xxxzz_1, g_yyy_0_xxyy_1, g_yyy_0_xxyyy_1, g_yyy_0_xxyyz_1, g_yyy_0_xxyz_1, g_yyy_0_xxyzz_1, g_yyy_0_xxzz_1, g_yyy_0_xxzzz_1, g_yyy_0_xyyy_1, g_yyy_0_xyyyy_1, g_yyy_0_xyyyz_1, g_yyy_0_xyyz_1, g_yyy_0_xyyzz_1, g_yyy_0_xyzz_1, g_yyy_0_xyzzz_1, g_yyy_0_xzzz_1, g_yyy_0_xzzzz_1, g_yyy_0_yyyy_1, g_yyy_0_yyyyy_1, g_yyy_0_yyyyz_1, g_yyy_0_yyyz_1, g_yyy_0_yyyzz_1, g_yyy_0_yyzz_1, g_yyy_0_yyzzz_1, g_yyy_0_yzzz_1, g_yyy_0_yzzzz_1, g_yyy_0_zzzz_1, g_yyy_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxxxx_0[i] = 5.0 * g_yyy_0_xxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxx_1[i] * wa_x[i];

        g_xyyy_0_xxxxy_0[i] = 4.0 * g_yyy_0_xxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxy_1[i] * wa_x[i];

        g_xyyy_0_xxxxz_0[i] = 4.0 * g_yyy_0_xxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxz_1[i] * wa_x[i];

        g_xyyy_0_xxxyy_0[i] = 3.0 * g_yyy_0_xxyy_1[i] * fi_acd_0 + g_yyy_0_xxxyy_1[i] * wa_x[i];

        g_xyyy_0_xxxyz_0[i] = 3.0 * g_yyy_0_xxyz_1[i] * fi_acd_0 + g_yyy_0_xxxyz_1[i] * wa_x[i];

        g_xyyy_0_xxxzz_0[i] = 3.0 * g_yyy_0_xxzz_1[i] * fi_acd_0 + g_yyy_0_xxxzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyy_0[i] = 2.0 * g_yyy_0_xyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyy_1[i] * wa_x[i];

        g_xyyy_0_xxyyz_0[i] = 2.0 * g_yyy_0_xyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyz_1[i] * wa_x[i];

        g_xyyy_0_xxyzz_0[i] = 2.0 * g_yyy_0_xyzz_1[i] * fi_acd_0 + g_yyy_0_xxyzz_1[i] * wa_x[i];

        g_xyyy_0_xxzzz_0[i] = 2.0 * g_yyy_0_xzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyy_0[i] = g_yyy_0_yyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyy_1[i] * wa_x[i];

        g_xyyy_0_xyyyz_0[i] = g_yyy_0_yyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyz_1[i] * wa_x[i];

        g_xyyy_0_xyyzz_0[i] = g_yyy_0_yyzz_1[i] * fi_acd_0 + g_yyy_0_xyyzz_1[i] * wa_x[i];

        g_xyyy_0_xyzzz_0[i] = g_yyy_0_yzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzz_1[i] * wa_x[i];

        g_xyyy_0_xzzzz_0[i] = g_yyy_0_zzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyy_0[i] = g_yyy_0_yyyyy_1[i] * wa_x[i];

        g_xyyy_0_yyyyz_0[i] = g_yyy_0_yyyyz_1[i] * wa_x[i];

        g_xyyy_0_yyyzz_0[i] = g_yyy_0_yyyzz_1[i] * wa_x[i];

        g_xyyy_0_yyzzz_0[i] = g_yyy_0_yyzzz_1[i] * wa_x[i];

        g_xyyy_0_yzzzz_0[i] = g_yyy_0_yzzzz_1[i] * wa_x[i];

        g_xyyy_0_zzzzz_0[i] = g_yyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 147-168 components of targeted buffer : GSH

    auto g_xyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 147);

    auto g_xyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 148);

    auto g_xyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 149);

    auto g_xyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 150);

    auto g_xyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 151);

    auto g_xyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 152);

    auto g_xyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 153);

    auto g_xyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 154);

    auto g_xyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 155);

    auto g_xyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 156);

    auto g_xyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 157);

    auto g_xyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 158);

    auto g_xyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 159);

    auto g_xyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 160);

    auto g_xyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 161);

    auto g_xyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 162);

    auto g_xyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 163);

    auto g_xyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 164);

    auto g_xyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 165);

    auto g_xyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 166);

    auto g_xyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 167);

    #pragma omp simd aligned(g_xyy_0_xxxxx_1, g_xyy_0_xxxxy_1, g_xyy_0_xxxyy_1, g_xyy_0_xxyyy_1, g_xyy_0_xyyyy_1, g_xyyz_0_xxxxx_0, g_xyyz_0_xxxxy_0, g_xyyz_0_xxxxz_0, g_xyyz_0_xxxyy_0, g_xyyz_0_xxxyz_0, g_xyyz_0_xxxzz_0, g_xyyz_0_xxyyy_0, g_xyyz_0_xxyyz_0, g_xyyz_0_xxyzz_0, g_xyyz_0_xxzzz_0, g_xyyz_0_xyyyy_0, g_xyyz_0_xyyyz_0, g_xyyz_0_xyyzz_0, g_xyyz_0_xyzzz_0, g_xyyz_0_xzzzz_0, g_xyyz_0_yyyyy_0, g_xyyz_0_yyyyz_0, g_xyyz_0_yyyzz_0, g_xyyz_0_yyzzz_0, g_xyyz_0_yzzzz_0, g_xyyz_0_zzzzz_0, g_yyz_0_xxxxz_1, g_yyz_0_xxxyz_1, g_yyz_0_xxxz_1, g_yyz_0_xxxzz_1, g_yyz_0_xxyyz_1, g_yyz_0_xxyz_1, g_yyz_0_xxyzz_1, g_yyz_0_xxzz_1, g_yyz_0_xxzzz_1, g_yyz_0_xyyyz_1, g_yyz_0_xyyz_1, g_yyz_0_xyyzz_1, g_yyz_0_xyzz_1, g_yyz_0_xyzzz_1, g_yyz_0_xzzz_1, g_yyz_0_xzzzz_1, g_yyz_0_yyyyy_1, g_yyz_0_yyyyz_1, g_yyz_0_yyyz_1, g_yyz_0_yyyzz_1, g_yyz_0_yyzz_1, g_yyz_0_yyzzz_1, g_yyz_0_yzzz_1, g_yyz_0_yzzzz_1, g_yyz_0_zzzz_1, g_yyz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxxxx_0[i] = g_xyy_0_xxxxx_1[i] * wa_z[i];

        g_xyyz_0_xxxxy_0[i] = g_xyy_0_xxxxy_1[i] * wa_z[i];

        g_xyyz_0_xxxxz_0[i] = 4.0 * g_yyz_0_xxxz_1[i] * fi_acd_0 + g_yyz_0_xxxxz_1[i] * wa_x[i];

        g_xyyz_0_xxxyy_0[i] = g_xyy_0_xxxyy_1[i] * wa_z[i];

        g_xyyz_0_xxxyz_0[i] = 3.0 * g_yyz_0_xxyz_1[i] * fi_acd_0 + g_yyz_0_xxxyz_1[i] * wa_x[i];

        g_xyyz_0_xxxzz_0[i] = 3.0 * g_yyz_0_xxzz_1[i] * fi_acd_0 + g_yyz_0_xxxzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyy_0[i] = g_xyy_0_xxyyy_1[i] * wa_z[i];

        g_xyyz_0_xxyyz_0[i] = 2.0 * g_yyz_0_xyyz_1[i] * fi_acd_0 + g_yyz_0_xxyyz_1[i] * wa_x[i];

        g_xyyz_0_xxyzz_0[i] = 2.0 * g_yyz_0_xyzz_1[i] * fi_acd_0 + g_yyz_0_xxyzz_1[i] * wa_x[i];

        g_xyyz_0_xxzzz_0[i] = 2.0 * g_yyz_0_xzzz_1[i] * fi_acd_0 + g_yyz_0_xxzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyy_0[i] = g_xyy_0_xyyyy_1[i] * wa_z[i];

        g_xyyz_0_xyyyz_0[i] = g_yyz_0_yyyz_1[i] * fi_acd_0 + g_yyz_0_xyyyz_1[i] * wa_x[i];

        g_xyyz_0_xyyzz_0[i] = g_yyz_0_yyzz_1[i] * fi_acd_0 + g_yyz_0_xyyzz_1[i] * wa_x[i];

        g_xyyz_0_xyzzz_0[i] = g_yyz_0_yzzz_1[i] * fi_acd_0 + g_yyz_0_xyzzz_1[i] * wa_x[i];

        g_xyyz_0_xzzzz_0[i] = g_yyz_0_zzzz_1[i] * fi_acd_0 + g_yyz_0_xzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyy_0[i] = g_yyz_0_yyyyy_1[i] * wa_x[i];

        g_xyyz_0_yyyyz_0[i] = g_yyz_0_yyyyz_1[i] * wa_x[i];

        g_xyyz_0_yyyzz_0[i] = g_yyz_0_yyyzz_1[i] * wa_x[i];

        g_xyyz_0_yyzzz_0[i] = g_yyz_0_yyzzz_1[i] * wa_x[i];

        g_xyyz_0_yzzzz_0[i] = g_yyz_0_yzzzz_1[i] * wa_x[i];

        g_xyyz_0_zzzzz_0[i] = g_yyz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 168-189 components of targeted buffer : GSH

    auto g_xyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 168);

    auto g_xyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 169);

    auto g_xyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 170);

    auto g_xyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 171);

    auto g_xyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 172);

    auto g_xyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 173);

    auto g_xyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 174);

    auto g_xyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 175);

    auto g_xyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 176);

    auto g_xyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 177);

    auto g_xyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 178);

    auto g_xyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 179);

    auto g_xyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 180);

    auto g_xyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 181);

    auto g_xyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 182);

    auto g_xyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 183);

    auto g_xyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 184);

    auto g_xyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 185);

    auto g_xyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 186);

    auto g_xyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 187);

    auto g_xyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 188);

    #pragma omp simd aligned(g_xyzz_0_xxxxx_0, g_xyzz_0_xxxxy_0, g_xyzz_0_xxxxz_0, g_xyzz_0_xxxyy_0, g_xyzz_0_xxxyz_0, g_xyzz_0_xxxzz_0, g_xyzz_0_xxyyy_0, g_xyzz_0_xxyyz_0, g_xyzz_0_xxyzz_0, g_xyzz_0_xxzzz_0, g_xyzz_0_xyyyy_0, g_xyzz_0_xyyyz_0, g_xyzz_0_xyyzz_0, g_xyzz_0_xyzzz_0, g_xyzz_0_xzzzz_0, g_xyzz_0_yyyyy_0, g_xyzz_0_yyyyz_0, g_xyzz_0_yyyzz_0, g_xyzz_0_yyzzz_0, g_xyzz_0_yzzzz_0, g_xyzz_0_zzzzz_0, g_xzz_0_xxxxx_1, g_xzz_0_xxxxz_1, g_xzz_0_xxxzz_1, g_xzz_0_xxzzz_1, g_xzz_0_xzzzz_1, g_yzz_0_xxxxy_1, g_yzz_0_xxxy_1, g_yzz_0_xxxyy_1, g_yzz_0_xxxyz_1, g_yzz_0_xxyy_1, g_yzz_0_xxyyy_1, g_yzz_0_xxyyz_1, g_yzz_0_xxyz_1, g_yzz_0_xxyzz_1, g_yzz_0_xyyy_1, g_yzz_0_xyyyy_1, g_yzz_0_xyyyz_1, g_yzz_0_xyyz_1, g_yzz_0_xyyzz_1, g_yzz_0_xyzz_1, g_yzz_0_xyzzz_1, g_yzz_0_yyyy_1, g_yzz_0_yyyyy_1, g_yzz_0_yyyyz_1, g_yzz_0_yyyz_1, g_yzz_0_yyyzz_1, g_yzz_0_yyzz_1, g_yzz_0_yyzzz_1, g_yzz_0_yzzz_1, g_yzz_0_yzzzz_1, g_yzz_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxxxx_0[i] = g_xzz_0_xxxxx_1[i] * wa_y[i];

        g_xyzz_0_xxxxy_0[i] = 4.0 * g_yzz_0_xxxy_1[i] * fi_acd_0 + g_yzz_0_xxxxy_1[i] * wa_x[i];

        g_xyzz_0_xxxxz_0[i] = g_xzz_0_xxxxz_1[i] * wa_y[i];

        g_xyzz_0_xxxyy_0[i] = 3.0 * g_yzz_0_xxyy_1[i] * fi_acd_0 + g_yzz_0_xxxyy_1[i] * wa_x[i];

        g_xyzz_0_xxxyz_0[i] = 3.0 * g_yzz_0_xxyz_1[i] * fi_acd_0 + g_yzz_0_xxxyz_1[i] * wa_x[i];

        g_xyzz_0_xxxzz_0[i] = g_xzz_0_xxxzz_1[i] * wa_y[i];

        g_xyzz_0_xxyyy_0[i] = 2.0 * g_yzz_0_xyyy_1[i] * fi_acd_0 + g_yzz_0_xxyyy_1[i] * wa_x[i];

        g_xyzz_0_xxyyz_0[i] = 2.0 * g_yzz_0_xyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyz_1[i] * wa_x[i];

        g_xyzz_0_xxyzz_0[i] = 2.0 * g_yzz_0_xyzz_1[i] * fi_acd_0 + g_yzz_0_xxyzz_1[i] * wa_x[i];

        g_xyzz_0_xxzzz_0[i] = g_xzz_0_xxzzz_1[i] * wa_y[i];

        g_xyzz_0_xyyyy_0[i] = g_yzz_0_yyyy_1[i] * fi_acd_0 + g_yzz_0_xyyyy_1[i] * wa_x[i];

        g_xyzz_0_xyyyz_0[i] = g_yzz_0_yyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyz_1[i] * wa_x[i];

        g_xyzz_0_xyyzz_0[i] = g_yzz_0_yyzz_1[i] * fi_acd_0 + g_yzz_0_xyyzz_1[i] * wa_x[i];

        g_xyzz_0_xyzzz_0[i] = g_yzz_0_yzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzz_1[i] * wa_x[i];

        g_xyzz_0_xzzzz_0[i] = g_xzz_0_xzzzz_1[i] * wa_y[i];

        g_xyzz_0_yyyyy_0[i] = g_yzz_0_yyyyy_1[i] * wa_x[i];

        g_xyzz_0_yyyyz_0[i] = g_yzz_0_yyyyz_1[i] * wa_x[i];

        g_xyzz_0_yyyzz_0[i] = g_yzz_0_yyyzz_1[i] * wa_x[i];

        g_xyzz_0_yyzzz_0[i] = g_yzz_0_yyzzz_1[i] * wa_x[i];

        g_xyzz_0_yzzzz_0[i] = g_yzz_0_yzzzz_1[i] * wa_x[i];

        g_xyzz_0_zzzzz_0[i] = g_yzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 189-210 components of targeted buffer : GSH

    auto g_xzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 189);

    auto g_xzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 190);

    auto g_xzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 191);

    auto g_xzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 192);

    auto g_xzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 193);

    auto g_xzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 194);

    auto g_xzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 195);

    auto g_xzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 196);

    auto g_xzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 197);

    auto g_xzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 198);

    auto g_xzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 199);

    auto g_xzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 200);

    auto g_xzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 201);

    auto g_xzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 202);

    auto g_xzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 203);

    auto g_xzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 204);

    auto g_xzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 205);

    auto g_xzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 206);

    auto g_xzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 207);

    auto g_xzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 208);

    auto g_xzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 209);

    #pragma omp simd aligned(g_xzzz_0_xxxxx_0, g_xzzz_0_xxxxy_0, g_xzzz_0_xxxxz_0, g_xzzz_0_xxxyy_0, g_xzzz_0_xxxyz_0, g_xzzz_0_xxxzz_0, g_xzzz_0_xxyyy_0, g_xzzz_0_xxyyz_0, g_xzzz_0_xxyzz_0, g_xzzz_0_xxzzz_0, g_xzzz_0_xyyyy_0, g_xzzz_0_xyyyz_0, g_xzzz_0_xyyzz_0, g_xzzz_0_xyzzz_0, g_xzzz_0_xzzzz_0, g_xzzz_0_yyyyy_0, g_xzzz_0_yyyyz_0, g_xzzz_0_yyyzz_0, g_xzzz_0_yyzzz_0, g_xzzz_0_yzzzz_0, g_xzzz_0_zzzzz_0, g_zzz_0_xxxx_1, g_zzz_0_xxxxx_1, g_zzz_0_xxxxy_1, g_zzz_0_xxxxz_1, g_zzz_0_xxxy_1, g_zzz_0_xxxyy_1, g_zzz_0_xxxyz_1, g_zzz_0_xxxz_1, g_zzz_0_xxxzz_1, g_zzz_0_xxyy_1, g_zzz_0_xxyyy_1, g_zzz_0_xxyyz_1, g_zzz_0_xxyz_1, g_zzz_0_xxyzz_1, g_zzz_0_xxzz_1, g_zzz_0_xxzzz_1, g_zzz_0_xyyy_1, g_zzz_0_xyyyy_1, g_zzz_0_xyyyz_1, g_zzz_0_xyyz_1, g_zzz_0_xyyzz_1, g_zzz_0_xyzz_1, g_zzz_0_xyzzz_1, g_zzz_0_xzzz_1, g_zzz_0_xzzzz_1, g_zzz_0_yyyy_1, g_zzz_0_yyyyy_1, g_zzz_0_yyyyz_1, g_zzz_0_yyyz_1, g_zzz_0_yyyzz_1, g_zzz_0_yyzz_1, g_zzz_0_yyzzz_1, g_zzz_0_yzzz_1, g_zzz_0_yzzzz_1, g_zzz_0_zzzz_1, g_zzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxxxx_0[i] = 5.0 * g_zzz_0_xxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxx_1[i] * wa_x[i];

        g_xzzz_0_xxxxy_0[i] = 4.0 * g_zzz_0_xxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxy_1[i] * wa_x[i];

        g_xzzz_0_xxxxz_0[i] = 4.0 * g_zzz_0_xxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxz_1[i] * wa_x[i];

        g_xzzz_0_xxxyy_0[i] = 3.0 * g_zzz_0_xxyy_1[i] * fi_acd_0 + g_zzz_0_xxxyy_1[i] * wa_x[i];

        g_xzzz_0_xxxyz_0[i] = 3.0 * g_zzz_0_xxyz_1[i] * fi_acd_0 + g_zzz_0_xxxyz_1[i] * wa_x[i];

        g_xzzz_0_xxxzz_0[i] = 3.0 * g_zzz_0_xxzz_1[i] * fi_acd_0 + g_zzz_0_xxxzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyy_0[i] = 2.0 * g_zzz_0_xyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyy_1[i] * wa_x[i];

        g_xzzz_0_xxyyz_0[i] = 2.0 * g_zzz_0_xyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyz_1[i] * wa_x[i];

        g_xzzz_0_xxyzz_0[i] = 2.0 * g_zzz_0_xyzz_1[i] * fi_acd_0 + g_zzz_0_xxyzz_1[i] * wa_x[i];

        g_xzzz_0_xxzzz_0[i] = 2.0 * g_zzz_0_xzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyy_0[i] = g_zzz_0_yyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyy_1[i] * wa_x[i];

        g_xzzz_0_xyyyz_0[i] = g_zzz_0_yyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyz_1[i] * wa_x[i];

        g_xzzz_0_xyyzz_0[i] = g_zzz_0_yyzz_1[i] * fi_acd_0 + g_zzz_0_xyyzz_1[i] * wa_x[i];

        g_xzzz_0_xyzzz_0[i] = g_zzz_0_yzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzz_1[i] * wa_x[i];

        g_xzzz_0_xzzzz_0[i] = g_zzz_0_zzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyy_0[i] = g_zzz_0_yyyyy_1[i] * wa_x[i];

        g_xzzz_0_yyyyz_0[i] = g_zzz_0_yyyyz_1[i] * wa_x[i];

        g_xzzz_0_yyyzz_0[i] = g_zzz_0_yyyzz_1[i] * wa_x[i];

        g_xzzz_0_yyzzz_0[i] = g_zzz_0_yyzzz_1[i] * wa_x[i];

        g_xzzz_0_yzzzz_0[i] = g_zzz_0_yzzzz_1[i] * wa_x[i];

        g_xzzz_0_zzzzz_0[i] = g_zzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 210-231 components of targeted buffer : GSH

    auto g_yyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 210);

    auto g_yyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 211);

    auto g_yyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 212);

    auto g_yyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 213);

    auto g_yyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 214);

    auto g_yyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 215);

    auto g_yyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 216);

    auto g_yyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 217);

    auto g_yyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 218);

    auto g_yyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 219);

    auto g_yyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 220);

    auto g_yyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 221);

    auto g_yyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 222);

    auto g_yyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 223);

    auto g_yyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 224);

    auto g_yyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 225);

    auto g_yyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 226);

    auto g_yyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 227);

    auto g_yyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 228);

    auto g_yyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 229);

    auto g_yyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 230);

    #pragma omp simd aligned(g_yy_0_xxxxx_0, g_yy_0_xxxxx_1, g_yy_0_xxxxy_0, g_yy_0_xxxxy_1, g_yy_0_xxxxz_0, g_yy_0_xxxxz_1, g_yy_0_xxxyy_0, g_yy_0_xxxyy_1, g_yy_0_xxxyz_0, g_yy_0_xxxyz_1, g_yy_0_xxxzz_0, g_yy_0_xxxzz_1, g_yy_0_xxyyy_0, g_yy_0_xxyyy_1, g_yy_0_xxyyz_0, g_yy_0_xxyyz_1, g_yy_0_xxyzz_0, g_yy_0_xxyzz_1, g_yy_0_xxzzz_0, g_yy_0_xxzzz_1, g_yy_0_xyyyy_0, g_yy_0_xyyyy_1, g_yy_0_xyyyz_0, g_yy_0_xyyyz_1, g_yy_0_xyyzz_0, g_yy_0_xyyzz_1, g_yy_0_xyzzz_0, g_yy_0_xyzzz_1, g_yy_0_xzzzz_0, g_yy_0_xzzzz_1, g_yy_0_yyyyy_0, g_yy_0_yyyyy_1, g_yy_0_yyyyz_0, g_yy_0_yyyyz_1, g_yy_0_yyyzz_0, g_yy_0_yyyzz_1, g_yy_0_yyzzz_0, g_yy_0_yyzzz_1, g_yy_0_yzzzz_0, g_yy_0_yzzzz_1, g_yy_0_zzzzz_0, g_yy_0_zzzzz_1, g_yyy_0_xxxx_1, g_yyy_0_xxxxx_1, g_yyy_0_xxxxy_1, g_yyy_0_xxxxz_1, g_yyy_0_xxxy_1, g_yyy_0_xxxyy_1, g_yyy_0_xxxyz_1, g_yyy_0_xxxz_1, g_yyy_0_xxxzz_1, g_yyy_0_xxyy_1, g_yyy_0_xxyyy_1, g_yyy_0_xxyyz_1, g_yyy_0_xxyz_1, g_yyy_0_xxyzz_1, g_yyy_0_xxzz_1, g_yyy_0_xxzzz_1, g_yyy_0_xyyy_1, g_yyy_0_xyyyy_1, g_yyy_0_xyyyz_1, g_yyy_0_xyyz_1, g_yyy_0_xyyzz_1, g_yyy_0_xyzz_1, g_yyy_0_xyzzz_1, g_yyy_0_xzzz_1, g_yyy_0_xzzzz_1, g_yyy_0_yyyy_1, g_yyy_0_yyyyy_1, g_yyy_0_yyyyz_1, g_yyy_0_yyyz_1, g_yyy_0_yyyzz_1, g_yyy_0_yyzz_1, g_yyy_0_yyzzz_1, g_yyy_0_yzzz_1, g_yyy_0_yzzzz_1, g_yyy_0_zzzz_1, g_yyy_0_zzzzz_1, g_yyyy_0_xxxxx_0, g_yyyy_0_xxxxy_0, g_yyyy_0_xxxxz_0, g_yyyy_0_xxxyy_0, g_yyyy_0_xxxyz_0, g_yyyy_0_xxxzz_0, g_yyyy_0_xxyyy_0, g_yyyy_0_xxyyz_0, g_yyyy_0_xxyzz_0, g_yyyy_0_xxzzz_0, g_yyyy_0_xyyyy_0, g_yyyy_0_xyyyz_0, g_yyyy_0_xyyzz_0, g_yyyy_0_xyzzz_0, g_yyyy_0_xzzzz_0, g_yyyy_0_yyyyy_0, g_yyyy_0_yyyyz_0, g_yyyy_0_yyyzz_0, g_yyyy_0_yyzzz_0, g_yyyy_0_yzzzz_0, g_yyyy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxxxx_0[i] = 3.0 * g_yy_0_xxxxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxx_1[i] * fz_be_0 + g_yyy_0_xxxxx_1[i] * wa_y[i];

        g_yyyy_0_xxxxy_0[i] = 3.0 * g_yy_0_xxxxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxy_1[i] * fz_be_0 + g_yyy_0_xxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxy_1[i] * wa_y[i];

        g_yyyy_0_xxxxz_0[i] = 3.0 * g_yy_0_xxxxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxz_1[i] * fz_be_0 + g_yyy_0_xxxxz_1[i] * wa_y[i];

        g_yyyy_0_xxxyy_0[i] = 3.0 * g_yy_0_xxxyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxy_1[i] * fi_acd_0 + g_yyy_0_xxxyy_1[i] * wa_y[i];

        g_yyyy_0_xxxyz_0[i] = 3.0 * g_yy_0_xxxyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyz_1[i] * fz_be_0 + g_yyy_0_xxxz_1[i] * fi_acd_0 + g_yyy_0_xxxyz_1[i] * wa_y[i];

        g_yyyy_0_xxxzz_0[i] = 3.0 * g_yy_0_xxxzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxzz_1[i] * fz_be_0 + g_yyy_0_xxxzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyy_0[i] = 3.0 * g_yy_0_xxyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxyy_1[i] * fi_acd_0 + g_yyy_0_xxyyy_1[i] * wa_y[i];

        g_yyyy_0_xxyyz_0[i] = 3.0 * g_yy_0_xxyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxyz_1[i] * fi_acd_0 + g_yyy_0_xxyyz_1[i] * wa_y[i];

        g_yyyy_0_xxyzz_0[i] = 3.0 * g_yy_0_xxyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyzz_1[i] * fz_be_0 + g_yyy_0_xxzz_1[i] * fi_acd_0 + g_yyy_0_xxyzz_1[i] * wa_y[i];

        g_yyyy_0_xxzzz_0[i] = 3.0 * g_yy_0_xxzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxzzz_1[i] * fz_be_0 + g_yyy_0_xxzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyy_0[i] = 3.0 * g_yy_0_xyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_0_xyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyy_1[i] * wa_y[i];

        g_yyyy_0_xyyyz_0[i] = 3.0 * g_yy_0_xyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyz_1[i] * wa_y[i];

        g_yyyy_0_xyyzz_0[i] = 3.0 * g_yy_0_xyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xyzz_1[i] * fi_acd_0 + g_yyy_0_xyyzz_1[i] * wa_y[i];

        g_yyyy_0_xyzzz_0[i] = 3.0 * g_yy_0_xyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyzzz_1[i] * fz_be_0 + g_yyy_0_xzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzz_1[i] * wa_y[i];

        g_yyyy_0_xzzzz_0[i] = 3.0 * g_yy_0_xzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyy_0[i] = 3.0 * g_yy_0_yyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_0_yyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyy_1[i] * wa_y[i];

        g_yyyy_0_yyyyz_0[i] = 3.0 * g_yy_0_yyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_0_yyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyz_1[i] * wa_y[i];

        g_yyyy_0_yyyzz_0[i] = 3.0 * g_yy_0_yyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_yyzz_1[i] * fi_acd_0 + g_yyy_0_yyyzz_1[i] * wa_y[i];

        g_yyyy_0_yyzzz_0[i] = 3.0 * g_yy_0_yyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzz_1[i] * wa_y[i];

        g_yyyy_0_yzzzz_0[i] = 3.0 * g_yy_0_yzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzzzz_1[i] * fz_be_0 + g_yyy_0_zzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzz_1[i] * wa_y[i];

        g_yyyy_0_zzzzz_0[i] = 3.0 * g_yy_0_zzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 231-252 components of targeted buffer : GSH

    auto g_yyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 231);

    auto g_yyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 232);

    auto g_yyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 233);

    auto g_yyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 234);

    auto g_yyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 235);

    auto g_yyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 236);

    auto g_yyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 237);

    auto g_yyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 238);

    auto g_yyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 239);

    auto g_yyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 240);

    auto g_yyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 241);

    auto g_yyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 242);

    auto g_yyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 243);

    auto g_yyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 244);

    auto g_yyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 245);

    auto g_yyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 246);

    auto g_yyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 247);

    auto g_yyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 248);

    auto g_yyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 249);

    auto g_yyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 250);

    auto g_yyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 251);

    #pragma omp simd aligned(g_yyy_0_xxxx_1, g_yyy_0_xxxxx_1, g_yyy_0_xxxxy_1, g_yyy_0_xxxxz_1, g_yyy_0_xxxy_1, g_yyy_0_xxxyy_1, g_yyy_0_xxxyz_1, g_yyy_0_xxxz_1, g_yyy_0_xxxzz_1, g_yyy_0_xxyy_1, g_yyy_0_xxyyy_1, g_yyy_0_xxyyz_1, g_yyy_0_xxyz_1, g_yyy_0_xxyzz_1, g_yyy_0_xxzz_1, g_yyy_0_xxzzz_1, g_yyy_0_xyyy_1, g_yyy_0_xyyyy_1, g_yyy_0_xyyyz_1, g_yyy_0_xyyz_1, g_yyy_0_xyyzz_1, g_yyy_0_xyzz_1, g_yyy_0_xyzzz_1, g_yyy_0_xzzz_1, g_yyy_0_xzzzz_1, g_yyy_0_yyyy_1, g_yyy_0_yyyyy_1, g_yyy_0_yyyyz_1, g_yyy_0_yyyz_1, g_yyy_0_yyyzz_1, g_yyy_0_yyzz_1, g_yyy_0_yyzzz_1, g_yyy_0_yzzz_1, g_yyy_0_yzzzz_1, g_yyy_0_zzzz_1, g_yyy_0_zzzzz_1, g_yyyz_0_xxxxx_0, g_yyyz_0_xxxxy_0, g_yyyz_0_xxxxz_0, g_yyyz_0_xxxyy_0, g_yyyz_0_xxxyz_0, g_yyyz_0_xxxzz_0, g_yyyz_0_xxyyy_0, g_yyyz_0_xxyyz_0, g_yyyz_0_xxyzz_0, g_yyyz_0_xxzzz_0, g_yyyz_0_xyyyy_0, g_yyyz_0_xyyyz_0, g_yyyz_0_xyyzz_0, g_yyyz_0_xyzzz_0, g_yyyz_0_xzzzz_0, g_yyyz_0_yyyyy_0, g_yyyz_0_yyyyz_0, g_yyyz_0_yyyzz_0, g_yyyz_0_yyzzz_0, g_yyyz_0_yzzzz_0, g_yyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxxxx_0[i] = g_yyy_0_xxxxx_1[i] * wa_z[i];

        g_yyyz_0_xxxxy_0[i] = g_yyy_0_xxxxy_1[i] * wa_z[i];

        g_yyyz_0_xxxxz_0[i] = g_yyy_0_xxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxz_1[i] * wa_z[i];

        g_yyyz_0_xxxyy_0[i] = g_yyy_0_xxxyy_1[i] * wa_z[i];

        g_yyyz_0_xxxyz_0[i] = g_yyy_0_xxxy_1[i] * fi_acd_0 + g_yyy_0_xxxyz_1[i] * wa_z[i];

        g_yyyz_0_xxxzz_0[i] = 2.0 * g_yyy_0_xxxz_1[i] * fi_acd_0 + g_yyy_0_xxxzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyy_0[i] = g_yyy_0_xxyyy_1[i] * wa_z[i];

        g_yyyz_0_xxyyz_0[i] = g_yyy_0_xxyy_1[i] * fi_acd_0 + g_yyy_0_xxyyz_1[i] * wa_z[i];

        g_yyyz_0_xxyzz_0[i] = 2.0 * g_yyy_0_xxyz_1[i] * fi_acd_0 + g_yyy_0_xxyzz_1[i] * wa_z[i];

        g_yyyz_0_xxzzz_0[i] = 3.0 * g_yyy_0_xxzz_1[i] * fi_acd_0 + g_yyy_0_xxzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyy_0[i] = g_yyy_0_xyyyy_1[i] * wa_z[i];

        g_yyyz_0_xyyyz_0[i] = g_yyy_0_xyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyz_1[i] * wa_z[i];

        g_yyyz_0_xyyzz_0[i] = 2.0 * g_yyy_0_xyyz_1[i] * fi_acd_0 + g_yyy_0_xyyzz_1[i] * wa_z[i];

        g_yyyz_0_xyzzz_0[i] = 3.0 * g_yyy_0_xyzz_1[i] * fi_acd_0 + g_yyy_0_xyzzz_1[i] * wa_z[i];

        g_yyyz_0_xzzzz_0[i] = 4.0 * g_yyy_0_xzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyy_0[i] = g_yyy_0_yyyyy_1[i] * wa_z[i];

        g_yyyz_0_yyyyz_0[i] = g_yyy_0_yyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyz_1[i] * wa_z[i];

        g_yyyz_0_yyyzz_0[i] = 2.0 * g_yyy_0_yyyz_1[i] * fi_acd_0 + g_yyy_0_yyyzz_1[i] * wa_z[i];

        g_yyyz_0_yyzzz_0[i] = 3.0 * g_yyy_0_yyzz_1[i] * fi_acd_0 + g_yyy_0_yyzzz_1[i] * wa_z[i];

        g_yyyz_0_yzzzz_0[i] = 4.0 * g_yyy_0_yzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzz_1[i] * wa_z[i];

        g_yyyz_0_zzzzz_0[i] = 5.0 * g_yyy_0_zzzz_1[i] * fi_acd_0 + g_yyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 252-273 components of targeted buffer : GSH

    auto g_yyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 252);

    auto g_yyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 253);

    auto g_yyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 254);

    auto g_yyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 255);

    auto g_yyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 256);

    auto g_yyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 257);

    auto g_yyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 258);

    auto g_yyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 259);

    auto g_yyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 260);

    auto g_yyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 261);

    auto g_yyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 262);

    auto g_yyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 263);

    auto g_yyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 264);

    auto g_yyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 265);

    auto g_yyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 266);

    auto g_yyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 267);

    auto g_yyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 268);

    auto g_yyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 269);

    auto g_yyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 270);

    auto g_yyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 271);

    auto g_yyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 272);

    #pragma omp simd aligned(g_yy_0_xxxxy_0, g_yy_0_xxxxy_1, g_yy_0_xxxyy_0, g_yy_0_xxxyy_1, g_yy_0_xxyyy_0, g_yy_0_xxyyy_1, g_yy_0_xyyyy_0, g_yy_0_xyyyy_1, g_yy_0_yyyyy_0, g_yy_0_yyyyy_1, g_yyz_0_xxxxy_1, g_yyz_0_xxxyy_1, g_yyz_0_xxyyy_1, g_yyz_0_xyyyy_1, g_yyz_0_yyyyy_1, g_yyzz_0_xxxxx_0, g_yyzz_0_xxxxy_0, g_yyzz_0_xxxxz_0, g_yyzz_0_xxxyy_0, g_yyzz_0_xxxyz_0, g_yyzz_0_xxxzz_0, g_yyzz_0_xxyyy_0, g_yyzz_0_xxyyz_0, g_yyzz_0_xxyzz_0, g_yyzz_0_xxzzz_0, g_yyzz_0_xyyyy_0, g_yyzz_0_xyyyz_0, g_yyzz_0_xyyzz_0, g_yyzz_0_xyzzz_0, g_yyzz_0_xzzzz_0, g_yyzz_0_yyyyy_0, g_yyzz_0_yyyyz_0, g_yyzz_0_yyyzz_0, g_yyzz_0_yyzzz_0, g_yyzz_0_yzzzz_0, g_yyzz_0_zzzzz_0, g_yzz_0_xxxxx_1, g_yzz_0_xxxxz_1, g_yzz_0_xxxyz_1, g_yzz_0_xxxz_1, g_yzz_0_xxxzz_1, g_yzz_0_xxyyz_1, g_yzz_0_xxyz_1, g_yzz_0_xxyzz_1, g_yzz_0_xxzz_1, g_yzz_0_xxzzz_1, g_yzz_0_xyyyz_1, g_yzz_0_xyyz_1, g_yzz_0_xyyzz_1, g_yzz_0_xyzz_1, g_yzz_0_xyzzz_1, g_yzz_0_xzzz_1, g_yzz_0_xzzzz_1, g_yzz_0_yyyyz_1, g_yzz_0_yyyz_1, g_yzz_0_yyyzz_1, g_yzz_0_yyzz_1, g_yzz_0_yyzzz_1, g_yzz_0_yzzz_1, g_yzz_0_yzzzz_1, g_yzz_0_zzzz_1, g_yzz_0_zzzzz_1, g_zz_0_xxxxx_0, g_zz_0_xxxxx_1, g_zz_0_xxxxz_0, g_zz_0_xxxxz_1, g_zz_0_xxxyz_0, g_zz_0_xxxyz_1, g_zz_0_xxxzz_0, g_zz_0_xxxzz_1, g_zz_0_xxyyz_0, g_zz_0_xxyyz_1, g_zz_0_xxyzz_0, g_zz_0_xxyzz_1, g_zz_0_xxzzz_0, g_zz_0_xxzzz_1, g_zz_0_xyyyz_0, g_zz_0_xyyyz_1, g_zz_0_xyyzz_0, g_zz_0_xyyzz_1, g_zz_0_xyzzz_0, g_zz_0_xyzzz_1, g_zz_0_xzzzz_0, g_zz_0_xzzzz_1, g_zz_0_yyyyz_0, g_zz_0_yyyyz_1, g_zz_0_yyyzz_0, g_zz_0_yyyzz_1, g_zz_0_yyzzz_0, g_zz_0_yyzzz_1, g_zz_0_yzzzz_0, g_zz_0_yzzzz_1, g_zz_0_zzzzz_0, g_zz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxxxx_0[i] = g_zz_0_xxxxx_0[i] * fbe_0 - g_zz_0_xxxxx_1[i] * fz_be_0 + g_yzz_0_xxxxx_1[i] * wa_y[i];

        g_yyzz_0_xxxxy_0[i] = g_yy_0_xxxxy_0[i] * fbe_0 - g_yy_0_xxxxy_1[i] * fz_be_0 + g_yyz_0_xxxxy_1[i] * wa_z[i];

        g_yyzz_0_xxxxz_0[i] = g_zz_0_xxxxz_0[i] * fbe_0 - g_zz_0_xxxxz_1[i] * fz_be_0 + g_yzz_0_xxxxz_1[i] * wa_y[i];

        g_yyzz_0_xxxyy_0[i] = g_yy_0_xxxyy_0[i] * fbe_0 - g_yy_0_xxxyy_1[i] * fz_be_0 + g_yyz_0_xxxyy_1[i] * wa_z[i];

        g_yyzz_0_xxxyz_0[i] = g_zz_0_xxxyz_0[i] * fbe_0 - g_zz_0_xxxyz_1[i] * fz_be_0 + g_yzz_0_xxxz_1[i] * fi_acd_0 + g_yzz_0_xxxyz_1[i] * wa_y[i];

        g_yyzz_0_xxxzz_0[i] = g_zz_0_xxxzz_0[i] * fbe_0 - g_zz_0_xxxzz_1[i] * fz_be_0 + g_yzz_0_xxxzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyy_0[i] = g_yy_0_xxyyy_0[i] * fbe_0 - g_yy_0_xxyyy_1[i] * fz_be_0 + g_yyz_0_xxyyy_1[i] * wa_z[i];

        g_yyzz_0_xxyyz_0[i] = g_zz_0_xxyyz_0[i] * fbe_0 - g_zz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxyz_1[i] * fi_acd_0 + g_yzz_0_xxyyz_1[i] * wa_y[i];

        g_yyzz_0_xxyzz_0[i] = g_zz_0_xxyzz_0[i] * fbe_0 - g_zz_0_xxyzz_1[i] * fz_be_0 + g_yzz_0_xxzz_1[i] * fi_acd_0 + g_yzz_0_xxyzz_1[i] * wa_y[i];

        g_yyzz_0_xxzzz_0[i] = g_zz_0_xxzzz_0[i] * fbe_0 - g_zz_0_xxzzz_1[i] * fz_be_0 + g_yzz_0_xxzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyy_0[i] = g_yy_0_xyyyy_0[i] * fbe_0 - g_yy_0_xyyyy_1[i] * fz_be_0 + g_yyz_0_xyyyy_1[i] * wa_z[i];

        g_yyzz_0_xyyyz_0[i] = g_zz_0_xyyyz_0[i] * fbe_0 - g_zz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyz_1[i] * wa_y[i];

        g_yyzz_0_xyyzz_0[i] = g_zz_0_xyyzz_0[i] * fbe_0 - g_zz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xyzz_1[i] * fi_acd_0 + g_yzz_0_xyyzz_1[i] * wa_y[i];

        g_yyzz_0_xyzzz_0[i] = g_zz_0_xyzzz_0[i] * fbe_0 - g_zz_0_xyzzz_1[i] * fz_be_0 + g_yzz_0_xzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzz_1[i] * wa_y[i];

        g_yyzz_0_xzzzz_0[i] = g_zz_0_xzzzz_0[i] * fbe_0 - g_zz_0_xzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyy_0[i] = g_yy_0_yyyyy_0[i] * fbe_0 - g_yy_0_yyyyy_1[i] * fz_be_0 + g_yyz_0_yyyyy_1[i] * wa_z[i];

        g_yyzz_0_yyyyz_0[i] = g_zz_0_yyyyz_0[i] * fbe_0 - g_zz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_0_yyyz_1[i] * fi_acd_0 + g_yzz_0_yyyyz_1[i] * wa_y[i];

        g_yyzz_0_yyyzz_0[i] = g_zz_0_yyyzz_0[i] * fbe_0 - g_zz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_yyzz_1[i] * fi_acd_0 + g_yzz_0_yyyzz_1[i] * wa_y[i];

        g_yyzz_0_yyzzz_0[i] = g_zz_0_yyzzz_0[i] * fbe_0 - g_zz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yzzz_1[i] * fi_acd_0 + g_yzz_0_yyzzz_1[i] * wa_y[i];

        g_yyzz_0_yzzzz_0[i] = g_zz_0_yzzzz_0[i] * fbe_0 - g_zz_0_yzzzz_1[i] * fz_be_0 + g_yzz_0_zzzz_1[i] * fi_acd_0 + g_yzz_0_yzzzz_1[i] * wa_y[i];

        g_yyzz_0_zzzzz_0[i] = g_zz_0_zzzzz_0[i] * fbe_0 - g_zz_0_zzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 273-294 components of targeted buffer : GSH

    auto g_yzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 273);

    auto g_yzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 274);

    auto g_yzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 275);

    auto g_yzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 276);

    auto g_yzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 277);

    auto g_yzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 278);

    auto g_yzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 279);

    auto g_yzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 280);

    auto g_yzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 281);

    auto g_yzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 282);

    auto g_yzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 283);

    auto g_yzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 284);

    auto g_yzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 285);

    auto g_yzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 286);

    auto g_yzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 287);

    auto g_yzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 288);

    auto g_yzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 289);

    auto g_yzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 290);

    auto g_yzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 291);

    auto g_yzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 292);

    auto g_yzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 293);

    #pragma omp simd aligned(g_yzzz_0_xxxxx_0, g_yzzz_0_xxxxy_0, g_yzzz_0_xxxxz_0, g_yzzz_0_xxxyy_0, g_yzzz_0_xxxyz_0, g_yzzz_0_xxxzz_0, g_yzzz_0_xxyyy_0, g_yzzz_0_xxyyz_0, g_yzzz_0_xxyzz_0, g_yzzz_0_xxzzz_0, g_yzzz_0_xyyyy_0, g_yzzz_0_xyyyz_0, g_yzzz_0_xyyzz_0, g_yzzz_0_xyzzz_0, g_yzzz_0_xzzzz_0, g_yzzz_0_yyyyy_0, g_yzzz_0_yyyyz_0, g_yzzz_0_yyyzz_0, g_yzzz_0_yyzzz_0, g_yzzz_0_yzzzz_0, g_yzzz_0_zzzzz_0, g_zzz_0_xxxx_1, g_zzz_0_xxxxx_1, g_zzz_0_xxxxy_1, g_zzz_0_xxxxz_1, g_zzz_0_xxxy_1, g_zzz_0_xxxyy_1, g_zzz_0_xxxyz_1, g_zzz_0_xxxz_1, g_zzz_0_xxxzz_1, g_zzz_0_xxyy_1, g_zzz_0_xxyyy_1, g_zzz_0_xxyyz_1, g_zzz_0_xxyz_1, g_zzz_0_xxyzz_1, g_zzz_0_xxzz_1, g_zzz_0_xxzzz_1, g_zzz_0_xyyy_1, g_zzz_0_xyyyy_1, g_zzz_0_xyyyz_1, g_zzz_0_xyyz_1, g_zzz_0_xyyzz_1, g_zzz_0_xyzz_1, g_zzz_0_xyzzz_1, g_zzz_0_xzzz_1, g_zzz_0_xzzzz_1, g_zzz_0_yyyy_1, g_zzz_0_yyyyy_1, g_zzz_0_yyyyz_1, g_zzz_0_yyyz_1, g_zzz_0_yyyzz_1, g_zzz_0_yyzz_1, g_zzz_0_yyzzz_1, g_zzz_0_yzzz_1, g_zzz_0_yzzzz_1, g_zzz_0_zzzz_1, g_zzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxxxx_0[i] = g_zzz_0_xxxxx_1[i] * wa_y[i];

        g_yzzz_0_xxxxy_0[i] = g_zzz_0_xxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxy_1[i] * wa_y[i];

        g_yzzz_0_xxxxz_0[i] = g_zzz_0_xxxxz_1[i] * wa_y[i];

        g_yzzz_0_xxxyy_0[i] = 2.0 * g_zzz_0_xxxy_1[i] * fi_acd_0 + g_zzz_0_xxxyy_1[i] * wa_y[i];

        g_yzzz_0_xxxyz_0[i] = g_zzz_0_xxxz_1[i] * fi_acd_0 + g_zzz_0_xxxyz_1[i] * wa_y[i];

        g_yzzz_0_xxxzz_0[i] = g_zzz_0_xxxzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyy_0[i] = 3.0 * g_zzz_0_xxyy_1[i] * fi_acd_0 + g_zzz_0_xxyyy_1[i] * wa_y[i];

        g_yzzz_0_xxyyz_0[i] = 2.0 * g_zzz_0_xxyz_1[i] * fi_acd_0 + g_zzz_0_xxyyz_1[i] * wa_y[i];

        g_yzzz_0_xxyzz_0[i] = g_zzz_0_xxzz_1[i] * fi_acd_0 + g_zzz_0_xxyzz_1[i] * wa_y[i];

        g_yzzz_0_xxzzz_0[i] = g_zzz_0_xxzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyy_0[i] = 4.0 * g_zzz_0_xyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyy_1[i] * wa_y[i];

        g_yzzz_0_xyyyz_0[i] = 3.0 * g_zzz_0_xyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyz_1[i] * wa_y[i];

        g_yzzz_0_xyyzz_0[i] = 2.0 * g_zzz_0_xyzz_1[i] * fi_acd_0 + g_zzz_0_xyyzz_1[i] * wa_y[i];

        g_yzzz_0_xyzzz_0[i] = g_zzz_0_xzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzz_1[i] * wa_y[i];

        g_yzzz_0_xzzzz_0[i] = g_zzz_0_xzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyy_0[i] = 5.0 * g_zzz_0_yyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyy_1[i] * wa_y[i];

        g_yzzz_0_yyyyz_0[i] = 4.0 * g_zzz_0_yyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyz_1[i] * wa_y[i];

        g_yzzz_0_yyyzz_0[i] = 3.0 * g_zzz_0_yyzz_1[i] * fi_acd_0 + g_zzz_0_yyyzz_1[i] * wa_y[i];

        g_yzzz_0_yyzzz_0[i] = 2.0 * g_zzz_0_yzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzz_1[i] * wa_y[i];

        g_yzzz_0_yzzzz_0[i] = g_zzz_0_zzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzz_1[i] * wa_y[i];

        g_yzzz_0_zzzzz_0[i] = g_zzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 294-315 components of targeted buffer : GSH

    auto g_zzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_gsh + 294);

    auto g_zzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_gsh + 295);

    auto g_zzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_gsh + 296);

    auto g_zzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_gsh + 297);

    auto g_zzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_gsh + 298);

    auto g_zzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_gsh + 299);

    auto g_zzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_gsh + 300);

    auto g_zzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_gsh + 301);

    auto g_zzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_gsh + 302);

    auto g_zzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_gsh + 303);

    auto g_zzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_gsh + 304);

    auto g_zzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_gsh + 305);

    auto g_zzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_gsh + 306);

    auto g_zzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_gsh + 307);

    auto g_zzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_gsh + 308);

    auto g_zzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_gsh + 309);

    auto g_zzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_gsh + 310);

    auto g_zzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_gsh + 311);

    auto g_zzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_gsh + 312);

    auto g_zzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_gsh + 313);

    auto g_zzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_gsh + 314);

    #pragma omp simd aligned(g_zz_0_xxxxx_0, g_zz_0_xxxxx_1, g_zz_0_xxxxy_0, g_zz_0_xxxxy_1, g_zz_0_xxxxz_0, g_zz_0_xxxxz_1, g_zz_0_xxxyy_0, g_zz_0_xxxyy_1, g_zz_0_xxxyz_0, g_zz_0_xxxyz_1, g_zz_0_xxxzz_0, g_zz_0_xxxzz_1, g_zz_0_xxyyy_0, g_zz_0_xxyyy_1, g_zz_0_xxyyz_0, g_zz_0_xxyyz_1, g_zz_0_xxyzz_0, g_zz_0_xxyzz_1, g_zz_0_xxzzz_0, g_zz_0_xxzzz_1, g_zz_0_xyyyy_0, g_zz_0_xyyyy_1, g_zz_0_xyyyz_0, g_zz_0_xyyyz_1, g_zz_0_xyyzz_0, g_zz_0_xyyzz_1, g_zz_0_xyzzz_0, g_zz_0_xyzzz_1, g_zz_0_xzzzz_0, g_zz_0_xzzzz_1, g_zz_0_yyyyy_0, g_zz_0_yyyyy_1, g_zz_0_yyyyz_0, g_zz_0_yyyyz_1, g_zz_0_yyyzz_0, g_zz_0_yyyzz_1, g_zz_0_yyzzz_0, g_zz_0_yyzzz_1, g_zz_0_yzzzz_0, g_zz_0_yzzzz_1, g_zz_0_zzzzz_0, g_zz_0_zzzzz_1, g_zzz_0_xxxx_1, g_zzz_0_xxxxx_1, g_zzz_0_xxxxy_1, g_zzz_0_xxxxz_1, g_zzz_0_xxxy_1, g_zzz_0_xxxyy_1, g_zzz_0_xxxyz_1, g_zzz_0_xxxz_1, g_zzz_0_xxxzz_1, g_zzz_0_xxyy_1, g_zzz_0_xxyyy_1, g_zzz_0_xxyyz_1, g_zzz_0_xxyz_1, g_zzz_0_xxyzz_1, g_zzz_0_xxzz_1, g_zzz_0_xxzzz_1, g_zzz_0_xyyy_1, g_zzz_0_xyyyy_1, g_zzz_0_xyyyz_1, g_zzz_0_xyyz_1, g_zzz_0_xyyzz_1, g_zzz_0_xyzz_1, g_zzz_0_xyzzz_1, g_zzz_0_xzzz_1, g_zzz_0_xzzzz_1, g_zzz_0_yyyy_1, g_zzz_0_yyyyy_1, g_zzz_0_yyyyz_1, g_zzz_0_yyyz_1, g_zzz_0_yyyzz_1, g_zzz_0_yyzz_1, g_zzz_0_yyzzz_1, g_zzz_0_yzzz_1, g_zzz_0_yzzzz_1, g_zzz_0_zzzz_1, g_zzz_0_zzzzz_1, g_zzzz_0_xxxxx_0, g_zzzz_0_xxxxy_0, g_zzzz_0_xxxxz_0, g_zzzz_0_xxxyy_0, g_zzzz_0_xxxyz_0, g_zzzz_0_xxxzz_0, g_zzzz_0_xxyyy_0, g_zzzz_0_xxyyz_0, g_zzzz_0_xxyzz_0, g_zzzz_0_xxzzz_0, g_zzzz_0_xyyyy_0, g_zzzz_0_xyyyz_0, g_zzzz_0_xyyzz_0, g_zzzz_0_xyzzz_0, g_zzzz_0_xzzzz_0, g_zzzz_0_yyyyy_0, g_zzzz_0_yyyyz_0, g_zzzz_0_yyyzz_0, g_zzzz_0_yyzzz_0, g_zzzz_0_yzzzz_0, g_zzzz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxxxx_0[i] = 3.0 * g_zz_0_xxxxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxx_1[i] * fz_be_0 + g_zzz_0_xxxxx_1[i] * wa_z[i];

        g_zzzz_0_xxxxy_0[i] = 3.0 * g_zz_0_xxxxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxy_1[i] * fz_be_0 + g_zzz_0_xxxxy_1[i] * wa_z[i];

        g_zzzz_0_xxxxz_0[i] = 3.0 * g_zz_0_xxxxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxz_1[i] * fz_be_0 + g_zzz_0_xxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxz_1[i] * wa_z[i];

        g_zzzz_0_xxxyy_0[i] = 3.0 * g_zz_0_xxxyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyy_1[i] * fz_be_0 + g_zzz_0_xxxyy_1[i] * wa_z[i];

        g_zzzz_0_xxxyz_0[i] = 3.0 * g_zz_0_xxxyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyz_1[i] * fz_be_0 + g_zzz_0_xxxy_1[i] * fi_acd_0 + g_zzz_0_xxxyz_1[i] * wa_z[i];

        g_zzzz_0_xxxzz_0[i] = 3.0 * g_zz_0_xxxzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxz_1[i] * fi_acd_0 + g_zzz_0_xxxzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyy_0[i] = 3.0 * g_zz_0_xxyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyy_1[i] * fz_be_0 + g_zzz_0_xxyyy_1[i] * wa_z[i];

        g_zzzz_0_xxyyz_0[i] = 3.0 * g_zz_0_xxyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyz_1[i] * fz_be_0 + g_zzz_0_xxyy_1[i] * fi_acd_0 + g_zzz_0_xxyyz_1[i] * wa_z[i];

        g_zzzz_0_xxyzz_0[i] = 3.0 * g_zz_0_xxyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxyz_1[i] * fi_acd_0 + g_zzz_0_xxyzz_1[i] * wa_z[i];

        g_zzzz_0_xxzzz_0[i] = 3.0 * g_zz_0_xxzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxzz_1[i] * fi_acd_0 + g_zzz_0_xxzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyy_0[i] = 3.0 * g_zz_0_xyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyy_1[i] * fz_be_0 + g_zzz_0_xyyyy_1[i] * wa_z[i];

        g_zzzz_0_xyyyz_0[i] = 3.0 * g_zz_0_xyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyz_1[i] * fz_be_0 + g_zzz_0_xyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyz_1[i] * wa_z[i];

        g_zzzz_0_xyyzz_0[i] = 3.0 * g_zz_0_xyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xyyz_1[i] * fi_acd_0 + g_zzz_0_xyyzz_1[i] * wa_z[i];

        g_zzzz_0_xyzzz_0[i] = 3.0 * g_zz_0_xyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xyzz_1[i] * fi_acd_0 + g_zzz_0_xyzzz_1[i] * wa_z[i];

        g_zzzz_0_xzzzz_0[i] = 3.0 * g_zz_0_xzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyy_0[i] = 3.0 * g_zz_0_yyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyy_1[i] * fz_be_0 + g_zzz_0_yyyyy_1[i] * wa_z[i];

        g_zzzz_0_yyyyz_0[i] = 3.0 * g_zz_0_yyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyz_1[i] * fz_be_0 + g_zzz_0_yyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyz_1[i] * wa_z[i];

        g_zzzz_0_yyyzz_0[i] = 3.0 * g_zz_0_yyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yyyz_1[i] * fi_acd_0 + g_zzz_0_yyyzz_1[i] * wa_z[i];

        g_zzzz_0_yyzzz_0[i] = 3.0 * g_zz_0_yyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_yyzz_1[i] * fi_acd_0 + g_zzz_0_yyzzz_1[i] * wa_z[i];

        g_zzzz_0_yzzzz_0[i] = 3.0 * g_zz_0_yzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_yzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzz_1[i] * wa_z[i];

        g_zzzz_0_zzzzz_0[i] = 3.0 * g_zz_0_zzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_zzzz_1[i] * fi_acd_0 + g_zzz_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

