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

#include "ThreeCenterElectronRepulsionPrimRecISF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isf,
                                 size_t idx_eri_0_gsf,
                                 size_t idx_eri_1_gsf,
                                 size_t idx_eri_1_hsd,
                                 size_t idx_eri_1_hsf,
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

    /// Set up components of auxilary buffer : GSF

    auto g_xxxx_0_xxx_0 = pbuffer.data(idx_eri_0_gsf);

    auto g_xxxx_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 1);

    auto g_xxxx_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 2);

    auto g_xxxx_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 3);

    auto g_xxxx_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 4);

    auto g_xxxx_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 5);

    auto g_xxxx_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 6);

    auto g_xxxx_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 7);

    auto g_xxxx_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 8);

    auto g_xxxx_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 9);

    auto g_xxxy_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 10);

    auto g_xxxy_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 12);

    auto g_xxxy_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 15);

    auto g_xxxz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 20);

    auto g_xxxz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 21);

    auto g_xxxz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 23);

    auto g_xxyy_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 30);

    auto g_xxyy_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 31);

    auto g_xxyy_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 32);

    auto g_xxyy_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 33);

    auto g_xxyy_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 34);

    auto g_xxyy_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 35);

    auto g_xxyy_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 36);

    auto g_xxyy_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 37);

    auto g_xxyy_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 38);

    auto g_xxyy_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 39);

    auto g_xxzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 50);

    auto g_xxzz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 51);

    auto g_xxzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 52);

    auto g_xxzz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 53);

    auto g_xxzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 54);

    auto g_xxzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 55);

    auto g_xxzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 56);

    auto g_xxzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 57);

    auto g_xxzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 58);

    auto g_xxzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 59);

    auto g_xyyy_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 61);

    auto g_xyyy_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 63);

    auto g_xyyy_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 64);

    auto g_xyyy_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 66);

    auto g_xyyy_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 67);

    auto g_xyyy_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 68);

    auto g_xyyy_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 69);

    auto g_xzzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 92);

    auto g_xzzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 94);

    auto g_xzzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 95);

    auto g_xzzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 96);

    auto g_xzzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 97);

    auto g_xzzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 98);

    auto g_xzzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 99);

    auto g_yyyy_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 100);

    auto g_yyyy_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 101);

    auto g_yyyy_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 102);

    auto g_yyyy_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 103);

    auto g_yyyy_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 104);

    auto g_yyyy_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 105);

    auto g_yyyy_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 106);

    auto g_yyyy_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 107);

    auto g_yyyy_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 108);

    auto g_yyyy_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 109);

    auto g_yyyz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 111);

    auto g_yyyz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 113);

    auto g_yyyz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 116);

    auto g_yyzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 120);

    auto g_yyzz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 121);

    auto g_yyzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 122);

    auto g_yyzz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 123);

    auto g_yyzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 124);

    auto g_yyzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 125);

    auto g_yyzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 126);

    auto g_yyzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 127);

    auto g_yyzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 128);

    auto g_yyzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 129);

    auto g_yzzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 130);

    auto g_yzzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 132);

    auto g_yzzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 134);

    auto g_yzzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 135);

    auto g_yzzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 137);

    auto g_yzzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 138);

    auto g_yzzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 139);

    auto g_zzzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 140);

    auto g_zzzz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 141);

    auto g_zzzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 142);

    auto g_zzzz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 143);

    auto g_zzzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 144);

    auto g_zzzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 145);

    auto g_zzzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 146);

    auto g_zzzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 147);

    auto g_zzzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 148);

    auto g_zzzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 149);

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

    auto g_xxxy_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 10);

    auto g_xxxy_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 12);

    auto g_xxxy_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 15);

    auto g_xxxz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 20);

    auto g_xxxz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 21);

    auto g_xxxz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 23);

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

    auto g_xyyy_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 69);

    auto g_xzzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 92);

    auto g_xzzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 94);

    auto g_xzzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 95);

    auto g_xzzz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 96);

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

    auto g_yyyz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 111);

    auto g_yyyz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 113);

    auto g_yyyz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 116);

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

    auto g_yzzz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 130);

    auto g_yzzz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 132);

    auto g_yzzz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 134);

    auto g_yzzz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 135);

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

    /// Set up components of auxilary buffer : HSD

    auto g_xxxxx_0_xx_1 = pbuffer.data(idx_eri_1_hsd);

    auto g_xxxxx_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 1);

    auto g_xxxxx_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 2);

    auto g_xxxxx_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 3);

    auto g_xxxxx_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 4);

    auto g_xxxxx_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 5);

    auto g_xxxxz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 14);

    auto g_xxxxz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 16);

    auto g_xxxxz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 17);

    auto g_xxxyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 18);

    auto g_xxxyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 19);

    auto g_xxxyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 20);

    auto g_xxxyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 21);

    auto g_xxxyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 22);

    auto g_xxxyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 23);

    auto g_xxxzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 30);

    auto g_xxxzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 31);

    auto g_xxxzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 32);

    auto g_xxxzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 33);

    auto g_xxxzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 34);

    auto g_xxxzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 35);

    auto g_xxyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 36);

    auto g_xxyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 37);

    auto g_xxyyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 38);

    auto g_xxyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 39);

    auto g_xxyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 40);

    auto g_xxyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 41);

    auto g_xxzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 54);

    auto g_xxzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 55);

    auto g_xxzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 56);

    auto g_xxzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 57);

    auto g_xxzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 58);

    auto g_xxzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 59);

    auto g_xyyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 61);

    auto g_xyyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 63);

    auto g_xyyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 64);

    auto g_xyyzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 76);

    auto g_xzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 86);

    auto g_xzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 88);

    auto g_xzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 89);

    auto g_yyyyy_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 90);

    auto g_yyyyy_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 91);

    auto g_yyyyy_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 92);

    auto g_yyyyy_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 93);

    auto g_yyyyy_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 94);

    auto g_yyyyy_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 95);

    auto g_yyyyz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 98);

    auto g_yyyyz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 100);

    auto g_yyyyz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 101);

    auto g_yyyzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 102);

    auto g_yyyzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 103);

    auto g_yyyzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 104);

    auto g_yyyzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 105);

    auto g_yyyzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 106);

    auto g_yyyzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 107);

    auto g_yyzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 108);

    auto g_yyzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 109);

    auto g_yyzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 110);

    auto g_yyzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 111);

    auto g_yyzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 112);

    auto g_yyzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 113);

    auto g_yzzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 115);

    auto g_yzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 116);

    auto g_yzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 117);

    auto g_yzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 118);

    auto g_yzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 119);

    auto g_zzzzz_0_xx_1 = pbuffer.data(idx_eri_1_hsd + 120);

    auto g_zzzzz_0_xy_1 = pbuffer.data(idx_eri_1_hsd + 121);

    auto g_zzzzz_0_xz_1 = pbuffer.data(idx_eri_1_hsd + 122);

    auto g_zzzzz_0_yy_1 = pbuffer.data(idx_eri_1_hsd + 123);

    auto g_zzzzz_0_yz_1 = pbuffer.data(idx_eri_1_hsd + 124);

    auto g_zzzzz_0_zz_1 = pbuffer.data(idx_eri_1_hsd + 125);

    /// Set up components of auxilary buffer : HSF

    auto g_xxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_hsf);

    auto g_xxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 1);

    auto g_xxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 2);

    auto g_xxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 3);

    auto g_xxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 4);

    auto g_xxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 5);

    auto g_xxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 6);

    auto g_xxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 7);

    auto g_xxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 8);

    auto g_xxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 9);

    auto g_xxxxy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 10);

    auto g_xxxxy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 11);

    auto g_xxxxy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 12);

    auto g_xxxxy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 13);

    auto g_xxxxy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 15);

    auto g_xxxxy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 16);

    auto g_xxxxz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 20);

    auto g_xxxxz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 21);

    auto g_xxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 22);

    auto g_xxxxz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 23);

    auto g_xxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 24);

    auto g_xxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 25);

    auto g_xxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 27);

    auto g_xxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 28);

    auto g_xxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 29);

    auto g_xxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 30);

    auto g_xxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 31);

    auto g_xxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 32);

    auto g_xxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 33);

    auto g_xxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 34);

    auto g_xxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 35);

    auto g_xxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 36);

    auto g_xxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 37);

    auto g_xxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 38);

    auto g_xxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 39);

    auto g_xxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 50);

    auto g_xxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 51);

    auto g_xxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 52);

    auto g_xxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 53);

    auto g_xxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 54);

    auto g_xxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 55);

    auto g_xxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 56);

    auto g_xxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 57);

    auto g_xxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 58);

    auto g_xxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 59);

    auto g_xxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 60);

    auto g_xxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 61);

    auto g_xxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 62);

    auto g_xxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 63);

    auto g_xxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 64);

    auto g_xxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 65);

    auto g_xxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 66);

    auto g_xxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 67);

    auto g_xxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 68);

    auto g_xxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 69);

    auto g_xxyyz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 71);

    auto g_xxyyz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 73);

    auto g_xxyzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 80);

    auto g_xxyzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 82);

    auto g_xxyzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 85);

    auto g_xxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 90);

    auto g_xxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 91);

    auto g_xxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 92);

    auto g_xxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 93);

    auto g_xxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 94);

    auto g_xxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 95);

    auto g_xxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 96);

    auto g_xxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 97);

    auto g_xxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 98);

    auto g_xxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 99);

    auto g_xyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 100);

    auto g_xyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 101);

    auto g_xyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 103);

    auto g_xyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 104);

    auto g_xyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 106);

    auto g_xyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 107);

    auto g_xyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 108);

    auto g_xyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 109);

    auto g_xyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 124);

    auto g_xyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 126);

    auto g_xyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 127);

    auto g_xyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 128);

    auto g_xyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 129);

    auto g_xzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 140);

    auto g_xzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 142);

    auto g_xzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 144);

    auto g_xzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 145);

    auto g_xzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 146);

    auto g_xzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 147);

    auto g_xzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 148);

    auto g_xzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 149);

    auto g_yyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 150);

    auto g_yyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 151);

    auto g_yyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 152);

    auto g_yyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 153);

    auto g_yyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 154);

    auto g_yyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 155);

    auto g_yyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 156);

    auto g_yyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 157);

    auto g_yyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 158);

    auto g_yyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 159);

    auto g_yyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 161);

    auto g_yyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 162);

    auto g_yyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 163);

    auto g_yyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 164);

    auto g_yyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 165);

    auto g_yyyyz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 166);

    auto g_yyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 167);

    auto g_yyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 168);

    auto g_yyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 169);

    auto g_yyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 170);

    auto g_yyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 171);

    auto g_yyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 172);

    auto g_yyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 173);

    auto g_yyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 174);

    auto g_yyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 175);

    auto g_yyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 176);

    auto g_yyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 177);

    auto g_yyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 178);

    auto g_yyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 179);

    auto g_yyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 180);

    auto g_yyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 181);

    auto g_yyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 182);

    auto g_yyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 183);

    auto g_yyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 184);

    auto g_yyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 185);

    auto g_yyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 186);

    auto g_yyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 187);

    auto g_yyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 188);

    auto g_yyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 189);

    auto g_yzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 190);

    auto g_yzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 191);

    auto g_yzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 192);

    auto g_yzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 193);

    auto g_yzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 194);

    auto g_yzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 195);

    auto g_yzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 196);

    auto g_yzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 197);

    auto g_yzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 198);

    auto g_yzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 199);

    auto g_zzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_hsf + 200);

    auto g_zzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_hsf + 201);

    auto g_zzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_hsf + 202);

    auto g_zzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_hsf + 203);

    auto g_zzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_hsf + 204);

    auto g_zzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_hsf + 205);

    auto g_zzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_hsf + 206);

    auto g_zzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_hsf + 207);

    auto g_zzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_hsf + 208);

    auto g_zzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_hsf + 209);

    /// Set up 0-10 components of targeted buffer : ISF

    auto g_xxxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_isf);

    auto g_xxxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 1);

    auto g_xxxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 2);

    auto g_xxxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 3);

    auto g_xxxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 4);

    auto g_xxxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 5);

    auto g_xxxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 6);

    auto g_xxxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 7);

    auto g_xxxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 8);

    auto g_xxxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 9);

    #pragma omp simd aligned(g_xxxx_0_xxx_0, g_xxxx_0_xxx_1, g_xxxx_0_xxy_0, g_xxxx_0_xxy_1, g_xxxx_0_xxz_0, g_xxxx_0_xxz_1, g_xxxx_0_xyy_0, g_xxxx_0_xyy_1, g_xxxx_0_xyz_0, g_xxxx_0_xyz_1, g_xxxx_0_xzz_0, g_xxxx_0_xzz_1, g_xxxx_0_yyy_0, g_xxxx_0_yyy_1, g_xxxx_0_yyz_0, g_xxxx_0_yyz_1, g_xxxx_0_yzz_0, g_xxxx_0_yzz_1, g_xxxx_0_zzz_0, g_xxxx_0_zzz_1, g_xxxxx_0_xx_1, g_xxxxx_0_xxx_1, g_xxxxx_0_xxy_1, g_xxxxx_0_xxz_1, g_xxxxx_0_xy_1, g_xxxxx_0_xyy_1, g_xxxxx_0_xyz_1, g_xxxxx_0_xz_1, g_xxxxx_0_xzz_1, g_xxxxx_0_yy_1, g_xxxxx_0_yyy_1, g_xxxxx_0_yyz_1, g_xxxxx_0_yz_1, g_xxxxx_0_yzz_1, g_xxxxx_0_zz_1, g_xxxxx_0_zzz_1, g_xxxxxx_0_xxx_0, g_xxxxxx_0_xxy_0, g_xxxxxx_0_xxz_0, g_xxxxxx_0_xyy_0, g_xxxxxx_0_xyz_0, g_xxxxxx_0_xzz_0, g_xxxxxx_0_yyy_0, g_xxxxxx_0_yyz_0, g_xxxxxx_0_yzz_0, g_xxxxxx_0_zzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xxx_0[i] = 5.0 * g_xxxx_0_xxx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxx_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xx_1[i] * fi_acd_0 + g_xxxxx_0_xxx_1[i] * wa_x[i];

        g_xxxxxx_0_xxy_0[i] = 5.0 * g_xxxx_0_xxy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xy_1[i] * fi_acd_0 + g_xxxxx_0_xxy_1[i] * wa_x[i];

        g_xxxxxx_0_xxz_0[i] = 5.0 * g_xxxx_0_xxz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xz_1[i] * fi_acd_0 + g_xxxxx_0_xxz_1[i] * wa_x[i];

        g_xxxxxx_0_xyy_0[i] = 5.0 * g_xxxx_0_xyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyy_1[i] * fz_be_0 + g_xxxxx_0_yy_1[i] * fi_acd_0 + g_xxxxx_0_xyy_1[i] * wa_x[i];

        g_xxxxxx_0_xyz_0[i] = 5.0 * g_xxxx_0_xyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyz_1[i] * fz_be_0 + g_xxxxx_0_yz_1[i] * fi_acd_0 + g_xxxxx_0_xyz_1[i] * wa_x[i];

        g_xxxxxx_0_xzz_0[i] = 5.0 * g_xxxx_0_xzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xzz_1[i] * fz_be_0 + g_xxxxx_0_zz_1[i] * fi_acd_0 + g_xxxxx_0_xzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyy_0[i] = 5.0 * g_xxxx_0_yyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyy_1[i] * fz_be_0 + g_xxxxx_0_yyy_1[i] * wa_x[i];

        g_xxxxxx_0_yyz_0[i] = 5.0 * g_xxxx_0_yyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyz_1[i] * fz_be_0 + g_xxxxx_0_yyz_1[i] * wa_x[i];

        g_xxxxxx_0_yzz_0[i] = 5.0 * g_xxxx_0_yzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yzz_1[i] * fz_be_0 + g_xxxxx_0_yzz_1[i] * wa_x[i];

        g_xxxxxx_0_zzz_0[i] = 5.0 * g_xxxx_0_zzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zzz_1[i] * fz_be_0 + g_xxxxx_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : ISF

    auto g_xxxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 10);

    auto g_xxxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 11);

    auto g_xxxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 12);

    auto g_xxxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 13);

    auto g_xxxxxy_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 14);

    auto g_xxxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 15);

    auto g_xxxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 16);

    auto g_xxxxxy_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 17);

    auto g_xxxxxy_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 18);

    auto g_xxxxxy_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 19);

    #pragma omp simd aligned(g_xxxxx_0_xx_1, g_xxxxx_0_xxx_1, g_xxxxx_0_xxy_1, g_xxxxx_0_xxz_1, g_xxxxx_0_xy_1, g_xxxxx_0_xyy_1, g_xxxxx_0_xyz_1, g_xxxxx_0_xz_1, g_xxxxx_0_xzz_1, g_xxxxx_0_yy_1, g_xxxxx_0_yyy_1, g_xxxxx_0_yyz_1, g_xxxxx_0_yz_1, g_xxxxx_0_yzz_1, g_xxxxx_0_zz_1, g_xxxxx_0_zzz_1, g_xxxxxy_0_xxx_0, g_xxxxxy_0_xxy_0, g_xxxxxy_0_xxz_0, g_xxxxxy_0_xyy_0, g_xxxxxy_0_xyz_0, g_xxxxxy_0_xzz_0, g_xxxxxy_0_yyy_0, g_xxxxxy_0_yyz_0, g_xxxxxy_0_yzz_0, g_xxxxxy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xxx_0[i] = g_xxxxx_0_xxx_1[i] * wa_y[i];

        g_xxxxxy_0_xxy_0[i] = g_xxxxx_0_xx_1[i] * fi_acd_0 + g_xxxxx_0_xxy_1[i] * wa_y[i];

        g_xxxxxy_0_xxz_0[i] = g_xxxxx_0_xxz_1[i] * wa_y[i];

        g_xxxxxy_0_xyy_0[i] = 2.0 * g_xxxxx_0_xy_1[i] * fi_acd_0 + g_xxxxx_0_xyy_1[i] * wa_y[i];

        g_xxxxxy_0_xyz_0[i] = g_xxxxx_0_xz_1[i] * fi_acd_0 + g_xxxxx_0_xyz_1[i] * wa_y[i];

        g_xxxxxy_0_xzz_0[i] = g_xxxxx_0_xzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyy_0[i] = 3.0 * g_xxxxx_0_yy_1[i] * fi_acd_0 + g_xxxxx_0_yyy_1[i] * wa_y[i];

        g_xxxxxy_0_yyz_0[i] = 2.0 * g_xxxxx_0_yz_1[i] * fi_acd_0 + g_xxxxx_0_yyz_1[i] * wa_y[i];

        g_xxxxxy_0_yzz_0[i] = g_xxxxx_0_zz_1[i] * fi_acd_0 + g_xxxxx_0_yzz_1[i] * wa_y[i];

        g_xxxxxy_0_zzz_0[i] = g_xxxxx_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 20-30 components of targeted buffer : ISF

    auto g_xxxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 20);

    auto g_xxxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 21);

    auto g_xxxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 22);

    auto g_xxxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 23);

    auto g_xxxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 24);

    auto g_xxxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 25);

    auto g_xxxxxz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 26);

    auto g_xxxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 27);

    auto g_xxxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 28);

    auto g_xxxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 29);

    #pragma omp simd aligned(g_xxxxx_0_xx_1, g_xxxxx_0_xxx_1, g_xxxxx_0_xxy_1, g_xxxxx_0_xxz_1, g_xxxxx_0_xy_1, g_xxxxx_0_xyy_1, g_xxxxx_0_xyz_1, g_xxxxx_0_xz_1, g_xxxxx_0_xzz_1, g_xxxxx_0_yy_1, g_xxxxx_0_yyy_1, g_xxxxx_0_yyz_1, g_xxxxx_0_yz_1, g_xxxxx_0_yzz_1, g_xxxxx_0_zz_1, g_xxxxx_0_zzz_1, g_xxxxxz_0_xxx_0, g_xxxxxz_0_xxy_0, g_xxxxxz_0_xxz_0, g_xxxxxz_0_xyy_0, g_xxxxxz_0_xyz_0, g_xxxxxz_0_xzz_0, g_xxxxxz_0_yyy_0, g_xxxxxz_0_yyz_0, g_xxxxxz_0_yzz_0, g_xxxxxz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xxx_0[i] = g_xxxxx_0_xxx_1[i] * wa_z[i];

        g_xxxxxz_0_xxy_0[i] = g_xxxxx_0_xxy_1[i] * wa_z[i];

        g_xxxxxz_0_xxz_0[i] = g_xxxxx_0_xx_1[i] * fi_acd_0 + g_xxxxx_0_xxz_1[i] * wa_z[i];

        g_xxxxxz_0_xyy_0[i] = g_xxxxx_0_xyy_1[i] * wa_z[i];

        g_xxxxxz_0_xyz_0[i] = g_xxxxx_0_xy_1[i] * fi_acd_0 + g_xxxxx_0_xyz_1[i] * wa_z[i];

        g_xxxxxz_0_xzz_0[i] = 2.0 * g_xxxxx_0_xz_1[i] * fi_acd_0 + g_xxxxx_0_xzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyy_0[i] = g_xxxxx_0_yyy_1[i] * wa_z[i];

        g_xxxxxz_0_yyz_0[i] = g_xxxxx_0_yy_1[i] * fi_acd_0 + g_xxxxx_0_yyz_1[i] * wa_z[i];

        g_xxxxxz_0_yzz_0[i] = 2.0 * g_xxxxx_0_yz_1[i] * fi_acd_0 + g_xxxxx_0_yzz_1[i] * wa_z[i];

        g_xxxxxz_0_zzz_0[i] = 3.0 * g_xxxxx_0_zz_1[i] * fi_acd_0 + g_xxxxx_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 30-40 components of targeted buffer : ISF

    auto g_xxxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 30);

    auto g_xxxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 31);

    auto g_xxxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 32);

    auto g_xxxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 33);

    auto g_xxxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 34);

    auto g_xxxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 35);

    auto g_xxxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 36);

    auto g_xxxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 37);

    auto g_xxxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 38);

    auto g_xxxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 39);

    #pragma omp simd aligned(g_xxxx_0_xxx_0, g_xxxx_0_xxx_1, g_xxxx_0_xxz_0, g_xxxx_0_xxz_1, g_xxxx_0_xzz_0, g_xxxx_0_xzz_1, g_xxxxy_0_xxx_1, g_xxxxy_0_xxz_1, g_xxxxy_0_xzz_1, g_xxxxyy_0_xxx_0, g_xxxxyy_0_xxy_0, g_xxxxyy_0_xxz_0, g_xxxxyy_0_xyy_0, g_xxxxyy_0_xyz_0, g_xxxxyy_0_xzz_0, g_xxxxyy_0_yyy_0, g_xxxxyy_0_yyz_0, g_xxxxyy_0_yzz_0, g_xxxxyy_0_zzz_0, g_xxxyy_0_xxy_1, g_xxxyy_0_xy_1, g_xxxyy_0_xyy_1, g_xxxyy_0_xyz_1, g_xxxyy_0_yy_1, g_xxxyy_0_yyy_1, g_xxxyy_0_yyz_1, g_xxxyy_0_yz_1, g_xxxyy_0_yzz_1, g_xxxyy_0_zzz_1, g_xxyy_0_xxy_0, g_xxyy_0_xxy_1, g_xxyy_0_xyy_0, g_xxyy_0_xyy_1, g_xxyy_0_xyz_0, g_xxyy_0_xyz_1, g_xxyy_0_yyy_0, g_xxyy_0_yyy_1, g_xxyy_0_yyz_0, g_xxyy_0_yyz_1, g_xxyy_0_yzz_0, g_xxyy_0_yzz_1, g_xxyy_0_zzz_0, g_xxyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xxx_0[i] = g_xxxx_0_xxx_0[i] * fbe_0 - g_xxxx_0_xxx_1[i] * fz_be_0 + g_xxxxy_0_xxx_1[i] * wa_y[i];

        g_xxxxyy_0_xxy_0[i] = 3.0 * g_xxyy_0_xxy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xy_1[i] * fi_acd_0 + g_xxxyy_0_xxy_1[i] * wa_x[i];

        g_xxxxyy_0_xxz_0[i] = g_xxxx_0_xxz_0[i] * fbe_0 - g_xxxx_0_xxz_1[i] * fz_be_0 + g_xxxxy_0_xxz_1[i] * wa_y[i];

        g_xxxxyy_0_xyy_0[i] = 3.0 * g_xxyy_0_xyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyy_1[i] * fz_be_0 + g_xxxyy_0_yy_1[i] * fi_acd_0 + g_xxxyy_0_xyy_1[i] * wa_x[i];

        g_xxxxyy_0_xyz_0[i] = 3.0 * g_xxyy_0_xyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyz_1[i] * fz_be_0 + g_xxxyy_0_yz_1[i] * fi_acd_0 + g_xxxyy_0_xyz_1[i] * wa_x[i];

        g_xxxxyy_0_xzz_0[i] = g_xxxx_0_xzz_0[i] * fbe_0 - g_xxxx_0_xzz_1[i] * fz_be_0 + g_xxxxy_0_xzz_1[i] * wa_y[i];

        g_xxxxyy_0_yyy_0[i] = 3.0 * g_xxyy_0_yyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyy_1[i] * fz_be_0 + g_xxxyy_0_yyy_1[i] * wa_x[i];

        g_xxxxyy_0_yyz_0[i] = 3.0 * g_xxyy_0_yyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyz_1[i] * fz_be_0 + g_xxxyy_0_yyz_1[i] * wa_x[i];

        g_xxxxyy_0_yzz_0[i] = 3.0 * g_xxyy_0_yzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yzz_1[i] * fz_be_0 + g_xxxyy_0_yzz_1[i] * wa_x[i];

        g_xxxxyy_0_zzz_0[i] = 3.0 * g_xxyy_0_zzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zzz_1[i] * fz_be_0 + g_xxxyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 40-50 components of targeted buffer : ISF

    auto g_xxxxyz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 40);

    auto g_xxxxyz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 41);

    auto g_xxxxyz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 42);

    auto g_xxxxyz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 43);

    auto g_xxxxyz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 44);

    auto g_xxxxyz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 45);

    auto g_xxxxyz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 46);

    auto g_xxxxyz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 47);

    auto g_xxxxyz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 48);

    auto g_xxxxyz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 49);

    #pragma omp simd aligned(g_xxxxy_0_xxy_1, g_xxxxy_0_xyy_1, g_xxxxy_0_yyy_1, g_xxxxyz_0_xxx_0, g_xxxxyz_0_xxy_0, g_xxxxyz_0_xxz_0, g_xxxxyz_0_xyy_0, g_xxxxyz_0_xyz_0, g_xxxxyz_0_xzz_0, g_xxxxyz_0_yyy_0, g_xxxxyz_0_yyz_0, g_xxxxyz_0_yzz_0, g_xxxxyz_0_zzz_0, g_xxxxz_0_xxx_1, g_xxxxz_0_xxz_1, g_xxxxz_0_xyz_1, g_xxxxz_0_xz_1, g_xxxxz_0_xzz_1, g_xxxxz_0_yyz_1, g_xxxxz_0_yz_1, g_xxxxz_0_yzz_1, g_xxxxz_0_zz_1, g_xxxxz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xxx_0[i] = g_xxxxz_0_xxx_1[i] * wa_y[i];

        g_xxxxyz_0_xxy_0[i] = g_xxxxy_0_xxy_1[i] * wa_z[i];

        g_xxxxyz_0_xxz_0[i] = g_xxxxz_0_xxz_1[i] * wa_y[i];

        g_xxxxyz_0_xyy_0[i] = g_xxxxy_0_xyy_1[i] * wa_z[i];

        g_xxxxyz_0_xyz_0[i] = g_xxxxz_0_xz_1[i] * fi_acd_0 + g_xxxxz_0_xyz_1[i] * wa_y[i];

        g_xxxxyz_0_xzz_0[i] = g_xxxxz_0_xzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyy_0[i] = g_xxxxy_0_yyy_1[i] * wa_z[i];

        g_xxxxyz_0_yyz_0[i] = 2.0 * g_xxxxz_0_yz_1[i] * fi_acd_0 + g_xxxxz_0_yyz_1[i] * wa_y[i];

        g_xxxxyz_0_yzz_0[i] = g_xxxxz_0_zz_1[i] * fi_acd_0 + g_xxxxz_0_yzz_1[i] * wa_y[i];

        g_xxxxyz_0_zzz_0[i] = g_xxxxz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 50-60 components of targeted buffer : ISF

    auto g_xxxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 50);

    auto g_xxxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 51);

    auto g_xxxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 52);

    auto g_xxxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 53);

    auto g_xxxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 54);

    auto g_xxxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 55);

    auto g_xxxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 56);

    auto g_xxxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 57);

    auto g_xxxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 58);

    auto g_xxxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 59);

    #pragma omp simd aligned(g_xxxx_0_xxx_0, g_xxxx_0_xxx_1, g_xxxx_0_xxy_0, g_xxxx_0_xxy_1, g_xxxx_0_xyy_0, g_xxxx_0_xyy_1, g_xxxxz_0_xxx_1, g_xxxxz_0_xxy_1, g_xxxxz_0_xyy_1, g_xxxxzz_0_xxx_0, g_xxxxzz_0_xxy_0, g_xxxxzz_0_xxz_0, g_xxxxzz_0_xyy_0, g_xxxxzz_0_xyz_0, g_xxxxzz_0_xzz_0, g_xxxxzz_0_yyy_0, g_xxxxzz_0_yyz_0, g_xxxxzz_0_yzz_0, g_xxxxzz_0_zzz_0, g_xxxzz_0_xxz_1, g_xxxzz_0_xyz_1, g_xxxzz_0_xz_1, g_xxxzz_0_xzz_1, g_xxxzz_0_yyy_1, g_xxxzz_0_yyz_1, g_xxxzz_0_yz_1, g_xxxzz_0_yzz_1, g_xxxzz_0_zz_1, g_xxxzz_0_zzz_1, g_xxzz_0_xxz_0, g_xxzz_0_xxz_1, g_xxzz_0_xyz_0, g_xxzz_0_xyz_1, g_xxzz_0_xzz_0, g_xxzz_0_xzz_1, g_xxzz_0_yyy_0, g_xxzz_0_yyy_1, g_xxzz_0_yyz_0, g_xxzz_0_yyz_1, g_xxzz_0_yzz_0, g_xxzz_0_yzz_1, g_xxzz_0_zzz_0, g_xxzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xxx_0[i] = g_xxxx_0_xxx_0[i] * fbe_0 - g_xxxx_0_xxx_1[i] * fz_be_0 + g_xxxxz_0_xxx_1[i] * wa_z[i];

        g_xxxxzz_0_xxy_0[i] = g_xxxx_0_xxy_0[i] * fbe_0 - g_xxxx_0_xxy_1[i] * fz_be_0 + g_xxxxz_0_xxy_1[i] * wa_z[i];

        g_xxxxzz_0_xxz_0[i] = 3.0 * g_xxzz_0_xxz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xz_1[i] * fi_acd_0 + g_xxxzz_0_xxz_1[i] * wa_x[i];

        g_xxxxzz_0_xyy_0[i] = g_xxxx_0_xyy_0[i] * fbe_0 - g_xxxx_0_xyy_1[i] * fz_be_0 + g_xxxxz_0_xyy_1[i] * wa_z[i];

        g_xxxxzz_0_xyz_0[i] = 3.0 * g_xxzz_0_xyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyz_1[i] * fz_be_0 + g_xxxzz_0_yz_1[i] * fi_acd_0 + g_xxxzz_0_xyz_1[i] * wa_x[i];

        g_xxxxzz_0_xzz_0[i] = 3.0 * g_xxzz_0_xzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xzz_1[i] * fz_be_0 + g_xxxzz_0_zz_1[i] * fi_acd_0 + g_xxxzz_0_xzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyy_0[i] = 3.0 * g_xxzz_0_yyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyy_1[i] * fz_be_0 + g_xxxzz_0_yyy_1[i] * wa_x[i];

        g_xxxxzz_0_yyz_0[i] = 3.0 * g_xxzz_0_yyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyz_1[i] * fz_be_0 + g_xxxzz_0_yyz_1[i] * wa_x[i];

        g_xxxxzz_0_yzz_0[i] = 3.0 * g_xxzz_0_yzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yzz_1[i] * fz_be_0 + g_xxxzz_0_yzz_1[i] * wa_x[i];

        g_xxxxzz_0_zzz_0[i] = 3.0 * g_xxzz_0_zzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zzz_1[i] * fz_be_0 + g_xxxzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 60-70 components of targeted buffer : ISF

    auto g_xxxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 60);

    auto g_xxxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 61);

    auto g_xxxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 62);

    auto g_xxxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 63);

    auto g_xxxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 64);

    auto g_xxxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 65);

    auto g_xxxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 66);

    auto g_xxxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 67);

    auto g_xxxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 68);

    auto g_xxxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 69);

    #pragma omp simd aligned(g_xxxy_0_xxx_0, g_xxxy_0_xxx_1, g_xxxy_0_xxz_0, g_xxxy_0_xxz_1, g_xxxy_0_xzz_0, g_xxxy_0_xzz_1, g_xxxyy_0_xxx_1, g_xxxyy_0_xxz_1, g_xxxyy_0_xzz_1, g_xxxyyy_0_xxx_0, g_xxxyyy_0_xxy_0, g_xxxyyy_0_xxz_0, g_xxxyyy_0_xyy_0, g_xxxyyy_0_xyz_0, g_xxxyyy_0_xzz_0, g_xxxyyy_0_yyy_0, g_xxxyyy_0_yyz_0, g_xxxyyy_0_yzz_0, g_xxxyyy_0_zzz_0, g_xxyyy_0_xxy_1, g_xxyyy_0_xy_1, g_xxyyy_0_xyy_1, g_xxyyy_0_xyz_1, g_xxyyy_0_yy_1, g_xxyyy_0_yyy_1, g_xxyyy_0_yyz_1, g_xxyyy_0_yz_1, g_xxyyy_0_yzz_1, g_xxyyy_0_zzz_1, g_xyyy_0_xxy_0, g_xyyy_0_xxy_1, g_xyyy_0_xyy_0, g_xyyy_0_xyy_1, g_xyyy_0_xyz_0, g_xyyy_0_xyz_1, g_xyyy_0_yyy_0, g_xyyy_0_yyy_1, g_xyyy_0_yyz_0, g_xyyy_0_yyz_1, g_xyyy_0_yzz_0, g_xyyy_0_yzz_1, g_xyyy_0_zzz_0, g_xyyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xxx_0[i] = 2.0 * g_xxxy_0_xxx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxx_1[i] * fz_be_0 + g_xxxyy_0_xxx_1[i] * wa_y[i];

        g_xxxyyy_0_xxy_0[i] = 2.0 * g_xyyy_0_xxy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xy_1[i] * fi_acd_0 + g_xxyyy_0_xxy_1[i] * wa_x[i];

        g_xxxyyy_0_xxz_0[i] = 2.0 * g_xxxy_0_xxz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxz_1[i] * fz_be_0 + g_xxxyy_0_xxz_1[i] * wa_y[i];

        g_xxxyyy_0_xyy_0[i] = 2.0 * g_xyyy_0_xyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyy_1[i] * fz_be_0 + g_xxyyy_0_yy_1[i] * fi_acd_0 + g_xxyyy_0_xyy_1[i] * wa_x[i];

        g_xxxyyy_0_xyz_0[i] = 2.0 * g_xyyy_0_xyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyz_1[i] * fz_be_0 + g_xxyyy_0_yz_1[i] * fi_acd_0 + g_xxyyy_0_xyz_1[i] * wa_x[i];

        g_xxxyyy_0_xzz_0[i] = 2.0 * g_xxxy_0_xzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xzz_1[i] * fz_be_0 + g_xxxyy_0_xzz_1[i] * wa_y[i];

        g_xxxyyy_0_yyy_0[i] = 2.0 * g_xyyy_0_yyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyy_1[i] * fz_be_0 + g_xxyyy_0_yyy_1[i] * wa_x[i];

        g_xxxyyy_0_yyz_0[i] = 2.0 * g_xyyy_0_yyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyz_1[i] * fz_be_0 + g_xxyyy_0_yyz_1[i] * wa_x[i];

        g_xxxyyy_0_yzz_0[i] = 2.0 * g_xyyy_0_yzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yzz_1[i] * fz_be_0 + g_xxyyy_0_yzz_1[i] * wa_x[i];

        g_xxxyyy_0_zzz_0[i] = 2.0 * g_xyyy_0_zzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zzz_1[i] * fz_be_0 + g_xxyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 70-80 components of targeted buffer : ISF

    auto g_xxxyyz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 70);

    auto g_xxxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 71);

    auto g_xxxyyz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 72);

    auto g_xxxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 73);

    auto g_xxxyyz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 74);

    auto g_xxxyyz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 75);

    auto g_xxxyyz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 76);

    auto g_xxxyyz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 77);

    auto g_xxxyyz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 78);

    auto g_xxxyyz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 79);

    #pragma omp simd aligned(g_xxxyy_0_xx_1, g_xxxyy_0_xxx_1, g_xxxyy_0_xxy_1, g_xxxyy_0_xxz_1, g_xxxyy_0_xy_1, g_xxxyy_0_xyy_1, g_xxxyy_0_xyz_1, g_xxxyy_0_xz_1, g_xxxyy_0_xzz_1, g_xxxyy_0_yy_1, g_xxxyy_0_yyy_1, g_xxxyy_0_yyz_1, g_xxxyy_0_yz_1, g_xxxyy_0_yzz_1, g_xxxyy_0_zz_1, g_xxxyy_0_zzz_1, g_xxxyyz_0_xxx_0, g_xxxyyz_0_xxy_0, g_xxxyyz_0_xxz_0, g_xxxyyz_0_xyy_0, g_xxxyyz_0_xyz_0, g_xxxyyz_0_xzz_0, g_xxxyyz_0_yyy_0, g_xxxyyz_0_yyz_0, g_xxxyyz_0_yzz_0, g_xxxyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xxx_0[i] = g_xxxyy_0_xxx_1[i] * wa_z[i];

        g_xxxyyz_0_xxy_0[i] = g_xxxyy_0_xxy_1[i] * wa_z[i];

        g_xxxyyz_0_xxz_0[i] = g_xxxyy_0_xx_1[i] * fi_acd_0 + g_xxxyy_0_xxz_1[i] * wa_z[i];

        g_xxxyyz_0_xyy_0[i] = g_xxxyy_0_xyy_1[i] * wa_z[i];

        g_xxxyyz_0_xyz_0[i] = g_xxxyy_0_xy_1[i] * fi_acd_0 + g_xxxyy_0_xyz_1[i] * wa_z[i];

        g_xxxyyz_0_xzz_0[i] = 2.0 * g_xxxyy_0_xz_1[i] * fi_acd_0 + g_xxxyy_0_xzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyy_0[i] = g_xxxyy_0_yyy_1[i] * wa_z[i];

        g_xxxyyz_0_yyz_0[i] = g_xxxyy_0_yy_1[i] * fi_acd_0 + g_xxxyy_0_yyz_1[i] * wa_z[i];

        g_xxxyyz_0_yzz_0[i] = 2.0 * g_xxxyy_0_yz_1[i] * fi_acd_0 + g_xxxyy_0_yzz_1[i] * wa_z[i];

        g_xxxyyz_0_zzz_0[i] = 3.0 * g_xxxyy_0_zz_1[i] * fi_acd_0 + g_xxxyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 80-90 components of targeted buffer : ISF

    auto g_xxxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 80);

    auto g_xxxyzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 81);

    auto g_xxxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 82);

    auto g_xxxyzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 83);

    auto g_xxxyzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 84);

    auto g_xxxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 85);

    auto g_xxxyzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 86);

    auto g_xxxyzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 87);

    auto g_xxxyzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 88);

    auto g_xxxyzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 89);

    #pragma omp simd aligned(g_xxxyzz_0_xxx_0, g_xxxyzz_0_xxy_0, g_xxxyzz_0_xxz_0, g_xxxyzz_0_xyy_0, g_xxxyzz_0_xyz_0, g_xxxyzz_0_xzz_0, g_xxxyzz_0_yyy_0, g_xxxyzz_0_yyz_0, g_xxxyzz_0_yzz_0, g_xxxyzz_0_zzz_0, g_xxxzz_0_xx_1, g_xxxzz_0_xxx_1, g_xxxzz_0_xxy_1, g_xxxzz_0_xxz_1, g_xxxzz_0_xy_1, g_xxxzz_0_xyy_1, g_xxxzz_0_xyz_1, g_xxxzz_0_xz_1, g_xxxzz_0_xzz_1, g_xxxzz_0_yy_1, g_xxxzz_0_yyy_1, g_xxxzz_0_yyz_1, g_xxxzz_0_yz_1, g_xxxzz_0_yzz_1, g_xxxzz_0_zz_1, g_xxxzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xxx_0[i] = g_xxxzz_0_xxx_1[i] * wa_y[i];

        g_xxxyzz_0_xxy_0[i] = g_xxxzz_0_xx_1[i] * fi_acd_0 + g_xxxzz_0_xxy_1[i] * wa_y[i];

        g_xxxyzz_0_xxz_0[i] = g_xxxzz_0_xxz_1[i] * wa_y[i];

        g_xxxyzz_0_xyy_0[i] = 2.0 * g_xxxzz_0_xy_1[i] * fi_acd_0 + g_xxxzz_0_xyy_1[i] * wa_y[i];

        g_xxxyzz_0_xyz_0[i] = g_xxxzz_0_xz_1[i] * fi_acd_0 + g_xxxzz_0_xyz_1[i] * wa_y[i];

        g_xxxyzz_0_xzz_0[i] = g_xxxzz_0_xzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyy_0[i] = 3.0 * g_xxxzz_0_yy_1[i] * fi_acd_0 + g_xxxzz_0_yyy_1[i] * wa_y[i];

        g_xxxyzz_0_yyz_0[i] = 2.0 * g_xxxzz_0_yz_1[i] * fi_acd_0 + g_xxxzz_0_yyz_1[i] * wa_y[i];

        g_xxxyzz_0_yzz_0[i] = g_xxxzz_0_zz_1[i] * fi_acd_0 + g_xxxzz_0_yzz_1[i] * wa_y[i];

        g_xxxyzz_0_zzz_0[i] = g_xxxzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 90-100 components of targeted buffer : ISF

    auto g_xxxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 90);

    auto g_xxxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 91);

    auto g_xxxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 92);

    auto g_xxxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 93);

    auto g_xxxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 94);

    auto g_xxxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 95);

    auto g_xxxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 96);

    auto g_xxxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 97);

    auto g_xxxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 98);

    auto g_xxxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 99);

    #pragma omp simd aligned(g_xxxz_0_xxx_0, g_xxxz_0_xxx_1, g_xxxz_0_xxy_0, g_xxxz_0_xxy_1, g_xxxz_0_xyy_0, g_xxxz_0_xyy_1, g_xxxzz_0_xxx_1, g_xxxzz_0_xxy_1, g_xxxzz_0_xyy_1, g_xxxzzz_0_xxx_0, g_xxxzzz_0_xxy_0, g_xxxzzz_0_xxz_0, g_xxxzzz_0_xyy_0, g_xxxzzz_0_xyz_0, g_xxxzzz_0_xzz_0, g_xxxzzz_0_yyy_0, g_xxxzzz_0_yyz_0, g_xxxzzz_0_yzz_0, g_xxxzzz_0_zzz_0, g_xxzzz_0_xxz_1, g_xxzzz_0_xyz_1, g_xxzzz_0_xz_1, g_xxzzz_0_xzz_1, g_xxzzz_0_yyy_1, g_xxzzz_0_yyz_1, g_xxzzz_0_yz_1, g_xxzzz_0_yzz_1, g_xxzzz_0_zz_1, g_xxzzz_0_zzz_1, g_xzzz_0_xxz_0, g_xzzz_0_xxz_1, g_xzzz_0_xyz_0, g_xzzz_0_xyz_1, g_xzzz_0_xzz_0, g_xzzz_0_xzz_1, g_xzzz_0_yyy_0, g_xzzz_0_yyy_1, g_xzzz_0_yyz_0, g_xzzz_0_yyz_1, g_xzzz_0_yzz_0, g_xzzz_0_yzz_1, g_xzzz_0_zzz_0, g_xzzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xxx_0[i] = 2.0 * g_xxxz_0_xxx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxx_1[i] * fz_be_0 + g_xxxzz_0_xxx_1[i] * wa_z[i];

        g_xxxzzz_0_xxy_0[i] = 2.0 * g_xxxz_0_xxy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxy_1[i] * fz_be_0 + g_xxxzz_0_xxy_1[i] * wa_z[i];

        g_xxxzzz_0_xxz_0[i] = 2.0 * g_xzzz_0_xxz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xz_1[i] * fi_acd_0 + g_xxzzz_0_xxz_1[i] * wa_x[i];

        g_xxxzzz_0_xyy_0[i] = 2.0 * g_xxxz_0_xyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xyy_1[i] * fz_be_0 + g_xxxzz_0_xyy_1[i] * wa_z[i];

        g_xxxzzz_0_xyz_0[i] = 2.0 * g_xzzz_0_xyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyz_1[i] * fz_be_0 + g_xxzzz_0_yz_1[i] * fi_acd_0 + g_xxzzz_0_xyz_1[i] * wa_x[i];

        g_xxxzzz_0_xzz_0[i] = 2.0 * g_xzzz_0_xzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xzz_1[i] * fz_be_0 + g_xxzzz_0_zz_1[i] * fi_acd_0 + g_xxzzz_0_xzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyy_0[i] = 2.0 * g_xzzz_0_yyy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyy_1[i] * fz_be_0 + g_xxzzz_0_yyy_1[i] * wa_x[i];

        g_xxxzzz_0_yyz_0[i] = 2.0 * g_xzzz_0_yyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyz_1[i] * fz_be_0 + g_xxzzz_0_yyz_1[i] * wa_x[i];

        g_xxxzzz_0_yzz_0[i] = 2.0 * g_xzzz_0_yzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yzz_1[i] * fz_be_0 + g_xxzzz_0_yzz_1[i] * wa_x[i];

        g_xxxzzz_0_zzz_0[i] = 2.0 * g_xzzz_0_zzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zzz_1[i] * fz_be_0 + g_xxzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 100-110 components of targeted buffer : ISF

    auto g_xxyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 100);

    auto g_xxyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 101);

    auto g_xxyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 102);

    auto g_xxyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 103);

    auto g_xxyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 104);

    auto g_xxyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 105);

    auto g_xxyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 106);

    auto g_xxyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 107);

    auto g_xxyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 108);

    auto g_xxyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 109);

    #pragma omp simd aligned(g_xxyy_0_xxx_0, g_xxyy_0_xxx_1, g_xxyy_0_xxz_0, g_xxyy_0_xxz_1, g_xxyy_0_xzz_0, g_xxyy_0_xzz_1, g_xxyyy_0_xxx_1, g_xxyyy_0_xxz_1, g_xxyyy_0_xzz_1, g_xxyyyy_0_xxx_0, g_xxyyyy_0_xxy_0, g_xxyyyy_0_xxz_0, g_xxyyyy_0_xyy_0, g_xxyyyy_0_xyz_0, g_xxyyyy_0_xzz_0, g_xxyyyy_0_yyy_0, g_xxyyyy_0_yyz_0, g_xxyyyy_0_yzz_0, g_xxyyyy_0_zzz_0, g_xyyyy_0_xxy_1, g_xyyyy_0_xy_1, g_xyyyy_0_xyy_1, g_xyyyy_0_xyz_1, g_xyyyy_0_yy_1, g_xyyyy_0_yyy_1, g_xyyyy_0_yyz_1, g_xyyyy_0_yz_1, g_xyyyy_0_yzz_1, g_xyyyy_0_zzz_1, g_yyyy_0_xxy_0, g_yyyy_0_xxy_1, g_yyyy_0_xyy_0, g_yyyy_0_xyy_1, g_yyyy_0_xyz_0, g_yyyy_0_xyz_1, g_yyyy_0_yyy_0, g_yyyy_0_yyy_1, g_yyyy_0_yyz_0, g_yyyy_0_yyz_1, g_yyyy_0_yzz_0, g_yyyy_0_yzz_1, g_yyyy_0_zzz_0, g_yyyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xxx_0[i] = 3.0 * g_xxyy_0_xxx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxx_1[i] * fz_be_0 + g_xxyyy_0_xxx_1[i] * wa_y[i];

        g_xxyyyy_0_xxy_0[i] = g_yyyy_0_xxy_0[i] * fbe_0 - g_yyyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xy_1[i] * fi_acd_0 + g_xyyyy_0_xxy_1[i] * wa_x[i];

        g_xxyyyy_0_xxz_0[i] = 3.0 * g_xxyy_0_xxz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxz_1[i] * fz_be_0 + g_xxyyy_0_xxz_1[i] * wa_y[i];

        g_xxyyyy_0_xyy_0[i] = g_yyyy_0_xyy_0[i] * fbe_0 - g_yyyy_0_xyy_1[i] * fz_be_0 + g_xyyyy_0_yy_1[i] * fi_acd_0 + g_xyyyy_0_xyy_1[i] * wa_x[i];

        g_xxyyyy_0_xyz_0[i] = g_yyyy_0_xyz_0[i] * fbe_0 - g_yyyy_0_xyz_1[i] * fz_be_0 + g_xyyyy_0_yz_1[i] * fi_acd_0 + g_xyyyy_0_xyz_1[i] * wa_x[i];

        g_xxyyyy_0_xzz_0[i] = 3.0 * g_xxyy_0_xzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xzz_1[i] * fz_be_0 + g_xxyyy_0_xzz_1[i] * wa_y[i];

        g_xxyyyy_0_yyy_0[i] = g_yyyy_0_yyy_0[i] * fbe_0 - g_yyyy_0_yyy_1[i] * fz_be_0 + g_xyyyy_0_yyy_1[i] * wa_x[i];

        g_xxyyyy_0_yyz_0[i] = g_yyyy_0_yyz_0[i] * fbe_0 - g_yyyy_0_yyz_1[i] * fz_be_0 + g_xyyyy_0_yyz_1[i] * wa_x[i];

        g_xxyyyy_0_yzz_0[i] = g_yyyy_0_yzz_0[i] * fbe_0 - g_yyyy_0_yzz_1[i] * fz_be_0 + g_xyyyy_0_yzz_1[i] * wa_x[i];

        g_xxyyyy_0_zzz_0[i] = g_yyyy_0_zzz_0[i] * fbe_0 - g_yyyy_0_zzz_1[i] * fz_be_0 + g_xyyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 110-120 components of targeted buffer : ISF

    auto g_xxyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 110);

    auto g_xxyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 111);

    auto g_xxyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 112);

    auto g_xxyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 113);

    auto g_xxyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 114);

    auto g_xxyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 115);

    auto g_xxyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 116);

    auto g_xxyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 117);

    auto g_xxyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 118);

    auto g_xxyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 119);

    #pragma omp simd aligned(g_xxyyy_0_xx_1, g_xxyyy_0_xxx_1, g_xxyyy_0_xxy_1, g_xxyyy_0_xxz_1, g_xxyyy_0_xy_1, g_xxyyy_0_xyy_1, g_xxyyy_0_xyz_1, g_xxyyy_0_xz_1, g_xxyyy_0_xzz_1, g_xxyyy_0_yy_1, g_xxyyy_0_yyy_1, g_xxyyy_0_yyz_1, g_xxyyy_0_yz_1, g_xxyyy_0_yzz_1, g_xxyyy_0_zz_1, g_xxyyy_0_zzz_1, g_xxyyyz_0_xxx_0, g_xxyyyz_0_xxy_0, g_xxyyyz_0_xxz_0, g_xxyyyz_0_xyy_0, g_xxyyyz_0_xyz_0, g_xxyyyz_0_xzz_0, g_xxyyyz_0_yyy_0, g_xxyyyz_0_yyz_0, g_xxyyyz_0_yzz_0, g_xxyyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xxx_0[i] = g_xxyyy_0_xxx_1[i] * wa_z[i];

        g_xxyyyz_0_xxy_0[i] = g_xxyyy_0_xxy_1[i] * wa_z[i];

        g_xxyyyz_0_xxz_0[i] = g_xxyyy_0_xx_1[i] * fi_acd_0 + g_xxyyy_0_xxz_1[i] * wa_z[i];

        g_xxyyyz_0_xyy_0[i] = g_xxyyy_0_xyy_1[i] * wa_z[i];

        g_xxyyyz_0_xyz_0[i] = g_xxyyy_0_xy_1[i] * fi_acd_0 + g_xxyyy_0_xyz_1[i] * wa_z[i];

        g_xxyyyz_0_xzz_0[i] = 2.0 * g_xxyyy_0_xz_1[i] * fi_acd_0 + g_xxyyy_0_xzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyy_0[i] = g_xxyyy_0_yyy_1[i] * wa_z[i];

        g_xxyyyz_0_yyz_0[i] = g_xxyyy_0_yy_1[i] * fi_acd_0 + g_xxyyy_0_yyz_1[i] * wa_z[i];

        g_xxyyyz_0_yzz_0[i] = 2.0 * g_xxyyy_0_yz_1[i] * fi_acd_0 + g_xxyyy_0_yzz_1[i] * wa_z[i];

        g_xxyyyz_0_zzz_0[i] = 3.0 * g_xxyyy_0_zz_1[i] * fi_acd_0 + g_xxyyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 120-130 components of targeted buffer : ISF

    auto g_xxyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 120);

    auto g_xxyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 121);

    auto g_xxyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 122);

    auto g_xxyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 123);

    auto g_xxyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 124);

    auto g_xxyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 125);

    auto g_xxyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 126);

    auto g_xxyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 127);

    auto g_xxyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 128);

    auto g_xxyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 129);

    #pragma omp simd aligned(g_xxyy_0_xxy_0, g_xxyy_0_xxy_1, g_xxyy_0_xyy_0, g_xxyy_0_xyy_1, g_xxyyz_0_xxy_1, g_xxyyz_0_xyy_1, g_xxyyzz_0_xxx_0, g_xxyyzz_0_xxy_0, g_xxyyzz_0_xxz_0, g_xxyyzz_0_xyy_0, g_xxyyzz_0_xyz_0, g_xxyyzz_0_xzz_0, g_xxyyzz_0_yyy_0, g_xxyyzz_0_yyz_0, g_xxyyzz_0_yzz_0, g_xxyyzz_0_zzz_0, g_xxyzz_0_xxx_1, g_xxyzz_0_xxz_1, g_xxyzz_0_xzz_1, g_xxzz_0_xxx_0, g_xxzz_0_xxx_1, g_xxzz_0_xxz_0, g_xxzz_0_xxz_1, g_xxzz_0_xzz_0, g_xxzz_0_xzz_1, g_xyyzz_0_xyz_1, g_xyyzz_0_yyy_1, g_xyyzz_0_yyz_1, g_xyyzz_0_yz_1, g_xyyzz_0_yzz_1, g_xyyzz_0_zzz_1, g_yyzz_0_xyz_0, g_yyzz_0_xyz_1, g_yyzz_0_yyy_0, g_yyzz_0_yyy_1, g_yyzz_0_yyz_0, g_yyzz_0_yyz_1, g_yyzz_0_yzz_0, g_yyzz_0_yzz_1, g_yyzz_0_zzz_0, g_yyzz_0_zzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzz_0_xxx_0[i] = g_xxzz_0_xxx_0[i] * fbe_0 - g_xxzz_0_xxx_1[i] * fz_be_0 + g_xxyzz_0_xxx_1[i] * wa_y[i];

        g_xxyyzz_0_xxy_0[i] = g_xxyy_0_xxy_0[i] * fbe_0 - g_xxyy_0_xxy_1[i] * fz_be_0 + g_xxyyz_0_xxy_1[i] * wa_z[i];

        g_xxyyzz_0_xxz_0[i] = g_xxzz_0_xxz_0[i] * fbe_0 - g_xxzz_0_xxz_1[i] * fz_be_0 + g_xxyzz_0_xxz_1[i] * wa_y[i];

        g_xxyyzz_0_xyy_0[i] = g_xxyy_0_xyy_0[i] * fbe_0 - g_xxyy_0_xyy_1[i] * fz_be_0 + g_xxyyz_0_xyy_1[i] * wa_z[i];

        g_xxyyzz_0_xyz_0[i] = g_yyzz_0_xyz_0[i] * fbe_0 - g_yyzz_0_xyz_1[i] * fz_be_0 + g_xyyzz_0_yz_1[i] * fi_acd_0 + g_xyyzz_0_xyz_1[i] * wa_x[i];

        g_xxyyzz_0_xzz_0[i] = g_xxzz_0_xzz_0[i] * fbe_0 - g_xxzz_0_xzz_1[i] * fz_be_0 + g_xxyzz_0_xzz_1[i] * wa_y[i];

        g_xxyyzz_0_yyy_0[i] = g_yyzz_0_yyy_0[i] * fbe_0 - g_yyzz_0_yyy_1[i] * fz_be_0 + g_xyyzz_0_yyy_1[i] * wa_x[i];

        g_xxyyzz_0_yyz_0[i] = g_yyzz_0_yyz_0[i] * fbe_0 - g_yyzz_0_yyz_1[i] * fz_be_0 + g_xyyzz_0_yyz_1[i] * wa_x[i];

        g_xxyyzz_0_yzz_0[i] = g_yyzz_0_yzz_0[i] * fbe_0 - g_yyzz_0_yzz_1[i] * fz_be_0 + g_xyyzz_0_yzz_1[i] * wa_x[i];

        g_xxyyzz_0_zzz_0[i] = g_yyzz_0_zzz_0[i] * fbe_0 - g_yyzz_0_zzz_1[i] * fz_be_0 + g_xyyzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 130-140 components of targeted buffer : ISF

    auto g_xxyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 130);

    auto g_xxyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 131);

    auto g_xxyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 132);

    auto g_xxyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 133);

    auto g_xxyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 134);

    auto g_xxyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 135);

    auto g_xxyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 136);

    auto g_xxyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 137);

    auto g_xxyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 138);

    auto g_xxyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 139);

    #pragma omp simd aligned(g_xxyzzz_0_xxx_0, g_xxyzzz_0_xxy_0, g_xxyzzz_0_xxz_0, g_xxyzzz_0_xyy_0, g_xxyzzz_0_xyz_0, g_xxyzzz_0_xzz_0, g_xxyzzz_0_yyy_0, g_xxyzzz_0_yyz_0, g_xxyzzz_0_yzz_0, g_xxyzzz_0_zzz_0, g_xxzzz_0_xx_1, g_xxzzz_0_xxx_1, g_xxzzz_0_xxy_1, g_xxzzz_0_xxz_1, g_xxzzz_0_xy_1, g_xxzzz_0_xyy_1, g_xxzzz_0_xyz_1, g_xxzzz_0_xz_1, g_xxzzz_0_xzz_1, g_xxzzz_0_yy_1, g_xxzzz_0_yyy_1, g_xxzzz_0_yyz_1, g_xxzzz_0_yz_1, g_xxzzz_0_yzz_1, g_xxzzz_0_zz_1, g_xxzzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xxx_0[i] = g_xxzzz_0_xxx_1[i] * wa_y[i];

        g_xxyzzz_0_xxy_0[i] = g_xxzzz_0_xx_1[i] * fi_acd_0 + g_xxzzz_0_xxy_1[i] * wa_y[i];

        g_xxyzzz_0_xxz_0[i] = g_xxzzz_0_xxz_1[i] * wa_y[i];

        g_xxyzzz_0_xyy_0[i] = 2.0 * g_xxzzz_0_xy_1[i] * fi_acd_0 + g_xxzzz_0_xyy_1[i] * wa_y[i];

        g_xxyzzz_0_xyz_0[i] = g_xxzzz_0_xz_1[i] * fi_acd_0 + g_xxzzz_0_xyz_1[i] * wa_y[i];

        g_xxyzzz_0_xzz_0[i] = g_xxzzz_0_xzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyy_0[i] = 3.0 * g_xxzzz_0_yy_1[i] * fi_acd_0 + g_xxzzz_0_yyy_1[i] * wa_y[i];

        g_xxyzzz_0_yyz_0[i] = 2.0 * g_xxzzz_0_yz_1[i] * fi_acd_0 + g_xxzzz_0_yyz_1[i] * wa_y[i];

        g_xxyzzz_0_yzz_0[i] = g_xxzzz_0_zz_1[i] * fi_acd_0 + g_xxzzz_0_yzz_1[i] * wa_y[i];

        g_xxyzzz_0_zzz_0[i] = g_xxzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 140-150 components of targeted buffer : ISF

    auto g_xxzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 140);

    auto g_xxzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 141);

    auto g_xxzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 142);

    auto g_xxzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 143);

    auto g_xxzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 144);

    auto g_xxzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 145);

    auto g_xxzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 146);

    auto g_xxzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 147);

    auto g_xxzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 148);

    auto g_xxzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 149);

    #pragma omp simd aligned(g_xxzz_0_xxx_0, g_xxzz_0_xxx_1, g_xxzz_0_xxy_0, g_xxzz_0_xxy_1, g_xxzz_0_xyy_0, g_xxzz_0_xyy_1, g_xxzzz_0_xxx_1, g_xxzzz_0_xxy_1, g_xxzzz_0_xyy_1, g_xxzzzz_0_xxx_0, g_xxzzzz_0_xxy_0, g_xxzzzz_0_xxz_0, g_xxzzzz_0_xyy_0, g_xxzzzz_0_xyz_0, g_xxzzzz_0_xzz_0, g_xxzzzz_0_yyy_0, g_xxzzzz_0_yyz_0, g_xxzzzz_0_yzz_0, g_xxzzzz_0_zzz_0, g_xzzzz_0_xxz_1, g_xzzzz_0_xyz_1, g_xzzzz_0_xz_1, g_xzzzz_0_xzz_1, g_xzzzz_0_yyy_1, g_xzzzz_0_yyz_1, g_xzzzz_0_yz_1, g_xzzzz_0_yzz_1, g_xzzzz_0_zz_1, g_xzzzz_0_zzz_1, g_zzzz_0_xxz_0, g_zzzz_0_xxz_1, g_zzzz_0_xyz_0, g_zzzz_0_xyz_1, g_zzzz_0_xzz_0, g_zzzz_0_xzz_1, g_zzzz_0_yyy_0, g_zzzz_0_yyy_1, g_zzzz_0_yyz_0, g_zzzz_0_yyz_1, g_zzzz_0_yzz_0, g_zzzz_0_yzz_1, g_zzzz_0_zzz_0, g_zzzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xxx_0[i] = 3.0 * g_xxzz_0_xxx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxx_1[i] * fz_be_0 + g_xxzzz_0_xxx_1[i] * wa_z[i];

        g_xxzzzz_0_xxy_0[i] = 3.0 * g_xxzz_0_xxy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxy_1[i] * fz_be_0 + g_xxzzz_0_xxy_1[i] * wa_z[i];

        g_xxzzzz_0_xxz_0[i] = g_zzzz_0_xxz_0[i] * fbe_0 - g_zzzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xz_1[i] * fi_acd_0 + g_xzzzz_0_xxz_1[i] * wa_x[i];

        g_xxzzzz_0_xyy_0[i] = 3.0 * g_xxzz_0_xyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyy_1[i] * fz_be_0 + g_xxzzz_0_xyy_1[i] * wa_z[i];

        g_xxzzzz_0_xyz_0[i] = g_zzzz_0_xyz_0[i] * fbe_0 - g_zzzz_0_xyz_1[i] * fz_be_0 + g_xzzzz_0_yz_1[i] * fi_acd_0 + g_xzzzz_0_xyz_1[i] * wa_x[i];

        g_xxzzzz_0_xzz_0[i] = g_zzzz_0_xzz_0[i] * fbe_0 - g_zzzz_0_xzz_1[i] * fz_be_0 + g_xzzzz_0_zz_1[i] * fi_acd_0 + g_xzzzz_0_xzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyy_0[i] = g_zzzz_0_yyy_0[i] * fbe_0 - g_zzzz_0_yyy_1[i] * fz_be_0 + g_xzzzz_0_yyy_1[i] * wa_x[i];

        g_xxzzzz_0_yyz_0[i] = g_zzzz_0_yyz_0[i] * fbe_0 - g_zzzz_0_yyz_1[i] * fz_be_0 + g_xzzzz_0_yyz_1[i] * wa_x[i];

        g_xxzzzz_0_yzz_0[i] = g_zzzz_0_yzz_0[i] * fbe_0 - g_zzzz_0_yzz_1[i] * fz_be_0 + g_xzzzz_0_yzz_1[i] * wa_x[i];

        g_xxzzzz_0_zzz_0[i] = g_zzzz_0_zzz_0[i] * fbe_0 - g_zzzz_0_zzz_1[i] * fz_be_0 + g_xzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 150-160 components of targeted buffer : ISF

    auto g_xyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 150);

    auto g_xyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 151);

    auto g_xyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 152);

    auto g_xyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 153);

    auto g_xyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 154);

    auto g_xyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 155);

    auto g_xyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 156);

    auto g_xyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 157);

    auto g_xyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 158);

    auto g_xyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 159);

    #pragma omp simd aligned(g_xyyyyy_0_xxx_0, g_xyyyyy_0_xxy_0, g_xyyyyy_0_xxz_0, g_xyyyyy_0_xyy_0, g_xyyyyy_0_xyz_0, g_xyyyyy_0_xzz_0, g_xyyyyy_0_yyy_0, g_xyyyyy_0_yyz_0, g_xyyyyy_0_yzz_0, g_xyyyyy_0_zzz_0, g_yyyyy_0_xx_1, g_yyyyy_0_xxx_1, g_yyyyy_0_xxy_1, g_yyyyy_0_xxz_1, g_yyyyy_0_xy_1, g_yyyyy_0_xyy_1, g_yyyyy_0_xyz_1, g_yyyyy_0_xz_1, g_yyyyy_0_xzz_1, g_yyyyy_0_yy_1, g_yyyyy_0_yyy_1, g_yyyyy_0_yyz_1, g_yyyyy_0_yz_1, g_yyyyy_0_yzz_1, g_yyyyy_0_zz_1, g_yyyyy_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xxx_0[i] = 3.0 * g_yyyyy_0_xx_1[i] * fi_acd_0 + g_yyyyy_0_xxx_1[i] * wa_x[i];

        g_xyyyyy_0_xxy_0[i] = 2.0 * g_yyyyy_0_xy_1[i] * fi_acd_0 + g_yyyyy_0_xxy_1[i] * wa_x[i];

        g_xyyyyy_0_xxz_0[i] = 2.0 * g_yyyyy_0_xz_1[i] * fi_acd_0 + g_yyyyy_0_xxz_1[i] * wa_x[i];

        g_xyyyyy_0_xyy_0[i] = g_yyyyy_0_yy_1[i] * fi_acd_0 + g_yyyyy_0_xyy_1[i] * wa_x[i];

        g_xyyyyy_0_xyz_0[i] = g_yyyyy_0_yz_1[i] * fi_acd_0 + g_yyyyy_0_xyz_1[i] * wa_x[i];

        g_xyyyyy_0_xzz_0[i] = g_yyyyy_0_zz_1[i] * fi_acd_0 + g_yyyyy_0_xzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyy_0[i] = g_yyyyy_0_yyy_1[i] * wa_x[i];

        g_xyyyyy_0_yyz_0[i] = g_yyyyy_0_yyz_1[i] * wa_x[i];

        g_xyyyyy_0_yzz_0[i] = g_yyyyy_0_yzz_1[i] * wa_x[i];

        g_xyyyyy_0_zzz_0[i] = g_yyyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 160-170 components of targeted buffer : ISF

    auto g_xyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 160);

    auto g_xyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 161);

    auto g_xyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 162);

    auto g_xyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 163);

    auto g_xyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 164);

    auto g_xyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 165);

    auto g_xyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 166);

    auto g_xyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 167);

    auto g_xyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 168);

    auto g_xyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 169);

    #pragma omp simd aligned(g_xyyyy_0_xxx_1, g_xyyyy_0_xxy_1, g_xyyyy_0_xyy_1, g_xyyyyz_0_xxx_0, g_xyyyyz_0_xxy_0, g_xyyyyz_0_xxz_0, g_xyyyyz_0_xyy_0, g_xyyyyz_0_xyz_0, g_xyyyyz_0_xzz_0, g_xyyyyz_0_yyy_0, g_xyyyyz_0_yyz_0, g_xyyyyz_0_yzz_0, g_xyyyyz_0_zzz_0, g_yyyyz_0_xxz_1, g_yyyyz_0_xyz_1, g_yyyyz_0_xz_1, g_yyyyz_0_xzz_1, g_yyyyz_0_yyy_1, g_yyyyz_0_yyz_1, g_yyyyz_0_yz_1, g_yyyyz_0_yzz_1, g_yyyyz_0_zz_1, g_yyyyz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xxx_0[i] = g_xyyyy_0_xxx_1[i] * wa_z[i];

        g_xyyyyz_0_xxy_0[i] = g_xyyyy_0_xxy_1[i] * wa_z[i];

        g_xyyyyz_0_xxz_0[i] = 2.0 * g_yyyyz_0_xz_1[i] * fi_acd_0 + g_yyyyz_0_xxz_1[i] * wa_x[i];

        g_xyyyyz_0_xyy_0[i] = g_xyyyy_0_xyy_1[i] * wa_z[i];

        g_xyyyyz_0_xyz_0[i] = g_yyyyz_0_yz_1[i] * fi_acd_0 + g_yyyyz_0_xyz_1[i] * wa_x[i];

        g_xyyyyz_0_xzz_0[i] = g_yyyyz_0_zz_1[i] * fi_acd_0 + g_yyyyz_0_xzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyy_0[i] = g_yyyyz_0_yyy_1[i] * wa_x[i];

        g_xyyyyz_0_yyz_0[i] = g_yyyyz_0_yyz_1[i] * wa_x[i];

        g_xyyyyz_0_yzz_0[i] = g_yyyyz_0_yzz_1[i] * wa_x[i];

        g_xyyyyz_0_zzz_0[i] = g_yyyyz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 170-180 components of targeted buffer : ISF

    auto g_xyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 170);

    auto g_xyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 171);

    auto g_xyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 172);

    auto g_xyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 173);

    auto g_xyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 174);

    auto g_xyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 175);

    auto g_xyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 176);

    auto g_xyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 177);

    auto g_xyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 178);

    auto g_xyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 179);

    #pragma omp simd aligned(g_xyyyzz_0_xxx_0, g_xyyyzz_0_xxy_0, g_xyyyzz_0_xxz_0, g_xyyyzz_0_xyy_0, g_xyyyzz_0_xyz_0, g_xyyyzz_0_xzz_0, g_xyyyzz_0_yyy_0, g_xyyyzz_0_yyz_0, g_xyyyzz_0_yzz_0, g_xyyyzz_0_zzz_0, g_yyyzz_0_xx_1, g_yyyzz_0_xxx_1, g_yyyzz_0_xxy_1, g_yyyzz_0_xxz_1, g_yyyzz_0_xy_1, g_yyyzz_0_xyy_1, g_yyyzz_0_xyz_1, g_yyyzz_0_xz_1, g_yyyzz_0_xzz_1, g_yyyzz_0_yy_1, g_yyyzz_0_yyy_1, g_yyyzz_0_yyz_1, g_yyyzz_0_yz_1, g_yyyzz_0_yzz_1, g_yyyzz_0_zz_1, g_yyyzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xxx_0[i] = 3.0 * g_yyyzz_0_xx_1[i] * fi_acd_0 + g_yyyzz_0_xxx_1[i] * wa_x[i];

        g_xyyyzz_0_xxy_0[i] = 2.0 * g_yyyzz_0_xy_1[i] * fi_acd_0 + g_yyyzz_0_xxy_1[i] * wa_x[i];

        g_xyyyzz_0_xxz_0[i] = 2.0 * g_yyyzz_0_xz_1[i] * fi_acd_0 + g_yyyzz_0_xxz_1[i] * wa_x[i];

        g_xyyyzz_0_xyy_0[i] = g_yyyzz_0_yy_1[i] * fi_acd_0 + g_yyyzz_0_xyy_1[i] * wa_x[i];

        g_xyyyzz_0_xyz_0[i] = g_yyyzz_0_yz_1[i] * fi_acd_0 + g_yyyzz_0_xyz_1[i] * wa_x[i];

        g_xyyyzz_0_xzz_0[i] = g_yyyzz_0_zz_1[i] * fi_acd_0 + g_yyyzz_0_xzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyy_0[i] = g_yyyzz_0_yyy_1[i] * wa_x[i];

        g_xyyyzz_0_yyz_0[i] = g_yyyzz_0_yyz_1[i] * wa_x[i];

        g_xyyyzz_0_yzz_0[i] = g_yyyzz_0_yzz_1[i] * wa_x[i];

        g_xyyyzz_0_zzz_0[i] = g_yyyzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 180-190 components of targeted buffer : ISF

    auto g_xyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 180);

    auto g_xyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 181);

    auto g_xyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 182);

    auto g_xyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 183);

    auto g_xyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 184);

    auto g_xyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 185);

    auto g_xyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 186);

    auto g_xyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 187);

    auto g_xyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 188);

    auto g_xyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 189);

    #pragma omp simd aligned(g_xyyzzz_0_xxx_0, g_xyyzzz_0_xxy_0, g_xyyzzz_0_xxz_0, g_xyyzzz_0_xyy_0, g_xyyzzz_0_xyz_0, g_xyyzzz_0_xzz_0, g_xyyzzz_0_yyy_0, g_xyyzzz_0_yyz_0, g_xyyzzz_0_yzz_0, g_xyyzzz_0_zzz_0, g_yyzzz_0_xx_1, g_yyzzz_0_xxx_1, g_yyzzz_0_xxy_1, g_yyzzz_0_xxz_1, g_yyzzz_0_xy_1, g_yyzzz_0_xyy_1, g_yyzzz_0_xyz_1, g_yyzzz_0_xz_1, g_yyzzz_0_xzz_1, g_yyzzz_0_yy_1, g_yyzzz_0_yyy_1, g_yyzzz_0_yyz_1, g_yyzzz_0_yz_1, g_yyzzz_0_yzz_1, g_yyzzz_0_zz_1, g_yyzzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xxx_0[i] = 3.0 * g_yyzzz_0_xx_1[i] * fi_acd_0 + g_yyzzz_0_xxx_1[i] * wa_x[i];

        g_xyyzzz_0_xxy_0[i] = 2.0 * g_yyzzz_0_xy_1[i] * fi_acd_0 + g_yyzzz_0_xxy_1[i] * wa_x[i];

        g_xyyzzz_0_xxz_0[i] = 2.0 * g_yyzzz_0_xz_1[i] * fi_acd_0 + g_yyzzz_0_xxz_1[i] * wa_x[i];

        g_xyyzzz_0_xyy_0[i] = g_yyzzz_0_yy_1[i] * fi_acd_0 + g_yyzzz_0_xyy_1[i] * wa_x[i];

        g_xyyzzz_0_xyz_0[i] = g_yyzzz_0_yz_1[i] * fi_acd_0 + g_yyzzz_0_xyz_1[i] * wa_x[i];

        g_xyyzzz_0_xzz_0[i] = g_yyzzz_0_zz_1[i] * fi_acd_0 + g_yyzzz_0_xzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyy_0[i] = g_yyzzz_0_yyy_1[i] * wa_x[i];

        g_xyyzzz_0_yyz_0[i] = g_yyzzz_0_yyz_1[i] * wa_x[i];

        g_xyyzzz_0_yzz_0[i] = g_yyzzz_0_yzz_1[i] * wa_x[i];

        g_xyyzzz_0_zzz_0[i] = g_yyzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 190-200 components of targeted buffer : ISF

    auto g_xyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 190);

    auto g_xyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 191);

    auto g_xyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 192);

    auto g_xyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 193);

    auto g_xyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 194);

    auto g_xyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 195);

    auto g_xyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 196);

    auto g_xyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 197);

    auto g_xyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 198);

    auto g_xyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 199);

    #pragma omp simd aligned(g_xyzzzz_0_xxx_0, g_xyzzzz_0_xxy_0, g_xyzzzz_0_xxz_0, g_xyzzzz_0_xyy_0, g_xyzzzz_0_xyz_0, g_xyzzzz_0_xzz_0, g_xyzzzz_0_yyy_0, g_xyzzzz_0_yyz_0, g_xyzzzz_0_yzz_0, g_xyzzzz_0_zzz_0, g_xzzzz_0_xxx_1, g_xzzzz_0_xxz_1, g_xzzzz_0_xzz_1, g_yzzzz_0_xxy_1, g_yzzzz_0_xy_1, g_yzzzz_0_xyy_1, g_yzzzz_0_xyz_1, g_yzzzz_0_yy_1, g_yzzzz_0_yyy_1, g_yzzzz_0_yyz_1, g_yzzzz_0_yz_1, g_yzzzz_0_yzz_1, g_yzzzz_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xxx_0[i] = g_xzzzz_0_xxx_1[i] * wa_y[i];

        g_xyzzzz_0_xxy_0[i] = 2.0 * g_yzzzz_0_xy_1[i] * fi_acd_0 + g_yzzzz_0_xxy_1[i] * wa_x[i];

        g_xyzzzz_0_xxz_0[i] = g_xzzzz_0_xxz_1[i] * wa_y[i];

        g_xyzzzz_0_xyy_0[i] = g_yzzzz_0_yy_1[i] * fi_acd_0 + g_yzzzz_0_xyy_1[i] * wa_x[i];

        g_xyzzzz_0_xyz_0[i] = g_yzzzz_0_yz_1[i] * fi_acd_0 + g_yzzzz_0_xyz_1[i] * wa_x[i];

        g_xyzzzz_0_xzz_0[i] = g_xzzzz_0_xzz_1[i] * wa_y[i];

        g_xyzzzz_0_yyy_0[i] = g_yzzzz_0_yyy_1[i] * wa_x[i];

        g_xyzzzz_0_yyz_0[i] = g_yzzzz_0_yyz_1[i] * wa_x[i];

        g_xyzzzz_0_yzz_0[i] = g_yzzzz_0_yzz_1[i] * wa_x[i];

        g_xyzzzz_0_zzz_0[i] = g_yzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 200-210 components of targeted buffer : ISF

    auto g_xzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 200);

    auto g_xzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 201);

    auto g_xzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 202);

    auto g_xzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 203);

    auto g_xzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 204);

    auto g_xzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 205);

    auto g_xzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 206);

    auto g_xzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 207);

    auto g_xzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 208);

    auto g_xzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 209);

    #pragma omp simd aligned(g_xzzzzz_0_xxx_0, g_xzzzzz_0_xxy_0, g_xzzzzz_0_xxz_0, g_xzzzzz_0_xyy_0, g_xzzzzz_0_xyz_0, g_xzzzzz_0_xzz_0, g_xzzzzz_0_yyy_0, g_xzzzzz_0_yyz_0, g_xzzzzz_0_yzz_0, g_xzzzzz_0_zzz_0, g_zzzzz_0_xx_1, g_zzzzz_0_xxx_1, g_zzzzz_0_xxy_1, g_zzzzz_0_xxz_1, g_zzzzz_0_xy_1, g_zzzzz_0_xyy_1, g_zzzzz_0_xyz_1, g_zzzzz_0_xz_1, g_zzzzz_0_xzz_1, g_zzzzz_0_yy_1, g_zzzzz_0_yyy_1, g_zzzzz_0_yyz_1, g_zzzzz_0_yz_1, g_zzzzz_0_yzz_1, g_zzzzz_0_zz_1, g_zzzzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xxx_0[i] = 3.0 * g_zzzzz_0_xx_1[i] * fi_acd_0 + g_zzzzz_0_xxx_1[i] * wa_x[i];

        g_xzzzzz_0_xxy_0[i] = 2.0 * g_zzzzz_0_xy_1[i] * fi_acd_0 + g_zzzzz_0_xxy_1[i] * wa_x[i];

        g_xzzzzz_0_xxz_0[i] = 2.0 * g_zzzzz_0_xz_1[i] * fi_acd_0 + g_zzzzz_0_xxz_1[i] * wa_x[i];

        g_xzzzzz_0_xyy_0[i] = g_zzzzz_0_yy_1[i] * fi_acd_0 + g_zzzzz_0_xyy_1[i] * wa_x[i];

        g_xzzzzz_0_xyz_0[i] = g_zzzzz_0_yz_1[i] * fi_acd_0 + g_zzzzz_0_xyz_1[i] * wa_x[i];

        g_xzzzzz_0_xzz_0[i] = g_zzzzz_0_zz_1[i] * fi_acd_0 + g_zzzzz_0_xzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyy_0[i] = g_zzzzz_0_yyy_1[i] * wa_x[i];

        g_xzzzzz_0_yyz_0[i] = g_zzzzz_0_yyz_1[i] * wa_x[i];

        g_xzzzzz_0_yzz_0[i] = g_zzzzz_0_yzz_1[i] * wa_x[i];

        g_xzzzzz_0_zzz_0[i] = g_zzzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 210-220 components of targeted buffer : ISF

    auto g_yyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 210);

    auto g_yyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 211);

    auto g_yyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 212);

    auto g_yyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 213);

    auto g_yyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 214);

    auto g_yyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 215);

    auto g_yyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 216);

    auto g_yyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 217);

    auto g_yyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 218);

    auto g_yyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 219);

    #pragma omp simd aligned(g_yyyy_0_xxx_0, g_yyyy_0_xxx_1, g_yyyy_0_xxy_0, g_yyyy_0_xxy_1, g_yyyy_0_xxz_0, g_yyyy_0_xxz_1, g_yyyy_0_xyy_0, g_yyyy_0_xyy_1, g_yyyy_0_xyz_0, g_yyyy_0_xyz_1, g_yyyy_0_xzz_0, g_yyyy_0_xzz_1, g_yyyy_0_yyy_0, g_yyyy_0_yyy_1, g_yyyy_0_yyz_0, g_yyyy_0_yyz_1, g_yyyy_0_yzz_0, g_yyyy_0_yzz_1, g_yyyy_0_zzz_0, g_yyyy_0_zzz_1, g_yyyyy_0_xx_1, g_yyyyy_0_xxx_1, g_yyyyy_0_xxy_1, g_yyyyy_0_xxz_1, g_yyyyy_0_xy_1, g_yyyyy_0_xyy_1, g_yyyyy_0_xyz_1, g_yyyyy_0_xz_1, g_yyyyy_0_xzz_1, g_yyyyy_0_yy_1, g_yyyyy_0_yyy_1, g_yyyyy_0_yyz_1, g_yyyyy_0_yz_1, g_yyyyy_0_yzz_1, g_yyyyy_0_zz_1, g_yyyyy_0_zzz_1, g_yyyyyy_0_xxx_0, g_yyyyyy_0_xxy_0, g_yyyyyy_0_xxz_0, g_yyyyyy_0_xyy_0, g_yyyyyy_0_xyz_0, g_yyyyyy_0_xzz_0, g_yyyyyy_0_yyy_0, g_yyyyyy_0_yyz_0, g_yyyyyy_0_yzz_0, g_yyyyyy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xxx_0[i] = 5.0 * g_yyyy_0_xxx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxx_1[i] * fz_be_0 + g_yyyyy_0_xxx_1[i] * wa_y[i];

        g_yyyyyy_0_xxy_0[i] = 5.0 * g_yyyy_0_xxy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxy_1[i] * fz_be_0 + g_yyyyy_0_xx_1[i] * fi_acd_0 + g_yyyyy_0_xxy_1[i] * wa_y[i];

        g_yyyyyy_0_xxz_0[i] = 5.0 * g_yyyy_0_xxz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxz_1[i] * fz_be_0 + g_yyyyy_0_xxz_1[i] * wa_y[i];

        g_yyyyyy_0_xyy_0[i] = 5.0 * g_yyyy_0_xyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xy_1[i] * fi_acd_0 + g_yyyyy_0_xyy_1[i] * wa_y[i];

        g_yyyyyy_0_xyz_0[i] = 5.0 * g_yyyy_0_xyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyz_1[i] * fz_be_0 + g_yyyyy_0_xz_1[i] * fi_acd_0 + g_yyyyy_0_xyz_1[i] * wa_y[i];

        g_yyyyyy_0_xzz_0[i] = 5.0 * g_yyyy_0_xzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xzz_1[i] * fz_be_0 + g_yyyyy_0_xzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyy_0[i] = 5.0 * g_yyyy_0_yyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_yy_1[i] * fi_acd_0 + g_yyyyy_0_yyy_1[i] * wa_y[i];

        g_yyyyyy_0_yyz_0[i] = 5.0 * g_yyyy_0_yyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_yz_1[i] * fi_acd_0 + g_yyyyy_0_yyz_1[i] * wa_y[i];

        g_yyyyyy_0_yzz_0[i] = 5.0 * g_yyyy_0_yzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yzz_1[i] * fz_be_0 + g_yyyyy_0_zz_1[i] * fi_acd_0 + g_yyyyy_0_yzz_1[i] * wa_y[i];

        g_yyyyyy_0_zzz_0[i] = 5.0 * g_yyyy_0_zzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zzz_1[i] * fz_be_0 + g_yyyyy_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 220-230 components of targeted buffer : ISF

    auto g_yyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 220);

    auto g_yyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 221);

    auto g_yyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 222);

    auto g_yyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 223);

    auto g_yyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 224);

    auto g_yyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 225);

    auto g_yyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 226);

    auto g_yyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 227);

    auto g_yyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 228);

    auto g_yyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 229);

    #pragma omp simd aligned(g_yyyyy_0_xx_1, g_yyyyy_0_xxx_1, g_yyyyy_0_xxy_1, g_yyyyy_0_xxz_1, g_yyyyy_0_xy_1, g_yyyyy_0_xyy_1, g_yyyyy_0_xyz_1, g_yyyyy_0_xz_1, g_yyyyy_0_xzz_1, g_yyyyy_0_yy_1, g_yyyyy_0_yyy_1, g_yyyyy_0_yyz_1, g_yyyyy_0_yz_1, g_yyyyy_0_yzz_1, g_yyyyy_0_zz_1, g_yyyyy_0_zzz_1, g_yyyyyz_0_xxx_0, g_yyyyyz_0_xxy_0, g_yyyyyz_0_xxz_0, g_yyyyyz_0_xyy_0, g_yyyyyz_0_xyz_0, g_yyyyyz_0_xzz_0, g_yyyyyz_0_yyy_0, g_yyyyyz_0_yyz_0, g_yyyyyz_0_yzz_0, g_yyyyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xxx_0[i] = g_yyyyy_0_xxx_1[i] * wa_z[i];

        g_yyyyyz_0_xxy_0[i] = g_yyyyy_0_xxy_1[i] * wa_z[i];

        g_yyyyyz_0_xxz_0[i] = g_yyyyy_0_xx_1[i] * fi_acd_0 + g_yyyyy_0_xxz_1[i] * wa_z[i];

        g_yyyyyz_0_xyy_0[i] = g_yyyyy_0_xyy_1[i] * wa_z[i];

        g_yyyyyz_0_xyz_0[i] = g_yyyyy_0_xy_1[i] * fi_acd_0 + g_yyyyy_0_xyz_1[i] * wa_z[i];

        g_yyyyyz_0_xzz_0[i] = 2.0 * g_yyyyy_0_xz_1[i] * fi_acd_0 + g_yyyyy_0_xzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyy_0[i] = g_yyyyy_0_yyy_1[i] * wa_z[i];

        g_yyyyyz_0_yyz_0[i] = g_yyyyy_0_yy_1[i] * fi_acd_0 + g_yyyyy_0_yyz_1[i] * wa_z[i];

        g_yyyyyz_0_yzz_0[i] = 2.0 * g_yyyyy_0_yz_1[i] * fi_acd_0 + g_yyyyy_0_yzz_1[i] * wa_z[i];

        g_yyyyyz_0_zzz_0[i] = 3.0 * g_yyyyy_0_zz_1[i] * fi_acd_0 + g_yyyyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 230-240 components of targeted buffer : ISF

    auto g_yyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 230);

    auto g_yyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 231);

    auto g_yyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 232);

    auto g_yyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 233);

    auto g_yyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 234);

    auto g_yyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 235);

    auto g_yyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 236);

    auto g_yyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 237);

    auto g_yyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 238);

    auto g_yyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 239);

    #pragma omp simd aligned(g_yyyy_0_xxy_0, g_yyyy_0_xxy_1, g_yyyy_0_xyy_0, g_yyyy_0_xyy_1, g_yyyy_0_yyy_0, g_yyyy_0_yyy_1, g_yyyyz_0_xxy_1, g_yyyyz_0_xyy_1, g_yyyyz_0_yyy_1, g_yyyyzz_0_xxx_0, g_yyyyzz_0_xxy_0, g_yyyyzz_0_xxz_0, g_yyyyzz_0_xyy_0, g_yyyyzz_0_xyz_0, g_yyyyzz_0_xzz_0, g_yyyyzz_0_yyy_0, g_yyyyzz_0_yyz_0, g_yyyyzz_0_yzz_0, g_yyyyzz_0_zzz_0, g_yyyzz_0_xxx_1, g_yyyzz_0_xxz_1, g_yyyzz_0_xyz_1, g_yyyzz_0_xz_1, g_yyyzz_0_xzz_1, g_yyyzz_0_yyz_1, g_yyyzz_0_yz_1, g_yyyzz_0_yzz_1, g_yyyzz_0_zz_1, g_yyyzz_0_zzz_1, g_yyzz_0_xxx_0, g_yyzz_0_xxx_1, g_yyzz_0_xxz_0, g_yyzz_0_xxz_1, g_yyzz_0_xyz_0, g_yyzz_0_xyz_1, g_yyzz_0_xzz_0, g_yyzz_0_xzz_1, g_yyzz_0_yyz_0, g_yyzz_0_yyz_1, g_yyzz_0_yzz_0, g_yyzz_0_yzz_1, g_yyzz_0_zzz_0, g_yyzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xxx_0[i] = 3.0 * g_yyzz_0_xxx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxx_1[i] * fz_be_0 + g_yyyzz_0_xxx_1[i] * wa_y[i];

        g_yyyyzz_0_xxy_0[i] = g_yyyy_0_xxy_0[i] * fbe_0 - g_yyyy_0_xxy_1[i] * fz_be_0 + g_yyyyz_0_xxy_1[i] * wa_z[i];

        g_yyyyzz_0_xxz_0[i] = 3.0 * g_yyzz_0_xxz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxz_1[i] * fz_be_0 + g_yyyzz_0_xxz_1[i] * wa_y[i];

        g_yyyyzz_0_xyy_0[i] = g_yyyy_0_xyy_0[i] * fbe_0 - g_yyyy_0_xyy_1[i] * fz_be_0 + g_yyyyz_0_xyy_1[i] * wa_z[i];

        g_yyyyzz_0_xyz_0[i] = 3.0 * g_yyzz_0_xyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyz_1[i] * fz_be_0 + g_yyyzz_0_xz_1[i] * fi_acd_0 + g_yyyzz_0_xyz_1[i] * wa_y[i];

        g_yyyyzz_0_xzz_0[i] = 3.0 * g_yyzz_0_xzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xzz_1[i] * fz_be_0 + g_yyyzz_0_xzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyy_0[i] = g_yyyy_0_yyy_0[i] * fbe_0 - g_yyyy_0_yyy_1[i] * fz_be_0 + g_yyyyz_0_yyy_1[i] * wa_z[i];

        g_yyyyzz_0_yyz_0[i] = 3.0 * g_yyzz_0_yyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_yz_1[i] * fi_acd_0 + g_yyyzz_0_yyz_1[i] * wa_y[i];

        g_yyyyzz_0_yzz_0[i] = 3.0 * g_yyzz_0_yzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yzz_1[i] * fz_be_0 + g_yyyzz_0_zz_1[i] * fi_acd_0 + g_yyyzz_0_yzz_1[i] * wa_y[i];

        g_yyyyzz_0_zzz_0[i] = 3.0 * g_yyzz_0_zzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zzz_1[i] * fz_be_0 + g_yyyzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 240-250 components of targeted buffer : ISF

    auto g_yyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 240);

    auto g_yyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 241);

    auto g_yyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 242);

    auto g_yyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 243);

    auto g_yyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 244);

    auto g_yyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 245);

    auto g_yyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 246);

    auto g_yyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 247);

    auto g_yyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 248);

    auto g_yyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 249);

    #pragma omp simd aligned(g_yyyz_0_xxy_0, g_yyyz_0_xxy_1, g_yyyz_0_xyy_0, g_yyyz_0_xyy_1, g_yyyz_0_yyy_0, g_yyyz_0_yyy_1, g_yyyzz_0_xxy_1, g_yyyzz_0_xyy_1, g_yyyzz_0_yyy_1, g_yyyzzz_0_xxx_0, g_yyyzzz_0_xxy_0, g_yyyzzz_0_xxz_0, g_yyyzzz_0_xyy_0, g_yyyzzz_0_xyz_0, g_yyyzzz_0_xzz_0, g_yyyzzz_0_yyy_0, g_yyyzzz_0_yyz_0, g_yyyzzz_0_yzz_0, g_yyyzzz_0_zzz_0, g_yyzzz_0_xxx_1, g_yyzzz_0_xxz_1, g_yyzzz_0_xyz_1, g_yyzzz_0_xz_1, g_yyzzz_0_xzz_1, g_yyzzz_0_yyz_1, g_yyzzz_0_yz_1, g_yyzzz_0_yzz_1, g_yyzzz_0_zz_1, g_yyzzz_0_zzz_1, g_yzzz_0_xxx_0, g_yzzz_0_xxx_1, g_yzzz_0_xxz_0, g_yzzz_0_xxz_1, g_yzzz_0_xyz_0, g_yzzz_0_xyz_1, g_yzzz_0_xzz_0, g_yzzz_0_xzz_1, g_yzzz_0_yyz_0, g_yzzz_0_yyz_1, g_yzzz_0_yzz_0, g_yzzz_0_yzz_1, g_yzzz_0_zzz_0, g_yzzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xxx_0[i] = 2.0 * g_yzzz_0_xxx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxx_1[i] * fz_be_0 + g_yyzzz_0_xxx_1[i] * wa_y[i];

        g_yyyzzz_0_xxy_0[i] = 2.0 * g_yyyz_0_xxy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxy_1[i] * fz_be_0 + g_yyyzz_0_xxy_1[i] * wa_z[i];

        g_yyyzzz_0_xxz_0[i] = 2.0 * g_yzzz_0_xxz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxz_1[i] * fz_be_0 + g_yyzzz_0_xxz_1[i] * wa_y[i];

        g_yyyzzz_0_xyy_0[i] = 2.0 * g_yyyz_0_xyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xyy_1[i] * fz_be_0 + g_yyyzz_0_xyy_1[i] * wa_z[i];

        g_yyyzzz_0_xyz_0[i] = 2.0 * g_yzzz_0_xyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyz_1[i] * fz_be_0 + g_yyzzz_0_xz_1[i] * fi_acd_0 + g_yyzzz_0_xyz_1[i] * wa_y[i];

        g_yyyzzz_0_xzz_0[i] = 2.0 * g_yzzz_0_xzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xzz_1[i] * fz_be_0 + g_yyzzz_0_xzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyy_0[i] = 2.0 * g_yyyz_0_yyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yyy_1[i] * fz_be_0 + g_yyyzz_0_yyy_1[i] * wa_z[i];

        g_yyyzzz_0_yyz_0[i] = 2.0 * g_yzzz_0_yyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_yz_1[i] * fi_acd_0 + g_yyzzz_0_yyz_1[i] * wa_y[i];

        g_yyyzzz_0_yzz_0[i] = 2.0 * g_yzzz_0_yzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yzz_1[i] * fz_be_0 + g_yyzzz_0_zz_1[i] * fi_acd_0 + g_yyzzz_0_yzz_1[i] * wa_y[i];

        g_yyyzzz_0_zzz_0[i] = 2.0 * g_yzzz_0_zzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zzz_1[i] * fz_be_0 + g_yyzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 250-260 components of targeted buffer : ISF

    auto g_yyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 250);

    auto g_yyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 251);

    auto g_yyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 252);

    auto g_yyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 253);

    auto g_yyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 254);

    auto g_yyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 255);

    auto g_yyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 256);

    auto g_yyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 257);

    auto g_yyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 258);

    auto g_yyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 259);

    #pragma omp simd aligned(g_yyzz_0_xxy_0, g_yyzz_0_xxy_1, g_yyzz_0_xyy_0, g_yyzz_0_xyy_1, g_yyzz_0_yyy_0, g_yyzz_0_yyy_1, g_yyzzz_0_xxy_1, g_yyzzz_0_xyy_1, g_yyzzz_0_yyy_1, g_yyzzzz_0_xxx_0, g_yyzzzz_0_xxy_0, g_yyzzzz_0_xxz_0, g_yyzzzz_0_xyy_0, g_yyzzzz_0_xyz_0, g_yyzzzz_0_xzz_0, g_yyzzzz_0_yyy_0, g_yyzzzz_0_yyz_0, g_yyzzzz_0_yzz_0, g_yyzzzz_0_zzz_0, g_yzzzz_0_xxx_1, g_yzzzz_0_xxz_1, g_yzzzz_0_xyz_1, g_yzzzz_0_xz_1, g_yzzzz_0_xzz_1, g_yzzzz_0_yyz_1, g_yzzzz_0_yz_1, g_yzzzz_0_yzz_1, g_yzzzz_0_zz_1, g_yzzzz_0_zzz_1, g_zzzz_0_xxx_0, g_zzzz_0_xxx_1, g_zzzz_0_xxz_0, g_zzzz_0_xxz_1, g_zzzz_0_xyz_0, g_zzzz_0_xyz_1, g_zzzz_0_xzz_0, g_zzzz_0_xzz_1, g_zzzz_0_yyz_0, g_zzzz_0_yyz_1, g_zzzz_0_yzz_0, g_zzzz_0_yzz_1, g_zzzz_0_zzz_0, g_zzzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xxx_0[i] = g_zzzz_0_xxx_0[i] * fbe_0 - g_zzzz_0_xxx_1[i] * fz_be_0 + g_yzzzz_0_xxx_1[i] * wa_y[i];

        g_yyzzzz_0_xxy_0[i] = 3.0 * g_yyzz_0_xxy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxy_1[i] * fz_be_0 + g_yyzzz_0_xxy_1[i] * wa_z[i];

        g_yyzzzz_0_xxz_0[i] = g_zzzz_0_xxz_0[i] * fbe_0 - g_zzzz_0_xxz_1[i] * fz_be_0 + g_yzzzz_0_xxz_1[i] * wa_y[i];

        g_yyzzzz_0_xyy_0[i] = 3.0 * g_yyzz_0_xyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyy_1[i] * fz_be_0 + g_yyzzz_0_xyy_1[i] * wa_z[i];

        g_yyzzzz_0_xyz_0[i] = g_zzzz_0_xyz_0[i] * fbe_0 - g_zzzz_0_xyz_1[i] * fz_be_0 + g_yzzzz_0_xz_1[i] * fi_acd_0 + g_yzzzz_0_xyz_1[i] * wa_y[i];

        g_yyzzzz_0_xzz_0[i] = g_zzzz_0_xzz_0[i] * fbe_0 - g_zzzz_0_xzz_1[i] * fz_be_0 + g_yzzzz_0_xzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyy_0[i] = 3.0 * g_yyzz_0_yyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyy_1[i] * fz_be_0 + g_yyzzz_0_yyy_1[i] * wa_z[i];

        g_yyzzzz_0_yyz_0[i] = g_zzzz_0_yyz_0[i] * fbe_0 - g_zzzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_yz_1[i] * fi_acd_0 + g_yzzzz_0_yyz_1[i] * wa_y[i];

        g_yyzzzz_0_yzz_0[i] = g_zzzz_0_yzz_0[i] * fbe_0 - g_zzzz_0_yzz_1[i] * fz_be_0 + g_yzzzz_0_zz_1[i] * fi_acd_0 + g_yzzzz_0_yzz_1[i] * wa_y[i];

        g_yyzzzz_0_zzz_0[i] = g_zzzz_0_zzz_0[i] * fbe_0 - g_zzzz_0_zzz_1[i] * fz_be_0 + g_yzzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 260-270 components of targeted buffer : ISF

    auto g_yzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 260);

    auto g_yzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 261);

    auto g_yzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 262);

    auto g_yzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 263);

    auto g_yzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 264);

    auto g_yzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 265);

    auto g_yzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 266);

    auto g_yzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 267);

    auto g_yzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 268);

    auto g_yzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 269);

    #pragma omp simd aligned(g_yzzzzz_0_xxx_0, g_yzzzzz_0_xxy_0, g_yzzzzz_0_xxz_0, g_yzzzzz_0_xyy_0, g_yzzzzz_0_xyz_0, g_yzzzzz_0_xzz_0, g_yzzzzz_0_yyy_0, g_yzzzzz_0_yyz_0, g_yzzzzz_0_yzz_0, g_yzzzzz_0_zzz_0, g_zzzzz_0_xx_1, g_zzzzz_0_xxx_1, g_zzzzz_0_xxy_1, g_zzzzz_0_xxz_1, g_zzzzz_0_xy_1, g_zzzzz_0_xyy_1, g_zzzzz_0_xyz_1, g_zzzzz_0_xz_1, g_zzzzz_0_xzz_1, g_zzzzz_0_yy_1, g_zzzzz_0_yyy_1, g_zzzzz_0_yyz_1, g_zzzzz_0_yz_1, g_zzzzz_0_yzz_1, g_zzzzz_0_zz_1, g_zzzzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xxx_0[i] = g_zzzzz_0_xxx_1[i] * wa_y[i];

        g_yzzzzz_0_xxy_0[i] = g_zzzzz_0_xx_1[i] * fi_acd_0 + g_zzzzz_0_xxy_1[i] * wa_y[i];

        g_yzzzzz_0_xxz_0[i] = g_zzzzz_0_xxz_1[i] * wa_y[i];

        g_yzzzzz_0_xyy_0[i] = 2.0 * g_zzzzz_0_xy_1[i] * fi_acd_0 + g_zzzzz_0_xyy_1[i] * wa_y[i];

        g_yzzzzz_0_xyz_0[i] = g_zzzzz_0_xz_1[i] * fi_acd_0 + g_zzzzz_0_xyz_1[i] * wa_y[i];

        g_yzzzzz_0_xzz_0[i] = g_zzzzz_0_xzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyy_0[i] = 3.0 * g_zzzzz_0_yy_1[i] * fi_acd_0 + g_zzzzz_0_yyy_1[i] * wa_y[i];

        g_yzzzzz_0_yyz_0[i] = 2.0 * g_zzzzz_0_yz_1[i] * fi_acd_0 + g_zzzzz_0_yyz_1[i] * wa_y[i];

        g_yzzzzz_0_yzz_0[i] = g_zzzzz_0_zz_1[i] * fi_acd_0 + g_zzzzz_0_yzz_1[i] * wa_y[i];

        g_yzzzzz_0_zzz_0[i] = g_zzzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 270-280 components of targeted buffer : ISF

    auto g_zzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_isf + 270);

    auto g_zzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_isf + 271);

    auto g_zzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_isf + 272);

    auto g_zzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_isf + 273);

    auto g_zzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_isf + 274);

    auto g_zzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_isf + 275);

    auto g_zzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_isf + 276);

    auto g_zzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_isf + 277);

    auto g_zzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_isf + 278);

    auto g_zzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_isf + 279);

    #pragma omp simd aligned(g_zzzz_0_xxx_0, g_zzzz_0_xxx_1, g_zzzz_0_xxy_0, g_zzzz_0_xxy_1, g_zzzz_0_xxz_0, g_zzzz_0_xxz_1, g_zzzz_0_xyy_0, g_zzzz_0_xyy_1, g_zzzz_0_xyz_0, g_zzzz_0_xyz_1, g_zzzz_0_xzz_0, g_zzzz_0_xzz_1, g_zzzz_0_yyy_0, g_zzzz_0_yyy_1, g_zzzz_0_yyz_0, g_zzzz_0_yyz_1, g_zzzz_0_yzz_0, g_zzzz_0_yzz_1, g_zzzz_0_zzz_0, g_zzzz_0_zzz_1, g_zzzzz_0_xx_1, g_zzzzz_0_xxx_1, g_zzzzz_0_xxy_1, g_zzzzz_0_xxz_1, g_zzzzz_0_xy_1, g_zzzzz_0_xyy_1, g_zzzzz_0_xyz_1, g_zzzzz_0_xz_1, g_zzzzz_0_xzz_1, g_zzzzz_0_yy_1, g_zzzzz_0_yyy_1, g_zzzzz_0_yyz_1, g_zzzzz_0_yz_1, g_zzzzz_0_yzz_1, g_zzzzz_0_zz_1, g_zzzzz_0_zzz_1, g_zzzzzz_0_xxx_0, g_zzzzzz_0_xxy_0, g_zzzzzz_0_xxz_0, g_zzzzzz_0_xyy_0, g_zzzzzz_0_xyz_0, g_zzzzzz_0_xzz_0, g_zzzzzz_0_yyy_0, g_zzzzzz_0_yyz_0, g_zzzzzz_0_yzz_0, g_zzzzzz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xxx_0[i] = 5.0 * g_zzzz_0_xxx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxx_1[i] * fz_be_0 + g_zzzzz_0_xxx_1[i] * wa_z[i];

        g_zzzzzz_0_xxy_0[i] = 5.0 * g_zzzz_0_xxy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxy_1[i] * fz_be_0 + g_zzzzz_0_xxy_1[i] * wa_z[i];

        g_zzzzzz_0_xxz_0[i] = 5.0 * g_zzzz_0_xxz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxz_1[i] * fz_be_0 + g_zzzzz_0_xx_1[i] * fi_acd_0 + g_zzzzz_0_xxz_1[i] * wa_z[i];

        g_zzzzzz_0_xyy_0[i] = 5.0 * g_zzzz_0_xyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyy_1[i] * fz_be_0 + g_zzzzz_0_xyy_1[i] * wa_z[i];

        g_zzzzzz_0_xyz_0[i] = 5.0 * g_zzzz_0_xyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyz_1[i] * fz_be_0 + g_zzzzz_0_xy_1[i] * fi_acd_0 + g_zzzzz_0_xyz_1[i] * wa_z[i];

        g_zzzzzz_0_xzz_0[i] = 5.0 * g_zzzz_0_xzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xz_1[i] * fi_acd_0 + g_zzzzz_0_xzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyy_0[i] = 5.0 * g_zzzz_0_yyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyy_1[i] * fz_be_0 + g_zzzzz_0_yyy_1[i] * wa_z[i];

        g_zzzzzz_0_yyz_0[i] = 5.0 * g_zzzz_0_yyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyz_1[i] * fz_be_0 + g_zzzzz_0_yy_1[i] * fi_acd_0 + g_zzzzz_0_yyz_1[i] * wa_z[i];

        g_zzzzzz_0_yzz_0[i] = 5.0 * g_zzzz_0_yzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_yz_1[i] * fi_acd_0 + g_zzzzz_0_yzz_1[i] * wa_z[i];

        g_zzzzzz_0_zzz_0[i] = 5.0 * g_zzzz_0_zzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_zz_1[i] * fi_acd_0 + g_zzzzz_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

