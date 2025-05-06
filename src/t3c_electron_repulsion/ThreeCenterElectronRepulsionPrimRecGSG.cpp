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

#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsg,
                                 size_t idx_eri_0_dsg,
                                 size_t idx_eri_1_dsg,
                                 size_t idx_eri_1_fsf,
                                 size_t idx_eri_1_fsg,
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

    /// Set up components of auxilary buffer : DSG

    auto g_xx_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg);

    auto g_xx_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 1);

    auto g_xx_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 2);

    auto g_xx_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 3);

    auto g_xx_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 4);

    auto g_xx_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 5);

    auto g_xx_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 6);

    auto g_xx_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 7);

    auto g_xx_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 8);

    auto g_xx_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 9);

    auto g_xx_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 10);

    auto g_xx_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 11);

    auto g_xx_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 12);

    auto g_xx_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 13);

    auto g_xx_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 14);

    auto g_yy_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 45);

    auto g_yy_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 46);

    auto g_yy_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 47);

    auto g_yy_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 48);

    auto g_yy_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 49);

    auto g_yy_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 50);

    auto g_yy_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 51);

    auto g_yy_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 52);

    auto g_yy_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 53);

    auto g_yy_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 54);

    auto g_yy_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 55);

    auto g_yy_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 56);

    auto g_yy_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 57);

    auto g_yy_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 58);

    auto g_yy_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 59);

    auto g_zz_0_xxxx_0 = pbuffer.data(idx_eri_0_dsg + 75);

    auto g_zz_0_xxxy_0 = pbuffer.data(idx_eri_0_dsg + 76);

    auto g_zz_0_xxxz_0 = pbuffer.data(idx_eri_0_dsg + 77);

    auto g_zz_0_xxyy_0 = pbuffer.data(idx_eri_0_dsg + 78);

    auto g_zz_0_xxyz_0 = pbuffer.data(idx_eri_0_dsg + 79);

    auto g_zz_0_xxzz_0 = pbuffer.data(idx_eri_0_dsg + 80);

    auto g_zz_0_xyyy_0 = pbuffer.data(idx_eri_0_dsg + 81);

    auto g_zz_0_xyyz_0 = pbuffer.data(idx_eri_0_dsg + 82);

    auto g_zz_0_xyzz_0 = pbuffer.data(idx_eri_0_dsg + 83);

    auto g_zz_0_xzzz_0 = pbuffer.data(idx_eri_0_dsg + 84);

    auto g_zz_0_yyyy_0 = pbuffer.data(idx_eri_0_dsg + 85);

    auto g_zz_0_yyyz_0 = pbuffer.data(idx_eri_0_dsg + 86);

    auto g_zz_0_yyzz_0 = pbuffer.data(idx_eri_0_dsg + 87);

    auto g_zz_0_yzzz_0 = pbuffer.data(idx_eri_0_dsg + 88);

    auto g_zz_0_zzzz_0 = pbuffer.data(idx_eri_0_dsg + 89);

    /// Set up components of auxilary buffer : DSG

    auto g_xx_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg);

    auto g_xx_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 1);

    auto g_xx_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 2);

    auto g_xx_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 3);

    auto g_xx_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 4);

    auto g_xx_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 5);

    auto g_xx_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 6);

    auto g_xx_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 7);

    auto g_xx_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 8);

    auto g_xx_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 9);

    auto g_xx_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 10);

    auto g_xx_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 11);

    auto g_xx_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 12);

    auto g_xx_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 13);

    auto g_xx_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 14);

    auto g_yy_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg + 45);

    auto g_yy_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 46);

    auto g_yy_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 47);

    auto g_yy_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 48);

    auto g_yy_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 49);

    auto g_yy_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 50);

    auto g_yy_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 51);

    auto g_yy_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 52);

    auto g_yy_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 53);

    auto g_yy_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 54);

    auto g_yy_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 55);

    auto g_yy_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 56);

    auto g_yy_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 57);

    auto g_yy_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 58);

    auto g_yy_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 59);

    auto g_zz_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg + 75);

    auto g_zz_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 76);

    auto g_zz_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 77);

    auto g_zz_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 78);

    auto g_zz_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 79);

    auto g_zz_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 80);

    auto g_zz_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 81);

    auto g_zz_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 82);

    auto g_zz_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 83);

    auto g_zz_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 84);

    auto g_zz_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 85);

    auto g_zz_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 86);

    auto g_zz_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 87);

    auto g_zz_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 88);

    auto g_zz_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 89);

    /// Set up components of auxilary buffer : FSF

    auto g_xxx_0_xxx_1 = pbuffer.data(idx_eri_1_fsf);

    auto g_xxx_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 1);

    auto g_xxx_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 2);

    auto g_xxx_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 3);

    auto g_xxx_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 4);

    auto g_xxx_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 5);

    auto g_xxx_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 6);

    auto g_xxx_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 7);

    auto g_xxx_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 8);

    auto g_xxx_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 9);

    auto g_xxz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 22);

    auto g_xxz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 24);

    auto g_xxz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 25);

    auto g_xxz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 27);

    auto g_xxz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 28);

    auto g_xxz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 29);

    auto g_xyy_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 31);

    auto g_xyy_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 33);

    auto g_xyy_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 34);

    auto g_xyy_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 36);

    auto g_xyy_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 37);

    auto g_xyy_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 38);

    auto g_xzz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 52);

    auto g_xzz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 54);

    auto g_xzz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 55);

    auto g_xzz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 57);

    auto g_xzz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 58);

    auto g_xzz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 59);

    auto g_yyy_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 60);

    auto g_yyy_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 61);

    auto g_yyy_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 62);

    auto g_yyy_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 63);

    auto g_yyy_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 64);

    auto g_yyy_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 65);

    auto g_yyy_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 66);

    auto g_yyy_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 67);

    auto g_yyy_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 68);

    auto g_yyy_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 69);

    auto g_yyz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 72);

    auto g_yyz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 74);

    auto g_yyz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 75);

    auto g_yyz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 77);

    auto g_yyz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 78);

    auto g_yyz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 79);

    auto g_yzz_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 81);

    auto g_yzz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 82);

    auto g_yzz_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 83);

    auto g_yzz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 84);

    auto g_yzz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 85);

    auto g_yzz_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 86);

    auto g_yzz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 87);

    auto g_yzz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 88);

    auto g_yzz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 89);

    auto g_zzz_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 90);

    auto g_zzz_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 91);

    auto g_zzz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 92);

    auto g_zzz_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 93);

    auto g_zzz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 94);

    auto g_zzz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 95);

    auto g_zzz_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 96);

    auto g_zzz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 97);

    auto g_zzz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 98);

    auto g_zzz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 99);

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

    auto g_xxy_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 16);

    auto g_xxy_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 17);

    auto g_xxy_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 18);

    auto g_xxy_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 20);

    auto g_xxy_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 21);

    auto g_xxy_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 24);

    auto g_xxy_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 25);

    auto g_xxz_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 30);

    auto g_xxz_0_xxxy_1 = pbuffer.data(idx_eri_1_fsg + 31);

    auto g_xxz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 32);

    auto g_xxz_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 33);

    auto g_xxz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 34);

    auto g_xxz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 35);

    auto g_xxz_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 36);

    auto g_xxz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 37);

    auto g_xxz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 38);

    auto g_xxz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 39);

    auto g_xxz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 41);

    auto g_xxz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 42);

    auto g_xxz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 43);

    auto g_xxz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 44);

    auto g_xyy_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 45);

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

    auto g_xzz_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 75);

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

    auto g_yyz_0_xxxz_1 = pbuffer.data(idx_eri_1_fsg + 107);

    auto g_yyz_0_xxyy_1 = pbuffer.data(idx_eri_1_fsg + 108);

    auto g_yyz_0_xxyz_1 = pbuffer.data(idx_eri_1_fsg + 109);

    auto g_yyz_0_xxzz_1 = pbuffer.data(idx_eri_1_fsg + 110);

    auto g_yyz_0_xyyy_1 = pbuffer.data(idx_eri_1_fsg + 111);

    auto g_yyz_0_xyyz_1 = pbuffer.data(idx_eri_1_fsg + 112);

    auto g_yyz_0_xyzz_1 = pbuffer.data(idx_eri_1_fsg + 113);

    auto g_yyz_0_xzzz_1 = pbuffer.data(idx_eri_1_fsg + 114);

    auto g_yyz_0_yyyy_1 = pbuffer.data(idx_eri_1_fsg + 115);

    auto g_yyz_0_yyyz_1 = pbuffer.data(idx_eri_1_fsg + 116);

    auto g_yyz_0_yyzz_1 = pbuffer.data(idx_eri_1_fsg + 117);

    auto g_yyz_0_yzzz_1 = pbuffer.data(idx_eri_1_fsg + 118);

    auto g_yyz_0_zzzz_1 = pbuffer.data(idx_eri_1_fsg + 119);

    auto g_yzz_0_xxxx_1 = pbuffer.data(idx_eri_1_fsg + 120);

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

    /// Set up 0-15 components of targeted buffer : GSG

    auto g_xxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg);

    auto g_xxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 1);

    auto g_xxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 2);

    auto g_xxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 3);

    auto g_xxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 4);

    auto g_xxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 5);

    auto g_xxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 6);

    auto g_xxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 7);

    auto g_xxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 8);

    auto g_xxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 9);

    auto g_xxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 10);

    auto g_xxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 11);

    auto g_xxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 12);

    auto g_xxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 13);

    auto g_xxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 14);

    #pragma omp simd aligned(g_xx_0_xxxx_0, g_xx_0_xxxx_1, g_xx_0_xxxy_0, g_xx_0_xxxy_1, g_xx_0_xxxz_0, g_xx_0_xxxz_1, g_xx_0_xxyy_0, g_xx_0_xxyy_1, g_xx_0_xxyz_0, g_xx_0_xxyz_1, g_xx_0_xxzz_0, g_xx_0_xxzz_1, g_xx_0_xyyy_0, g_xx_0_xyyy_1, g_xx_0_xyyz_0, g_xx_0_xyyz_1, g_xx_0_xyzz_0, g_xx_0_xyzz_1, g_xx_0_xzzz_0, g_xx_0_xzzz_1, g_xx_0_yyyy_0, g_xx_0_yyyy_1, g_xx_0_yyyz_0, g_xx_0_yyyz_1, g_xx_0_yyzz_0, g_xx_0_yyzz_1, g_xx_0_yzzz_0, g_xx_0_yzzz_1, g_xx_0_zzzz_0, g_xx_0_zzzz_1, g_xxx_0_xxx_1, g_xxx_0_xxxx_1, g_xxx_0_xxxy_1, g_xxx_0_xxxz_1, g_xxx_0_xxy_1, g_xxx_0_xxyy_1, g_xxx_0_xxyz_1, g_xxx_0_xxz_1, g_xxx_0_xxzz_1, g_xxx_0_xyy_1, g_xxx_0_xyyy_1, g_xxx_0_xyyz_1, g_xxx_0_xyz_1, g_xxx_0_xyzz_1, g_xxx_0_xzz_1, g_xxx_0_xzzz_1, g_xxx_0_yyy_1, g_xxx_0_yyyy_1, g_xxx_0_yyyz_1, g_xxx_0_yyz_1, g_xxx_0_yyzz_1, g_xxx_0_yzz_1, g_xxx_0_yzzz_1, g_xxx_0_zzz_1, g_xxx_0_zzzz_1, g_xxxx_0_xxxx_0, g_xxxx_0_xxxy_0, g_xxxx_0_xxxz_0, g_xxxx_0_xxyy_0, g_xxxx_0_xxyz_0, g_xxxx_0_xxzz_0, g_xxxx_0_xyyy_0, g_xxxx_0_xyyz_0, g_xxxx_0_xyzz_0, g_xxxx_0_xzzz_0, g_xxxx_0_yyyy_0, g_xxxx_0_yyyz_0, g_xxxx_0_yyzz_0, g_xxxx_0_yzzz_0, g_xxxx_0_zzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxxx_0[i] = 3.0 * g_xx_0_xxxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxxx_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxx_1[i] * fi_acd_0 + g_xxx_0_xxxx_1[i] * wa_x[i];

        g_xxxx_0_xxxy_0[i] = 3.0 * g_xx_0_xxxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxy_1[i] * fi_acd_0 + g_xxx_0_xxxy_1[i] * wa_x[i];

        g_xxxx_0_xxxz_0[i] = 3.0 * g_xx_0_xxxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxz_1[i] * fi_acd_0 + g_xxx_0_xxxz_1[i] * wa_x[i];

        g_xxxx_0_xxyy_0[i] = 3.0 * g_xx_0_xxyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyy_1[i] * fi_acd_0 + g_xxx_0_xxyy_1[i] * wa_x[i];

        g_xxxx_0_xxyz_0[i] = 3.0 * g_xx_0_xxyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyz_1[i] * fi_acd_0 + g_xxx_0_xxyz_1[i] * wa_x[i];

        g_xxxx_0_xxzz_0[i] = 3.0 * g_xx_0_xxzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xzz_1[i] * fi_acd_0 + g_xxx_0_xxzz_1[i] * wa_x[i];

        g_xxxx_0_xyyy_0[i] = 3.0 * g_xx_0_xyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyyy_1[i] * fz_be_0 + g_xxx_0_yyy_1[i] * fi_acd_0 + g_xxx_0_xyyy_1[i] * wa_x[i];

        g_xxxx_0_xyyz_0[i] = 3.0 * g_xx_0_xyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyz_1[i] * fz_be_0 + g_xxx_0_yyz_1[i] * fi_acd_0 + g_xxx_0_xyyz_1[i] * wa_x[i];

        g_xxxx_0_xyzz_0[i] = 3.0 * g_xx_0_xyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyzz_1[i] * fz_be_0 + g_xxx_0_yzz_1[i] * fi_acd_0 + g_xxx_0_xyzz_1[i] * wa_x[i];

        g_xxxx_0_xzzz_0[i] = 3.0 * g_xx_0_xzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzzz_1[i] * fz_be_0 + g_xxx_0_zzz_1[i] * fi_acd_0 + g_xxx_0_xzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyy_0[i] = 3.0 * g_xx_0_yyyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyyy_1[i] * fz_be_0 + g_xxx_0_yyyy_1[i] * wa_x[i];

        g_xxxx_0_yyyz_0[i] = 3.0 * g_xx_0_yyyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyz_1[i] * fz_be_0 + g_xxx_0_yyyz_1[i] * wa_x[i];

        g_xxxx_0_yyzz_0[i] = 3.0 * g_xx_0_yyzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyzz_1[i] * fz_be_0 + g_xxx_0_yyzz_1[i] * wa_x[i];

        g_xxxx_0_yzzz_0[i] = 3.0 * g_xx_0_yzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzzz_1[i] * fz_be_0 + g_xxx_0_yzzz_1[i] * wa_x[i];

        g_xxxx_0_zzzz_0[i] = 3.0 * g_xx_0_zzzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzzz_1[i] * fz_be_0 + g_xxx_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : GSG

    auto g_xxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 15);

    auto g_xxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 16);

    auto g_xxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 17);

    auto g_xxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 18);

    auto g_xxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 19);

    auto g_xxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 20);

    auto g_xxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 21);

    auto g_xxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 22);

    auto g_xxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 23);

    auto g_xxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 24);

    auto g_xxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 25);

    auto g_xxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 26);

    auto g_xxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 27);

    auto g_xxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 28);

    auto g_xxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 29);

    #pragma omp simd aligned(g_xxx_0_xxx_1, g_xxx_0_xxxx_1, g_xxx_0_xxxy_1, g_xxx_0_xxxz_1, g_xxx_0_xxy_1, g_xxx_0_xxyy_1, g_xxx_0_xxyz_1, g_xxx_0_xxz_1, g_xxx_0_xxzz_1, g_xxx_0_xyy_1, g_xxx_0_xyyy_1, g_xxx_0_xyyz_1, g_xxx_0_xyz_1, g_xxx_0_xyzz_1, g_xxx_0_xzz_1, g_xxx_0_xzzz_1, g_xxx_0_yyy_1, g_xxx_0_yyyy_1, g_xxx_0_yyyz_1, g_xxx_0_yyz_1, g_xxx_0_yyzz_1, g_xxx_0_yzz_1, g_xxx_0_yzzz_1, g_xxx_0_zzz_1, g_xxx_0_zzzz_1, g_xxxy_0_xxxx_0, g_xxxy_0_xxxy_0, g_xxxy_0_xxxz_0, g_xxxy_0_xxyy_0, g_xxxy_0_xxyz_0, g_xxxy_0_xxzz_0, g_xxxy_0_xyyy_0, g_xxxy_0_xyyz_0, g_xxxy_0_xyzz_0, g_xxxy_0_xzzz_0, g_xxxy_0_yyyy_0, g_xxxy_0_yyyz_0, g_xxxy_0_yyzz_0, g_xxxy_0_yzzz_0, g_xxxy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxxx_0[i] = g_xxx_0_xxxx_1[i] * wa_y[i];

        g_xxxy_0_xxxy_0[i] = g_xxx_0_xxx_1[i] * fi_acd_0 + g_xxx_0_xxxy_1[i] * wa_y[i];

        g_xxxy_0_xxxz_0[i] = g_xxx_0_xxxz_1[i] * wa_y[i];

        g_xxxy_0_xxyy_0[i] = 2.0 * g_xxx_0_xxy_1[i] * fi_acd_0 + g_xxx_0_xxyy_1[i] * wa_y[i];

        g_xxxy_0_xxyz_0[i] = g_xxx_0_xxz_1[i] * fi_acd_0 + g_xxx_0_xxyz_1[i] * wa_y[i];

        g_xxxy_0_xxzz_0[i] = g_xxx_0_xxzz_1[i] * wa_y[i];

        g_xxxy_0_xyyy_0[i] = 3.0 * g_xxx_0_xyy_1[i] * fi_acd_0 + g_xxx_0_xyyy_1[i] * wa_y[i];

        g_xxxy_0_xyyz_0[i] = 2.0 * g_xxx_0_xyz_1[i] * fi_acd_0 + g_xxx_0_xyyz_1[i] * wa_y[i];

        g_xxxy_0_xyzz_0[i] = g_xxx_0_xzz_1[i] * fi_acd_0 + g_xxx_0_xyzz_1[i] * wa_y[i];

        g_xxxy_0_xzzz_0[i] = g_xxx_0_xzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyy_0[i] = 4.0 * g_xxx_0_yyy_1[i] * fi_acd_0 + g_xxx_0_yyyy_1[i] * wa_y[i];

        g_xxxy_0_yyyz_0[i] = 3.0 * g_xxx_0_yyz_1[i] * fi_acd_0 + g_xxx_0_yyyz_1[i] * wa_y[i];

        g_xxxy_0_yyzz_0[i] = 2.0 * g_xxx_0_yzz_1[i] * fi_acd_0 + g_xxx_0_yyzz_1[i] * wa_y[i];

        g_xxxy_0_yzzz_0[i] = g_xxx_0_zzz_1[i] * fi_acd_0 + g_xxx_0_yzzz_1[i] * wa_y[i];

        g_xxxy_0_zzzz_0[i] = g_xxx_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 30-45 components of targeted buffer : GSG

    auto g_xxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 30);

    auto g_xxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 31);

    auto g_xxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 32);

    auto g_xxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 33);

    auto g_xxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 34);

    auto g_xxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 35);

    auto g_xxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 36);

    auto g_xxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 37);

    auto g_xxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 38);

    auto g_xxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 39);

    auto g_xxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 40);

    auto g_xxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 41);

    auto g_xxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 42);

    auto g_xxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 43);

    auto g_xxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 44);

    #pragma omp simd aligned(g_xxx_0_xxx_1, g_xxx_0_xxxx_1, g_xxx_0_xxxy_1, g_xxx_0_xxxz_1, g_xxx_0_xxy_1, g_xxx_0_xxyy_1, g_xxx_0_xxyz_1, g_xxx_0_xxz_1, g_xxx_0_xxzz_1, g_xxx_0_xyy_1, g_xxx_0_xyyy_1, g_xxx_0_xyyz_1, g_xxx_0_xyz_1, g_xxx_0_xyzz_1, g_xxx_0_xzz_1, g_xxx_0_xzzz_1, g_xxx_0_yyy_1, g_xxx_0_yyyy_1, g_xxx_0_yyyz_1, g_xxx_0_yyz_1, g_xxx_0_yyzz_1, g_xxx_0_yzz_1, g_xxx_0_yzzz_1, g_xxx_0_zzz_1, g_xxx_0_zzzz_1, g_xxxz_0_xxxx_0, g_xxxz_0_xxxy_0, g_xxxz_0_xxxz_0, g_xxxz_0_xxyy_0, g_xxxz_0_xxyz_0, g_xxxz_0_xxzz_0, g_xxxz_0_xyyy_0, g_xxxz_0_xyyz_0, g_xxxz_0_xyzz_0, g_xxxz_0_xzzz_0, g_xxxz_0_yyyy_0, g_xxxz_0_yyyz_0, g_xxxz_0_yyzz_0, g_xxxz_0_yzzz_0, g_xxxz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxxx_0[i] = g_xxx_0_xxxx_1[i] * wa_z[i];

        g_xxxz_0_xxxy_0[i] = g_xxx_0_xxxy_1[i] * wa_z[i];

        g_xxxz_0_xxxz_0[i] = g_xxx_0_xxx_1[i] * fi_acd_0 + g_xxx_0_xxxz_1[i] * wa_z[i];

        g_xxxz_0_xxyy_0[i] = g_xxx_0_xxyy_1[i] * wa_z[i];

        g_xxxz_0_xxyz_0[i] = g_xxx_0_xxy_1[i] * fi_acd_0 + g_xxx_0_xxyz_1[i] * wa_z[i];

        g_xxxz_0_xxzz_0[i] = 2.0 * g_xxx_0_xxz_1[i] * fi_acd_0 + g_xxx_0_xxzz_1[i] * wa_z[i];

        g_xxxz_0_xyyy_0[i] = g_xxx_0_xyyy_1[i] * wa_z[i];

        g_xxxz_0_xyyz_0[i] = g_xxx_0_xyy_1[i] * fi_acd_0 + g_xxx_0_xyyz_1[i] * wa_z[i];

        g_xxxz_0_xyzz_0[i] = 2.0 * g_xxx_0_xyz_1[i] * fi_acd_0 + g_xxx_0_xyzz_1[i] * wa_z[i];

        g_xxxz_0_xzzz_0[i] = 3.0 * g_xxx_0_xzz_1[i] * fi_acd_0 + g_xxx_0_xzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyy_0[i] = g_xxx_0_yyyy_1[i] * wa_z[i];

        g_xxxz_0_yyyz_0[i] = g_xxx_0_yyy_1[i] * fi_acd_0 + g_xxx_0_yyyz_1[i] * wa_z[i];

        g_xxxz_0_yyzz_0[i] = 2.0 * g_xxx_0_yyz_1[i] * fi_acd_0 + g_xxx_0_yyzz_1[i] * wa_z[i];

        g_xxxz_0_yzzz_0[i] = 3.0 * g_xxx_0_yzz_1[i] * fi_acd_0 + g_xxx_0_yzzz_1[i] * wa_z[i];

        g_xxxz_0_zzzz_0[i] = 4.0 * g_xxx_0_zzz_1[i] * fi_acd_0 + g_xxx_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 45-60 components of targeted buffer : GSG

    auto g_xxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 45);

    auto g_xxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 46);

    auto g_xxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 47);

    auto g_xxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 48);

    auto g_xxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 49);

    auto g_xxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 50);

    auto g_xxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 51);

    auto g_xxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 52);

    auto g_xxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 53);

    auto g_xxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 54);

    auto g_xxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 55);

    auto g_xxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 56);

    auto g_xxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 57);

    auto g_xxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 58);

    auto g_xxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 59);

    #pragma omp simd aligned(g_xx_0_xxxx_0, g_xx_0_xxxx_1, g_xx_0_xxxz_0, g_xx_0_xxxz_1, g_xx_0_xxzz_0, g_xx_0_xxzz_1, g_xx_0_xzzz_0, g_xx_0_xzzz_1, g_xxy_0_xxxx_1, g_xxy_0_xxxz_1, g_xxy_0_xxzz_1, g_xxy_0_xzzz_1, g_xxyy_0_xxxx_0, g_xxyy_0_xxxy_0, g_xxyy_0_xxxz_0, g_xxyy_0_xxyy_0, g_xxyy_0_xxyz_0, g_xxyy_0_xxzz_0, g_xxyy_0_xyyy_0, g_xxyy_0_xyyz_0, g_xxyy_0_xyzz_0, g_xxyy_0_xzzz_0, g_xxyy_0_yyyy_0, g_xxyy_0_yyyz_0, g_xxyy_0_yyzz_0, g_xxyy_0_yzzz_0, g_xxyy_0_zzzz_0, g_xyy_0_xxxy_1, g_xyy_0_xxy_1, g_xyy_0_xxyy_1, g_xyy_0_xxyz_1, g_xyy_0_xyy_1, g_xyy_0_xyyy_1, g_xyy_0_xyyz_1, g_xyy_0_xyz_1, g_xyy_0_xyzz_1, g_xyy_0_yyy_1, g_xyy_0_yyyy_1, g_xyy_0_yyyz_1, g_xyy_0_yyz_1, g_xyy_0_yyzz_1, g_xyy_0_yzz_1, g_xyy_0_yzzz_1, g_xyy_0_zzzz_1, g_yy_0_xxxy_0, g_yy_0_xxxy_1, g_yy_0_xxyy_0, g_yy_0_xxyy_1, g_yy_0_xxyz_0, g_yy_0_xxyz_1, g_yy_0_xyyy_0, g_yy_0_xyyy_1, g_yy_0_xyyz_0, g_yy_0_xyyz_1, g_yy_0_xyzz_0, g_yy_0_xyzz_1, g_yy_0_yyyy_0, g_yy_0_yyyy_1, g_yy_0_yyyz_0, g_yy_0_yyyz_1, g_yy_0_yyzz_0, g_yy_0_yyzz_1, g_yy_0_yzzz_0, g_yy_0_yzzz_1, g_yy_0_zzzz_0, g_yy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxxx_0[i] = g_xx_0_xxxx_0[i] * fbe_0 - g_xx_0_xxxx_1[i] * fz_be_0 + g_xxy_0_xxxx_1[i] * wa_y[i];

        g_xxyy_0_xxxy_0[i] = g_yy_0_xxxy_0[i] * fbe_0 - g_yy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxy_1[i] * fi_acd_0 + g_xyy_0_xxxy_1[i] * wa_x[i];

        g_xxyy_0_xxxz_0[i] = g_xx_0_xxxz_0[i] * fbe_0 - g_xx_0_xxxz_1[i] * fz_be_0 + g_xxy_0_xxxz_1[i] * wa_y[i];

        g_xxyy_0_xxyy_0[i] = g_yy_0_xxyy_0[i] * fbe_0 - g_yy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyy_1[i] * fi_acd_0 + g_xyy_0_xxyy_1[i] * wa_x[i];

        g_xxyy_0_xxyz_0[i] = g_yy_0_xxyz_0[i] * fbe_0 - g_yy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyz_1[i] * fi_acd_0 + g_xyy_0_xxyz_1[i] * wa_x[i];

        g_xxyy_0_xxzz_0[i] = g_xx_0_xxzz_0[i] * fbe_0 - g_xx_0_xxzz_1[i] * fz_be_0 + g_xxy_0_xxzz_1[i] * wa_y[i];

        g_xxyy_0_xyyy_0[i] = g_yy_0_xyyy_0[i] * fbe_0 - g_yy_0_xyyy_1[i] * fz_be_0 + g_xyy_0_yyy_1[i] * fi_acd_0 + g_xyy_0_xyyy_1[i] * wa_x[i];

        g_xxyy_0_xyyz_0[i] = g_yy_0_xyyz_0[i] * fbe_0 - g_yy_0_xyyz_1[i] * fz_be_0 + g_xyy_0_yyz_1[i] * fi_acd_0 + g_xyy_0_xyyz_1[i] * wa_x[i];

        g_xxyy_0_xyzz_0[i] = g_yy_0_xyzz_0[i] * fbe_0 - g_yy_0_xyzz_1[i] * fz_be_0 + g_xyy_0_yzz_1[i] * fi_acd_0 + g_xyy_0_xyzz_1[i] * wa_x[i];

        g_xxyy_0_xzzz_0[i] = g_xx_0_xzzz_0[i] * fbe_0 - g_xx_0_xzzz_1[i] * fz_be_0 + g_xxy_0_xzzz_1[i] * wa_y[i];

        g_xxyy_0_yyyy_0[i] = g_yy_0_yyyy_0[i] * fbe_0 - g_yy_0_yyyy_1[i] * fz_be_0 + g_xyy_0_yyyy_1[i] * wa_x[i];

        g_xxyy_0_yyyz_0[i] = g_yy_0_yyyz_0[i] * fbe_0 - g_yy_0_yyyz_1[i] * fz_be_0 + g_xyy_0_yyyz_1[i] * wa_x[i];

        g_xxyy_0_yyzz_0[i] = g_yy_0_yyzz_0[i] * fbe_0 - g_yy_0_yyzz_1[i] * fz_be_0 + g_xyy_0_yyzz_1[i] * wa_x[i];

        g_xxyy_0_yzzz_0[i] = g_yy_0_yzzz_0[i] * fbe_0 - g_yy_0_yzzz_1[i] * fz_be_0 + g_xyy_0_yzzz_1[i] * wa_x[i];

        g_xxyy_0_zzzz_0[i] = g_yy_0_zzzz_0[i] * fbe_0 - g_yy_0_zzzz_1[i] * fz_be_0 + g_xyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 60-75 components of targeted buffer : GSG

    auto g_xxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 60);

    auto g_xxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 61);

    auto g_xxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 62);

    auto g_xxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 63);

    auto g_xxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 64);

    auto g_xxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 65);

    auto g_xxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 66);

    auto g_xxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 67);

    auto g_xxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 68);

    auto g_xxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 69);

    auto g_xxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 70);

    auto g_xxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 71);

    auto g_xxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 72);

    auto g_xxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 73);

    auto g_xxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 74);

    #pragma omp simd aligned(g_xxy_0_xxxy_1, g_xxy_0_xxyy_1, g_xxy_0_xyyy_1, g_xxy_0_yyyy_1, g_xxyz_0_xxxx_0, g_xxyz_0_xxxy_0, g_xxyz_0_xxxz_0, g_xxyz_0_xxyy_0, g_xxyz_0_xxyz_0, g_xxyz_0_xxzz_0, g_xxyz_0_xyyy_0, g_xxyz_0_xyyz_0, g_xxyz_0_xyzz_0, g_xxyz_0_xzzz_0, g_xxyz_0_yyyy_0, g_xxyz_0_yyyz_0, g_xxyz_0_yyzz_0, g_xxyz_0_yzzz_0, g_xxyz_0_zzzz_0, g_xxz_0_xxxx_1, g_xxz_0_xxxz_1, g_xxz_0_xxyz_1, g_xxz_0_xxz_1, g_xxz_0_xxzz_1, g_xxz_0_xyyz_1, g_xxz_0_xyz_1, g_xxz_0_xyzz_1, g_xxz_0_xzz_1, g_xxz_0_xzzz_1, g_xxz_0_yyyz_1, g_xxz_0_yyz_1, g_xxz_0_yyzz_1, g_xxz_0_yzz_1, g_xxz_0_yzzz_1, g_xxz_0_zzz_1, g_xxz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxxx_0[i] = g_xxz_0_xxxx_1[i] * wa_y[i];

        g_xxyz_0_xxxy_0[i] = g_xxy_0_xxxy_1[i] * wa_z[i];

        g_xxyz_0_xxxz_0[i] = g_xxz_0_xxxz_1[i] * wa_y[i];

        g_xxyz_0_xxyy_0[i] = g_xxy_0_xxyy_1[i] * wa_z[i];

        g_xxyz_0_xxyz_0[i] = g_xxz_0_xxz_1[i] * fi_acd_0 + g_xxz_0_xxyz_1[i] * wa_y[i];

        g_xxyz_0_xxzz_0[i] = g_xxz_0_xxzz_1[i] * wa_y[i];

        g_xxyz_0_xyyy_0[i] = g_xxy_0_xyyy_1[i] * wa_z[i];

        g_xxyz_0_xyyz_0[i] = 2.0 * g_xxz_0_xyz_1[i] * fi_acd_0 + g_xxz_0_xyyz_1[i] * wa_y[i];

        g_xxyz_0_xyzz_0[i] = g_xxz_0_xzz_1[i] * fi_acd_0 + g_xxz_0_xyzz_1[i] * wa_y[i];

        g_xxyz_0_xzzz_0[i] = g_xxz_0_xzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyy_0[i] = g_xxy_0_yyyy_1[i] * wa_z[i];

        g_xxyz_0_yyyz_0[i] = 3.0 * g_xxz_0_yyz_1[i] * fi_acd_0 + g_xxz_0_yyyz_1[i] * wa_y[i];

        g_xxyz_0_yyzz_0[i] = 2.0 * g_xxz_0_yzz_1[i] * fi_acd_0 + g_xxz_0_yyzz_1[i] * wa_y[i];

        g_xxyz_0_yzzz_0[i] = g_xxz_0_zzz_1[i] * fi_acd_0 + g_xxz_0_yzzz_1[i] * wa_y[i];

        g_xxyz_0_zzzz_0[i] = g_xxz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 75-90 components of targeted buffer : GSG

    auto g_xxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 75);

    auto g_xxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 76);

    auto g_xxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 77);

    auto g_xxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 78);

    auto g_xxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 79);

    auto g_xxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 80);

    auto g_xxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 81);

    auto g_xxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 82);

    auto g_xxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 83);

    auto g_xxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 84);

    auto g_xxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 85);

    auto g_xxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 86);

    auto g_xxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 87);

    auto g_xxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 88);

    auto g_xxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 89);

    #pragma omp simd aligned(g_xx_0_xxxx_0, g_xx_0_xxxx_1, g_xx_0_xxxy_0, g_xx_0_xxxy_1, g_xx_0_xxyy_0, g_xx_0_xxyy_1, g_xx_0_xyyy_0, g_xx_0_xyyy_1, g_xxz_0_xxxx_1, g_xxz_0_xxxy_1, g_xxz_0_xxyy_1, g_xxz_0_xyyy_1, g_xxzz_0_xxxx_0, g_xxzz_0_xxxy_0, g_xxzz_0_xxxz_0, g_xxzz_0_xxyy_0, g_xxzz_0_xxyz_0, g_xxzz_0_xxzz_0, g_xxzz_0_xyyy_0, g_xxzz_0_xyyz_0, g_xxzz_0_xyzz_0, g_xxzz_0_xzzz_0, g_xxzz_0_yyyy_0, g_xxzz_0_yyyz_0, g_xxzz_0_yyzz_0, g_xxzz_0_yzzz_0, g_xxzz_0_zzzz_0, g_xzz_0_xxxz_1, g_xzz_0_xxyz_1, g_xzz_0_xxz_1, g_xzz_0_xxzz_1, g_xzz_0_xyyz_1, g_xzz_0_xyz_1, g_xzz_0_xyzz_1, g_xzz_0_xzz_1, g_xzz_0_xzzz_1, g_xzz_0_yyyy_1, g_xzz_0_yyyz_1, g_xzz_0_yyz_1, g_xzz_0_yyzz_1, g_xzz_0_yzz_1, g_xzz_0_yzzz_1, g_xzz_0_zzz_1, g_xzz_0_zzzz_1, g_zz_0_xxxz_0, g_zz_0_xxxz_1, g_zz_0_xxyz_0, g_zz_0_xxyz_1, g_zz_0_xxzz_0, g_zz_0_xxzz_1, g_zz_0_xyyz_0, g_zz_0_xyyz_1, g_zz_0_xyzz_0, g_zz_0_xyzz_1, g_zz_0_xzzz_0, g_zz_0_xzzz_1, g_zz_0_yyyy_0, g_zz_0_yyyy_1, g_zz_0_yyyz_0, g_zz_0_yyyz_1, g_zz_0_yyzz_0, g_zz_0_yyzz_1, g_zz_0_yzzz_0, g_zz_0_yzzz_1, g_zz_0_zzzz_0, g_zz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxxx_0[i] = g_xx_0_xxxx_0[i] * fbe_0 - g_xx_0_xxxx_1[i] * fz_be_0 + g_xxz_0_xxxx_1[i] * wa_z[i];

        g_xxzz_0_xxxy_0[i] = g_xx_0_xxxy_0[i] * fbe_0 - g_xx_0_xxxy_1[i] * fz_be_0 + g_xxz_0_xxxy_1[i] * wa_z[i];

        g_xxzz_0_xxxz_0[i] = g_zz_0_xxxz_0[i] * fbe_0 - g_zz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxz_1[i] * fi_acd_0 + g_xzz_0_xxxz_1[i] * wa_x[i];

        g_xxzz_0_xxyy_0[i] = g_xx_0_xxyy_0[i] * fbe_0 - g_xx_0_xxyy_1[i] * fz_be_0 + g_xxz_0_xxyy_1[i] * wa_z[i];

        g_xxzz_0_xxyz_0[i] = g_zz_0_xxyz_0[i] * fbe_0 - g_zz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyz_1[i] * fi_acd_0 + g_xzz_0_xxyz_1[i] * wa_x[i];

        g_xxzz_0_xxzz_0[i] = g_zz_0_xxzz_0[i] * fbe_0 - g_zz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xzz_1[i] * fi_acd_0 + g_xzz_0_xxzz_1[i] * wa_x[i];

        g_xxzz_0_xyyy_0[i] = g_xx_0_xyyy_0[i] * fbe_0 - g_xx_0_xyyy_1[i] * fz_be_0 + g_xxz_0_xyyy_1[i] * wa_z[i];

        g_xxzz_0_xyyz_0[i] = g_zz_0_xyyz_0[i] * fbe_0 - g_zz_0_xyyz_1[i] * fz_be_0 + g_xzz_0_yyz_1[i] * fi_acd_0 + g_xzz_0_xyyz_1[i] * wa_x[i];

        g_xxzz_0_xyzz_0[i] = g_zz_0_xyzz_0[i] * fbe_0 - g_zz_0_xyzz_1[i] * fz_be_0 + g_xzz_0_yzz_1[i] * fi_acd_0 + g_xzz_0_xyzz_1[i] * wa_x[i];

        g_xxzz_0_xzzz_0[i] = g_zz_0_xzzz_0[i] * fbe_0 - g_zz_0_xzzz_1[i] * fz_be_0 + g_xzz_0_zzz_1[i] * fi_acd_0 + g_xzz_0_xzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyy_0[i] = g_zz_0_yyyy_0[i] * fbe_0 - g_zz_0_yyyy_1[i] * fz_be_0 + g_xzz_0_yyyy_1[i] * wa_x[i];

        g_xxzz_0_yyyz_0[i] = g_zz_0_yyyz_0[i] * fbe_0 - g_zz_0_yyyz_1[i] * fz_be_0 + g_xzz_0_yyyz_1[i] * wa_x[i];

        g_xxzz_0_yyzz_0[i] = g_zz_0_yyzz_0[i] * fbe_0 - g_zz_0_yyzz_1[i] * fz_be_0 + g_xzz_0_yyzz_1[i] * wa_x[i];

        g_xxzz_0_yzzz_0[i] = g_zz_0_yzzz_0[i] * fbe_0 - g_zz_0_yzzz_1[i] * fz_be_0 + g_xzz_0_yzzz_1[i] * wa_x[i];

        g_xxzz_0_zzzz_0[i] = g_zz_0_zzzz_0[i] * fbe_0 - g_zz_0_zzzz_1[i] * fz_be_0 + g_xzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 90-105 components of targeted buffer : GSG

    auto g_xyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 90);

    auto g_xyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 91);

    auto g_xyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 92);

    auto g_xyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 93);

    auto g_xyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 94);

    auto g_xyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 95);

    auto g_xyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 96);

    auto g_xyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 97);

    auto g_xyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 98);

    auto g_xyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 99);

    auto g_xyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 100);

    auto g_xyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 101);

    auto g_xyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 102);

    auto g_xyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 103);

    auto g_xyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 104);

    #pragma omp simd aligned(g_xyyy_0_xxxx_0, g_xyyy_0_xxxy_0, g_xyyy_0_xxxz_0, g_xyyy_0_xxyy_0, g_xyyy_0_xxyz_0, g_xyyy_0_xxzz_0, g_xyyy_0_xyyy_0, g_xyyy_0_xyyz_0, g_xyyy_0_xyzz_0, g_xyyy_0_xzzz_0, g_xyyy_0_yyyy_0, g_xyyy_0_yyyz_0, g_xyyy_0_yyzz_0, g_xyyy_0_yzzz_0, g_xyyy_0_zzzz_0, g_yyy_0_xxx_1, g_yyy_0_xxxx_1, g_yyy_0_xxxy_1, g_yyy_0_xxxz_1, g_yyy_0_xxy_1, g_yyy_0_xxyy_1, g_yyy_0_xxyz_1, g_yyy_0_xxz_1, g_yyy_0_xxzz_1, g_yyy_0_xyy_1, g_yyy_0_xyyy_1, g_yyy_0_xyyz_1, g_yyy_0_xyz_1, g_yyy_0_xyzz_1, g_yyy_0_xzz_1, g_yyy_0_xzzz_1, g_yyy_0_yyy_1, g_yyy_0_yyyy_1, g_yyy_0_yyyz_1, g_yyy_0_yyz_1, g_yyy_0_yyzz_1, g_yyy_0_yzz_1, g_yyy_0_yzzz_1, g_yyy_0_zzz_1, g_yyy_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxxx_0[i] = 4.0 * g_yyy_0_xxx_1[i] * fi_acd_0 + g_yyy_0_xxxx_1[i] * wa_x[i];

        g_xyyy_0_xxxy_0[i] = 3.0 * g_yyy_0_xxy_1[i] * fi_acd_0 + g_yyy_0_xxxy_1[i] * wa_x[i];

        g_xyyy_0_xxxz_0[i] = 3.0 * g_yyy_0_xxz_1[i] * fi_acd_0 + g_yyy_0_xxxz_1[i] * wa_x[i];

        g_xyyy_0_xxyy_0[i] = 2.0 * g_yyy_0_xyy_1[i] * fi_acd_0 + g_yyy_0_xxyy_1[i] * wa_x[i];

        g_xyyy_0_xxyz_0[i] = 2.0 * g_yyy_0_xyz_1[i] * fi_acd_0 + g_yyy_0_xxyz_1[i] * wa_x[i];

        g_xyyy_0_xxzz_0[i] = 2.0 * g_yyy_0_xzz_1[i] * fi_acd_0 + g_yyy_0_xxzz_1[i] * wa_x[i];

        g_xyyy_0_xyyy_0[i] = g_yyy_0_yyy_1[i] * fi_acd_0 + g_yyy_0_xyyy_1[i] * wa_x[i];

        g_xyyy_0_xyyz_0[i] = g_yyy_0_yyz_1[i] * fi_acd_0 + g_yyy_0_xyyz_1[i] * wa_x[i];

        g_xyyy_0_xyzz_0[i] = g_yyy_0_yzz_1[i] * fi_acd_0 + g_yyy_0_xyzz_1[i] * wa_x[i];

        g_xyyy_0_xzzz_0[i] = g_yyy_0_zzz_1[i] * fi_acd_0 + g_yyy_0_xzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyy_0[i] = g_yyy_0_yyyy_1[i] * wa_x[i];

        g_xyyy_0_yyyz_0[i] = g_yyy_0_yyyz_1[i] * wa_x[i];

        g_xyyy_0_yyzz_0[i] = g_yyy_0_yyzz_1[i] * wa_x[i];

        g_xyyy_0_yzzz_0[i] = g_yyy_0_yzzz_1[i] * wa_x[i];

        g_xyyy_0_zzzz_0[i] = g_yyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 105-120 components of targeted buffer : GSG

    auto g_xyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 105);

    auto g_xyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 106);

    auto g_xyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 107);

    auto g_xyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 108);

    auto g_xyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 109);

    auto g_xyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 110);

    auto g_xyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 111);

    auto g_xyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 112);

    auto g_xyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 113);

    auto g_xyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 114);

    auto g_xyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 115);

    auto g_xyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 116);

    auto g_xyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 117);

    auto g_xyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 118);

    auto g_xyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 119);

    #pragma omp simd aligned(g_xyy_0_xxxx_1, g_xyy_0_xxxy_1, g_xyy_0_xxyy_1, g_xyy_0_xyyy_1, g_xyyz_0_xxxx_0, g_xyyz_0_xxxy_0, g_xyyz_0_xxxz_0, g_xyyz_0_xxyy_0, g_xyyz_0_xxyz_0, g_xyyz_0_xxzz_0, g_xyyz_0_xyyy_0, g_xyyz_0_xyyz_0, g_xyyz_0_xyzz_0, g_xyyz_0_xzzz_0, g_xyyz_0_yyyy_0, g_xyyz_0_yyyz_0, g_xyyz_0_yyzz_0, g_xyyz_0_yzzz_0, g_xyyz_0_zzzz_0, g_yyz_0_xxxz_1, g_yyz_0_xxyz_1, g_yyz_0_xxz_1, g_yyz_0_xxzz_1, g_yyz_0_xyyz_1, g_yyz_0_xyz_1, g_yyz_0_xyzz_1, g_yyz_0_xzz_1, g_yyz_0_xzzz_1, g_yyz_0_yyyy_1, g_yyz_0_yyyz_1, g_yyz_0_yyz_1, g_yyz_0_yyzz_1, g_yyz_0_yzz_1, g_yyz_0_yzzz_1, g_yyz_0_zzz_1, g_yyz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxxx_0[i] = g_xyy_0_xxxx_1[i] * wa_z[i];

        g_xyyz_0_xxxy_0[i] = g_xyy_0_xxxy_1[i] * wa_z[i];

        g_xyyz_0_xxxz_0[i] = 3.0 * g_yyz_0_xxz_1[i] * fi_acd_0 + g_yyz_0_xxxz_1[i] * wa_x[i];

        g_xyyz_0_xxyy_0[i] = g_xyy_0_xxyy_1[i] * wa_z[i];

        g_xyyz_0_xxyz_0[i] = 2.0 * g_yyz_0_xyz_1[i] * fi_acd_0 + g_yyz_0_xxyz_1[i] * wa_x[i];

        g_xyyz_0_xxzz_0[i] = 2.0 * g_yyz_0_xzz_1[i] * fi_acd_0 + g_yyz_0_xxzz_1[i] * wa_x[i];

        g_xyyz_0_xyyy_0[i] = g_xyy_0_xyyy_1[i] * wa_z[i];

        g_xyyz_0_xyyz_0[i] = g_yyz_0_yyz_1[i] * fi_acd_0 + g_yyz_0_xyyz_1[i] * wa_x[i];

        g_xyyz_0_xyzz_0[i] = g_yyz_0_yzz_1[i] * fi_acd_0 + g_yyz_0_xyzz_1[i] * wa_x[i];

        g_xyyz_0_xzzz_0[i] = g_yyz_0_zzz_1[i] * fi_acd_0 + g_yyz_0_xzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyy_0[i] = g_yyz_0_yyyy_1[i] * wa_x[i];

        g_xyyz_0_yyyz_0[i] = g_yyz_0_yyyz_1[i] * wa_x[i];

        g_xyyz_0_yyzz_0[i] = g_yyz_0_yyzz_1[i] * wa_x[i];

        g_xyyz_0_yzzz_0[i] = g_yyz_0_yzzz_1[i] * wa_x[i];

        g_xyyz_0_zzzz_0[i] = g_yyz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 120-135 components of targeted buffer : GSG

    auto g_xyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 120);

    auto g_xyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 121);

    auto g_xyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 122);

    auto g_xyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 123);

    auto g_xyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 124);

    auto g_xyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 125);

    auto g_xyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 126);

    auto g_xyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 127);

    auto g_xyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 128);

    auto g_xyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 129);

    auto g_xyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 130);

    auto g_xyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 131);

    auto g_xyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 132);

    auto g_xyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 133);

    auto g_xyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 134);

    #pragma omp simd aligned(g_xyzz_0_xxxx_0, g_xyzz_0_xxxy_0, g_xyzz_0_xxxz_0, g_xyzz_0_xxyy_0, g_xyzz_0_xxyz_0, g_xyzz_0_xxzz_0, g_xyzz_0_xyyy_0, g_xyzz_0_xyyz_0, g_xyzz_0_xyzz_0, g_xyzz_0_xzzz_0, g_xyzz_0_yyyy_0, g_xyzz_0_yyyz_0, g_xyzz_0_yyzz_0, g_xyzz_0_yzzz_0, g_xyzz_0_zzzz_0, g_xzz_0_xxxx_1, g_xzz_0_xxxz_1, g_xzz_0_xxzz_1, g_xzz_0_xzzz_1, g_yzz_0_xxxy_1, g_yzz_0_xxy_1, g_yzz_0_xxyy_1, g_yzz_0_xxyz_1, g_yzz_0_xyy_1, g_yzz_0_xyyy_1, g_yzz_0_xyyz_1, g_yzz_0_xyz_1, g_yzz_0_xyzz_1, g_yzz_0_yyy_1, g_yzz_0_yyyy_1, g_yzz_0_yyyz_1, g_yzz_0_yyz_1, g_yzz_0_yyzz_1, g_yzz_0_yzz_1, g_yzz_0_yzzz_1, g_yzz_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxxx_0[i] = g_xzz_0_xxxx_1[i] * wa_y[i];

        g_xyzz_0_xxxy_0[i] = 3.0 * g_yzz_0_xxy_1[i] * fi_acd_0 + g_yzz_0_xxxy_1[i] * wa_x[i];

        g_xyzz_0_xxxz_0[i] = g_xzz_0_xxxz_1[i] * wa_y[i];

        g_xyzz_0_xxyy_0[i] = 2.0 * g_yzz_0_xyy_1[i] * fi_acd_0 + g_yzz_0_xxyy_1[i] * wa_x[i];

        g_xyzz_0_xxyz_0[i] = 2.0 * g_yzz_0_xyz_1[i] * fi_acd_0 + g_yzz_0_xxyz_1[i] * wa_x[i];

        g_xyzz_0_xxzz_0[i] = g_xzz_0_xxzz_1[i] * wa_y[i];

        g_xyzz_0_xyyy_0[i] = g_yzz_0_yyy_1[i] * fi_acd_0 + g_yzz_0_xyyy_1[i] * wa_x[i];

        g_xyzz_0_xyyz_0[i] = g_yzz_0_yyz_1[i] * fi_acd_0 + g_yzz_0_xyyz_1[i] * wa_x[i];

        g_xyzz_0_xyzz_0[i] = g_yzz_0_yzz_1[i] * fi_acd_0 + g_yzz_0_xyzz_1[i] * wa_x[i];

        g_xyzz_0_xzzz_0[i] = g_xzz_0_xzzz_1[i] * wa_y[i];

        g_xyzz_0_yyyy_0[i] = g_yzz_0_yyyy_1[i] * wa_x[i];

        g_xyzz_0_yyyz_0[i] = g_yzz_0_yyyz_1[i] * wa_x[i];

        g_xyzz_0_yyzz_0[i] = g_yzz_0_yyzz_1[i] * wa_x[i];

        g_xyzz_0_yzzz_0[i] = g_yzz_0_yzzz_1[i] * wa_x[i];

        g_xyzz_0_zzzz_0[i] = g_yzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 135-150 components of targeted buffer : GSG

    auto g_xzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 135);

    auto g_xzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 136);

    auto g_xzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 137);

    auto g_xzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 138);

    auto g_xzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 139);

    auto g_xzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 140);

    auto g_xzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 141);

    auto g_xzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 142);

    auto g_xzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 143);

    auto g_xzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 144);

    auto g_xzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 145);

    auto g_xzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 146);

    auto g_xzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 147);

    auto g_xzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 148);

    auto g_xzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 149);

    #pragma omp simd aligned(g_xzzz_0_xxxx_0, g_xzzz_0_xxxy_0, g_xzzz_0_xxxz_0, g_xzzz_0_xxyy_0, g_xzzz_0_xxyz_0, g_xzzz_0_xxzz_0, g_xzzz_0_xyyy_0, g_xzzz_0_xyyz_0, g_xzzz_0_xyzz_0, g_xzzz_0_xzzz_0, g_xzzz_0_yyyy_0, g_xzzz_0_yyyz_0, g_xzzz_0_yyzz_0, g_xzzz_0_yzzz_0, g_xzzz_0_zzzz_0, g_zzz_0_xxx_1, g_zzz_0_xxxx_1, g_zzz_0_xxxy_1, g_zzz_0_xxxz_1, g_zzz_0_xxy_1, g_zzz_0_xxyy_1, g_zzz_0_xxyz_1, g_zzz_0_xxz_1, g_zzz_0_xxzz_1, g_zzz_0_xyy_1, g_zzz_0_xyyy_1, g_zzz_0_xyyz_1, g_zzz_0_xyz_1, g_zzz_0_xyzz_1, g_zzz_0_xzz_1, g_zzz_0_xzzz_1, g_zzz_0_yyy_1, g_zzz_0_yyyy_1, g_zzz_0_yyyz_1, g_zzz_0_yyz_1, g_zzz_0_yyzz_1, g_zzz_0_yzz_1, g_zzz_0_yzzz_1, g_zzz_0_zzz_1, g_zzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxxx_0[i] = 4.0 * g_zzz_0_xxx_1[i] * fi_acd_0 + g_zzz_0_xxxx_1[i] * wa_x[i];

        g_xzzz_0_xxxy_0[i] = 3.0 * g_zzz_0_xxy_1[i] * fi_acd_0 + g_zzz_0_xxxy_1[i] * wa_x[i];

        g_xzzz_0_xxxz_0[i] = 3.0 * g_zzz_0_xxz_1[i] * fi_acd_0 + g_zzz_0_xxxz_1[i] * wa_x[i];

        g_xzzz_0_xxyy_0[i] = 2.0 * g_zzz_0_xyy_1[i] * fi_acd_0 + g_zzz_0_xxyy_1[i] * wa_x[i];

        g_xzzz_0_xxyz_0[i] = 2.0 * g_zzz_0_xyz_1[i] * fi_acd_0 + g_zzz_0_xxyz_1[i] * wa_x[i];

        g_xzzz_0_xxzz_0[i] = 2.0 * g_zzz_0_xzz_1[i] * fi_acd_0 + g_zzz_0_xxzz_1[i] * wa_x[i];

        g_xzzz_0_xyyy_0[i] = g_zzz_0_yyy_1[i] * fi_acd_0 + g_zzz_0_xyyy_1[i] * wa_x[i];

        g_xzzz_0_xyyz_0[i] = g_zzz_0_yyz_1[i] * fi_acd_0 + g_zzz_0_xyyz_1[i] * wa_x[i];

        g_xzzz_0_xyzz_0[i] = g_zzz_0_yzz_1[i] * fi_acd_0 + g_zzz_0_xyzz_1[i] * wa_x[i];

        g_xzzz_0_xzzz_0[i] = g_zzz_0_zzz_1[i] * fi_acd_0 + g_zzz_0_xzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyy_0[i] = g_zzz_0_yyyy_1[i] * wa_x[i];

        g_xzzz_0_yyyz_0[i] = g_zzz_0_yyyz_1[i] * wa_x[i];

        g_xzzz_0_yyzz_0[i] = g_zzz_0_yyzz_1[i] * wa_x[i];

        g_xzzz_0_yzzz_0[i] = g_zzz_0_yzzz_1[i] * wa_x[i];

        g_xzzz_0_zzzz_0[i] = g_zzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 150-165 components of targeted buffer : GSG

    auto g_yyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 150);

    auto g_yyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 151);

    auto g_yyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 152);

    auto g_yyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 153);

    auto g_yyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 154);

    auto g_yyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 155);

    auto g_yyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 156);

    auto g_yyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 157);

    auto g_yyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 158);

    auto g_yyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 159);

    auto g_yyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 160);

    auto g_yyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 161);

    auto g_yyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 162);

    auto g_yyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 163);

    auto g_yyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 164);

    #pragma omp simd aligned(g_yy_0_xxxx_0, g_yy_0_xxxx_1, g_yy_0_xxxy_0, g_yy_0_xxxy_1, g_yy_0_xxxz_0, g_yy_0_xxxz_1, g_yy_0_xxyy_0, g_yy_0_xxyy_1, g_yy_0_xxyz_0, g_yy_0_xxyz_1, g_yy_0_xxzz_0, g_yy_0_xxzz_1, g_yy_0_xyyy_0, g_yy_0_xyyy_1, g_yy_0_xyyz_0, g_yy_0_xyyz_1, g_yy_0_xyzz_0, g_yy_0_xyzz_1, g_yy_0_xzzz_0, g_yy_0_xzzz_1, g_yy_0_yyyy_0, g_yy_0_yyyy_1, g_yy_0_yyyz_0, g_yy_0_yyyz_1, g_yy_0_yyzz_0, g_yy_0_yyzz_1, g_yy_0_yzzz_0, g_yy_0_yzzz_1, g_yy_0_zzzz_0, g_yy_0_zzzz_1, g_yyy_0_xxx_1, g_yyy_0_xxxx_1, g_yyy_0_xxxy_1, g_yyy_0_xxxz_1, g_yyy_0_xxy_1, g_yyy_0_xxyy_1, g_yyy_0_xxyz_1, g_yyy_0_xxz_1, g_yyy_0_xxzz_1, g_yyy_0_xyy_1, g_yyy_0_xyyy_1, g_yyy_0_xyyz_1, g_yyy_0_xyz_1, g_yyy_0_xyzz_1, g_yyy_0_xzz_1, g_yyy_0_xzzz_1, g_yyy_0_yyy_1, g_yyy_0_yyyy_1, g_yyy_0_yyyz_1, g_yyy_0_yyz_1, g_yyy_0_yyzz_1, g_yyy_0_yzz_1, g_yyy_0_yzzz_1, g_yyy_0_zzz_1, g_yyy_0_zzzz_1, g_yyyy_0_xxxx_0, g_yyyy_0_xxxy_0, g_yyyy_0_xxxz_0, g_yyyy_0_xxyy_0, g_yyyy_0_xxyz_0, g_yyyy_0_xxzz_0, g_yyyy_0_xyyy_0, g_yyyy_0_xyyz_0, g_yyyy_0_xyzz_0, g_yyyy_0_xzzz_0, g_yyyy_0_yyyy_0, g_yyyy_0_yyyz_0, g_yyyy_0_yyzz_0, g_yyyy_0_yzzz_0, g_yyyy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxxx_0[i] = 3.0 * g_yy_0_xxxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxxx_1[i] * fz_be_0 + g_yyy_0_xxxx_1[i] * wa_y[i];

        g_yyyy_0_xxxy_0[i] = 3.0 * g_yy_0_xxxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxy_1[i] * fz_be_0 + g_yyy_0_xxx_1[i] * fi_acd_0 + g_yyy_0_xxxy_1[i] * wa_y[i];

        g_yyyy_0_xxxz_0[i] = 3.0 * g_yy_0_xxxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxz_1[i] * fz_be_0 + g_yyy_0_xxxz_1[i] * wa_y[i];

        g_yyyy_0_xxyy_0[i] = 3.0 * g_yy_0_xxyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxy_1[i] * fi_acd_0 + g_yyy_0_xxyy_1[i] * wa_y[i];

        g_yyyy_0_xxyz_0[i] = 3.0 * g_yy_0_xxyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyz_1[i] * fz_be_0 + g_yyy_0_xxz_1[i] * fi_acd_0 + g_yyy_0_xxyz_1[i] * wa_y[i];

        g_yyyy_0_xxzz_0[i] = 3.0 * g_yy_0_xxzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxzz_1[i] * fz_be_0 + g_yyy_0_xxzz_1[i] * wa_y[i];

        g_yyyy_0_xyyy_0[i] = 3.0 * g_yy_0_xyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_xyy_1[i] * fi_acd_0 + g_yyy_0_xyyy_1[i] * wa_y[i];

        g_yyyy_0_xyyz_0[i] = 3.0 * g_yy_0_xyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xyz_1[i] * fi_acd_0 + g_yyy_0_xyyz_1[i] * wa_y[i];

        g_yyyy_0_xyzz_0[i] = 3.0 * g_yy_0_xyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyzz_1[i] * fz_be_0 + g_yyy_0_xzz_1[i] * fi_acd_0 + g_yyy_0_xyzz_1[i] * wa_y[i];

        g_yyyy_0_xzzz_0[i] = 3.0 * g_yy_0_xzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzzz_1[i] * fz_be_0 + g_yyy_0_xzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyy_0[i] = 3.0 * g_yy_0_yyyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyyy_1[i] * fz_be_0 + 4.0 * g_yyy_0_yyy_1[i] * fi_acd_0 + g_yyy_0_yyyy_1[i] * wa_y[i];

        g_yyyy_0_yyyz_0[i] = 3.0 * g_yy_0_yyyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyy_0_yyz_1[i] * fi_acd_0 + g_yyy_0_yyyz_1[i] * wa_y[i];

        g_yyyy_0_yyzz_0[i] = 3.0 * g_yy_0_yyzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yzz_1[i] * fi_acd_0 + g_yyy_0_yyzz_1[i] * wa_y[i];

        g_yyyy_0_yzzz_0[i] = 3.0 * g_yy_0_yzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzzz_1[i] * fz_be_0 + g_yyy_0_zzz_1[i] * fi_acd_0 + g_yyy_0_yzzz_1[i] * wa_y[i];

        g_yyyy_0_zzzz_0[i] = 3.0 * g_yy_0_zzzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzzz_1[i] * fz_be_0 + g_yyy_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 165-180 components of targeted buffer : GSG

    auto g_yyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 165);

    auto g_yyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 166);

    auto g_yyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 167);

    auto g_yyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 168);

    auto g_yyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 169);

    auto g_yyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 170);

    auto g_yyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 171);

    auto g_yyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 172);

    auto g_yyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 173);

    auto g_yyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 174);

    auto g_yyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 175);

    auto g_yyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 176);

    auto g_yyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 177);

    auto g_yyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 178);

    auto g_yyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 179);

    #pragma omp simd aligned(g_yyy_0_xxx_1, g_yyy_0_xxxx_1, g_yyy_0_xxxy_1, g_yyy_0_xxxz_1, g_yyy_0_xxy_1, g_yyy_0_xxyy_1, g_yyy_0_xxyz_1, g_yyy_0_xxz_1, g_yyy_0_xxzz_1, g_yyy_0_xyy_1, g_yyy_0_xyyy_1, g_yyy_0_xyyz_1, g_yyy_0_xyz_1, g_yyy_0_xyzz_1, g_yyy_0_xzz_1, g_yyy_0_xzzz_1, g_yyy_0_yyy_1, g_yyy_0_yyyy_1, g_yyy_0_yyyz_1, g_yyy_0_yyz_1, g_yyy_0_yyzz_1, g_yyy_0_yzz_1, g_yyy_0_yzzz_1, g_yyy_0_zzz_1, g_yyy_0_zzzz_1, g_yyyz_0_xxxx_0, g_yyyz_0_xxxy_0, g_yyyz_0_xxxz_0, g_yyyz_0_xxyy_0, g_yyyz_0_xxyz_0, g_yyyz_0_xxzz_0, g_yyyz_0_xyyy_0, g_yyyz_0_xyyz_0, g_yyyz_0_xyzz_0, g_yyyz_0_xzzz_0, g_yyyz_0_yyyy_0, g_yyyz_0_yyyz_0, g_yyyz_0_yyzz_0, g_yyyz_0_yzzz_0, g_yyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxxx_0[i] = g_yyy_0_xxxx_1[i] * wa_z[i];

        g_yyyz_0_xxxy_0[i] = g_yyy_0_xxxy_1[i] * wa_z[i];

        g_yyyz_0_xxxz_0[i] = g_yyy_0_xxx_1[i] * fi_acd_0 + g_yyy_0_xxxz_1[i] * wa_z[i];

        g_yyyz_0_xxyy_0[i] = g_yyy_0_xxyy_1[i] * wa_z[i];

        g_yyyz_0_xxyz_0[i] = g_yyy_0_xxy_1[i] * fi_acd_0 + g_yyy_0_xxyz_1[i] * wa_z[i];

        g_yyyz_0_xxzz_0[i] = 2.0 * g_yyy_0_xxz_1[i] * fi_acd_0 + g_yyy_0_xxzz_1[i] * wa_z[i];

        g_yyyz_0_xyyy_0[i] = g_yyy_0_xyyy_1[i] * wa_z[i];

        g_yyyz_0_xyyz_0[i] = g_yyy_0_xyy_1[i] * fi_acd_0 + g_yyy_0_xyyz_1[i] * wa_z[i];

        g_yyyz_0_xyzz_0[i] = 2.0 * g_yyy_0_xyz_1[i] * fi_acd_0 + g_yyy_0_xyzz_1[i] * wa_z[i];

        g_yyyz_0_xzzz_0[i] = 3.0 * g_yyy_0_xzz_1[i] * fi_acd_0 + g_yyy_0_xzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyy_0[i] = g_yyy_0_yyyy_1[i] * wa_z[i];

        g_yyyz_0_yyyz_0[i] = g_yyy_0_yyy_1[i] * fi_acd_0 + g_yyy_0_yyyz_1[i] * wa_z[i];

        g_yyyz_0_yyzz_0[i] = 2.0 * g_yyy_0_yyz_1[i] * fi_acd_0 + g_yyy_0_yyzz_1[i] * wa_z[i];

        g_yyyz_0_yzzz_0[i] = 3.0 * g_yyy_0_yzz_1[i] * fi_acd_0 + g_yyy_0_yzzz_1[i] * wa_z[i];

        g_yyyz_0_zzzz_0[i] = 4.0 * g_yyy_0_zzz_1[i] * fi_acd_0 + g_yyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 180-195 components of targeted buffer : GSG

    auto g_yyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 180);

    auto g_yyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 181);

    auto g_yyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 182);

    auto g_yyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 183);

    auto g_yyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 184);

    auto g_yyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 185);

    auto g_yyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 186);

    auto g_yyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 187);

    auto g_yyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 188);

    auto g_yyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 189);

    auto g_yyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 190);

    auto g_yyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 191);

    auto g_yyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 192);

    auto g_yyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 193);

    auto g_yyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 194);

    #pragma omp simd aligned(g_yy_0_xxxy_0, g_yy_0_xxxy_1, g_yy_0_xxyy_0, g_yy_0_xxyy_1, g_yy_0_xyyy_0, g_yy_0_xyyy_1, g_yy_0_yyyy_0, g_yy_0_yyyy_1, g_yyz_0_xxxy_1, g_yyz_0_xxyy_1, g_yyz_0_xyyy_1, g_yyz_0_yyyy_1, g_yyzz_0_xxxx_0, g_yyzz_0_xxxy_0, g_yyzz_0_xxxz_0, g_yyzz_0_xxyy_0, g_yyzz_0_xxyz_0, g_yyzz_0_xxzz_0, g_yyzz_0_xyyy_0, g_yyzz_0_xyyz_0, g_yyzz_0_xyzz_0, g_yyzz_0_xzzz_0, g_yyzz_0_yyyy_0, g_yyzz_0_yyyz_0, g_yyzz_0_yyzz_0, g_yyzz_0_yzzz_0, g_yyzz_0_zzzz_0, g_yzz_0_xxxx_1, g_yzz_0_xxxz_1, g_yzz_0_xxyz_1, g_yzz_0_xxz_1, g_yzz_0_xxzz_1, g_yzz_0_xyyz_1, g_yzz_0_xyz_1, g_yzz_0_xyzz_1, g_yzz_0_xzz_1, g_yzz_0_xzzz_1, g_yzz_0_yyyz_1, g_yzz_0_yyz_1, g_yzz_0_yyzz_1, g_yzz_0_yzz_1, g_yzz_0_yzzz_1, g_yzz_0_zzz_1, g_yzz_0_zzzz_1, g_zz_0_xxxx_0, g_zz_0_xxxx_1, g_zz_0_xxxz_0, g_zz_0_xxxz_1, g_zz_0_xxyz_0, g_zz_0_xxyz_1, g_zz_0_xxzz_0, g_zz_0_xxzz_1, g_zz_0_xyyz_0, g_zz_0_xyyz_1, g_zz_0_xyzz_0, g_zz_0_xyzz_1, g_zz_0_xzzz_0, g_zz_0_xzzz_1, g_zz_0_yyyz_0, g_zz_0_yyyz_1, g_zz_0_yyzz_0, g_zz_0_yyzz_1, g_zz_0_yzzz_0, g_zz_0_yzzz_1, g_zz_0_zzzz_0, g_zz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxxx_0[i] = g_zz_0_xxxx_0[i] * fbe_0 - g_zz_0_xxxx_1[i] * fz_be_0 + g_yzz_0_xxxx_1[i] * wa_y[i];

        g_yyzz_0_xxxy_0[i] = g_yy_0_xxxy_0[i] * fbe_0 - g_yy_0_xxxy_1[i] * fz_be_0 + g_yyz_0_xxxy_1[i] * wa_z[i];

        g_yyzz_0_xxxz_0[i] = g_zz_0_xxxz_0[i] * fbe_0 - g_zz_0_xxxz_1[i] * fz_be_0 + g_yzz_0_xxxz_1[i] * wa_y[i];

        g_yyzz_0_xxyy_0[i] = g_yy_0_xxyy_0[i] * fbe_0 - g_yy_0_xxyy_1[i] * fz_be_0 + g_yyz_0_xxyy_1[i] * wa_z[i];

        g_yyzz_0_xxyz_0[i] = g_zz_0_xxyz_0[i] * fbe_0 - g_zz_0_xxyz_1[i] * fz_be_0 + g_yzz_0_xxz_1[i] * fi_acd_0 + g_yzz_0_xxyz_1[i] * wa_y[i];

        g_yyzz_0_xxzz_0[i] = g_zz_0_xxzz_0[i] * fbe_0 - g_zz_0_xxzz_1[i] * fz_be_0 + g_yzz_0_xxzz_1[i] * wa_y[i];

        g_yyzz_0_xyyy_0[i] = g_yy_0_xyyy_0[i] * fbe_0 - g_yy_0_xyyy_1[i] * fz_be_0 + g_yyz_0_xyyy_1[i] * wa_z[i];

        g_yyzz_0_xyyz_0[i] = g_zz_0_xyyz_0[i] * fbe_0 - g_zz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xyz_1[i] * fi_acd_0 + g_yzz_0_xyyz_1[i] * wa_y[i];

        g_yyzz_0_xyzz_0[i] = g_zz_0_xyzz_0[i] * fbe_0 - g_zz_0_xyzz_1[i] * fz_be_0 + g_yzz_0_xzz_1[i] * fi_acd_0 + g_yzz_0_xyzz_1[i] * wa_y[i];

        g_yyzz_0_xzzz_0[i] = g_zz_0_xzzz_0[i] * fbe_0 - g_zz_0_xzzz_1[i] * fz_be_0 + g_yzz_0_xzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyy_0[i] = g_yy_0_yyyy_0[i] * fbe_0 - g_yy_0_yyyy_1[i] * fz_be_0 + g_yyz_0_yyyy_1[i] * wa_z[i];

        g_yyzz_0_yyyz_0[i] = g_zz_0_yyyz_0[i] * fbe_0 - g_zz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yzz_0_yyz_1[i] * fi_acd_0 + g_yzz_0_yyyz_1[i] * wa_y[i];

        g_yyzz_0_yyzz_0[i] = g_zz_0_yyzz_0[i] * fbe_0 - g_zz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yzz_1[i] * fi_acd_0 + g_yzz_0_yyzz_1[i] * wa_y[i];

        g_yyzz_0_yzzz_0[i] = g_zz_0_yzzz_0[i] * fbe_0 - g_zz_0_yzzz_1[i] * fz_be_0 + g_yzz_0_zzz_1[i] * fi_acd_0 + g_yzz_0_yzzz_1[i] * wa_y[i];

        g_yyzz_0_zzzz_0[i] = g_zz_0_zzzz_0[i] * fbe_0 - g_zz_0_zzzz_1[i] * fz_be_0 + g_yzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 195-210 components of targeted buffer : GSG

    auto g_yzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 195);

    auto g_yzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 196);

    auto g_yzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 197);

    auto g_yzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 198);

    auto g_yzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 199);

    auto g_yzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 200);

    auto g_yzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 201);

    auto g_yzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 202);

    auto g_yzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 203);

    auto g_yzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 204);

    auto g_yzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 205);

    auto g_yzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 206);

    auto g_yzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 207);

    auto g_yzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 208);

    auto g_yzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 209);

    #pragma omp simd aligned(g_yzzz_0_xxxx_0, g_yzzz_0_xxxy_0, g_yzzz_0_xxxz_0, g_yzzz_0_xxyy_0, g_yzzz_0_xxyz_0, g_yzzz_0_xxzz_0, g_yzzz_0_xyyy_0, g_yzzz_0_xyyz_0, g_yzzz_0_xyzz_0, g_yzzz_0_xzzz_0, g_yzzz_0_yyyy_0, g_yzzz_0_yyyz_0, g_yzzz_0_yyzz_0, g_yzzz_0_yzzz_0, g_yzzz_0_zzzz_0, g_zzz_0_xxx_1, g_zzz_0_xxxx_1, g_zzz_0_xxxy_1, g_zzz_0_xxxz_1, g_zzz_0_xxy_1, g_zzz_0_xxyy_1, g_zzz_0_xxyz_1, g_zzz_0_xxz_1, g_zzz_0_xxzz_1, g_zzz_0_xyy_1, g_zzz_0_xyyy_1, g_zzz_0_xyyz_1, g_zzz_0_xyz_1, g_zzz_0_xyzz_1, g_zzz_0_xzz_1, g_zzz_0_xzzz_1, g_zzz_0_yyy_1, g_zzz_0_yyyy_1, g_zzz_0_yyyz_1, g_zzz_0_yyz_1, g_zzz_0_yyzz_1, g_zzz_0_yzz_1, g_zzz_0_yzzz_1, g_zzz_0_zzz_1, g_zzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxxx_0[i] = g_zzz_0_xxxx_1[i] * wa_y[i];

        g_yzzz_0_xxxy_0[i] = g_zzz_0_xxx_1[i] * fi_acd_0 + g_zzz_0_xxxy_1[i] * wa_y[i];

        g_yzzz_0_xxxz_0[i] = g_zzz_0_xxxz_1[i] * wa_y[i];

        g_yzzz_0_xxyy_0[i] = 2.0 * g_zzz_0_xxy_1[i] * fi_acd_0 + g_zzz_0_xxyy_1[i] * wa_y[i];

        g_yzzz_0_xxyz_0[i] = g_zzz_0_xxz_1[i] * fi_acd_0 + g_zzz_0_xxyz_1[i] * wa_y[i];

        g_yzzz_0_xxzz_0[i] = g_zzz_0_xxzz_1[i] * wa_y[i];

        g_yzzz_0_xyyy_0[i] = 3.0 * g_zzz_0_xyy_1[i] * fi_acd_0 + g_zzz_0_xyyy_1[i] * wa_y[i];

        g_yzzz_0_xyyz_0[i] = 2.0 * g_zzz_0_xyz_1[i] * fi_acd_0 + g_zzz_0_xyyz_1[i] * wa_y[i];

        g_yzzz_0_xyzz_0[i] = g_zzz_0_xzz_1[i] * fi_acd_0 + g_zzz_0_xyzz_1[i] * wa_y[i];

        g_yzzz_0_xzzz_0[i] = g_zzz_0_xzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyy_0[i] = 4.0 * g_zzz_0_yyy_1[i] * fi_acd_0 + g_zzz_0_yyyy_1[i] * wa_y[i];

        g_yzzz_0_yyyz_0[i] = 3.0 * g_zzz_0_yyz_1[i] * fi_acd_0 + g_zzz_0_yyyz_1[i] * wa_y[i];

        g_yzzz_0_yyzz_0[i] = 2.0 * g_zzz_0_yzz_1[i] * fi_acd_0 + g_zzz_0_yyzz_1[i] * wa_y[i];

        g_yzzz_0_yzzz_0[i] = g_zzz_0_zzz_1[i] * fi_acd_0 + g_zzz_0_yzzz_1[i] * wa_y[i];

        g_yzzz_0_zzzz_0[i] = g_zzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 210-225 components of targeted buffer : GSG

    auto g_zzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_gsg + 210);

    auto g_zzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_gsg + 211);

    auto g_zzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_gsg + 212);

    auto g_zzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_gsg + 213);

    auto g_zzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_gsg + 214);

    auto g_zzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_gsg + 215);

    auto g_zzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_gsg + 216);

    auto g_zzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_gsg + 217);

    auto g_zzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_gsg + 218);

    auto g_zzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_gsg + 219);

    auto g_zzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_gsg + 220);

    auto g_zzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_gsg + 221);

    auto g_zzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_gsg + 222);

    auto g_zzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_gsg + 223);

    auto g_zzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_gsg + 224);

    #pragma omp simd aligned(g_zz_0_xxxx_0, g_zz_0_xxxx_1, g_zz_0_xxxy_0, g_zz_0_xxxy_1, g_zz_0_xxxz_0, g_zz_0_xxxz_1, g_zz_0_xxyy_0, g_zz_0_xxyy_1, g_zz_0_xxyz_0, g_zz_0_xxyz_1, g_zz_0_xxzz_0, g_zz_0_xxzz_1, g_zz_0_xyyy_0, g_zz_0_xyyy_1, g_zz_0_xyyz_0, g_zz_0_xyyz_1, g_zz_0_xyzz_0, g_zz_0_xyzz_1, g_zz_0_xzzz_0, g_zz_0_xzzz_1, g_zz_0_yyyy_0, g_zz_0_yyyy_1, g_zz_0_yyyz_0, g_zz_0_yyyz_1, g_zz_0_yyzz_0, g_zz_0_yyzz_1, g_zz_0_yzzz_0, g_zz_0_yzzz_1, g_zz_0_zzzz_0, g_zz_0_zzzz_1, g_zzz_0_xxx_1, g_zzz_0_xxxx_1, g_zzz_0_xxxy_1, g_zzz_0_xxxz_1, g_zzz_0_xxy_1, g_zzz_0_xxyy_1, g_zzz_0_xxyz_1, g_zzz_0_xxz_1, g_zzz_0_xxzz_1, g_zzz_0_xyy_1, g_zzz_0_xyyy_1, g_zzz_0_xyyz_1, g_zzz_0_xyz_1, g_zzz_0_xyzz_1, g_zzz_0_xzz_1, g_zzz_0_xzzz_1, g_zzz_0_yyy_1, g_zzz_0_yyyy_1, g_zzz_0_yyyz_1, g_zzz_0_yyz_1, g_zzz_0_yyzz_1, g_zzz_0_yzz_1, g_zzz_0_yzzz_1, g_zzz_0_zzz_1, g_zzz_0_zzzz_1, g_zzzz_0_xxxx_0, g_zzzz_0_xxxy_0, g_zzzz_0_xxxz_0, g_zzzz_0_xxyy_0, g_zzzz_0_xxyz_0, g_zzzz_0_xxzz_0, g_zzzz_0_xyyy_0, g_zzzz_0_xyyz_0, g_zzzz_0_xyzz_0, g_zzzz_0_xzzz_0, g_zzzz_0_yyyy_0, g_zzzz_0_yyyz_0, g_zzzz_0_yyzz_0, g_zzzz_0_yzzz_0, g_zzzz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxxx_0[i] = 3.0 * g_zz_0_xxxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxxx_1[i] * fz_be_0 + g_zzz_0_xxxx_1[i] * wa_z[i];

        g_zzzz_0_xxxy_0[i] = 3.0 * g_zz_0_xxxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxy_1[i] * fz_be_0 + g_zzz_0_xxxy_1[i] * wa_z[i];

        g_zzzz_0_xxxz_0[i] = 3.0 * g_zz_0_xxxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxz_1[i] * fz_be_0 + g_zzz_0_xxx_1[i] * fi_acd_0 + g_zzz_0_xxxz_1[i] * wa_z[i];

        g_zzzz_0_xxyy_0[i] = 3.0 * g_zz_0_xxyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxyy_1[i] * fz_be_0 + g_zzz_0_xxyy_1[i] * wa_z[i];

        g_zzzz_0_xxyz_0[i] = 3.0 * g_zz_0_xxyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyz_1[i] * fz_be_0 + g_zzz_0_xxy_1[i] * fi_acd_0 + g_zzz_0_xxyz_1[i] * wa_z[i];

        g_zzzz_0_xxzz_0[i] = 3.0 * g_zz_0_xxzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxz_1[i] * fi_acd_0 + g_zzz_0_xxzz_1[i] * wa_z[i];

        g_zzzz_0_xyyy_0[i] = 3.0 * g_zz_0_xyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyyy_1[i] * fz_be_0 + g_zzz_0_xyyy_1[i] * wa_z[i];

        g_zzzz_0_xyyz_0[i] = 3.0 * g_zz_0_xyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyz_1[i] * fz_be_0 + g_zzz_0_xyy_1[i] * fi_acd_0 + g_zzz_0_xyyz_1[i] * wa_z[i];

        g_zzzz_0_xyzz_0[i] = 3.0 * g_zz_0_xyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xyz_1[i] * fi_acd_0 + g_zzz_0_xyzz_1[i] * wa_z[i];

        g_zzzz_0_xzzz_0[i] = 3.0 * g_zz_0_xzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xzz_1[i] * fi_acd_0 + g_zzz_0_xzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyy_0[i] = 3.0 * g_zz_0_yyyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyyy_1[i] * fz_be_0 + g_zzz_0_yyyy_1[i] * wa_z[i];

        g_zzzz_0_yyyz_0[i] = 3.0 * g_zz_0_yyyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyz_1[i] * fz_be_0 + g_zzz_0_yyy_1[i] * fi_acd_0 + g_zzz_0_yyyz_1[i] * wa_z[i];

        g_zzzz_0_yyzz_0[i] = 3.0 * g_zz_0_yyzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yyz_1[i] * fi_acd_0 + g_zzz_0_yyzz_1[i] * wa_z[i];

        g_zzzz_0_yzzz_0[i] = 3.0 * g_zz_0_yzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_yzz_1[i] * fi_acd_0 + g_zzz_0_yzzz_1[i] * wa_z[i];

        g_zzzz_0_zzzz_0[i] = 3.0 * g_zz_0_zzzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_zzz_1[i] * fi_acd_0 + g_zzz_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

