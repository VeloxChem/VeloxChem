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

#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsf,
                                 size_t idx_eri_0_fsf,
                                 size_t idx_eri_1_fsf,
                                 size_t idx_eri_1_gsd,
                                 size_t idx_eri_1_gsf,
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

    /// Set up components of auxilary buffer : FSF

    auto g_xxx_0_xxx_0 = pbuffer.data(idx_eri_0_fsf);

    auto g_xxx_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 1);

    auto g_xxx_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 2);

    auto g_xxx_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 3);

    auto g_xxx_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 4);

    auto g_xxx_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 5);

    auto g_xxx_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 6);

    auto g_xxx_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 7);

    auto g_xxx_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 8);

    auto g_xxx_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 9);

    auto g_xxy_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 10);

    auto g_xxy_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 12);

    auto g_xxy_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 15);

    auto g_xxz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 20);

    auto g_xxz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 21);

    auto g_xxz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 23);

    auto g_xyy_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 31);

    auto g_xyy_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 33);

    auto g_xyy_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 34);

    auto g_xyy_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 36);

    auto g_xyy_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 37);

    auto g_xyy_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 38);

    auto g_xyy_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 39);

    auto g_xzz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 52);

    auto g_xzz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 54);

    auto g_xzz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 55);

    auto g_xzz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 56);

    auto g_xzz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 57);

    auto g_xzz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 58);

    auto g_xzz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 59);

    auto g_yyy_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 60);

    auto g_yyy_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 61);

    auto g_yyy_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 62);

    auto g_yyy_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 63);

    auto g_yyy_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 64);

    auto g_yyy_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 65);

    auto g_yyy_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 66);

    auto g_yyy_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 67);

    auto g_yyy_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 68);

    auto g_yyy_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 69);

    auto g_yyz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 71);

    auto g_yyz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 73);

    auto g_yyz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 76);

    auto g_yzz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 80);

    auto g_yzz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 82);

    auto g_yzz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 84);

    auto g_yzz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 85);

    auto g_yzz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 87);

    auto g_yzz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 88);

    auto g_yzz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 89);

    auto g_zzz_0_xxx_0 = pbuffer.data(idx_eri_0_fsf + 90);

    auto g_zzz_0_xxy_0 = pbuffer.data(idx_eri_0_fsf + 91);

    auto g_zzz_0_xxz_0 = pbuffer.data(idx_eri_0_fsf + 92);

    auto g_zzz_0_xyy_0 = pbuffer.data(idx_eri_0_fsf + 93);

    auto g_zzz_0_xyz_0 = pbuffer.data(idx_eri_0_fsf + 94);

    auto g_zzz_0_xzz_0 = pbuffer.data(idx_eri_0_fsf + 95);

    auto g_zzz_0_yyy_0 = pbuffer.data(idx_eri_0_fsf + 96);

    auto g_zzz_0_yyz_0 = pbuffer.data(idx_eri_0_fsf + 97);

    auto g_zzz_0_yzz_0 = pbuffer.data(idx_eri_0_fsf + 98);

    auto g_zzz_0_zzz_0 = pbuffer.data(idx_eri_0_fsf + 99);

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

    auto g_xxy_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 10);

    auto g_xxy_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 12);

    auto g_xxy_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 15);

    auto g_xxz_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 20);

    auto g_xxz_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 21);

    auto g_xxz_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 23);

    auto g_xyy_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 31);

    auto g_xyy_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 33);

    auto g_xyy_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 34);

    auto g_xyy_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 36);

    auto g_xyy_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 37);

    auto g_xyy_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 38);

    auto g_xyy_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 39);

    auto g_xzz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 52);

    auto g_xzz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 54);

    auto g_xzz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 55);

    auto g_xzz_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 56);

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

    auto g_yyz_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 71);

    auto g_yyz_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 73);

    auto g_yyz_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 76);

    auto g_yzz_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 80);

    auto g_yzz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 82);

    auto g_yzz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 84);

    auto g_yzz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 85);

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

    /// Set up components of auxilary buffer : GSD

    auto g_xxxx_0_xx_1 = pbuffer.data(idx_eri_1_gsd);

    auto g_xxxx_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 1);

    auto g_xxxx_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 2);

    auto g_xxxx_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 3);

    auto g_xxxx_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 4);

    auto g_xxxx_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 5);

    auto g_xxxz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 14);

    auto g_xxxz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 16);

    auto g_xxxz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 17);

    auto g_xxyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 18);

    auto g_xxyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 19);

    auto g_xxyy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 20);

    auto g_xxyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 21);

    auto g_xxyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 22);

    auto g_xxyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 23);

    auto g_xxzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 30);

    auto g_xxzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 31);

    auto g_xxzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 32);

    auto g_xxzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 33);

    auto g_xxzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 34);

    auto g_xxzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 35);

    auto g_xyyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 37);

    auto g_xyyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 39);

    auto g_xyyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 40);

    auto g_xzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 56);

    auto g_xzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 58);

    auto g_xzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 59);

    auto g_yyyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 60);

    auto g_yyyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 61);

    auto g_yyyy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 62);

    auto g_yyyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 63);

    auto g_yyyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 64);

    auto g_yyyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 65);

    auto g_yyyz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 68);

    auto g_yyyz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 70);

    auto g_yyyz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 71);

    auto g_yyzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 72);

    auto g_yyzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 73);

    auto g_yyzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 74);

    auto g_yyzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 75);

    auto g_yyzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 76);

    auto g_yyzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 77);

    auto g_yzzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 79);

    auto g_yzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 80);

    auto g_yzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 81);

    auto g_yzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 82);

    auto g_yzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 83);

    auto g_zzzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 84);

    auto g_zzzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 85);

    auto g_zzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 86);

    auto g_zzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 87);

    auto g_zzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 88);

    auto g_zzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 89);

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

    auto g_xxxy_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 11);

    auto g_xxxy_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 12);

    auto g_xxxy_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 13);

    auto g_xxxy_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 15);

    auto g_xxxy_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 16);

    auto g_xxxz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 20);

    auto g_xxxz_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 21);

    auto g_xxxz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 22);

    auto g_xxxz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 23);

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

    auto g_xyyy_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 60);

    auto g_xyyy_0_xxy_1 = pbuffer.data(idx_eri_1_gsf + 61);

    auto g_xyyy_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 63);

    auto g_xyyy_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 64);

    auto g_xyyy_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 66);

    auto g_xyyy_0_yyz_1 = pbuffer.data(idx_eri_1_gsf + 67);

    auto g_xyyy_0_yzz_1 = pbuffer.data(idx_eri_1_gsf + 68);

    auto g_xyyy_0_zzz_1 = pbuffer.data(idx_eri_1_gsf + 69);

    auto g_xzzz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 90);

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

    auto g_yyyz_0_xxz_1 = pbuffer.data(idx_eri_1_gsf + 112);

    auto g_yyyz_0_xyy_1 = pbuffer.data(idx_eri_1_gsf + 113);

    auto g_yyyz_0_xyz_1 = pbuffer.data(idx_eri_1_gsf + 114);

    auto g_yyyz_0_xzz_1 = pbuffer.data(idx_eri_1_gsf + 115);

    auto g_yyyz_0_yyy_1 = pbuffer.data(idx_eri_1_gsf + 116);

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

    auto g_yzzz_0_xxx_1 = pbuffer.data(idx_eri_1_gsf + 130);

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

    /// Set up 0-10 components of targeted buffer : HSF

    auto g_xxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_hsf);

    auto g_xxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 1);

    auto g_xxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 2);

    auto g_xxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 3);

    auto g_xxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 4);

    auto g_xxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 5);

    auto g_xxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 6);

    auto g_xxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 7);

    auto g_xxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 8);

    auto g_xxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 9);

    #pragma omp simd aligned(g_xxx_0_xxx_0, g_xxx_0_xxx_1, g_xxx_0_xxy_0, g_xxx_0_xxy_1, g_xxx_0_xxz_0, g_xxx_0_xxz_1, g_xxx_0_xyy_0, g_xxx_0_xyy_1, g_xxx_0_xyz_0, g_xxx_0_xyz_1, g_xxx_0_xzz_0, g_xxx_0_xzz_1, g_xxx_0_yyy_0, g_xxx_0_yyy_1, g_xxx_0_yyz_0, g_xxx_0_yyz_1, g_xxx_0_yzz_0, g_xxx_0_yzz_1, g_xxx_0_zzz_0, g_xxx_0_zzz_1, g_xxxx_0_xx_1, g_xxxx_0_xxx_1, g_xxxx_0_xxy_1, g_xxxx_0_xxz_1, g_xxxx_0_xy_1, g_xxxx_0_xyy_1, g_xxxx_0_xyz_1, g_xxxx_0_xz_1, g_xxxx_0_xzz_1, g_xxxx_0_yy_1, g_xxxx_0_yyy_1, g_xxxx_0_yyz_1, g_xxxx_0_yz_1, g_xxxx_0_yzz_1, g_xxxx_0_zz_1, g_xxxx_0_zzz_1, g_xxxxx_0_xxx_0, g_xxxxx_0_xxy_0, g_xxxxx_0_xxz_0, g_xxxxx_0_xyy_0, g_xxxxx_0_xyz_0, g_xxxxx_0_xzz_0, g_xxxxx_0_yyy_0, g_xxxxx_0_yyz_0, g_xxxxx_0_yzz_0, g_xxxxx_0_zzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xxx_0[i] = 4.0 * g_xxx_0_xxx_0[i] * fbe_0 - 4.0 * g_xxx_0_xxx_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xx_1[i] * fi_acd_0 + g_xxxx_0_xxx_1[i] * wa_x[i];

        g_xxxxx_0_xxy_0[i] = 4.0 * g_xxx_0_xxy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xy_1[i] * fi_acd_0 + g_xxxx_0_xxy_1[i] * wa_x[i];

        g_xxxxx_0_xxz_0[i] = 4.0 * g_xxx_0_xxz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xz_1[i] * fi_acd_0 + g_xxxx_0_xxz_1[i] * wa_x[i];

        g_xxxxx_0_xyy_0[i] = 4.0 * g_xxx_0_xyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xyy_1[i] * fz_be_0 + g_xxxx_0_yy_1[i] * fi_acd_0 + g_xxxx_0_xyy_1[i] * wa_x[i];

        g_xxxxx_0_xyz_0[i] = 4.0 * g_xxx_0_xyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyz_1[i] * fz_be_0 + g_xxxx_0_yz_1[i] * fi_acd_0 + g_xxxx_0_xyz_1[i] * wa_x[i];

        g_xxxxx_0_xzz_0[i] = 4.0 * g_xxx_0_xzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xzz_1[i] * fz_be_0 + g_xxxx_0_zz_1[i] * fi_acd_0 + g_xxxx_0_xzz_1[i] * wa_x[i];

        g_xxxxx_0_yyy_0[i] = 4.0 * g_xxx_0_yyy_0[i] * fbe_0 - 4.0 * g_xxx_0_yyy_1[i] * fz_be_0 + g_xxxx_0_yyy_1[i] * wa_x[i];

        g_xxxxx_0_yyz_0[i] = 4.0 * g_xxx_0_yyz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyz_1[i] * fz_be_0 + g_xxxx_0_yyz_1[i] * wa_x[i];

        g_xxxxx_0_yzz_0[i] = 4.0 * g_xxx_0_yzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yzz_1[i] * fz_be_0 + g_xxxx_0_yzz_1[i] * wa_x[i];

        g_xxxxx_0_zzz_0[i] = 4.0 * g_xxx_0_zzz_0[i] * fbe_0 - 4.0 * g_xxx_0_zzz_1[i] * fz_be_0 + g_xxxx_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : HSF

    auto g_xxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 10);

    auto g_xxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 11);

    auto g_xxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 12);

    auto g_xxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 13);

    auto g_xxxxy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 14);

    auto g_xxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 15);

    auto g_xxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 16);

    auto g_xxxxy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 17);

    auto g_xxxxy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 18);

    auto g_xxxxy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 19);

    #pragma omp simd aligned(g_xxxx_0_xx_1, g_xxxx_0_xxx_1, g_xxxx_0_xxy_1, g_xxxx_0_xxz_1, g_xxxx_0_xy_1, g_xxxx_0_xyy_1, g_xxxx_0_xyz_1, g_xxxx_0_xz_1, g_xxxx_0_xzz_1, g_xxxx_0_yy_1, g_xxxx_0_yyy_1, g_xxxx_0_yyz_1, g_xxxx_0_yz_1, g_xxxx_0_yzz_1, g_xxxx_0_zz_1, g_xxxx_0_zzz_1, g_xxxxy_0_xxx_0, g_xxxxy_0_xxy_0, g_xxxxy_0_xxz_0, g_xxxxy_0_xyy_0, g_xxxxy_0_xyz_0, g_xxxxy_0_xzz_0, g_xxxxy_0_yyy_0, g_xxxxy_0_yyz_0, g_xxxxy_0_yzz_0, g_xxxxy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xxx_0[i] = g_xxxx_0_xxx_1[i] * wa_y[i];

        g_xxxxy_0_xxy_0[i] = g_xxxx_0_xx_1[i] * fi_acd_0 + g_xxxx_0_xxy_1[i] * wa_y[i];

        g_xxxxy_0_xxz_0[i] = g_xxxx_0_xxz_1[i] * wa_y[i];

        g_xxxxy_0_xyy_0[i] = 2.0 * g_xxxx_0_xy_1[i] * fi_acd_0 + g_xxxx_0_xyy_1[i] * wa_y[i];

        g_xxxxy_0_xyz_0[i] = g_xxxx_0_xz_1[i] * fi_acd_0 + g_xxxx_0_xyz_1[i] * wa_y[i];

        g_xxxxy_0_xzz_0[i] = g_xxxx_0_xzz_1[i] * wa_y[i];

        g_xxxxy_0_yyy_0[i] = 3.0 * g_xxxx_0_yy_1[i] * fi_acd_0 + g_xxxx_0_yyy_1[i] * wa_y[i];

        g_xxxxy_0_yyz_0[i] = 2.0 * g_xxxx_0_yz_1[i] * fi_acd_0 + g_xxxx_0_yyz_1[i] * wa_y[i];

        g_xxxxy_0_yzz_0[i] = g_xxxx_0_zz_1[i] * fi_acd_0 + g_xxxx_0_yzz_1[i] * wa_y[i];

        g_xxxxy_0_zzz_0[i] = g_xxxx_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 20-30 components of targeted buffer : HSF

    auto g_xxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 20);

    auto g_xxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 21);

    auto g_xxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 22);

    auto g_xxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 23);

    auto g_xxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 24);

    auto g_xxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 25);

    auto g_xxxxz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 26);

    auto g_xxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 27);

    auto g_xxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 28);

    auto g_xxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 29);

    #pragma omp simd aligned(g_xxxx_0_xx_1, g_xxxx_0_xxx_1, g_xxxx_0_xxy_1, g_xxxx_0_xxz_1, g_xxxx_0_xy_1, g_xxxx_0_xyy_1, g_xxxx_0_xyz_1, g_xxxx_0_xz_1, g_xxxx_0_xzz_1, g_xxxx_0_yy_1, g_xxxx_0_yyy_1, g_xxxx_0_yyz_1, g_xxxx_0_yz_1, g_xxxx_0_yzz_1, g_xxxx_0_zz_1, g_xxxx_0_zzz_1, g_xxxxz_0_xxx_0, g_xxxxz_0_xxy_0, g_xxxxz_0_xxz_0, g_xxxxz_0_xyy_0, g_xxxxz_0_xyz_0, g_xxxxz_0_xzz_0, g_xxxxz_0_yyy_0, g_xxxxz_0_yyz_0, g_xxxxz_0_yzz_0, g_xxxxz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xxx_0[i] = g_xxxx_0_xxx_1[i] * wa_z[i];

        g_xxxxz_0_xxy_0[i] = g_xxxx_0_xxy_1[i] * wa_z[i];

        g_xxxxz_0_xxz_0[i] = g_xxxx_0_xx_1[i] * fi_acd_0 + g_xxxx_0_xxz_1[i] * wa_z[i];

        g_xxxxz_0_xyy_0[i] = g_xxxx_0_xyy_1[i] * wa_z[i];

        g_xxxxz_0_xyz_0[i] = g_xxxx_0_xy_1[i] * fi_acd_0 + g_xxxx_0_xyz_1[i] * wa_z[i];

        g_xxxxz_0_xzz_0[i] = 2.0 * g_xxxx_0_xz_1[i] * fi_acd_0 + g_xxxx_0_xzz_1[i] * wa_z[i];

        g_xxxxz_0_yyy_0[i] = g_xxxx_0_yyy_1[i] * wa_z[i];

        g_xxxxz_0_yyz_0[i] = g_xxxx_0_yy_1[i] * fi_acd_0 + g_xxxx_0_yyz_1[i] * wa_z[i];

        g_xxxxz_0_yzz_0[i] = 2.0 * g_xxxx_0_yz_1[i] * fi_acd_0 + g_xxxx_0_yzz_1[i] * wa_z[i];

        g_xxxxz_0_zzz_0[i] = 3.0 * g_xxxx_0_zz_1[i] * fi_acd_0 + g_xxxx_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 30-40 components of targeted buffer : HSF

    auto g_xxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 30);

    auto g_xxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 31);

    auto g_xxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 32);

    auto g_xxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 33);

    auto g_xxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 34);

    auto g_xxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 35);

    auto g_xxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 36);

    auto g_xxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 37);

    auto g_xxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 38);

    auto g_xxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 39);

    #pragma omp simd aligned(g_xxx_0_xxx_0, g_xxx_0_xxx_1, g_xxx_0_xxz_0, g_xxx_0_xxz_1, g_xxx_0_xzz_0, g_xxx_0_xzz_1, g_xxxy_0_xxx_1, g_xxxy_0_xxz_1, g_xxxy_0_xzz_1, g_xxxyy_0_xxx_0, g_xxxyy_0_xxy_0, g_xxxyy_0_xxz_0, g_xxxyy_0_xyy_0, g_xxxyy_0_xyz_0, g_xxxyy_0_xzz_0, g_xxxyy_0_yyy_0, g_xxxyy_0_yyz_0, g_xxxyy_0_yzz_0, g_xxxyy_0_zzz_0, g_xxyy_0_xxy_1, g_xxyy_0_xy_1, g_xxyy_0_xyy_1, g_xxyy_0_xyz_1, g_xxyy_0_yy_1, g_xxyy_0_yyy_1, g_xxyy_0_yyz_1, g_xxyy_0_yz_1, g_xxyy_0_yzz_1, g_xxyy_0_zzz_1, g_xyy_0_xxy_0, g_xyy_0_xxy_1, g_xyy_0_xyy_0, g_xyy_0_xyy_1, g_xyy_0_xyz_0, g_xyy_0_xyz_1, g_xyy_0_yyy_0, g_xyy_0_yyy_1, g_xyy_0_yyz_0, g_xyy_0_yyz_1, g_xyy_0_yzz_0, g_xyy_0_yzz_1, g_xyy_0_zzz_0, g_xyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xxx_0[i] = g_xxx_0_xxx_0[i] * fbe_0 - g_xxx_0_xxx_1[i] * fz_be_0 + g_xxxy_0_xxx_1[i] * wa_y[i];

        g_xxxyy_0_xxy_0[i] = 2.0 * g_xyy_0_xxy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xy_1[i] * fi_acd_0 + g_xxyy_0_xxy_1[i] * wa_x[i];

        g_xxxyy_0_xxz_0[i] = g_xxx_0_xxz_0[i] * fbe_0 - g_xxx_0_xxz_1[i] * fz_be_0 + g_xxxy_0_xxz_1[i] * wa_y[i];

        g_xxxyy_0_xyy_0[i] = 2.0 * g_xyy_0_xyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xyy_1[i] * fz_be_0 + g_xxyy_0_yy_1[i] * fi_acd_0 + g_xxyy_0_xyy_1[i] * wa_x[i];

        g_xxxyy_0_xyz_0[i] = 2.0 * g_xyy_0_xyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyz_1[i] * fz_be_0 + g_xxyy_0_yz_1[i] * fi_acd_0 + g_xxyy_0_xyz_1[i] * wa_x[i];

        g_xxxyy_0_xzz_0[i] = g_xxx_0_xzz_0[i] * fbe_0 - g_xxx_0_xzz_1[i] * fz_be_0 + g_xxxy_0_xzz_1[i] * wa_y[i];

        g_xxxyy_0_yyy_0[i] = 2.0 * g_xyy_0_yyy_0[i] * fbe_0 - 2.0 * g_xyy_0_yyy_1[i] * fz_be_0 + g_xxyy_0_yyy_1[i] * wa_x[i];

        g_xxxyy_0_yyz_0[i] = 2.0 * g_xyy_0_yyz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyz_1[i] * fz_be_0 + g_xxyy_0_yyz_1[i] * wa_x[i];

        g_xxxyy_0_yzz_0[i] = 2.0 * g_xyy_0_yzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yzz_1[i] * fz_be_0 + g_xxyy_0_yzz_1[i] * wa_x[i];

        g_xxxyy_0_zzz_0[i] = 2.0 * g_xyy_0_zzz_0[i] * fbe_0 - 2.0 * g_xyy_0_zzz_1[i] * fz_be_0 + g_xxyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 40-50 components of targeted buffer : HSF

    auto g_xxxyz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 40);

    auto g_xxxyz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 41);

    auto g_xxxyz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 42);

    auto g_xxxyz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 43);

    auto g_xxxyz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 44);

    auto g_xxxyz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 45);

    auto g_xxxyz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 46);

    auto g_xxxyz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 47);

    auto g_xxxyz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 48);

    auto g_xxxyz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 49);

    #pragma omp simd aligned(g_xxxy_0_xxy_1, g_xxxy_0_xyy_1, g_xxxy_0_yyy_1, g_xxxyz_0_xxx_0, g_xxxyz_0_xxy_0, g_xxxyz_0_xxz_0, g_xxxyz_0_xyy_0, g_xxxyz_0_xyz_0, g_xxxyz_0_xzz_0, g_xxxyz_0_yyy_0, g_xxxyz_0_yyz_0, g_xxxyz_0_yzz_0, g_xxxyz_0_zzz_0, g_xxxz_0_xxx_1, g_xxxz_0_xxz_1, g_xxxz_0_xyz_1, g_xxxz_0_xz_1, g_xxxz_0_xzz_1, g_xxxz_0_yyz_1, g_xxxz_0_yz_1, g_xxxz_0_yzz_1, g_xxxz_0_zz_1, g_xxxz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xxx_0[i] = g_xxxz_0_xxx_1[i] * wa_y[i];

        g_xxxyz_0_xxy_0[i] = g_xxxy_0_xxy_1[i] * wa_z[i];

        g_xxxyz_0_xxz_0[i] = g_xxxz_0_xxz_1[i] * wa_y[i];

        g_xxxyz_0_xyy_0[i] = g_xxxy_0_xyy_1[i] * wa_z[i];

        g_xxxyz_0_xyz_0[i] = g_xxxz_0_xz_1[i] * fi_acd_0 + g_xxxz_0_xyz_1[i] * wa_y[i];

        g_xxxyz_0_xzz_0[i] = g_xxxz_0_xzz_1[i] * wa_y[i];

        g_xxxyz_0_yyy_0[i] = g_xxxy_0_yyy_1[i] * wa_z[i];

        g_xxxyz_0_yyz_0[i] = 2.0 * g_xxxz_0_yz_1[i] * fi_acd_0 + g_xxxz_0_yyz_1[i] * wa_y[i];

        g_xxxyz_0_yzz_0[i] = g_xxxz_0_zz_1[i] * fi_acd_0 + g_xxxz_0_yzz_1[i] * wa_y[i];

        g_xxxyz_0_zzz_0[i] = g_xxxz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 50-60 components of targeted buffer : HSF

    auto g_xxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 50);

    auto g_xxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 51);

    auto g_xxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 52);

    auto g_xxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 53);

    auto g_xxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 54);

    auto g_xxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 55);

    auto g_xxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 56);

    auto g_xxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 57);

    auto g_xxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 58);

    auto g_xxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 59);

    #pragma omp simd aligned(g_xxx_0_xxx_0, g_xxx_0_xxx_1, g_xxx_0_xxy_0, g_xxx_0_xxy_1, g_xxx_0_xyy_0, g_xxx_0_xyy_1, g_xxxz_0_xxx_1, g_xxxz_0_xxy_1, g_xxxz_0_xyy_1, g_xxxzz_0_xxx_0, g_xxxzz_0_xxy_0, g_xxxzz_0_xxz_0, g_xxxzz_0_xyy_0, g_xxxzz_0_xyz_0, g_xxxzz_0_xzz_0, g_xxxzz_0_yyy_0, g_xxxzz_0_yyz_0, g_xxxzz_0_yzz_0, g_xxxzz_0_zzz_0, g_xxzz_0_xxz_1, g_xxzz_0_xyz_1, g_xxzz_0_xz_1, g_xxzz_0_xzz_1, g_xxzz_0_yyy_1, g_xxzz_0_yyz_1, g_xxzz_0_yz_1, g_xxzz_0_yzz_1, g_xxzz_0_zz_1, g_xxzz_0_zzz_1, g_xzz_0_xxz_0, g_xzz_0_xxz_1, g_xzz_0_xyz_0, g_xzz_0_xyz_1, g_xzz_0_xzz_0, g_xzz_0_xzz_1, g_xzz_0_yyy_0, g_xzz_0_yyy_1, g_xzz_0_yyz_0, g_xzz_0_yyz_1, g_xzz_0_yzz_0, g_xzz_0_yzz_1, g_xzz_0_zzz_0, g_xzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xxx_0[i] = g_xxx_0_xxx_0[i] * fbe_0 - g_xxx_0_xxx_1[i] * fz_be_0 + g_xxxz_0_xxx_1[i] * wa_z[i];

        g_xxxzz_0_xxy_0[i] = g_xxx_0_xxy_0[i] * fbe_0 - g_xxx_0_xxy_1[i] * fz_be_0 + g_xxxz_0_xxy_1[i] * wa_z[i];

        g_xxxzz_0_xxz_0[i] = 2.0 * g_xzz_0_xxz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xz_1[i] * fi_acd_0 + g_xxzz_0_xxz_1[i] * wa_x[i];

        g_xxxzz_0_xyy_0[i] = g_xxx_0_xyy_0[i] * fbe_0 - g_xxx_0_xyy_1[i] * fz_be_0 + g_xxxz_0_xyy_1[i] * wa_z[i];

        g_xxxzz_0_xyz_0[i] = 2.0 * g_xzz_0_xyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyz_1[i] * fz_be_0 + g_xxzz_0_yz_1[i] * fi_acd_0 + g_xxzz_0_xyz_1[i] * wa_x[i];

        g_xxxzz_0_xzz_0[i] = 2.0 * g_xzz_0_xzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xzz_1[i] * fz_be_0 + g_xxzz_0_zz_1[i] * fi_acd_0 + g_xxzz_0_xzz_1[i] * wa_x[i];

        g_xxxzz_0_yyy_0[i] = 2.0 * g_xzz_0_yyy_0[i] * fbe_0 - 2.0 * g_xzz_0_yyy_1[i] * fz_be_0 + g_xxzz_0_yyy_1[i] * wa_x[i];

        g_xxxzz_0_yyz_0[i] = 2.0 * g_xzz_0_yyz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyz_1[i] * fz_be_0 + g_xxzz_0_yyz_1[i] * wa_x[i];

        g_xxxzz_0_yzz_0[i] = 2.0 * g_xzz_0_yzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yzz_1[i] * fz_be_0 + g_xxzz_0_yzz_1[i] * wa_x[i];

        g_xxxzz_0_zzz_0[i] = 2.0 * g_xzz_0_zzz_0[i] * fbe_0 - 2.0 * g_xzz_0_zzz_1[i] * fz_be_0 + g_xxzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 60-70 components of targeted buffer : HSF

    auto g_xxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 60);

    auto g_xxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 61);

    auto g_xxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 62);

    auto g_xxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 63);

    auto g_xxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 64);

    auto g_xxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 65);

    auto g_xxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 66);

    auto g_xxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 67);

    auto g_xxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 68);

    auto g_xxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 69);

    #pragma omp simd aligned(g_xxy_0_xxx_0, g_xxy_0_xxx_1, g_xxy_0_xxz_0, g_xxy_0_xxz_1, g_xxy_0_xzz_0, g_xxy_0_xzz_1, g_xxyy_0_xxx_1, g_xxyy_0_xxz_1, g_xxyy_0_xzz_1, g_xxyyy_0_xxx_0, g_xxyyy_0_xxy_0, g_xxyyy_0_xxz_0, g_xxyyy_0_xyy_0, g_xxyyy_0_xyz_0, g_xxyyy_0_xzz_0, g_xxyyy_0_yyy_0, g_xxyyy_0_yyz_0, g_xxyyy_0_yzz_0, g_xxyyy_0_zzz_0, g_xyyy_0_xxy_1, g_xyyy_0_xy_1, g_xyyy_0_xyy_1, g_xyyy_0_xyz_1, g_xyyy_0_yy_1, g_xyyy_0_yyy_1, g_xyyy_0_yyz_1, g_xyyy_0_yz_1, g_xyyy_0_yzz_1, g_xyyy_0_zzz_1, g_yyy_0_xxy_0, g_yyy_0_xxy_1, g_yyy_0_xyy_0, g_yyy_0_xyy_1, g_yyy_0_xyz_0, g_yyy_0_xyz_1, g_yyy_0_yyy_0, g_yyy_0_yyy_1, g_yyy_0_yyz_0, g_yyy_0_yyz_1, g_yyy_0_yzz_0, g_yyy_0_yzz_1, g_yyy_0_zzz_0, g_yyy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xxx_0[i] = 2.0 * g_xxy_0_xxx_0[i] * fbe_0 - 2.0 * g_xxy_0_xxx_1[i] * fz_be_0 + g_xxyy_0_xxx_1[i] * wa_y[i];

        g_xxyyy_0_xxy_0[i] = g_yyy_0_xxy_0[i] * fbe_0 - g_yyy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xy_1[i] * fi_acd_0 + g_xyyy_0_xxy_1[i] * wa_x[i];

        g_xxyyy_0_xxz_0[i] = 2.0 * g_xxy_0_xxz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxz_1[i] * fz_be_0 + g_xxyy_0_xxz_1[i] * wa_y[i];

        g_xxyyy_0_xyy_0[i] = g_yyy_0_xyy_0[i] * fbe_0 - g_yyy_0_xyy_1[i] * fz_be_0 + g_xyyy_0_yy_1[i] * fi_acd_0 + g_xyyy_0_xyy_1[i] * wa_x[i];

        g_xxyyy_0_xyz_0[i] = g_yyy_0_xyz_0[i] * fbe_0 - g_yyy_0_xyz_1[i] * fz_be_0 + g_xyyy_0_yz_1[i] * fi_acd_0 + g_xyyy_0_xyz_1[i] * wa_x[i];

        g_xxyyy_0_xzz_0[i] = 2.0 * g_xxy_0_xzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xzz_1[i] * fz_be_0 + g_xxyy_0_xzz_1[i] * wa_y[i];

        g_xxyyy_0_yyy_0[i] = g_yyy_0_yyy_0[i] * fbe_0 - g_yyy_0_yyy_1[i] * fz_be_0 + g_xyyy_0_yyy_1[i] * wa_x[i];

        g_xxyyy_0_yyz_0[i] = g_yyy_0_yyz_0[i] * fbe_0 - g_yyy_0_yyz_1[i] * fz_be_0 + g_xyyy_0_yyz_1[i] * wa_x[i];

        g_xxyyy_0_yzz_0[i] = g_yyy_0_yzz_0[i] * fbe_0 - g_yyy_0_yzz_1[i] * fz_be_0 + g_xyyy_0_yzz_1[i] * wa_x[i];

        g_xxyyy_0_zzz_0[i] = g_yyy_0_zzz_0[i] * fbe_0 - g_yyy_0_zzz_1[i] * fz_be_0 + g_xyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 70-80 components of targeted buffer : HSF

    auto g_xxyyz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 70);

    auto g_xxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 71);

    auto g_xxyyz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 72);

    auto g_xxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 73);

    auto g_xxyyz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 74);

    auto g_xxyyz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 75);

    auto g_xxyyz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 76);

    auto g_xxyyz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 77);

    auto g_xxyyz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 78);

    auto g_xxyyz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 79);

    #pragma omp simd aligned(g_xxyy_0_xx_1, g_xxyy_0_xxx_1, g_xxyy_0_xxy_1, g_xxyy_0_xxz_1, g_xxyy_0_xy_1, g_xxyy_0_xyy_1, g_xxyy_0_xyz_1, g_xxyy_0_xz_1, g_xxyy_0_xzz_1, g_xxyy_0_yy_1, g_xxyy_0_yyy_1, g_xxyy_0_yyz_1, g_xxyy_0_yz_1, g_xxyy_0_yzz_1, g_xxyy_0_zz_1, g_xxyy_0_zzz_1, g_xxyyz_0_xxx_0, g_xxyyz_0_xxy_0, g_xxyyz_0_xxz_0, g_xxyyz_0_xyy_0, g_xxyyz_0_xyz_0, g_xxyyz_0_xzz_0, g_xxyyz_0_yyy_0, g_xxyyz_0_yyz_0, g_xxyyz_0_yzz_0, g_xxyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xxx_0[i] = g_xxyy_0_xxx_1[i] * wa_z[i];

        g_xxyyz_0_xxy_0[i] = g_xxyy_0_xxy_1[i] * wa_z[i];

        g_xxyyz_0_xxz_0[i] = g_xxyy_0_xx_1[i] * fi_acd_0 + g_xxyy_0_xxz_1[i] * wa_z[i];

        g_xxyyz_0_xyy_0[i] = g_xxyy_0_xyy_1[i] * wa_z[i];

        g_xxyyz_0_xyz_0[i] = g_xxyy_0_xy_1[i] * fi_acd_0 + g_xxyy_0_xyz_1[i] * wa_z[i];

        g_xxyyz_0_xzz_0[i] = 2.0 * g_xxyy_0_xz_1[i] * fi_acd_0 + g_xxyy_0_xzz_1[i] * wa_z[i];

        g_xxyyz_0_yyy_0[i] = g_xxyy_0_yyy_1[i] * wa_z[i];

        g_xxyyz_0_yyz_0[i] = g_xxyy_0_yy_1[i] * fi_acd_0 + g_xxyy_0_yyz_1[i] * wa_z[i];

        g_xxyyz_0_yzz_0[i] = 2.0 * g_xxyy_0_yz_1[i] * fi_acd_0 + g_xxyy_0_yzz_1[i] * wa_z[i];

        g_xxyyz_0_zzz_0[i] = 3.0 * g_xxyy_0_zz_1[i] * fi_acd_0 + g_xxyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 80-90 components of targeted buffer : HSF

    auto g_xxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 80);

    auto g_xxyzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 81);

    auto g_xxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 82);

    auto g_xxyzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 83);

    auto g_xxyzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 84);

    auto g_xxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 85);

    auto g_xxyzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 86);

    auto g_xxyzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 87);

    auto g_xxyzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 88);

    auto g_xxyzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 89);

    #pragma omp simd aligned(g_xxyzz_0_xxx_0, g_xxyzz_0_xxy_0, g_xxyzz_0_xxz_0, g_xxyzz_0_xyy_0, g_xxyzz_0_xyz_0, g_xxyzz_0_xzz_0, g_xxyzz_0_yyy_0, g_xxyzz_0_yyz_0, g_xxyzz_0_yzz_0, g_xxyzz_0_zzz_0, g_xxzz_0_xx_1, g_xxzz_0_xxx_1, g_xxzz_0_xxy_1, g_xxzz_0_xxz_1, g_xxzz_0_xy_1, g_xxzz_0_xyy_1, g_xxzz_0_xyz_1, g_xxzz_0_xz_1, g_xxzz_0_xzz_1, g_xxzz_0_yy_1, g_xxzz_0_yyy_1, g_xxzz_0_yyz_1, g_xxzz_0_yz_1, g_xxzz_0_yzz_1, g_xxzz_0_zz_1, g_xxzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xxx_0[i] = g_xxzz_0_xxx_1[i] * wa_y[i];

        g_xxyzz_0_xxy_0[i] = g_xxzz_0_xx_1[i] * fi_acd_0 + g_xxzz_0_xxy_1[i] * wa_y[i];

        g_xxyzz_0_xxz_0[i] = g_xxzz_0_xxz_1[i] * wa_y[i];

        g_xxyzz_0_xyy_0[i] = 2.0 * g_xxzz_0_xy_1[i] * fi_acd_0 + g_xxzz_0_xyy_1[i] * wa_y[i];

        g_xxyzz_0_xyz_0[i] = g_xxzz_0_xz_1[i] * fi_acd_0 + g_xxzz_0_xyz_1[i] * wa_y[i];

        g_xxyzz_0_xzz_0[i] = g_xxzz_0_xzz_1[i] * wa_y[i];

        g_xxyzz_0_yyy_0[i] = 3.0 * g_xxzz_0_yy_1[i] * fi_acd_0 + g_xxzz_0_yyy_1[i] * wa_y[i];

        g_xxyzz_0_yyz_0[i] = 2.0 * g_xxzz_0_yz_1[i] * fi_acd_0 + g_xxzz_0_yyz_1[i] * wa_y[i];

        g_xxyzz_0_yzz_0[i] = g_xxzz_0_zz_1[i] * fi_acd_0 + g_xxzz_0_yzz_1[i] * wa_y[i];

        g_xxyzz_0_zzz_0[i] = g_xxzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 90-100 components of targeted buffer : HSF

    auto g_xxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 90);

    auto g_xxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 91);

    auto g_xxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 92);

    auto g_xxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 93);

    auto g_xxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 94);

    auto g_xxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 95);

    auto g_xxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 96);

    auto g_xxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 97);

    auto g_xxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 98);

    auto g_xxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 99);

    #pragma omp simd aligned(g_xxz_0_xxx_0, g_xxz_0_xxx_1, g_xxz_0_xxy_0, g_xxz_0_xxy_1, g_xxz_0_xyy_0, g_xxz_0_xyy_1, g_xxzz_0_xxx_1, g_xxzz_0_xxy_1, g_xxzz_0_xyy_1, g_xxzzz_0_xxx_0, g_xxzzz_0_xxy_0, g_xxzzz_0_xxz_0, g_xxzzz_0_xyy_0, g_xxzzz_0_xyz_0, g_xxzzz_0_xzz_0, g_xxzzz_0_yyy_0, g_xxzzz_0_yyz_0, g_xxzzz_0_yzz_0, g_xxzzz_0_zzz_0, g_xzzz_0_xxz_1, g_xzzz_0_xyz_1, g_xzzz_0_xz_1, g_xzzz_0_xzz_1, g_xzzz_0_yyy_1, g_xzzz_0_yyz_1, g_xzzz_0_yz_1, g_xzzz_0_yzz_1, g_xzzz_0_zz_1, g_xzzz_0_zzz_1, g_zzz_0_xxz_0, g_zzz_0_xxz_1, g_zzz_0_xyz_0, g_zzz_0_xyz_1, g_zzz_0_xzz_0, g_zzz_0_xzz_1, g_zzz_0_yyy_0, g_zzz_0_yyy_1, g_zzz_0_yyz_0, g_zzz_0_yyz_1, g_zzz_0_yzz_0, g_zzz_0_yzz_1, g_zzz_0_zzz_0, g_zzz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xxx_0[i] = 2.0 * g_xxz_0_xxx_0[i] * fbe_0 - 2.0 * g_xxz_0_xxx_1[i] * fz_be_0 + g_xxzz_0_xxx_1[i] * wa_z[i];

        g_xxzzz_0_xxy_0[i] = 2.0 * g_xxz_0_xxy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxy_1[i] * fz_be_0 + g_xxzz_0_xxy_1[i] * wa_z[i];

        g_xxzzz_0_xxz_0[i] = g_zzz_0_xxz_0[i] * fbe_0 - g_zzz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xz_1[i] * fi_acd_0 + g_xzzz_0_xxz_1[i] * wa_x[i];

        g_xxzzz_0_xyy_0[i] = 2.0 * g_xxz_0_xyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xyy_1[i] * fz_be_0 + g_xxzz_0_xyy_1[i] * wa_z[i];

        g_xxzzz_0_xyz_0[i] = g_zzz_0_xyz_0[i] * fbe_0 - g_zzz_0_xyz_1[i] * fz_be_0 + g_xzzz_0_yz_1[i] * fi_acd_0 + g_xzzz_0_xyz_1[i] * wa_x[i];

        g_xxzzz_0_xzz_0[i] = g_zzz_0_xzz_0[i] * fbe_0 - g_zzz_0_xzz_1[i] * fz_be_0 + g_xzzz_0_zz_1[i] * fi_acd_0 + g_xzzz_0_xzz_1[i] * wa_x[i];

        g_xxzzz_0_yyy_0[i] = g_zzz_0_yyy_0[i] * fbe_0 - g_zzz_0_yyy_1[i] * fz_be_0 + g_xzzz_0_yyy_1[i] * wa_x[i];

        g_xxzzz_0_yyz_0[i] = g_zzz_0_yyz_0[i] * fbe_0 - g_zzz_0_yyz_1[i] * fz_be_0 + g_xzzz_0_yyz_1[i] * wa_x[i];

        g_xxzzz_0_yzz_0[i] = g_zzz_0_yzz_0[i] * fbe_0 - g_zzz_0_yzz_1[i] * fz_be_0 + g_xzzz_0_yzz_1[i] * wa_x[i];

        g_xxzzz_0_zzz_0[i] = g_zzz_0_zzz_0[i] * fbe_0 - g_zzz_0_zzz_1[i] * fz_be_0 + g_xzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 100-110 components of targeted buffer : HSF

    auto g_xyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 100);

    auto g_xyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 101);

    auto g_xyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 102);

    auto g_xyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 103);

    auto g_xyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 104);

    auto g_xyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 105);

    auto g_xyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 106);

    auto g_xyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 107);

    auto g_xyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 108);

    auto g_xyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 109);

    #pragma omp simd aligned(g_xyyyy_0_xxx_0, g_xyyyy_0_xxy_0, g_xyyyy_0_xxz_0, g_xyyyy_0_xyy_0, g_xyyyy_0_xyz_0, g_xyyyy_0_xzz_0, g_xyyyy_0_yyy_0, g_xyyyy_0_yyz_0, g_xyyyy_0_yzz_0, g_xyyyy_0_zzz_0, g_yyyy_0_xx_1, g_yyyy_0_xxx_1, g_yyyy_0_xxy_1, g_yyyy_0_xxz_1, g_yyyy_0_xy_1, g_yyyy_0_xyy_1, g_yyyy_0_xyz_1, g_yyyy_0_xz_1, g_yyyy_0_xzz_1, g_yyyy_0_yy_1, g_yyyy_0_yyy_1, g_yyyy_0_yyz_1, g_yyyy_0_yz_1, g_yyyy_0_yzz_1, g_yyyy_0_zz_1, g_yyyy_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xxx_0[i] = 3.0 * g_yyyy_0_xx_1[i] * fi_acd_0 + g_yyyy_0_xxx_1[i] * wa_x[i];

        g_xyyyy_0_xxy_0[i] = 2.0 * g_yyyy_0_xy_1[i] * fi_acd_0 + g_yyyy_0_xxy_1[i] * wa_x[i];

        g_xyyyy_0_xxz_0[i] = 2.0 * g_yyyy_0_xz_1[i] * fi_acd_0 + g_yyyy_0_xxz_1[i] * wa_x[i];

        g_xyyyy_0_xyy_0[i] = g_yyyy_0_yy_1[i] * fi_acd_0 + g_yyyy_0_xyy_1[i] * wa_x[i];

        g_xyyyy_0_xyz_0[i] = g_yyyy_0_yz_1[i] * fi_acd_0 + g_yyyy_0_xyz_1[i] * wa_x[i];

        g_xyyyy_0_xzz_0[i] = g_yyyy_0_zz_1[i] * fi_acd_0 + g_yyyy_0_xzz_1[i] * wa_x[i];

        g_xyyyy_0_yyy_0[i] = g_yyyy_0_yyy_1[i] * wa_x[i];

        g_xyyyy_0_yyz_0[i] = g_yyyy_0_yyz_1[i] * wa_x[i];

        g_xyyyy_0_yzz_0[i] = g_yyyy_0_yzz_1[i] * wa_x[i];

        g_xyyyy_0_zzz_0[i] = g_yyyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 110-120 components of targeted buffer : HSF

    auto g_xyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 110);

    auto g_xyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 111);

    auto g_xyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 112);

    auto g_xyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 113);

    auto g_xyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 114);

    auto g_xyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 115);

    auto g_xyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 116);

    auto g_xyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 117);

    auto g_xyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 118);

    auto g_xyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 119);

    #pragma omp simd aligned(g_xyyy_0_xxx_1, g_xyyy_0_xxy_1, g_xyyy_0_xyy_1, g_xyyyz_0_xxx_0, g_xyyyz_0_xxy_0, g_xyyyz_0_xxz_0, g_xyyyz_0_xyy_0, g_xyyyz_0_xyz_0, g_xyyyz_0_xzz_0, g_xyyyz_0_yyy_0, g_xyyyz_0_yyz_0, g_xyyyz_0_yzz_0, g_xyyyz_0_zzz_0, g_yyyz_0_xxz_1, g_yyyz_0_xyz_1, g_yyyz_0_xz_1, g_yyyz_0_xzz_1, g_yyyz_0_yyy_1, g_yyyz_0_yyz_1, g_yyyz_0_yz_1, g_yyyz_0_yzz_1, g_yyyz_0_zz_1, g_yyyz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xxx_0[i] = g_xyyy_0_xxx_1[i] * wa_z[i];

        g_xyyyz_0_xxy_0[i] = g_xyyy_0_xxy_1[i] * wa_z[i];

        g_xyyyz_0_xxz_0[i] = 2.0 * g_yyyz_0_xz_1[i] * fi_acd_0 + g_yyyz_0_xxz_1[i] * wa_x[i];

        g_xyyyz_0_xyy_0[i] = g_xyyy_0_xyy_1[i] * wa_z[i];

        g_xyyyz_0_xyz_0[i] = g_yyyz_0_yz_1[i] * fi_acd_0 + g_yyyz_0_xyz_1[i] * wa_x[i];

        g_xyyyz_0_xzz_0[i] = g_yyyz_0_zz_1[i] * fi_acd_0 + g_yyyz_0_xzz_1[i] * wa_x[i];

        g_xyyyz_0_yyy_0[i] = g_yyyz_0_yyy_1[i] * wa_x[i];

        g_xyyyz_0_yyz_0[i] = g_yyyz_0_yyz_1[i] * wa_x[i];

        g_xyyyz_0_yzz_0[i] = g_yyyz_0_yzz_1[i] * wa_x[i];

        g_xyyyz_0_zzz_0[i] = g_yyyz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 120-130 components of targeted buffer : HSF

    auto g_xyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 120);

    auto g_xyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 121);

    auto g_xyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 122);

    auto g_xyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 123);

    auto g_xyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 124);

    auto g_xyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 125);

    auto g_xyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 126);

    auto g_xyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 127);

    auto g_xyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 128);

    auto g_xyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 129);

    #pragma omp simd aligned(g_xyyzz_0_xxx_0, g_xyyzz_0_xxy_0, g_xyyzz_0_xxz_0, g_xyyzz_0_xyy_0, g_xyyzz_0_xyz_0, g_xyyzz_0_xzz_0, g_xyyzz_0_yyy_0, g_xyyzz_0_yyz_0, g_xyyzz_0_yzz_0, g_xyyzz_0_zzz_0, g_yyzz_0_xx_1, g_yyzz_0_xxx_1, g_yyzz_0_xxy_1, g_yyzz_0_xxz_1, g_yyzz_0_xy_1, g_yyzz_0_xyy_1, g_yyzz_0_xyz_1, g_yyzz_0_xz_1, g_yyzz_0_xzz_1, g_yyzz_0_yy_1, g_yyzz_0_yyy_1, g_yyzz_0_yyz_1, g_yyzz_0_yz_1, g_yyzz_0_yzz_1, g_yyzz_0_zz_1, g_yyzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xxx_0[i] = 3.0 * g_yyzz_0_xx_1[i] * fi_acd_0 + g_yyzz_0_xxx_1[i] * wa_x[i];

        g_xyyzz_0_xxy_0[i] = 2.0 * g_yyzz_0_xy_1[i] * fi_acd_0 + g_yyzz_0_xxy_1[i] * wa_x[i];

        g_xyyzz_0_xxz_0[i] = 2.0 * g_yyzz_0_xz_1[i] * fi_acd_0 + g_yyzz_0_xxz_1[i] * wa_x[i];

        g_xyyzz_0_xyy_0[i] = g_yyzz_0_yy_1[i] * fi_acd_0 + g_yyzz_0_xyy_1[i] * wa_x[i];

        g_xyyzz_0_xyz_0[i] = g_yyzz_0_yz_1[i] * fi_acd_0 + g_yyzz_0_xyz_1[i] * wa_x[i];

        g_xyyzz_0_xzz_0[i] = g_yyzz_0_zz_1[i] * fi_acd_0 + g_yyzz_0_xzz_1[i] * wa_x[i];

        g_xyyzz_0_yyy_0[i] = g_yyzz_0_yyy_1[i] * wa_x[i];

        g_xyyzz_0_yyz_0[i] = g_yyzz_0_yyz_1[i] * wa_x[i];

        g_xyyzz_0_yzz_0[i] = g_yyzz_0_yzz_1[i] * wa_x[i];

        g_xyyzz_0_zzz_0[i] = g_yyzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 130-140 components of targeted buffer : HSF

    auto g_xyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 130);

    auto g_xyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 131);

    auto g_xyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 132);

    auto g_xyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 133);

    auto g_xyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 134);

    auto g_xyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 135);

    auto g_xyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 136);

    auto g_xyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 137);

    auto g_xyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 138);

    auto g_xyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 139);

    #pragma omp simd aligned(g_xyzzz_0_xxx_0, g_xyzzz_0_xxy_0, g_xyzzz_0_xxz_0, g_xyzzz_0_xyy_0, g_xyzzz_0_xyz_0, g_xyzzz_0_xzz_0, g_xyzzz_0_yyy_0, g_xyzzz_0_yyz_0, g_xyzzz_0_yzz_0, g_xyzzz_0_zzz_0, g_xzzz_0_xxx_1, g_xzzz_0_xxz_1, g_xzzz_0_xzz_1, g_yzzz_0_xxy_1, g_yzzz_0_xy_1, g_yzzz_0_xyy_1, g_yzzz_0_xyz_1, g_yzzz_0_yy_1, g_yzzz_0_yyy_1, g_yzzz_0_yyz_1, g_yzzz_0_yz_1, g_yzzz_0_yzz_1, g_yzzz_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xxx_0[i] = g_xzzz_0_xxx_1[i] * wa_y[i];

        g_xyzzz_0_xxy_0[i] = 2.0 * g_yzzz_0_xy_1[i] * fi_acd_0 + g_yzzz_0_xxy_1[i] * wa_x[i];

        g_xyzzz_0_xxz_0[i] = g_xzzz_0_xxz_1[i] * wa_y[i];

        g_xyzzz_0_xyy_0[i] = g_yzzz_0_yy_1[i] * fi_acd_0 + g_yzzz_0_xyy_1[i] * wa_x[i];

        g_xyzzz_0_xyz_0[i] = g_yzzz_0_yz_1[i] * fi_acd_0 + g_yzzz_0_xyz_1[i] * wa_x[i];

        g_xyzzz_0_xzz_0[i] = g_xzzz_0_xzz_1[i] * wa_y[i];

        g_xyzzz_0_yyy_0[i] = g_yzzz_0_yyy_1[i] * wa_x[i];

        g_xyzzz_0_yyz_0[i] = g_yzzz_0_yyz_1[i] * wa_x[i];

        g_xyzzz_0_yzz_0[i] = g_yzzz_0_yzz_1[i] * wa_x[i];

        g_xyzzz_0_zzz_0[i] = g_yzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 140-150 components of targeted buffer : HSF

    auto g_xzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 140);

    auto g_xzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 141);

    auto g_xzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 142);

    auto g_xzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 143);

    auto g_xzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 144);

    auto g_xzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 145);

    auto g_xzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 146);

    auto g_xzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 147);

    auto g_xzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 148);

    auto g_xzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 149);

    #pragma omp simd aligned(g_xzzzz_0_xxx_0, g_xzzzz_0_xxy_0, g_xzzzz_0_xxz_0, g_xzzzz_0_xyy_0, g_xzzzz_0_xyz_0, g_xzzzz_0_xzz_0, g_xzzzz_0_yyy_0, g_xzzzz_0_yyz_0, g_xzzzz_0_yzz_0, g_xzzzz_0_zzz_0, g_zzzz_0_xx_1, g_zzzz_0_xxx_1, g_zzzz_0_xxy_1, g_zzzz_0_xxz_1, g_zzzz_0_xy_1, g_zzzz_0_xyy_1, g_zzzz_0_xyz_1, g_zzzz_0_xz_1, g_zzzz_0_xzz_1, g_zzzz_0_yy_1, g_zzzz_0_yyy_1, g_zzzz_0_yyz_1, g_zzzz_0_yz_1, g_zzzz_0_yzz_1, g_zzzz_0_zz_1, g_zzzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xxx_0[i] = 3.0 * g_zzzz_0_xx_1[i] * fi_acd_0 + g_zzzz_0_xxx_1[i] * wa_x[i];

        g_xzzzz_0_xxy_0[i] = 2.0 * g_zzzz_0_xy_1[i] * fi_acd_0 + g_zzzz_0_xxy_1[i] * wa_x[i];

        g_xzzzz_0_xxz_0[i] = 2.0 * g_zzzz_0_xz_1[i] * fi_acd_0 + g_zzzz_0_xxz_1[i] * wa_x[i];

        g_xzzzz_0_xyy_0[i] = g_zzzz_0_yy_1[i] * fi_acd_0 + g_zzzz_0_xyy_1[i] * wa_x[i];

        g_xzzzz_0_xyz_0[i] = g_zzzz_0_yz_1[i] * fi_acd_0 + g_zzzz_0_xyz_1[i] * wa_x[i];

        g_xzzzz_0_xzz_0[i] = g_zzzz_0_zz_1[i] * fi_acd_0 + g_zzzz_0_xzz_1[i] * wa_x[i];

        g_xzzzz_0_yyy_0[i] = g_zzzz_0_yyy_1[i] * wa_x[i];

        g_xzzzz_0_yyz_0[i] = g_zzzz_0_yyz_1[i] * wa_x[i];

        g_xzzzz_0_yzz_0[i] = g_zzzz_0_yzz_1[i] * wa_x[i];

        g_xzzzz_0_zzz_0[i] = g_zzzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 150-160 components of targeted buffer : HSF

    auto g_yyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 150);

    auto g_yyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 151);

    auto g_yyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 152);

    auto g_yyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 153);

    auto g_yyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 154);

    auto g_yyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 155);

    auto g_yyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 156);

    auto g_yyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 157);

    auto g_yyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 158);

    auto g_yyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 159);

    #pragma omp simd aligned(g_yyy_0_xxx_0, g_yyy_0_xxx_1, g_yyy_0_xxy_0, g_yyy_0_xxy_1, g_yyy_0_xxz_0, g_yyy_0_xxz_1, g_yyy_0_xyy_0, g_yyy_0_xyy_1, g_yyy_0_xyz_0, g_yyy_0_xyz_1, g_yyy_0_xzz_0, g_yyy_0_xzz_1, g_yyy_0_yyy_0, g_yyy_0_yyy_1, g_yyy_0_yyz_0, g_yyy_0_yyz_1, g_yyy_0_yzz_0, g_yyy_0_yzz_1, g_yyy_0_zzz_0, g_yyy_0_zzz_1, g_yyyy_0_xx_1, g_yyyy_0_xxx_1, g_yyyy_0_xxy_1, g_yyyy_0_xxz_1, g_yyyy_0_xy_1, g_yyyy_0_xyy_1, g_yyyy_0_xyz_1, g_yyyy_0_xz_1, g_yyyy_0_xzz_1, g_yyyy_0_yy_1, g_yyyy_0_yyy_1, g_yyyy_0_yyz_1, g_yyyy_0_yz_1, g_yyyy_0_yzz_1, g_yyyy_0_zz_1, g_yyyy_0_zzz_1, g_yyyyy_0_xxx_0, g_yyyyy_0_xxy_0, g_yyyyy_0_xxz_0, g_yyyyy_0_xyy_0, g_yyyyy_0_xyz_0, g_yyyyy_0_xzz_0, g_yyyyy_0_yyy_0, g_yyyyy_0_yyz_0, g_yyyyy_0_yzz_0, g_yyyyy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xxx_0[i] = 4.0 * g_yyy_0_xxx_0[i] * fbe_0 - 4.0 * g_yyy_0_xxx_1[i] * fz_be_0 + g_yyyy_0_xxx_1[i] * wa_y[i];

        g_yyyyy_0_xxy_0[i] = 4.0 * g_yyy_0_xxy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxy_1[i] * fz_be_0 + g_yyyy_0_xx_1[i] * fi_acd_0 + g_yyyy_0_xxy_1[i] * wa_y[i];

        g_yyyyy_0_xxz_0[i] = 4.0 * g_yyy_0_xxz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxz_1[i] * fz_be_0 + g_yyyy_0_xxz_1[i] * wa_y[i];

        g_yyyyy_0_xyy_0[i] = 4.0 * g_yyy_0_xyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xyy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xy_1[i] * fi_acd_0 + g_yyyy_0_xyy_1[i] * wa_y[i];

        g_yyyyy_0_xyz_0[i] = 4.0 * g_yyy_0_xyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyz_1[i] * fz_be_0 + g_yyyy_0_xz_1[i] * fi_acd_0 + g_yyyy_0_xyz_1[i] * wa_y[i];

        g_yyyyy_0_xzz_0[i] = 4.0 * g_yyy_0_xzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xzz_1[i] * fz_be_0 + g_yyyy_0_xzz_1[i] * wa_y[i];

        g_yyyyy_0_yyy_0[i] = 4.0 * g_yyy_0_yyy_0[i] * fbe_0 - 4.0 * g_yyy_0_yyy_1[i] * fz_be_0 + 3.0 * g_yyyy_0_yy_1[i] * fi_acd_0 + g_yyyy_0_yyy_1[i] * wa_y[i];

        g_yyyyy_0_yyz_0[i] = 4.0 * g_yyy_0_yyz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_yz_1[i] * fi_acd_0 + g_yyyy_0_yyz_1[i] * wa_y[i];

        g_yyyyy_0_yzz_0[i] = 4.0 * g_yyy_0_yzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yzz_1[i] * fz_be_0 + g_yyyy_0_zz_1[i] * fi_acd_0 + g_yyyy_0_yzz_1[i] * wa_y[i];

        g_yyyyy_0_zzz_0[i] = 4.0 * g_yyy_0_zzz_0[i] * fbe_0 - 4.0 * g_yyy_0_zzz_1[i] * fz_be_0 + g_yyyy_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 160-170 components of targeted buffer : HSF

    auto g_yyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 160);

    auto g_yyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 161);

    auto g_yyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 162);

    auto g_yyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 163);

    auto g_yyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 164);

    auto g_yyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 165);

    auto g_yyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 166);

    auto g_yyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 167);

    auto g_yyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 168);

    auto g_yyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 169);

    #pragma omp simd aligned(g_yyyy_0_xx_1, g_yyyy_0_xxx_1, g_yyyy_0_xxy_1, g_yyyy_0_xxz_1, g_yyyy_0_xy_1, g_yyyy_0_xyy_1, g_yyyy_0_xyz_1, g_yyyy_0_xz_1, g_yyyy_0_xzz_1, g_yyyy_0_yy_1, g_yyyy_0_yyy_1, g_yyyy_0_yyz_1, g_yyyy_0_yz_1, g_yyyy_0_yzz_1, g_yyyy_0_zz_1, g_yyyy_0_zzz_1, g_yyyyz_0_xxx_0, g_yyyyz_0_xxy_0, g_yyyyz_0_xxz_0, g_yyyyz_0_xyy_0, g_yyyyz_0_xyz_0, g_yyyyz_0_xzz_0, g_yyyyz_0_yyy_0, g_yyyyz_0_yyz_0, g_yyyyz_0_yzz_0, g_yyyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xxx_0[i] = g_yyyy_0_xxx_1[i] * wa_z[i];

        g_yyyyz_0_xxy_0[i] = g_yyyy_0_xxy_1[i] * wa_z[i];

        g_yyyyz_0_xxz_0[i] = g_yyyy_0_xx_1[i] * fi_acd_0 + g_yyyy_0_xxz_1[i] * wa_z[i];

        g_yyyyz_0_xyy_0[i] = g_yyyy_0_xyy_1[i] * wa_z[i];

        g_yyyyz_0_xyz_0[i] = g_yyyy_0_xy_1[i] * fi_acd_0 + g_yyyy_0_xyz_1[i] * wa_z[i];

        g_yyyyz_0_xzz_0[i] = 2.0 * g_yyyy_0_xz_1[i] * fi_acd_0 + g_yyyy_0_xzz_1[i] * wa_z[i];

        g_yyyyz_0_yyy_0[i] = g_yyyy_0_yyy_1[i] * wa_z[i];

        g_yyyyz_0_yyz_0[i] = g_yyyy_0_yy_1[i] * fi_acd_0 + g_yyyy_0_yyz_1[i] * wa_z[i];

        g_yyyyz_0_yzz_0[i] = 2.0 * g_yyyy_0_yz_1[i] * fi_acd_0 + g_yyyy_0_yzz_1[i] * wa_z[i];

        g_yyyyz_0_zzz_0[i] = 3.0 * g_yyyy_0_zz_1[i] * fi_acd_0 + g_yyyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 170-180 components of targeted buffer : HSF

    auto g_yyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 170);

    auto g_yyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 171);

    auto g_yyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 172);

    auto g_yyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 173);

    auto g_yyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 174);

    auto g_yyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 175);

    auto g_yyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 176);

    auto g_yyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 177);

    auto g_yyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 178);

    auto g_yyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 179);

    #pragma omp simd aligned(g_yyy_0_xxy_0, g_yyy_0_xxy_1, g_yyy_0_xyy_0, g_yyy_0_xyy_1, g_yyy_0_yyy_0, g_yyy_0_yyy_1, g_yyyz_0_xxy_1, g_yyyz_0_xyy_1, g_yyyz_0_yyy_1, g_yyyzz_0_xxx_0, g_yyyzz_0_xxy_0, g_yyyzz_0_xxz_0, g_yyyzz_0_xyy_0, g_yyyzz_0_xyz_0, g_yyyzz_0_xzz_0, g_yyyzz_0_yyy_0, g_yyyzz_0_yyz_0, g_yyyzz_0_yzz_0, g_yyyzz_0_zzz_0, g_yyzz_0_xxx_1, g_yyzz_0_xxz_1, g_yyzz_0_xyz_1, g_yyzz_0_xz_1, g_yyzz_0_xzz_1, g_yyzz_0_yyz_1, g_yyzz_0_yz_1, g_yyzz_0_yzz_1, g_yyzz_0_zz_1, g_yyzz_0_zzz_1, g_yzz_0_xxx_0, g_yzz_0_xxx_1, g_yzz_0_xxz_0, g_yzz_0_xxz_1, g_yzz_0_xyz_0, g_yzz_0_xyz_1, g_yzz_0_xzz_0, g_yzz_0_xzz_1, g_yzz_0_yyz_0, g_yzz_0_yyz_1, g_yzz_0_yzz_0, g_yzz_0_yzz_1, g_yzz_0_zzz_0, g_yzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xxx_0[i] = 2.0 * g_yzz_0_xxx_0[i] * fbe_0 - 2.0 * g_yzz_0_xxx_1[i] * fz_be_0 + g_yyzz_0_xxx_1[i] * wa_y[i];

        g_yyyzz_0_xxy_0[i] = g_yyy_0_xxy_0[i] * fbe_0 - g_yyy_0_xxy_1[i] * fz_be_0 + g_yyyz_0_xxy_1[i] * wa_z[i];

        g_yyyzz_0_xxz_0[i] = 2.0 * g_yzz_0_xxz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxz_1[i] * fz_be_0 + g_yyzz_0_xxz_1[i] * wa_y[i];

        g_yyyzz_0_xyy_0[i] = g_yyy_0_xyy_0[i] * fbe_0 - g_yyy_0_xyy_1[i] * fz_be_0 + g_yyyz_0_xyy_1[i] * wa_z[i];

        g_yyyzz_0_xyz_0[i] = 2.0 * g_yzz_0_xyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyz_1[i] * fz_be_0 + g_yyzz_0_xz_1[i] * fi_acd_0 + g_yyzz_0_xyz_1[i] * wa_y[i];

        g_yyyzz_0_xzz_0[i] = 2.0 * g_yzz_0_xzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xzz_1[i] * fz_be_0 + g_yyzz_0_xzz_1[i] * wa_y[i];

        g_yyyzz_0_yyy_0[i] = g_yyy_0_yyy_0[i] * fbe_0 - g_yyy_0_yyy_1[i] * fz_be_0 + g_yyyz_0_yyy_1[i] * wa_z[i];

        g_yyyzz_0_yyz_0[i] = 2.0 * g_yzz_0_yyz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_yz_1[i] * fi_acd_0 + g_yyzz_0_yyz_1[i] * wa_y[i];

        g_yyyzz_0_yzz_0[i] = 2.0 * g_yzz_0_yzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yzz_1[i] * fz_be_0 + g_yyzz_0_zz_1[i] * fi_acd_0 + g_yyzz_0_yzz_1[i] * wa_y[i];

        g_yyyzz_0_zzz_0[i] = 2.0 * g_yzz_0_zzz_0[i] * fbe_0 - 2.0 * g_yzz_0_zzz_1[i] * fz_be_0 + g_yyzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 180-190 components of targeted buffer : HSF

    auto g_yyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 180);

    auto g_yyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 181);

    auto g_yyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 182);

    auto g_yyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 183);

    auto g_yyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 184);

    auto g_yyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 185);

    auto g_yyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 186);

    auto g_yyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 187);

    auto g_yyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 188);

    auto g_yyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 189);

    #pragma omp simd aligned(g_yyz_0_xxy_0, g_yyz_0_xxy_1, g_yyz_0_xyy_0, g_yyz_0_xyy_1, g_yyz_0_yyy_0, g_yyz_0_yyy_1, g_yyzz_0_xxy_1, g_yyzz_0_xyy_1, g_yyzz_0_yyy_1, g_yyzzz_0_xxx_0, g_yyzzz_0_xxy_0, g_yyzzz_0_xxz_0, g_yyzzz_0_xyy_0, g_yyzzz_0_xyz_0, g_yyzzz_0_xzz_0, g_yyzzz_0_yyy_0, g_yyzzz_0_yyz_0, g_yyzzz_0_yzz_0, g_yyzzz_0_zzz_0, g_yzzz_0_xxx_1, g_yzzz_0_xxz_1, g_yzzz_0_xyz_1, g_yzzz_0_xz_1, g_yzzz_0_xzz_1, g_yzzz_0_yyz_1, g_yzzz_0_yz_1, g_yzzz_0_yzz_1, g_yzzz_0_zz_1, g_yzzz_0_zzz_1, g_zzz_0_xxx_0, g_zzz_0_xxx_1, g_zzz_0_xxz_0, g_zzz_0_xxz_1, g_zzz_0_xyz_0, g_zzz_0_xyz_1, g_zzz_0_xzz_0, g_zzz_0_xzz_1, g_zzz_0_yyz_0, g_zzz_0_yyz_1, g_zzz_0_yzz_0, g_zzz_0_yzz_1, g_zzz_0_zzz_0, g_zzz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xxx_0[i] = g_zzz_0_xxx_0[i] * fbe_0 - g_zzz_0_xxx_1[i] * fz_be_0 + g_yzzz_0_xxx_1[i] * wa_y[i];

        g_yyzzz_0_xxy_0[i] = 2.0 * g_yyz_0_xxy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxy_1[i] * fz_be_0 + g_yyzz_0_xxy_1[i] * wa_z[i];

        g_yyzzz_0_xxz_0[i] = g_zzz_0_xxz_0[i] * fbe_0 - g_zzz_0_xxz_1[i] * fz_be_0 + g_yzzz_0_xxz_1[i] * wa_y[i];

        g_yyzzz_0_xyy_0[i] = 2.0 * g_yyz_0_xyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xyy_1[i] * fz_be_0 + g_yyzz_0_xyy_1[i] * wa_z[i];

        g_yyzzz_0_xyz_0[i] = g_zzz_0_xyz_0[i] * fbe_0 - g_zzz_0_xyz_1[i] * fz_be_0 + g_yzzz_0_xz_1[i] * fi_acd_0 + g_yzzz_0_xyz_1[i] * wa_y[i];

        g_yyzzz_0_xzz_0[i] = g_zzz_0_xzz_0[i] * fbe_0 - g_zzz_0_xzz_1[i] * fz_be_0 + g_yzzz_0_xzz_1[i] * wa_y[i];

        g_yyzzz_0_yyy_0[i] = 2.0 * g_yyz_0_yyy_0[i] * fbe_0 - 2.0 * g_yyz_0_yyy_1[i] * fz_be_0 + g_yyzz_0_yyy_1[i] * wa_z[i];

        g_yyzzz_0_yyz_0[i] = g_zzz_0_yyz_0[i] * fbe_0 - g_zzz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_yz_1[i] * fi_acd_0 + g_yzzz_0_yyz_1[i] * wa_y[i];

        g_yyzzz_0_yzz_0[i] = g_zzz_0_yzz_0[i] * fbe_0 - g_zzz_0_yzz_1[i] * fz_be_0 + g_yzzz_0_zz_1[i] * fi_acd_0 + g_yzzz_0_yzz_1[i] * wa_y[i];

        g_yyzzz_0_zzz_0[i] = g_zzz_0_zzz_0[i] * fbe_0 - g_zzz_0_zzz_1[i] * fz_be_0 + g_yzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 190-200 components of targeted buffer : HSF

    auto g_yzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 190);

    auto g_yzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 191);

    auto g_yzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 192);

    auto g_yzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 193);

    auto g_yzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 194);

    auto g_yzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 195);

    auto g_yzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 196);

    auto g_yzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 197);

    auto g_yzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 198);

    auto g_yzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 199);

    #pragma omp simd aligned(g_yzzzz_0_xxx_0, g_yzzzz_0_xxy_0, g_yzzzz_0_xxz_0, g_yzzzz_0_xyy_0, g_yzzzz_0_xyz_0, g_yzzzz_0_xzz_0, g_yzzzz_0_yyy_0, g_yzzzz_0_yyz_0, g_yzzzz_0_yzz_0, g_yzzzz_0_zzz_0, g_zzzz_0_xx_1, g_zzzz_0_xxx_1, g_zzzz_0_xxy_1, g_zzzz_0_xxz_1, g_zzzz_0_xy_1, g_zzzz_0_xyy_1, g_zzzz_0_xyz_1, g_zzzz_0_xz_1, g_zzzz_0_xzz_1, g_zzzz_0_yy_1, g_zzzz_0_yyy_1, g_zzzz_0_yyz_1, g_zzzz_0_yz_1, g_zzzz_0_yzz_1, g_zzzz_0_zz_1, g_zzzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xxx_0[i] = g_zzzz_0_xxx_1[i] * wa_y[i];

        g_yzzzz_0_xxy_0[i] = g_zzzz_0_xx_1[i] * fi_acd_0 + g_zzzz_0_xxy_1[i] * wa_y[i];

        g_yzzzz_0_xxz_0[i] = g_zzzz_0_xxz_1[i] * wa_y[i];

        g_yzzzz_0_xyy_0[i] = 2.0 * g_zzzz_0_xy_1[i] * fi_acd_0 + g_zzzz_0_xyy_1[i] * wa_y[i];

        g_yzzzz_0_xyz_0[i] = g_zzzz_0_xz_1[i] * fi_acd_0 + g_zzzz_0_xyz_1[i] * wa_y[i];

        g_yzzzz_0_xzz_0[i] = g_zzzz_0_xzz_1[i] * wa_y[i];

        g_yzzzz_0_yyy_0[i] = 3.0 * g_zzzz_0_yy_1[i] * fi_acd_0 + g_zzzz_0_yyy_1[i] * wa_y[i];

        g_yzzzz_0_yyz_0[i] = 2.0 * g_zzzz_0_yz_1[i] * fi_acd_0 + g_zzzz_0_yyz_1[i] * wa_y[i];

        g_yzzzz_0_yzz_0[i] = g_zzzz_0_zz_1[i] * fi_acd_0 + g_zzzz_0_yzz_1[i] * wa_y[i];

        g_yzzzz_0_zzz_0[i] = g_zzzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 200-210 components of targeted buffer : HSF

    auto g_zzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_hsf + 200);

    auto g_zzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_hsf + 201);

    auto g_zzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_hsf + 202);

    auto g_zzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_hsf + 203);

    auto g_zzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_hsf + 204);

    auto g_zzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_hsf + 205);

    auto g_zzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_hsf + 206);

    auto g_zzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_hsf + 207);

    auto g_zzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_hsf + 208);

    auto g_zzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_hsf + 209);

    #pragma omp simd aligned(g_zzz_0_xxx_0, g_zzz_0_xxx_1, g_zzz_0_xxy_0, g_zzz_0_xxy_1, g_zzz_0_xxz_0, g_zzz_0_xxz_1, g_zzz_0_xyy_0, g_zzz_0_xyy_1, g_zzz_0_xyz_0, g_zzz_0_xyz_1, g_zzz_0_xzz_0, g_zzz_0_xzz_1, g_zzz_0_yyy_0, g_zzz_0_yyy_1, g_zzz_0_yyz_0, g_zzz_0_yyz_1, g_zzz_0_yzz_0, g_zzz_0_yzz_1, g_zzz_0_zzz_0, g_zzz_0_zzz_1, g_zzzz_0_xx_1, g_zzzz_0_xxx_1, g_zzzz_0_xxy_1, g_zzzz_0_xxz_1, g_zzzz_0_xy_1, g_zzzz_0_xyy_1, g_zzzz_0_xyz_1, g_zzzz_0_xz_1, g_zzzz_0_xzz_1, g_zzzz_0_yy_1, g_zzzz_0_yyy_1, g_zzzz_0_yyz_1, g_zzzz_0_yz_1, g_zzzz_0_yzz_1, g_zzzz_0_zz_1, g_zzzz_0_zzz_1, g_zzzzz_0_xxx_0, g_zzzzz_0_xxy_0, g_zzzzz_0_xxz_0, g_zzzzz_0_xyy_0, g_zzzzz_0_xyz_0, g_zzzzz_0_xzz_0, g_zzzzz_0_yyy_0, g_zzzzz_0_yyz_0, g_zzzzz_0_yzz_0, g_zzzzz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xxx_0[i] = 4.0 * g_zzz_0_xxx_0[i] * fbe_0 - 4.0 * g_zzz_0_xxx_1[i] * fz_be_0 + g_zzzz_0_xxx_1[i] * wa_z[i];

        g_zzzzz_0_xxy_0[i] = 4.0 * g_zzz_0_xxy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxy_1[i] * fz_be_0 + g_zzzz_0_xxy_1[i] * wa_z[i];

        g_zzzzz_0_xxz_0[i] = 4.0 * g_zzz_0_xxz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxz_1[i] * fz_be_0 + g_zzzz_0_xx_1[i] * fi_acd_0 + g_zzzz_0_xxz_1[i] * wa_z[i];

        g_zzzzz_0_xyy_0[i] = 4.0 * g_zzz_0_xyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xyy_1[i] * fz_be_0 + g_zzzz_0_xyy_1[i] * wa_z[i];

        g_zzzzz_0_xyz_0[i] = 4.0 * g_zzz_0_xyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyz_1[i] * fz_be_0 + g_zzzz_0_xy_1[i] * fi_acd_0 + g_zzzz_0_xyz_1[i] * wa_z[i];

        g_zzzzz_0_xzz_0[i] = 4.0 * g_zzz_0_xzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xz_1[i] * fi_acd_0 + g_zzzz_0_xzz_1[i] * wa_z[i];

        g_zzzzz_0_yyy_0[i] = 4.0 * g_zzz_0_yyy_0[i] * fbe_0 - 4.0 * g_zzz_0_yyy_1[i] * fz_be_0 + g_zzzz_0_yyy_1[i] * wa_z[i];

        g_zzzzz_0_yyz_0[i] = 4.0 * g_zzz_0_yyz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyz_1[i] * fz_be_0 + g_zzzz_0_yy_1[i] * fi_acd_0 + g_zzzz_0_yyz_1[i] * wa_z[i];

        g_zzzzz_0_yzz_0[i] = 4.0 * g_zzz_0_yzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_yz_1[i] * fi_acd_0 + g_zzzz_0_yzz_1[i] * wa_z[i];

        g_zzzzz_0_zzz_0[i] = 4.0 * g_zzz_0_zzz_0[i] * fbe_0 - 4.0 * g_zzz_0_zzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_zz_1[i] * fi_acd_0 + g_zzzz_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

