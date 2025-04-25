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

#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_fsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsg,
                                 size_t idx_eri_0_psg,
                                 size_t idx_eri_1_psg,
                                 size_t idx_eri_1_dsf,
                                 size_t idx_eri_1_dsg,
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

    /// Set up components of auxilary buffer : PSG

    auto g_x_0_xxxx_0 = pbuffer.data(idx_eri_0_psg);

    auto g_x_0_xxxy_0 = pbuffer.data(idx_eri_0_psg + 1);

    auto g_x_0_xxxz_0 = pbuffer.data(idx_eri_0_psg + 2);

    auto g_x_0_xxyy_0 = pbuffer.data(idx_eri_0_psg + 3);

    auto g_x_0_xxyz_0 = pbuffer.data(idx_eri_0_psg + 4);

    auto g_x_0_xxzz_0 = pbuffer.data(idx_eri_0_psg + 5);

    auto g_x_0_xyyy_0 = pbuffer.data(idx_eri_0_psg + 6);

    auto g_x_0_xyyz_0 = pbuffer.data(idx_eri_0_psg + 7);

    auto g_x_0_xyzz_0 = pbuffer.data(idx_eri_0_psg + 8);

    auto g_x_0_xzzz_0 = pbuffer.data(idx_eri_0_psg + 9);

    auto g_x_0_yyyy_0 = pbuffer.data(idx_eri_0_psg + 10);

    auto g_x_0_yyyz_0 = pbuffer.data(idx_eri_0_psg + 11);

    auto g_x_0_yyzz_0 = pbuffer.data(idx_eri_0_psg + 12);

    auto g_x_0_yzzz_0 = pbuffer.data(idx_eri_0_psg + 13);

    auto g_x_0_zzzz_0 = pbuffer.data(idx_eri_0_psg + 14);

    auto g_y_0_xxxx_0 = pbuffer.data(idx_eri_0_psg + 15);

    auto g_y_0_xxxy_0 = pbuffer.data(idx_eri_0_psg + 16);

    auto g_y_0_xxxz_0 = pbuffer.data(idx_eri_0_psg + 17);

    auto g_y_0_xxyy_0 = pbuffer.data(idx_eri_0_psg + 18);

    auto g_y_0_xxyz_0 = pbuffer.data(idx_eri_0_psg + 19);

    auto g_y_0_xxzz_0 = pbuffer.data(idx_eri_0_psg + 20);

    auto g_y_0_xyyy_0 = pbuffer.data(idx_eri_0_psg + 21);

    auto g_y_0_xyyz_0 = pbuffer.data(idx_eri_0_psg + 22);

    auto g_y_0_xyzz_0 = pbuffer.data(idx_eri_0_psg + 23);

    auto g_y_0_xzzz_0 = pbuffer.data(idx_eri_0_psg + 24);

    auto g_y_0_yyyy_0 = pbuffer.data(idx_eri_0_psg + 25);

    auto g_y_0_yyyz_0 = pbuffer.data(idx_eri_0_psg + 26);

    auto g_y_0_yyzz_0 = pbuffer.data(idx_eri_0_psg + 27);

    auto g_y_0_yzzz_0 = pbuffer.data(idx_eri_0_psg + 28);

    auto g_y_0_zzzz_0 = pbuffer.data(idx_eri_0_psg + 29);

    auto g_z_0_xxxx_0 = pbuffer.data(idx_eri_0_psg + 30);

    auto g_z_0_xxxy_0 = pbuffer.data(idx_eri_0_psg + 31);

    auto g_z_0_xxxz_0 = pbuffer.data(idx_eri_0_psg + 32);

    auto g_z_0_xxyy_0 = pbuffer.data(idx_eri_0_psg + 33);

    auto g_z_0_xxyz_0 = pbuffer.data(idx_eri_0_psg + 34);

    auto g_z_0_xxzz_0 = pbuffer.data(idx_eri_0_psg + 35);

    auto g_z_0_xyyy_0 = pbuffer.data(idx_eri_0_psg + 36);

    auto g_z_0_xyyz_0 = pbuffer.data(idx_eri_0_psg + 37);

    auto g_z_0_xyzz_0 = pbuffer.data(idx_eri_0_psg + 38);

    auto g_z_0_xzzz_0 = pbuffer.data(idx_eri_0_psg + 39);

    auto g_z_0_yyyy_0 = pbuffer.data(idx_eri_0_psg + 40);

    auto g_z_0_yyyz_0 = pbuffer.data(idx_eri_0_psg + 41);

    auto g_z_0_yyzz_0 = pbuffer.data(idx_eri_0_psg + 42);

    auto g_z_0_yzzz_0 = pbuffer.data(idx_eri_0_psg + 43);

    auto g_z_0_zzzz_0 = pbuffer.data(idx_eri_0_psg + 44);

    /// Set up components of auxilary buffer : PSG

    auto g_x_0_xxxx_1 = pbuffer.data(idx_eri_1_psg);

    auto g_x_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 1);

    auto g_x_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 2);

    auto g_x_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 3);

    auto g_x_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 4);

    auto g_x_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 5);

    auto g_x_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 6);

    auto g_x_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 7);

    auto g_x_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 8);

    auto g_x_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 9);

    auto g_x_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 10);

    auto g_x_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 11);

    auto g_x_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 12);

    auto g_x_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 13);

    auto g_x_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 14);

    auto g_y_0_xxxx_1 = pbuffer.data(idx_eri_1_psg + 15);

    auto g_y_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 16);

    auto g_y_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 17);

    auto g_y_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 18);

    auto g_y_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 19);

    auto g_y_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 20);

    auto g_y_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 21);

    auto g_y_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 22);

    auto g_y_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 23);

    auto g_y_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 24);

    auto g_y_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 25);

    auto g_y_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 26);

    auto g_y_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 27);

    auto g_y_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 28);

    auto g_y_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 29);

    auto g_z_0_xxxx_1 = pbuffer.data(idx_eri_1_psg + 30);

    auto g_z_0_xxxy_1 = pbuffer.data(idx_eri_1_psg + 31);

    auto g_z_0_xxxz_1 = pbuffer.data(idx_eri_1_psg + 32);

    auto g_z_0_xxyy_1 = pbuffer.data(idx_eri_1_psg + 33);

    auto g_z_0_xxyz_1 = pbuffer.data(idx_eri_1_psg + 34);

    auto g_z_0_xxzz_1 = pbuffer.data(idx_eri_1_psg + 35);

    auto g_z_0_xyyy_1 = pbuffer.data(idx_eri_1_psg + 36);

    auto g_z_0_xyyz_1 = pbuffer.data(idx_eri_1_psg + 37);

    auto g_z_0_xyzz_1 = pbuffer.data(idx_eri_1_psg + 38);

    auto g_z_0_xzzz_1 = pbuffer.data(idx_eri_1_psg + 39);

    auto g_z_0_yyyy_1 = pbuffer.data(idx_eri_1_psg + 40);

    auto g_z_0_yyyz_1 = pbuffer.data(idx_eri_1_psg + 41);

    auto g_z_0_yyzz_1 = pbuffer.data(idx_eri_1_psg + 42);

    auto g_z_0_yzzz_1 = pbuffer.data(idx_eri_1_psg + 43);

    auto g_z_0_zzzz_1 = pbuffer.data(idx_eri_1_psg + 44);

    /// Set up components of auxilary buffer : DSF

    auto g_xx_0_xxx_1 = pbuffer.data(idx_eri_1_dsf);

    auto g_xx_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 1);

    auto g_xx_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 2);

    auto g_xx_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 3);

    auto g_xx_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 4);

    auto g_xx_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 5);

    auto g_xx_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 6);

    auto g_xx_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 7);

    auto g_xx_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 8);

    auto g_xx_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 9);

    auto g_yy_0_xxx_1 = pbuffer.data(idx_eri_1_dsf + 30);

    auto g_yy_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 31);

    auto g_yy_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 32);

    auto g_yy_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 33);

    auto g_yy_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 34);

    auto g_yy_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 35);

    auto g_yy_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 36);

    auto g_yy_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 37);

    auto g_yy_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 38);

    auto g_yy_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 39);

    auto g_yz_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 44);

    auto g_yz_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 47);

    auto g_yz_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 48);

    auto g_zz_0_xxx_1 = pbuffer.data(idx_eri_1_dsf + 50);

    auto g_zz_0_xxy_1 = pbuffer.data(idx_eri_1_dsf + 51);

    auto g_zz_0_xxz_1 = pbuffer.data(idx_eri_1_dsf + 52);

    auto g_zz_0_xyy_1 = pbuffer.data(idx_eri_1_dsf + 53);

    auto g_zz_0_xyz_1 = pbuffer.data(idx_eri_1_dsf + 54);

    auto g_zz_0_xzz_1 = pbuffer.data(idx_eri_1_dsf + 55);

    auto g_zz_0_yyy_1 = pbuffer.data(idx_eri_1_dsf + 56);

    auto g_zz_0_yyz_1 = pbuffer.data(idx_eri_1_dsf + 57);

    auto g_zz_0_yzz_1 = pbuffer.data(idx_eri_1_dsf + 58);

    auto g_zz_0_zzz_1 = pbuffer.data(idx_eri_1_dsf + 59);

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

    auto g_xy_0_xxxy_1 = pbuffer.data(idx_eri_1_dsg + 16);

    auto g_xy_0_xxyy_1 = pbuffer.data(idx_eri_1_dsg + 18);

    auto g_xy_0_xyyy_1 = pbuffer.data(idx_eri_1_dsg + 21);

    auto g_xz_0_xxxx_1 = pbuffer.data(idx_eri_1_dsg + 30);

    auto g_xz_0_xxxz_1 = pbuffer.data(idx_eri_1_dsg + 32);

    auto g_xz_0_xxzz_1 = pbuffer.data(idx_eri_1_dsg + 35);

    auto g_xz_0_xzzz_1 = pbuffer.data(idx_eri_1_dsg + 39);

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

    auto g_yz_0_xxyz_1 = pbuffer.data(idx_eri_1_dsg + 64);

    auto g_yz_0_xyyz_1 = pbuffer.data(idx_eri_1_dsg + 67);

    auto g_yz_0_xyzz_1 = pbuffer.data(idx_eri_1_dsg + 68);

    auto g_yz_0_yyyy_1 = pbuffer.data(idx_eri_1_dsg + 70);

    auto g_yz_0_yyyz_1 = pbuffer.data(idx_eri_1_dsg + 71);

    auto g_yz_0_yyzz_1 = pbuffer.data(idx_eri_1_dsg + 72);

    auto g_yz_0_yzzz_1 = pbuffer.data(idx_eri_1_dsg + 73);

    auto g_yz_0_zzzz_1 = pbuffer.data(idx_eri_1_dsg + 74);

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

    /// Set up 0-15 components of targeted buffer : FSG

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

    #pragma omp simd aligned(g_x_0_xxxx_0, g_x_0_xxxx_1, g_x_0_xxxy_0, g_x_0_xxxy_1, g_x_0_xxxz_0, g_x_0_xxxz_1, g_x_0_xxyy_0, g_x_0_xxyy_1, g_x_0_xxyz_0, g_x_0_xxyz_1, g_x_0_xxzz_0, g_x_0_xxzz_1, g_x_0_xyyy_0, g_x_0_xyyy_1, g_x_0_xyyz_0, g_x_0_xyyz_1, g_x_0_xyzz_0, g_x_0_xyzz_1, g_x_0_xzzz_0, g_x_0_xzzz_1, g_x_0_yyyy_0, g_x_0_yyyy_1, g_x_0_yyyz_0, g_x_0_yyyz_1, g_x_0_yyzz_0, g_x_0_yyzz_1, g_x_0_yzzz_0, g_x_0_yzzz_1, g_x_0_zzzz_0, g_x_0_zzzz_1, g_xx_0_xxx_1, g_xx_0_xxxx_1, g_xx_0_xxxy_1, g_xx_0_xxxz_1, g_xx_0_xxy_1, g_xx_0_xxyy_1, g_xx_0_xxyz_1, g_xx_0_xxz_1, g_xx_0_xxzz_1, g_xx_0_xyy_1, g_xx_0_xyyy_1, g_xx_0_xyyz_1, g_xx_0_xyz_1, g_xx_0_xyzz_1, g_xx_0_xzz_1, g_xx_0_xzzz_1, g_xx_0_yyy_1, g_xx_0_yyyy_1, g_xx_0_yyyz_1, g_xx_0_yyz_1, g_xx_0_yyzz_1, g_xx_0_yzz_1, g_xx_0_yzzz_1, g_xx_0_zzz_1, g_xx_0_zzzz_1, g_xxx_0_xxxx_0, g_xxx_0_xxxy_0, g_xxx_0_xxxz_0, g_xxx_0_xxyy_0, g_xxx_0_xxyz_0, g_xxx_0_xxzz_0, g_xxx_0_xyyy_0, g_xxx_0_xyyz_0, g_xxx_0_xyzz_0, g_xxx_0_xzzz_0, g_xxx_0_yyyy_0, g_xxx_0_yyyz_0, g_xxx_0_yyzz_0, g_xxx_0_yzzz_0, g_xxx_0_zzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxx_0_xxxx_0[i] = 2.0 * g_x_0_xxxx_0[i] * fbe_0 - 2.0 * g_x_0_xxxx_1[i] * fz_be_0 + 4.0 * g_xx_0_xxx_1[i] * fi_acd_0 + g_xx_0_xxxx_1[i] * wa_x[i];

        g_xxx_0_xxxy_0[i] = 2.0 * g_x_0_xxxy_0[i] * fbe_0 - 2.0 * g_x_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xx_0_xxy_1[i] * fi_acd_0 + g_xx_0_xxxy_1[i] * wa_x[i];

        g_xxx_0_xxxz_0[i] = 2.0 * g_x_0_xxxz_0[i] * fbe_0 - 2.0 * g_x_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xx_0_xxz_1[i] * fi_acd_0 + g_xx_0_xxxz_1[i] * wa_x[i];

        g_xxx_0_xxyy_0[i] = 2.0 * g_x_0_xxyy_0[i] * fbe_0 - 2.0 * g_x_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xx_0_xyy_1[i] * fi_acd_0 + g_xx_0_xxyy_1[i] * wa_x[i];

        g_xxx_0_xxyz_0[i] = 2.0 * g_x_0_xxyz_0[i] * fbe_0 - 2.0 * g_x_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xx_0_xyz_1[i] * fi_acd_0 + g_xx_0_xxyz_1[i] * wa_x[i];

        g_xxx_0_xxzz_0[i] = 2.0 * g_x_0_xxzz_0[i] * fbe_0 - 2.0 * g_x_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xx_0_xzz_1[i] * fi_acd_0 + g_xx_0_xxzz_1[i] * wa_x[i];

        g_xxx_0_xyyy_0[i] = 2.0 * g_x_0_xyyy_0[i] * fbe_0 - 2.0 * g_x_0_xyyy_1[i] * fz_be_0 + g_xx_0_yyy_1[i] * fi_acd_0 + g_xx_0_xyyy_1[i] * wa_x[i];

        g_xxx_0_xyyz_0[i] = 2.0 * g_x_0_xyyz_0[i] * fbe_0 - 2.0 * g_x_0_xyyz_1[i] * fz_be_0 + g_xx_0_yyz_1[i] * fi_acd_0 + g_xx_0_xyyz_1[i] * wa_x[i];

        g_xxx_0_xyzz_0[i] = 2.0 * g_x_0_xyzz_0[i] * fbe_0 - 2.0 * g_x_0_xyzz_1[i] * fz_be_0 + g_xx_0_yzz_1[i] * fi_acd_0 + g_xx_0_xyzz_1[i] * wa_x[i];

        g_xxx_0_xzzz_0[i] = 2.0 * g_x_0_xzzz_0[i] * fbe_0 - 2.0 * g_x_0_xzzz_1[i] * fz_be_0 + g_xx_0_zzz_1[i] * fi_acd_0 + g_xx_0_xzzz_1[i] * wa_x[i];

        g_xxx_0_yyyy_0[i] = 2.0 * g_x_0_yyyy_0[i] * fbe_0 - 2.0 * g_x_0_yyyy_1[i] * fz_be_0 + g_xx_0_yyyy_1[i] * wa_x[i];

        g_xxx_0_yyyz_0[i] = 2.0 * g_x_0_yyyz_0[i] * fbe_0 - 2.0 * g_x_0_yyyz_1[i] * fz_be_0 + g_xx_0_yyyz_1[i] * wa_x[i];

        g_xxx_0_yyzz_0[i] = 2.0 * g_x_0_yyzz_0[i] * fbe_0 - 2.0 * g_x_0_yyzz_1[i] * fz_be_0 + g_xx_0_yyzz_1[i] * wa_x[i];

        g_xxx_0_yzzz_0[i] = 2.0 * g_x_0_yzzz_0[i] * fbe_0 - 2.0 * g_x_0_yzzz_1[i] * fz_be_0 + g_xx_0_yzzz_1[i] * wa_x[i];

        g_xxx_0_zzzz_0[i] = 2.0 * g_x_0_zzzz_0[i] * fbe_0 - 2.0 * g_x_0_zzzz_1[i] * fz_be_0 + g_xx_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : FSG

    auto g_xxy_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 15);

    auto g_xxy_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 16);

    auto g_xxy_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 17);

    auto g_xxy_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 18);

    auto g_xxy_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 19);

    auto g_xxy_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 20);

    auto g_xxy_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 21);

    auto g_xxy_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 22);

    auto g_xxy_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 23);

    auto g_xxy_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 24);

    auto g_xxy_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 25);

    auto g_xxy_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 26);

    auto g_xxy_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 27);

    auto g_xxy_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 28);

    auto g_xxy_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 29);

    #pragma omp simd aligned(g_xx_0_xxx_1, g_xx_0_xxxx_1, g_xx_0_xxxy_1, g_xx_0_xxxz_1, g_xx_0_xxy_1, g_xx_0_xxyy_1, g_xx_0_xxyz_1, g_xx_0_xxz_1, g_xx_0_xxzz_1, g_xx_0_xyy_1, g_xx_0_xyyy_1, g_xx_0_xyyz_1, g_xx_0_xyz_1, g_xx_0_xyzz_1, g_xx_0_xzz_1, g_xx_0_xzzz_1, g_xx_0_yyy_1, g_xx_0_yyyy_1, g_xx_0_yyyz_1, g_xx_0_yyz_1, g_xx_0_yyzz_1, g_xx_0_yzz_1, g_xx_0_yzzz_1, g_xx_0_zzz_1, g_xx_0_zzzz_1, g_xxy_0_xxxx_0, g_xxy_0_xxxy_0, g_xxy_0_xxxz_0, g_xxy_0_xxyy_0, g_xxy_0_xxyz_0, g_xxy_0_xxzz_0, g_xxy_0_xyyy_0, g_xxy_0_xyyz_0, g_xxy_0_xyzz_0, g_xxy_0_xzzz_0, g_xxy_0_yyyy_0, g_xxy_0_yyyz_0, g_xxy_0_yyzz_0, g_xxy_0_yzzz_0, g_xxy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxy_0_xxxx_0[i] = g_xx_0_xxxx_1[i] * wa_y[i];

        g_xxy_0_xxxy_0[i] = g_xx_0_xxx_1[i] * fi_acd_0 + g_xx_0_xxxy_1[i] * wa_y[i];

        g_xxy_0_xxxz_0[i] = g_xx_0_xxxz_1[i] * wa_y[i];

        g_xxy_0_xxyy_0[i] = 2.0 * g_xx_0_xxy_1[i] * fi_acd_0 + g_xx_0_xxyy_1[i] * wa_y[i];

        g_xxy_0_xxyz_0[i] = g_xx_0_xxz_1[i] * fi_acd_0 + g_xx_0_xxyz_1[i] * wa_y[i];

        g_xxy_0_xxzz_0[i] = g_xx_0_xxzz_1[i] * wa_y[i];

        g_xxy_0_xyyy_0[i] = 3.0 * g_xx_0_xyy_1[i] * fi_acd_0 + g_xx_0_xyyy_1[i] * wa_y[i];

        g_xxy_0_xyyz_0[i] = 2.0 * g_xx_0_xyz_1[i] * fi_acd_0 + g_xx_0_xyyz_1[i] * wa_y[i];

        g_xxy_0_xyzz_0[i] = g_xx_0_xzz_1[i] * fi_acd_0 + g_xx_0_xyzz_1[i] * wa_y[i];

        g_xxy_0_xzzz_0[i] = g_xx_0_xzzz_1[i] * wa_y[i];

        g_xxy_0_yyyy_0[i] = 4.0 * g_xx_0_yyy_1[i] * fi_acd_0 + g_xx_0_yyyy_1[i] * wa_y[i];

        g_xxy_0_yyyz_0[i] = 3.0 * g_xx_0_yyz_1[i] * fi_acd_0 + g_xx_0_yyyz_1[i] * wa_y[i];

        g_xxy_0_yyzz_0[i] = 2.0 * g_xx_0_yzz_1[i] * fi_acd_0 + g_xx_0_yyzz_1[i] * wa_y[i];

        g_xxy_0_yzzz_0[i] = g_xx_0_zzz_1[i] * fi_acd_0 + g_xx_0_yzzz_1[i] * wa_y[i];

        g_xxy_0_zzzz_0[i] = g_xx_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 30-45 components of targeted buffer : FSG

    auto g_xxz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 30);

    auto g_xxz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 31);

    auto g_xxz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 32);

    auto g_xxz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 33);

    auto g_xxz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 34);

    auto g_xxz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 35);

    auto g_xxz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 36);

    auto g_xxz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 37);

    auto g_xxz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 38);

    auto g_xxz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 39);

    auto g_xxz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 40);

    auto g_xxz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 41);

    auto g_xxz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 42);

    auto g_xxz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 43);

    auto g_xxz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 44);

    #pragma omp simd aligned(g_xx_0_xxx_1, g_xx_0_xxxx_1, g_xx_0_xxxy_1, g_xx_0_xxxz_1, g_xx_0_xxy_1, g_xx_0_xxyy_1, g_xx_0_xxyz_1, g_xx_0_xxz_1, g_xx_0_xxzz_1, g_xx_0_xyy_1, g_xx_0_xyyy_1, g_xx_0_xyyz_1, g_xx_0_xyz_1, g_xx_0_xyzz_1, g_xx_0_xzz_1, g_xx_0_xzzz_1, g_xx_0_yyy_1, g_xx_0_yyyy_1, g_xx_0_yyyz_1, g_xx_0_yyz_1, g_xx_0_yyzz_1, g_xx_0_yzz_1, g_xx_0_yzzz_1, g_xx_0_zzz_1, g_xx_0_zzzz_1, g_xxz_0_xxxx_0, g_xxz_0_xxxy_0, g_xxz_0_xxxz_0, g_xxz_0_xxyy_0, g_xxz_0_xxyz_0, g_xxz_0_xxzz_0, g_xxz_0_xyyy_0, g_xxz_0_xyyz_0, g_xxz_0_xyzz_0, g_xxz_0_xzzz_0, g_xxz_0_yyyy_0, g_xxz_0_yyyz_0, g_xxz_0_yyzz_0, g_xxz_0_yzzz_0, g_xxz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxz_0_xxxx_0[i] = g_xx_0_xxxx_1[i] * wa_z[i];

        g_xxz_0_xxxy_0[i] = g_xx_0_xxxy_1[i] * wa_z[i];

        g_xxz_0_xxxz_0[i] = g_xx_0_xxx_1[i] * fi_acd_0 + g_xx_0_xxxz_1[i] * wa_z[i];

        g_xxz_0_xxyy_0[i] = g_xx_0_xxyy_1[i] * wa_z[i];

        g_xxz_0_xxyz_0[i] = g_xx_0_xxy_1[i] * fi_acd_0 + g_xx_0_xxyz_1[i] * wa_z[i];

        g_xxz_0_xxzz_0[i] = 2.0 * g_xx_0_xxz_1[i] * fi_acd_0 + g_xx_0_xxzz_1[i] * wa_z[i];

        g_xxz_0_xyyy_0[i] = g_xx_0_xyyy_1[i] * wa_z[i];

        g_xxz_0_xyyz_0[i] = g_xx_0_xyy_1[i] * fi_acd_0 + g_xx_0_xyyz_1[i] * wa_z[i];

        g_xxz_0_xyzz_0[i] = 2.0 * g_xx_0_xyz_1[i] * fi_acd_0 + g_xx_0_xyzz_1[i] * wa_z[i];

        g_xxz_0_xzzz_0[i] = 3.0 * g_xx_0_xzz_1[i] * fi_acd_0 + g_xx_0_xzzz_1[i] * wa_z[i];

        g_xxz_0_yyyy_0[i] = g_xx_0_yyyy_1[i] * wa_z[i];

        g_xxz_0_yyyz_0[i] = g_xx_0_yyy_1[i] * fi_acd_0 + g_xx_0_yyyz_1[i] * wa_z[i];

        g_xxz_0_yyzz_0[i] = 2.0 * g_xx_0_yyz_1[i] * fi_acd_0 + g_xx_0_yyzz_1[i] * wa_z[i];

        g_xxz_0_yzzz_0[i] = 3.0 * g_xx_0_yzz_1[i] * fi_acd_0 + g_xx_0_yzzz_1[i] * wa_z[i];

        g_xxz_0_zzzz_0[i] = 4.0 * g_xx_0_zzz_1[i] * fi_acd_0 + g_xx_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 45-60 components of targeted buffer : FSG

    auto g_xyy_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 45);

    auto g_xyy_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 46);

    auto g_xyy_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 47);

    auto g_xyy_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 48);

    auto g_xyy_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 49);

    auto g_xyy_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 50);

    auto g_xyy_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 51);

    auto g_xyy_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 52);

    auto g_xyy_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 53);

    auto g_xyy_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 54);

    auto g_xyy_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 55);

    auto g_xyy_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 56);

    auto g_xyy_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 57);

    auto g_xyy_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 58);

    auto g_xyy_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 59);

    #pragma omp simd aligned(g_xyy_0_xxxx_0, g_xyy_0_xxxy_0, g_xyy_0_xxxz_0, g_xyy_0_xxyy_0, g_xyy_0_xxyz_0, g_xyy_0_xxzz_0, g_xyy_0_xyyy_0, g_xyy_0_xyyz_0, g_xyy_0_xyzz_0, g_xyy_0_xzzz_0, g_xyy_0_yyyy_0, g_xyy_0_yyyz_0, g_xyy_0_yyzz_0, g_xyy_0_yzzz_0, g_xyy_0_zzzz_0, g_yy_0_xxx_1, g_yy_0_xxxx_1, g_yy_0_xxxy_1, g_yy_0_xxxz_1, g_yy_0_xxy_1, g_yy_0_xxyy_1, g_yy_0_xxyz_1, g_yy_0_xxz_1, g_yy_0_xxzz_1, g_yy_0_xyy_1, g_yy_0_xyyy_1, g_yy_0_xyyz_1, g_yy_0_xyz_1, g_yy_0_xyzz_1, g_yy_0_xzz_1, g_yy_0_xzzz_1, g_yy_0_yyy_1, g_yy_0_yyyy_1, g_yy_0_yyyz_1, g_yy_0_yyz_1, g_yy_0_yyzz_1, g_yy_0_yzz_1, g_yy_0_yzzz_1, g_yy_0_zzz_1, g_yy_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyy_0_xxxx_0[i] = 4.0 * g_yy_0_xxx_1[i] * fi_acd_0 + g_yy_0_xxxx_1[i] * wa_x[i];

        g_xyy_0_xxxy_0[i] = 3.0 * g_yy_0_xxy_1[i] * fi_acd_0 + g_yy_0_xxxy_1[i] * wa_x[i];

        g_xyy_0_xxxz_0[i] = 3.0 * g_yy_0_xxz_1[i] * fi_acd_0 + g_yy_0_xxxz_1[i] * wa_x[i];

        g_xyy_0_xxyy_0[i] = 2.0 * g_yy_0_xyy_1[i] * fi_acd_0 + g_yy_0_xxyy_1[i] * wa_x[i];

        g_xyy_0_xxyz_0[i] = 2.0 * g_yy_0_xyz_1[i] * fi_acd_0 + g_yy_0_xxyz_1[i] * wa_x[i];

        g_xyy_0_xxzz_0[i] = 2.0 * g_yy_0_xzz_1[i] * fi_acd_0 + g_yy_0_xxzz_1[i] * wa_x[i];

        g_xyy_0_xyyy_0[i] = g_yy_0_yyy_1[i] * fi_acd_0 + g_yy_0_xyyy_1[i] * wa_x[i];

        g_xyy_0_xyyz_0[i] = g_yy_0_yyz_1[i] * fi_acd_0 + g_yy_0_xyyz_1[i] * wa_x[i];

        g_xyy_0_xyzz_0[i] = g_yy_0_yzz_1[i] * fi_acd_0 + g_yy_0_xyzz_1[i] * wa_x[i];

        g_xyy_0_xzzz_0[i] = g_yy_0_zzz_1[i] * fi_acd_0 + g_yy_0_xzzz_1[i] * wa_x[i];

        g_xyy_0_yyyy_0[i] = g_yy_0_yyyy_1[i] * wa_x[i];

        g_xyy_0_yyyz_0[i] = g_yy_0_yyyz_1[i] * wa_x[i];

        g_xyy_0_yyzz_0[i] = g_yy_0_yyzz_1[i] * wa_x[i];

        g_xyy_0_yzzz_0[i] = g_yy_0_yzzz_1[i] * wa_x[i];

        g_xyy_0_zzzz_0[i] = g_yy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 60-75 components of targeted buffer : FSG

    auto g_xyz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 60);

    auto g_xyz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 61);

    auto g_xyz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 62);

    auto g_xyz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 63);

    auto g_xyz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 64);

    auto g_xyz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 65);

    auto g_xyz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 66);

    auto g_xyz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 67);

    auto g_xyz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 68);

    auto g_xyz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 69);

    auto g_xyz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 70);

    auto g_xyz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 71);

    auto g_xyz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 72);

    auto g_xyz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 73);

    auto g_xyz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 74);

    #pragma omp simd aligned(g_xy_0_xxxy_1, g_xy_0_xxyy_1, g_xy_0_xyyy_1, g_xyz_0_xxxx_0, g_xyz_0_xxxy_0, g_xyz_0_xxxz_0, g_xyz_0_xxyy_0, g_xyz_0_xxyz_0, g_xyz_0_xxzz_0, g_xyz_0_xyyy_0, g_xyz_0_xyyz_0, g_xyz_0_xyzz_0, g_xyz_0_xzzz_0, g_xyz_0_yyyy_0, g_xyz_0_yyyz_0, g_xyz_0_yyzz_0, g_xyz_0_yzzz_0, g_xyz_0_zzzz_0, g_xz_0_xxxx_1, g_xz_0_xxxz_1, g_xz_0_xxzz_1, g_xz_0_xzzz_1, g_yz_0_xxyz_1, g_yz_0_xyyz_1, g_yz_0_xyz_1, g_yz_0_xyzz_1, g_yz_0_yyyy_1, g_yz_0_yyyz_1, g_yz_0_yyz_1, g_yz_0_yyzz_1, g_yz_0_yzz_1, g_yz_0_yzzz_1, g_yz_0_zzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyz_0_xxxx_0[i] = g_xz_0_xxxx_1[i] * wa_y[i];

        g_xyz_0_xxxy_0[i] = g_xy_0_xxxy_1[i] * wa_z[i];

        g_xyz_0_xxxz_0[i] = g_xz_0_xxxz_1[i] * wa_y[i];

        g_xyz_0_xxyy_0[i] = g_xy_0_xxyy_1[i] * wa_z[i];

        g_xyz_0_xxyz_0[i] = 2.0 * g_yz_0_xyz_1[i] * fi_acd_0 + g_yz_0_xxyz_1[i] * wa_x[i];

        g_xyz_0_xxzz_0[i] = g_xz_0_xxzz_1[i] * wa_y[i];

        g_xyz_0_xyyy_0[i] = g_xy_0_xyyy_1[i] * wa_z[i];

        g_xyz_0_xyyz_0[i] = g_yz_0_yyz_1[i] * fi_acd_0 + g_yz_0_xyyz_1[i] * wa_x[i];

        g_xyz_0_xyzz_0[i] = g_yz_0_yzz_1[i] * fi_acd_0 + g_yz_0_xyzz_1[i] * wa_x[i];

        g_xyz_0_xzzz_0[i] = g_xz_0_xzzz_1[i] * wa_y[i];

        g_xyz_0_yyyy_0[i] = g_yz_0_yyyy_1[i] * wa_x[i];

        g_xyz_0_yyyz_0[i] = g_yz_0_yyyz_1[i] * wa_x[i];

        g_xyz_0_yyzz_0[i] = g_yz_0_yyzz_1[i] * wa_x[i];

        g_xyz_0_yzzz_0[i] = g_yz_0_yzzz_1[i] * wa_x[i];

        g_xyz_0_zzzz_0[i] = g_yz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 75-90 components of targeted buffer : FSG

    auto g_xzz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 75);

    auto g_xzz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 76);

    auto g_xzz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 77);

    auto g_xzz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 78);

    auto g_xzz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 79);

    auto g_xzz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 80);

    auto g_xzz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 81);

    auto g_xzz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 82);

    auto g_xzz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 83);

    auto g_xzz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 84);

    auto g_xzz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 85);

    auto g_xzz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 86);

    auto g_xzz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 87);

    auto g_xzz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 88);

    auto g_xzz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 89);

    #pragma omp simd aligned(g_xzz_0_xxxx_0, g_xzz_0_xxxy_0, g_xzz_0_xxxz_0, g_xzz_0_xxyy_0, g_xzz_0_xxyz_0, g_xzz_0_xxzz_0, g_xzz_0_xyyy_0, g_xzz_0_xyyz_0, g_xzz_0_xyzz_0, g_xzz_0_xzzz_0, g_xzz_0_yyyy_0, g_xzz_0_yyyz_0, g_xzz_0_yyzz_0, g_xzz_0_yzzz_0, g_xzz_0_zzzz_0, g_zz_0_xxx_1, g_zz_0_xxxx_1, g_zz_0_xxxy_1, g_zz_0_xxxz_1, g_zz_0_xxy_1, g_zz_0_xxyy_1, g_zz_0_xxyz_1, g_zz_0_xxz_1, g_zz_0_xxzz_1, g_zz_0_xyy_1, g_zz_0_xyyy_1, g_zz_0_xyyz_1, g_zz_0_xyz_1, g_zz_0_xyzz_1, g_zz_0_xzz_1, g_zz_0_xzzz_1, g_zz_0_yyy_1, g_zz_0_yyyy_1, g_zz_0_yyyz_1, g_zz_0_yyz_1, g_zz_0_yyzz_1, g_zz_0_yzz_1, g_zz_0_yzzz_1, g_zz_0_zzz_1, g_zz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzz_0_xxxx_0[i] = 4.0 * g_zz_0_xxx_1[i] * fi_acd_0 + g_zz_0_xxxx_1[i] * wa_x[i];

        g_xzz_0_xxxy_0[i] = 3.0 * g_zz_0_xxy_1[i] * fi_acd_0 + g_zz_0_xxxy_1[i] * wa_x[i];

        g_xzz_0_xxxz_0[i] = 3.0 * g_zz_0_xxz_1[i] * fi_acd_0 + g_zz_0_xxxz_1[i] * wa_x[i];

        g_xzz_0_xxyy_0[i] = 2.0 * g_zz_0_xyy_1[i] * fi_acd_0 + g_zz_0_xxyy_1[i] * wa_x[i];

        g_xzz_0_xxyz_0[i] = 2.0 * g_zz_0_xyz_1[i] * fi_acd_0 + g_zz_0_xxyz_1[i] * wa_x[i];

        g_xzz_0_xxzz_0[i] = 2.0 * g_zz_0_xzz_1[i] * fi_acd_0 + g_zz_0_xxzz_1[i] * wa_x[i];

        g_xzz_0_xyyy_0[i] = g_zz_0_yyy_1[i] * fi_acd_0 + g_zz_0_xyyy_1[i] * wa_x[i];

        g_xzz_0_xyyz_0[i] = g_zz_0_yyz_1[i] * fi_acd_0 + g_zz_0_xyyz_1[i] * wa_x[i];

        g_xzz_0_xyzz_0[i] = g_zz_0_yzz_1[i] * fi_acd_0 + g_zz_0_xyzz_1[i] * wa_x[i];

        g_xzz_0_xzzz_0[i] = g_zz_0_zzz_1[i] * fi_acd_0 + g_zz_0_xzzz_1[i] * wa_x[i];

        g_xzz_0_yyyy_0[i] = g_zz_0_yyyy_1[i] * wa_x[i];

        g_xzz_0_yyyz_0[i] = g_zz_0_yyyz_1[i] * wa_x[i];

        g_xzz_0_yyzz_0[i] = g_zz_0_yyzz_1[i] * wa_x[i];

        g_xzz_0_yzzz_0[i] = g_zz_0_yzzz_1[i] * wa_x[i];

        g_xzz_0_zzzz_0[i] = g_zz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 90-105 components of targeted buffer : FSG

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

    #pragma omp simd aligned(g_y_0_xxxx_0, g_y_0_xxxx_1, g_y_0_xxxy_0, g_y_0_xxxy_1, g_y_0_xxxz_0, g_y_0_xxxz_1, g_y_0_xxyy_0, g_y_0_xxyy_1, g_y_0_xxyz_0, g_y_0_xxyz_1, g_y_0_xxzz_0, g_y_0_xxzz_1, g_y_0_xyyy_0, g_y_0_xyyy_1, g_y_0_xyyz_0, g_y_0_xyyz_1, g_y_0_xyzz_0, g_y_0_xyzz_1, g_y_0_xzzz_0, g_y_0_xzzz_1, g_y_0_yyyy_0, g_y_0_yyyy_1, g_y_0_yyyz_0, g_y_0_yyyz_1, g_y_0_yyzz_0, g_y_0_yyzz_1, g_y_0_yzzz_0, g_y_0_yzzz_1, g_y_0_zzzz_0, g_y_0_zzzz_1, g_yy_0_xxx_1, g_yy_0_xxxx_1, g_yy_0_xxxy_1, g_yy_0_xxxz_1, g_yy_0_xxy_1, g_yy_0_xxyy_1, g_yy_0_xxyz_1, g_yy_0_xxz_1, g_yy_0_xxzz_1, g_yy_0_xyy_1, g_yy_0_xyyy_1, g_yy_0_xyyz_1, g_yy_0_xyz_1, g_yy_0_xyzz_1, g_yy_0_xzz_1, g_yy_0_xzzz_1, g_yy_0_yyy_1, g_yy_0_yyyy_1, g_yy_0_yyyz_1, g_yy_0_yyz_1, g_yy_0_yyzz_1, g_yy_0_yzz_1, g_yy_0_yzzz_1, g_yy_0_zzz_1, g_yy_0_zzzz_1, g_yyy_0_xxxx_0, g_yyy_0_xxxy_0, g_yyy_0_xxxz_0, g_yyy_0_xxyy_0, g_yyy_0_xxyz_0, g_yyy_0_xxzz_0, g_yyy_0_xyyy_0, g_yyy_0_xyyz_0, g_yyy_0_xyzz_0, g_yyy_0_xzzz_0, g_yyy_0_yyyy_0, g_yyy_0_yyyz_0, g_yyy_0_yyzz_0, g_yyy_0_yzzz_0, g_yyy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyy_0_xxxx_0[i] = 2.0 * g_y_0_xxxx_0[i] * fbe_0 - 2.0 * g_y_0_xxxx_1[i] * fz_be_0 + g_yy_0_xxxx_1[i] * wa_y[i];

        g_yyy_0_xxxy_0[i] = 2.0 * g_y_0_xxxy_0[i] * fbe_0 - 2.0 * g_y_0_xxxy_1[i] * fz_be_0 + g_yy_0_xxx_1[i] * fi_acd_0 + g_yy_0_xxxy_1[i] * wa_y[i];

        g_yyy_0_xxxz_0[i] = 2.0 * g_y_0_xxxz_0[i] * fbe_0 - 2.0 * g_y_0_xxxz_1[i] * fz_be_0 + g_yy_0_xxxz_1[i] * wa_y[i];

        g_yyy_0_xxyy_0[i] = 2.0 * g_y_0_xxyy_0[i] * fbe_0 - 2.0 * g_y_0_xxyy_1[i] * fz_be_0 + 2.0 * g_yy_0_xxy_1[i] * fi_acd_0 + g_yy_0_xxyy_1[i] * wa_y[i];

        g_yyy_0_xxyz_0[i] = 2.0 * g_y_0_xxyz_0[i] * fbe_0 - 2.0 * g_y_0_xxyz_1[i] * fz_be_0 + g_yy_0_xxz_1[i] * fi_acd_0 + g_yy_0_xxyz_1[i] * wa_y[i];

        g_yyy_0_xxzz_0[i] = 2.0 * g_y_0_xxzz_0[i] * fbe_0 - 2.0 * g_y_0_xxzz_1[i] * fz_be_0 + g_yy_0_xxzz_1[i] * wa_y[i];

        g_yyy_0_xyyy_0[i] = 2.0 * g_y_0_xyyy_0[i] * fbe_0 - 2.0 * g_y_0_xyyy_1[i] * fz_be_0 + 3.0 * g_yy_0_xyy_1[i] * fi_acd_0 + g_yy_0_xyyy_1[i] * wa_y[i];

        g_yyy_0_xyyz_0[i] = 2.0 * g_y_0_xyyz_0[i] * fbe_0 - 2.0 * g_y_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yy_0_xyz_1[i] * fi_acd_0 + g_yy_0_xyyz_1[i] * wa_y[i];

        g_yyy_0_xyzz_0[i] = 2.0 * g_y_0_xyzz_0[i] * fbe_0 - 2.0 * g_y_0_xyzz_1[i] * fz_be_0 + g_yy_0_xzz_1[i] * fi_acd_0 + g_yy_0_xyzz_1[i] * wa_y[i];

        g_yyy_0_xzzz_0[i] = 2.0 * g_y_0_xzzz_0[i] * fbe_0 - 2.0 * g_y_0_xzzz_1[i] * fz_be_0 + g_yy_0_xzzz_1[i] * wa_y[i];

        g_yyy_0_yyyy_0[i] = 2.0 * g_y_0_yyyy_0[i] * fbe_0 - 2.0 * g_y_0_yyyy_1[i] * fz_be_0 + 4.0 * g_yy_0_yyy_1[i] * fi_acd_0 + g_yy_0_yyyy_1[i] * wa_y[i];

        g_yyy_0_yyyz_0[i] = 2.0 * g_y_0_yyyz_0[i] * fbe_0 - 2.0 * g_y_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yy_0_yyz_1[i] * fi_acd_0 + g_yy_0_yyyz_1[i] * wa_y[i];

        g_yyy_0_yyzz_0[i] = 2.0 * g_y_0_yyzz_0[i] * fbe_0 - 2.0 * g_y_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yy_0_yzz_1[i] * fi_acd_0 + g_yy_0_yyzz_1[i] * wa_y[i];

        g_yyy_0_yzzz_0[i] = 2.0 * g_y_0_yzzz_0[i] * fbe_0 - 2.0 * g_y_0_yzzz_1[i] * fz_be_0 + g_yy_0_zzz_1[i] * fi_acd_0 + g_yy_0_yzzz_1[i] * wa_y[i];

        g_yyy_0_zzzz_0[i] = 2.0 * g_y_0_zzzz_0[i] * fbe_0 - 2.0 * g_y_0_zzzz_1[i] * fz_be_0 + g_yy_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 105-120 components of targeted buffer : FSG

    auto g_yyz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 105);

    auto g_yyz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 106);

    auto g_yyz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 107);

    auto g_yyz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 108);

    auto g_yyz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 109);

    auto g_yyz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 110);

    auto g_yyz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 111);

    auto g_yyz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 112);

    auto g_yyz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 113);

    auto g_yyz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 114);

    auto g_yyz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 115);

    auto g_yyz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 116);

    auto g_yyz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 117);

    auto g_yyz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 118);

    auto g_yyz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 119);

    #pragma omp simd aligned(g_yy_0_xxx_1, g_yy_0_xxxx_1, g_yy_0_xxxy_1, g_yy_0_xxxz_1, g_yy_0_xxy_1, g_yy_0_xxyy_1, g_yy_0_xxyz_1, g_yy_0_xxz_1, g_yy_0_xxzz_1, g_yy_0_xyy_1, g_yy_0_xyyy_1, g_yy_0_xyyz_1, g_yy_0_xyz_1, g_yy_0_xyzz_1, g_yy_0_xzz_1, g_yy_0_xzzz_1, g_yy_0_yyy_1, g_yy_0_yyyy_1, g_yy_0_yyyz_1, g_yy_0_yyz_1, g_yy_0_yyzz_1, g_yy_0_yzz_1, g_yy_0_yzzz_1, g_yy_0_zzz_1, g_yy_0_zzzz_1, g_yyz_0_xxxx_0, g_yyz_0_xxxy_0, g_yyz_0_xxxz_0, g_yyz_0_xxyy_0, g_yyz_0_xxyz_0, g_yyz_0_xxzz_0, g_yyz_0_xyyy_0, g_yyz_0_xyyz_0, g_yyz_0_xyzz_0, g_yyz_0_xzzz_0, g_yyz_0_yyyy_0, g_yyz_0_yyyz_0, g_yyz_0_yyzz_0, g_yyz_0_yzzz_0, g_yyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyz_0_xxxx_0[i] = g_yy_0_xxxx_1[i] * wa_z[i];

        g_yyz_0_xxxy_0[i] = g_yy_0_xxxy_1[i] * wa_z[i];

        g_yyz_0_xxxz_0[i] = g_yy_0_xxx_1[i] * fi_acd_0 + g_yy_0_xxxz_1[i] * wa_z[i];

        g_yyz_0_xxyy_0[i] = g_yy_0_xxyy_1[i] * wa_z[i];

        g_yyz_0_xxyz_0[i] = g_yy_0_xxy_1[i] * fi_acd_0 + g_yy_0_xxyz_1[i] * wa_z[i];

        g_yyz_0_xxzz_0[i] = 2.0 * g_yy_0_xxz_1[i] * fi_acd_0 + g_yy_0_xxzz_1[i] * wa_z[i];

        g_yyz_0_xyyy_0[i] = g_yy_0_xyyy_1[i] * wa_z[i];

        g_yyz_0_xyyz_0[i] = g_yy_0_xyy_1[i] * fi_acd_0 + g_yy_0_xyyz_1[i] * wa_z[i];

        g_yyz_0_xyzz_0[i] = 2.0 * g_yy_0_xyz_1[i] * fi_acd_0 + g_yy_0_xyzz_1[i] * wa_z[i];

        g_yyz_0_xzzz_0[i] = 3.0 * g_yy_0_xzz_1[i] * fi_acd_0 + g_yy_0_xzzz_1[i] * wa_z[i];

        g_yyz_0_yyyy_0[i] = g_yy_0_yyyy_1[i] * wa_z[i];

        g_yyz_0_yyyz_0[i] = g_yy_0_yyy_1[i] * fi_acd_0 + g_yy_0_yyyz_1[i] * wa_z[i];

        g_yyz_0_yyzz_0[i] = 2.0 * g_yy_0_yyz_1[i] * fi_acd_0 + g_yy_0_yyzz_1[i] * wa_z[i];

        g_yyz_0_yzzz_0[i] = 3.0 * g_yy_0_yzz_1[i] * fi_acd_0 + g_yy_0_yzzz_1[i] * wa_z[i];

        g_yyz_0_zzzz_0[i] = 4.0 * g_yy_0_zzz_1[i] * fi_acd_0 + g_yy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 120-135 components of targeted buffer : FSG

    auto g_yzz_0_xxxx_0 = pbuffer.data(idx_eri_0_fsg + 120);

    auto g_yzz_0_xxxy_0 = pbuffer.data(idx_eri_0_fsg + 121);

    auto g_yzz_0_xxxz_0 = pbuffer.data(idx_eri_0_fsg + 122);

    auto g_yzz_0_xxyy_0 = pbuffer.data(idx_eri_0_fsg + 123);

    auto g_yzz_0_xxyz_0 = pbuffer.data(idx_eri_0_fsg + 124);

    auto g_yzz_0_xxzz_0 = pbuffer.data(idx_eri_0_fsg + 125);

    auto g_yzz_0_xyyy_0 = pbuffer.data(idx_eri_0_fsg + 126);

    auto g_yzz_0_xyyz_0 = pbuffer.data(idx_eri_0_fsg + 127);

    auto g_yzz_0_xyzz_0 = pbuffer.data(idx_eri_0_fsg + 128);

    auto g_yzz_0_xzzz_0 = pbuffer.data(idx_eri_0_fsg + 129);

    auto g_yzz_0_yyyy_0 = pbuffer.data(idx_eri_0_fsg + 130);

    auto g_yzz_0_yyyz_0 = pbuffer.data(idx_eri_0_fsg + 131);

    auto g_yzz_0_yyzz_0 = pbuffer.data(idx_eri_0_fsg + 132);

    auto g_yzz_0_yzzz_0 = pbuffer.data(idx_eri_0_fsg + 133);

    auto g_yzz_0_zzzz_0 = pbuffer.data(idx_eri_0_fsg + 134);

    #pragma omp simd aligned(g_yzz_0_xxxx_0, g_yzz_0_xxxy_0, g_yzz_0_xxxz_0, g_yzz_0_xxyy_0, g_yzz_0_xxyz_0, g_yzz_0_xxzz_0, g_yzz_0_xyyy_0, g_yzz_0_xyyz_0, g_yzz_0_xyzz_0, g_yzz_0_xzzz_0, g_yzz_0_yyyy_0, g_yzz_0_yyyz_0, g_yzz_0_yyzz_0, g_yzz_0_yzzz_0, g_yzz_0_zzzz_0, g_zz_0_xxx_1, g_zz_0_xxxx_1, g_zz_0_xxxy_1, g_zz_0_xxxz_1, g_zz_0_xxy_1, g_zz_0_xxyy_1, g_zz_0_xxyz_1, g_zz_0_xxz_1, g_zz_0_xxzz_1, g_zz_0_xyy_1, g_zz_0_xyyy_1, g_zz_0_xyyz_1, g_zz_0_xyz_1, g_zz_0_xyzz_1, g_zz_0_xzz_1, g_zz_0_xzzz_1, g_zz_0_yyy_1, g_zz_0_yyyy_1, g_zz_0_yyyz_1, g_zz_0_yyz_1, g_zz_0_yyzz_1, g_zz_0_yzz_1, g_zz_0_yzzz_1, g_zz_0_zzz_1, g_zz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzz_0_xxxx_0[i] = g_zz_0_xxxx_1[i] * wa_y[i];

        g_yzz_0_xxxy_0[i] = g_zz_0_xxx_1[i] * fi_acd_0 + g_zz_0_xxxy_1[i] * wa_y[i];

        g_yzz_0_xxxz_0[i] = g_zz_0_xxxz_1[i] * wa_y[i];

        g_yzz_0_xxyy_0[i] = 2.0 * g_zz_0_xxy_1[i] * fi_acd_0 + g_zz_0_xxyy_1[i] * wa_y[i];

        g_yzz_0_xxyz_0[i] = g_zz_0_xxz_1[i] * fi_acd_0 + g_zz_0_xxyz_1[i] * wa_y[i];

        g_yzz_0_xxzz_0[i] = g_zz_0_xxzz_1[i] * wa_y[i];

        g_yzz_0_xyyy_0[i] = 3.0 * g_zz_0_xyy_1[i] * fi_acd_0 + g_zz_0_xyyy_1[i] * wa_y[i];

        g_yzz_0_xyyz_0[i] = 2.0 * g_zz_0_xyz_1[i] * fi_acd_0 + g_zz_0_xyyz_1[i] * wa_y[i];

        g_yzz_0_xyzz_0[i] = g_zz_0_xzz_1[i] * fi_acd_0 + g_zz_0_xyzz_1[i] * wa_y[i];

        g_yzz_0_xzzz_0[i] = g_zz_0_xzzz_1[i] * wa_y[i];

        g_yzz_0_yyyy_0[i] = 4.0 * g_zz_0_yyy_1[i] * fi_acd_0 + g_zz_0_yyyy_1[i] * wa_y[i];

        g_yzz_0_yyyz_0[i] = 3.0 * g_zz_0_yyz_1[i] * fi_acd_0 + g_zz_0_yyyz_1[i] * wa_y[i];

        g_yzz_0_yyzz_0[i] = 2.0 * g_zz_0_yzz_1[i] * fi_acd_0 + g_zz_0_yyzz_1[i] * wa_y[i];

        g_yzz_0_yzzz_0[i] = g_zz_0_zzz_1[i] * fi_acd_0 + g_zz_0_yzzz_1[i] * wa_y[i];

        g_yzz_0_zzzz_0[i] = g_zz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 135-150 components of targeted buffer : FSG

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

    #pragma omp simd aligned(g_z_0_xxxx_0, g_z_0_xxxx_1, g_z_0_xxxy_0, g_z_0_xxxy_1, g_z_0_xxxz_0, g_z_0_xxxz_1, g_z_0_xxyy_0, g_z_0_xxyy_1, g_z_0_xxyz_0, g_z_0_xxyz_1, g_z_0_xxzz_0, g_z_0_xxzz_1, g_z_0_xyyy_0, g_z_0_xyyy_1, g_z_0_xyyz_0, g_z_0_xyyz_1, g_z_0_xyzz_0, g_z_0_xyzz_1, g_z_0_xzzz_0, g_z_0_xzzz_1, g_z_0_yyyy_0, g_z_0_yyyy_1, g_z_0_yyyz_0, g_z_0_yyyz_1, g_z_0_yyzz_0, g_z_0_yyzz_1, g_z_0_yzzz_0, g_z_0_yzzz_1, g_z_0_zzzz_0, g_z_0_zzzz_1, g_zz_0_xxx_1, g_zz_0_xxxx_1, g_zz_0_xxxy_1, g_zz_0_xxxz_1, g_zz_0_xxy_1, g_zz_0_xxyy_1, g_zz_0_xxyz_1, g_zz_0_xxz_1, g_zz_0_xxzz_1, g_zz_0_xyy_1, g_zz_0_xyyy_1, g_zz_0_xyyz_1, g_zz_0_xyz_1, g_zz_0_xyzz_1, g_zz_0_xzz_1, g_zz_0_xzzz_1, g_zz_0_yyy_1, g_zz_0_yyyy_1, g_zz_0_yyyz_1, g_zz_0_yyz_1, g_zz_0_yyzz_1, g_zz_0_yzz_1, g_zz_0_yzzz_1, g_zz_0_zzz_1, g_zz_0_zzzz_1, g_zzz_0_xxxx_0, g_zzz_0_xxxy_0, g_zzz_0_xxxz_0, g_zzz_0_xxyy_0, g_zzz_0_xxyz_0, g_zzz_0_xxzz_0, g_zzz_0_xyyy_0, g_zzz_0_xyyz_0, g_zzz_0_xyzz_0, g_zzz_0_xzzz_0, g_zzz_0_yyyy_0, g_zzz_0_yyyz_0, g_zzz_0_yyzz_0, g_zzz_0_yzzz_0, g_zzz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzz_0_xxxx_0[i] = 2.0 * g_z_0_xxxx_0[i] * fbe_0 - 2.0 * g_z_0_xxxx_1[i] * fz_be_0 + g_zz_0_xxxx_1[i] * wa_z[i];

        g_zzz_0_xxxy_0[i] = 2.0 * g_z_0_xxxy_0[i] * fbe_0 - 2.0 * g_z_0_xxxy_1[i] * fz_be_0 + g_zz_0_xxxy_1[i] * wa_z[i];

        g_zzz_0_xxxz_0[i] = 2.0 * g_z_0_xxxz_0[i] * fbe_0 - 2.0 * g_z_0_xxxz_1[i] * fz_be_0 + g_zz_0_xxx_1[i] * fi_acd_0 + g_zz_0_xxxz_1[i] * wa_z[i];

        g_zzz_0_xxyy_0[i] = 2.0 * g_z_0_xxyy_0[i] * fbe_0 - 2.0 * g_z_0_xxyy_1[i] * fz_be_0 + g_zz_0_xxyy_1[i] * wa_z[i];

        g_zzz_0_xxyz_0[i] = 2.0 * g_z_0_xxyz_0[i] * fbe_0 - 2.0 * g_z_0_xxyz_1[i] * fz_be_0 + g_zz_0_xxy_1[i] * fi_acd_0 + g_zz_0_xxyz_1[i] * wa_z[i];

        g_zzz_0_xxzz_0[i] = 2.0 * g_z_0_xxzz_0[i] * fbe_0 - 2.0 * g_z_0_xxzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xxz_1[i] * fi_acd_0 + g_zz_0_xxzz_1[i] * wa_z[i];

        g_zzz_0_xyyy_0[i] = 2.0 * g_z_0_xyyy_0[i] * fbe_0 - 2.0 * g_z_0_xyyy_1[i] * fz_be_0 + g_zz_0_xyyy_1[i] * wa_z[i];

        g_zzz_0_xyyz_0[i] = 2.0 * g_z_0_xyyz_0[i] * fbe_0 - 2.0 * g_z_0_xyyz_1[i] * fz_be_0 + g_zz_0_xyy_1[i] * fi_acd_0 + g_zz_0_xyyz_1[i] * wa_z[i];

        g_zzz_0_xyzz_0[i] = 2.0 * g_z_0_xyzz_0[i] * fbe_0 - 2.0 * g_z_0_xyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_xyz_1[i] * fi_acd_0 + g_zz_0_xyzz_1[i] * wa_z[i];

        g_zzz_0_xzzz_0[i] = 2.0 * g_z_0_xzzz_0[i] * fbe_0 - 2.0 * g_z_0_xzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_xzz_1[i] * fi_acd_0 + g_zz_0_xzzz_1[i] * wa_z[i];

        g_zzz_0_yyyy_0[i] = 2.0 * g_z_0_yyyy_0[i] * fbe_0 - 2.0 * g_z_0_yyyy_1[i] * fz_be_0 + g_zz_0_yyyy_1[i] * wa_z[i];

        g_zzz_0_yyyz_0[i] = 2.0 * g_z_0_yyyz_0[i] * fbe_0 - 2.0 * g_z_0_yyyz_1[i] * fz_be_0 + g_zz_0_yyy_1[i] * fi_acd_0 + g_zz_0_yyyz_1[i] * wa_z[i];

        g_zzz_0_yyzz_0[i] = 2.0 * g_z_0_yyzz_0[i] * fbe_0 - 2.0 * g_z_0_yyzz_1[i] * fz_be_0 + 2.0 * g_zz_0_yyz_1[i] * fi_acd_0 + g_zz_0_yyzz_1[i] * wa_z[i];

        g_zzz_0_yzzz_0[i] = 2.0 * g_z_0_yzzz_0[i] * fbe_0 - 2.0 * g_z_0_yzzz_1[i] * fz_be_0 + 3.0 * g_zz_0_yzz_1[i] * fi_acd_0 + g_zz_0_yzzz_1[i] * wa_z[i];

        g_zzz_0_zzzz_0[i] = 2.0 * g_z_0_zzzz_0[i] * fbe_0 - 2.0 * g_z_0_zzzz_1[i] * fz_be_0 + 4.0 * g_zz_0_zzz_1[i] * fi_acd_0 + g_zz_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

