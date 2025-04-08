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

#include "ElectronRepulsionPrimRecSHSF.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_shsf(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_shsf,
                                  size_t                idx_eri_0_sfsf,
                                  size_t                idx_eri_1_sfsf,
                                  size_t                idx_eri_1_sgsd,
                                  size_t                idx_eri_0_sgsf,
                                  size_t                idx_eri_1_sgsf,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SFSF

    auto g_0_xxx_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf);

    auto g_0_xxx_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 1);

    auto g_0_xxx_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 2);

    auto g_0_xxx_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 3);

    auto g_0_xxx_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 4);

    auto g_0_xxx_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 5);

    auto g_0_xxx_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 6);

    auto g_0_xxx_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 7);

    auto g_0_xxx_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 8);

    auto g_0_xxx_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 9);

    auto g_0_xxy_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 10);

    auto g_0_xxy_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 12);

    auto g_0_xxy_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 15);

    auto g_0_xxz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 20);

    auto g_0_xxz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 21);

    auto g_0_xxz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 23);

    auto g_0_xyy_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 31);

    auto g_0_xyy_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 33);

    auto g_0_xyy_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 34);

    auto g_0_xyy_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 36);

    auto g_0_xyy_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 37);

    auto g_0_xyy_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 38);

    auto g_0_xyy_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 39);

    auto g_0_xzz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 52);

    auto g_0_xzz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 54);

    auto g_0_xzz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 55);

    auto g_0_xzz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 56);

    auto g_0_xzz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 57);

    auto g_0_xzz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 58);

    auto g_0_xzz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 59);

    auto g_0_yyy_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 60);

    auto g_0_yyy_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 61);

    auto g_0_yyy_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 62);

    auto g_0_yyy_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 63);

    auto g_0_yyy_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 64);

    auto g_0_yyy_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 65);

    auto g_0_yyy_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 66);

    auto g_0_yyy_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 67);

    auto g_0_yyy_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 68);

    auto g_0_yyy_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 69);

    auto g_0_yyz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 71);

    auto g_0_yyz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 73);

    auto g_0_yyz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 76);

    auto g_0_yzz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 80);

    auto g_0_yzz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 82);

    auto g_0_yzz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 84);

    auto g_0_yzz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 85);

    auto g_0_yzz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 87);

    auto g_0_yzz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 88);

    auto g_0_yzz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 89);

    auto g_0_zzz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 90);

    auto g_0_zzz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 91);

    auto g_0_zzz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 92);

    auto g_0_zzz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 93);

    auto g_0_zzz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 94);

    auto g_0_zzz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 95);

    auto g_0_zzz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 96);

    auto g_0_zzz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 97);

    auto g_0_zzz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 98);

    auto g_0_zzz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 99);

    /// Set up components of auxilary buffer : SFSF

    auto g_0_xxx_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf);

    auto g_0_xxx_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 1);

    auto g_0_xxx_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 2);

    auto g_0_xxx_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 3);

    auto g_0_xxx_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 4);

    auto g_0_xxx_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 5);

    auto g_0_xxx_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 6);

    auto g_0_xxx_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 7);

    auto g_0_xxx_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 8);

    auto g_0_xxx_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 9);

    auto g_0_xxy_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 10);

    auto g_0_xxy_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 12);

    auto g_0_xxy_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 15);

    auto g_0_xxz_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 20);

    auto g_0_xxz_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 21);

    auto g_0_xxz_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 23);

    auto g_0_xyy_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 31);

    auto g_0_xyy_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 33);

    auto g_0_xyy_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 34);

    auto g_0_xyy_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 36);

    auto g_0_xyy_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 37);

    auto g_0_xyy_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 38);

    auto g_0_xyy_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 39);

    auto g_0_xzz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 52);

    auto g_0_xzz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 54);

    auto g_0_xzz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 55);

    auto g_0_xzz_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 56);

    auto g_0_xzz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 57);

    auto g_0_xzz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 58);

    auto g_0_xzz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 59);

    auto g_0_yyy_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 60);

    auto g_0_yyy_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 61);

    auto g_0_yyy_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 62);

    auto g_0_yyy_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 63);

    auto g_0_yyy_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 64);

    auto g_0_yyy_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 65);

    auto g_0_yyy_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 66);

    auto g_0_yyy_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 67);

    auto g_0_yyy_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 68);

    auto g_0_yyy_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 69);

    auto g_0_yyz_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 71);

    auto g_0_yyz_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 73);

    auto g_0_yyz_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 76);

    auto g_0_yzz_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 80);

    auto g_0_yzz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 82);

    auto g_0_yzz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 84);

    auto g_0_yzz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 85);

    auto g_0_yzz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 87);

    auto g_0_yzz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 88);

    auto g_0_yzz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 89);

    auto g_0_zzz_0_xxx_1 = pbuffer.data(idx_eri_1_sfsf + 90);

    auto g_0_zzz_0_xxy_1 = pbuffer.data(idx_eri_1_sfsf + 91);

    auto g_0_zzz_0_xxz_1 = pbuffer.data(idx_eri_1_sfsf + 92);

    auto g_0_zzz_0_xyy_1 = pbuffer.data(idx_eri_1_sfsf + 93);

    auto g_0_zzz_0_xyz_1 = pbuffer.data(idx_eri_1_sfsf + 94);

    auto g_0_zzz_0_xzz_1 = pbuffer.data(idx_eri_1_sfsf + 95);

    auto g_0_zzz_0_yyy_1 = pbuffer.data(idx_eri_1_sfsf + 96);

    auto g_0_zzz_0_yyz_1 = pbuffer.data(idx_eri_1_sfsf + 97);

    auto g_0_zzz_0_yzz_1 = pbuffer.data(idx_eri_1_sfsf + 98);

    auto g_0_zzz_0_zzz_1 = pbuffer.data(idx_eri_1_sfsf + 99);

    /// Set up components of auxilary buffer : SGSD

    auto g_0_xxxx_0_xx_1 = pbuffer.data(idx_eri_1_sgsd);

    auto g_0_xxxx_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 1);

    auto g_0_xxxx_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 2);

    auto g_0_xxxx_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 3);

    auto g_0_xxxx_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 4);

    auto g_0_xxxx_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 5);

    auto g_0_xxxz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 14);

    auto g_0_xxxz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 16);

    auto g_0_xxxz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 17);

    auto g_0_xxyy_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 18);

    auto g_0_xxyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 19);

    auto g_0_xxyy_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 20);

    auto g_0_xxyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 21);

    auto g_0_xxyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 22);

    auto g_0_xxyy_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 23);

    auto g_0_xxzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 30);

    auto g_0_xxzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 31);

    auto g_0_xxzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 32);

    auto g_0_xxzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 33);

    auto g_0_xxzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 34);

    auto g_0_xxzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 35);

    auto g_0_xyyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 37);

    auto g_0_xyyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 39);

    auto g_0_xyyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 40);

    auto g_0_xzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 56);

    auto g_0_xzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 58);

    auto g_0_xzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 59);

    auto g_0_yyyy_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 60);

    auto g_0_yyyy_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 61);

    auto g_0_yyyy_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 62);

    auto g_0_yyyy_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 63);

    auto g_0_yyyy_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 64);

    auto g_0_yyyy_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 65);

    auto g_0_yyyz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 68);

    auto g_0_yyyz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 70);

    auto g_0_yyyz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 71);

    auto g_0_yyzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 72);

    auto g_0_yyzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 73);

    auto g_0_yyzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 74);

    auto g_0_yyzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 75);

    auto g_0_yyzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 76);

    auto g_0_yyzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 77);

    auto g_0_yzzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 79);

    auto g_0_yzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 80);

    auto g_0_yzzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 81);

    auto g_0_yzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 82);

    auto g_0_yzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 83);

    auto g_0_zzzz_0_xx_1 = pbuffer.data(idx_eri_1_sgsd + 84);

    auto g_0_zzzz_0_xy_1 = pbuffer.data(idx_eri_1_sgsd + 85);

    auto g_0_zzzz_0_xz_1 = pbuffer.data(idx_eri_1_sgsd + 86);

    auto g_0_zzzz_0_yy_1 = pbuffer.data(idx_eri_1_sgsd + 87);

    auto g_0_zzzz_0_yz_1 = pbuffer.data(idx_eri_1_sgsd + 88);

    auto g_0_zzzz_0_zz_1 = pbuffer.data(idx_eri_1_sgsd + 89);

    /// Set up components of auxilary buffer : SGSF

    auto g_0_xxxx_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf);

    auto g_0_xxxx_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 1);

    auto g_0_xxxx_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 2);

    auto g_0_xxxx_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 3);

    auto g_0_xxxx_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 4);

    auto g_0_xxxx_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 5);

    auto g_0_xxxx_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 6);

    auto g_0_xxxx_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 7);

    auto g_0_xxxx_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 8);

    auto g_0_xxxx_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 9);

    auto g_0_xxxy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 10);

    auto g_0_xxxy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 11);

    auto g_0_xxxy_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 12);

    auto g_0_xxxy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 13);

    auto g_0_xxxy_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 15);

    auto g_0_xxxy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 16);

    auto g_0_xxxz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 20);

    auto g_0_xxxz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 21);

    auto g_0_xxxz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 22);

    auto g_0_xxxz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 23);

    auto g_0_xxxz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 24);

    auto g_0_xxxz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 25);

    auto g_0_xxxz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 27);

    auto g_0_xxxz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 28);

    auto g_0_xxxz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 29);

    auto g_0_xxyy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 30);

    auto g_0_xxyy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 31);

    auto g_0_xxyy_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 32);

    auto g_0_xxyy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 33);

    auto g_0_xxyy_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 34);

    auto g_0_xxyy_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 35);

    auto g_0_xxyy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 36);

    auto g_0_xxyy_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 37);

    auto g_0_xxyy_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 38);

    auto g_0_xxyy_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 39);

    auto g_0_xxzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 50);

    auto g_0_xxzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 51);

    auto g_0_xxzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 52);

    auto g_0_xxzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 53);

    auto g_0_xxzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 54);

    auto g_0_xxzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 55);

    auto g_0_xxzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 56);

    auto g_0_xxzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 57);

    auto g_0_xxzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 58);

    auto g_0_xxzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 59);

    auto g_0_xyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 60);

    auto g_0_xyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 61);

    auto g_0_xyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 63);

    auto g_0_xyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 64);

    auto g_0_xyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 66);

    auto g_0_xyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 67);

    auto g_0_xyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 68);

    auto g_0_xyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 69);

    auto g_0_xzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 90);

    auto g_0_xzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 92);

    auto g_0_xzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 94);

    auto g_0_xzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 95);

    auto g_0_xzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 96);

    auto g_0_xzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 97);

    auto g_0_xzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 98);

    auto g_0_xzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 99);

    auto g_0_yyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 100);

    auto g_0_yyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 101);

    auto g_0_yyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 102);

    auto g_0_yyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 103);

    auto g_0_yyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 104);

    auto g_0_yyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 105);

    auto g_0_yyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 106);

    auto g_0_yyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 107);

    auto g_0_yyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 108);

    auto g_0_yyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 109);

    auto g_0_yyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 111);

    auto g_0_yyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 112);

    auto g_0_yyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 113);

    auto g_0_yyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 114);

    auto g_0_yyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 115);

    auto g_0_yyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 116);

    auto g_0_yyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 117);

    auto g_0_yyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 118);

    auto g_0_yyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 119);

    auto g_0_yyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 120);

    auto g_0_yyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 121);

    auto g_0_yyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 122);

    auto g_0_yyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 123);

    auto g_0_yyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 124);

    auto g_0_yyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 125);

    auto g_0_yyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 126);

    auto g_0_yyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 127);

    auto g_0_yyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 128);

    auto g_0_yyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 129);

    auto g_0_yzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 130);

    auto g_0_yzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 131);

    auto g_0_yzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 132);

    auto g_0_yzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 133);

    auto g_0_yzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 134);

    auto g_0_yzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 135);

    auto g_0_yzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 136);

    auto g_0_yzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 137);

    auto g_0_yzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 138);

    auto g_0_yzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 139);

    auto g_0_zzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 140);

    auto g_0_zzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 141);

    auto g_0_zzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 142);

    auto g_0_zzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 143);

    auto g_0_zzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 144);

    auto g_0_zzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 145);

    auto g_0_zzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 146);

    auto g_0_zzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 147);

    auto g_0_zzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 148);

    auto g_0_zzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 149);

    /// Set up components of auxilary buffer : SGSF

    auto g_0_xxxx_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf);

    auto g_0_xxxx_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 1);

    auto g_0_xxxx_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 2);

    auto g_0_xxxx_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 3);

    auto g_0_xxxx_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 4);

    auto g_0_xxxx_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 5);

    auto g_0_xxxx_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 6);

    auto g_0_xxxx_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 7);

    auto g_0_xxxx_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 8);

    auto g_0_xxxx_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 9);

    auto g_0_xxxy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 10);

    auto g_0_xxxy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 11);

    auto g_0_xxxy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 12);

    auto g_0_xxxy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 13);

    auto g_0_xxxy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 15);

    auto g_0_xxxy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 16);

    auto g_0_xxxz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 20);

    auto g_0_xxxz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 21);

    auto g_0_xxxz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 22);

    auto g_0_xxxz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 23);

    auto g_0_xxxz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 24);

    auto g_0_xxxz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 25);

    auto g_0_xxxz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 27);

    auto g_0_xxxz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 28);

    auto g_0_xxxz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 29);

    auto g_0_xxyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 30);

    auto g_0_xxyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 31);

    auto g_0_xxyy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 32);

    auto g_0_xxyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 33);

    auto g_0_xxyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 34);

    auto g_0_xxyy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 35);

    auto g_0_xxyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 36);

    auto g_0_xxyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 37);

    auto g_0_xxyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 38);

    auto g_0_xxyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 39);

    auto g_0_xxzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 50);

    auto g_0_xxzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 51);

    auto g_0_xxzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 52);

    auto g_0_xxzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 53);

    auto g_0_xxzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 54);

    auto g_0_xxzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 55);

    auto g_0_xxzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 56);

    auto g_0_xxzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 57);

    auto g_0_xxzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 58);

    auto g_0_xxzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 59);

    auto g_0_xyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 60);

    auto g_0_xyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 61);

    auto g_0_xyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 63);

    auto g_0_xyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 64);

    auto g_0_xyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 66);

    auto g_0_xyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 67);

    auto g_0_xyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 68);

    auto g_0_xyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 69);

    auto g_0_xzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 90);

    auto g_0_xzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 92);

    auto g_0_xzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 94);

    auto g_0_xzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 95);

    auto g_0_xzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 96);

    auto g_0_xzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 97);

    auto g_0_xzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 98);

    auto g_0_xzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 99);

    auto g_0_yyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 100);

    auto g_0_yyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 101);

    auto g_0_yyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 102);

    auto g_0_yyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 103);

    auto g_0_yyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 104);

    auto g_0_yyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 105);

    auto g_0_yyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 106);

    auto g_0_yyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 107);

    auto g_0_yyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 108);

    auto g_0_yyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 109);

    auto g_0_yyyz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 111);

    auto g_0_yyyz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 112);

    auto g_0_yyyz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 113);

    auto g_0_yyyz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 114);

    auto g_0_yyyz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 115);

    auto g_0_yyyz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 116);

    auto g_0_yyyz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 117);

    auto g_0_yyyz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 118);

    auto g_0_yyyz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 119);

    auto g_0_yyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 120);

    auto g_0_yyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 121);

    auto g_0_yyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 122);

    auto g_0_yyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 123);

    auto g_0_yyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 124);

    auto g_0_yyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 125);

    auto g_0_yyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 126);

    auto g_0_yyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 127);

    auto g_0_yyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 128);

    auto g_0_yyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 129);

    auto g_0_yzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 130);

    auto g_0_yzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 131);

    auto g_0_yzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 132);

    auto g_0_yzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 133);

    auto g_0_yzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 134);

    auto g_0_yzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 135);

    auto g_0_yzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 136);

    auto g_0_yzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 137);

    auto g_0_yzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 138);

    auto g_0_yzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 139);

    auto g_0_zzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 140);

    auto g_0_zzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 141);

    auto g_0_zzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 142);

    auto g_0_zzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 143);

    auto g_0_zzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 144);

    auto g_0_zzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 145);

    auto g_0_zzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 146);

    auto g_0_zzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 147);

    auto g_0_zzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 148);

    auto g_0_zzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 149);

    /// Set up 0-10 components of targeted buffer : SHSF

    auto g_0_xxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_shsf);

    auto g_0_xxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 1);

    auto g_0_xxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 2);

    auto g_0_xxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 3);

    auto g_0_xxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 4);

    auto g_0_xxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 5);

    auto g_0_xxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 6);

    auto g_0_xxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 7);

    auto g_0_xxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 8);

    auto g_0_xxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 9);

#pragma omp simd aligned(g_0_xxx_0_xxx_0,       \
                             g_0_xxx_0_xxx_1,   \
                             g_0_xxx_0_xxy_0,   \
                             g_0_xxx_0_xxy_1,   \
                             g_0_xxx_0_xxz_0,   \
                             g_0_xxx_0_xxz_1,   \
                             g_0_xxx_0_xyy_0,   \
                             g_0_xxx_0_xyy_1,   \
                             g_0_xxx_0_xyz_0,   \
                             g_0_xxx_0_xyz_1,   \
                             g_0_xxx_0_xzz_0,   \
                             g_0_xxx_0_xzz_1,   \
                             g_0_xxx_0_yyy_0,   \
                             g_0_xxx_0_yyy_1,   \
                             g_0_xxx_0_yyz_0,   \
                             g_0_xxx_0_yyz_1,   \
                             g_0_xxx_0_yzz_0,   \
                             g_0_xxx_0_yzz_1,   \
                             g_0_xxx_0_zzz_0,   \
                             g_0_xxx_0_zzz_1,   \
                             g_0_xxxx_0_xx_1,   \
                             g_0_xxxx_0_xxx_0,  \
                             g_0_xxxx_0_xxx_1,  \
                             g_0_xxxx_0_xxy_0,  \
                             g_0_xxxx_0_xxy_1,  \
                             g_0_xxxx_0_xxz_0,  \
                             g_0_xxxx_0_xxz_1,  \
                             g_0_xxxx_0_xy_1,   \
                             g_0_xxxx_0_xyy_0,  \
                             g_0_xxxx_0_xyy_1,  \
                             g_0_xxxx_0_xyz_0,  \
                             g_0_xxxx_0_xyz_1,  \
                             g_0_xxxx_0_xz_1,   \
                             g_0_xxxx_0_xzz_0,  \
                             g_0_xxxx_0_xzz_1,  \
                             g_0_xxxx_0_yy_1,   \
                             g_0_xxxx_0_yyy_0,  \
                             g_0_xxxx_0_yyy_1,  \
                             g_0_xxxx_0_yyz_0,  \
                             g_0_xxxx_0_yyz_1,  \
                             g_0_xxxx_0_yz_1,   \
                             g_0_xxxx_0_yzz_0,  \
                             g_0_xxxx_0_yzz_1,  \
                             g_0_xxxx_0_zz_1,   \
                             g_0_xxxx_0_zzz_0,  \
                             g_0_xxxx_0_zzz_1,  \
                             g_0_xxxxx_0_xxx_0, \
                             g_0_xxxxx_0_xxy_0, \
                             g_0_xxxxx_0_xxz_0, \
                             g_0_xxxxx_0_xyy_0, \
                             g_0_xxxxx_0_xyz_0, \
                             g_0_xxxxx_0_xzz_0, \
                             g_0_xxxxx_0_yyy_0, \
                             g_0_xxxxx_0_yyz_0, \
                             g_0_xxxxx_0_yzz_0, \
                             g_0_xxxxx_0_zzz_0, \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xxx_0[i] = 4.0 * g_0_xxx_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xx_1[i] * fi_abcd_0 +
                               g_0_xxxx_0_xxx_0[i] * pb_x + g_0_xxxx_0_xxx_1[i] * wp_x[i];

        g_0_xxxxx_0_xxy_0[i] = 4.0 * g_0_xxx_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xy_1[i] * fi_abcd_0 +
                               g_0_xxxx_0_xxy_0[i] * pb_x + g_0_xxxx_0_xxy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxz_0[i] = 4.0 * g_0_xxx_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xz_1[i] * fi_abcd_0 +
                               g_0_xxxx_0_xxz_0[i] * pb_x + g_0_xxxx_0_xxz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyy_0[i] = 4.0 * g_0_xxx_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxx_0_yy_1[i] * fi_abcd_0 +
                               g_0_xxxx_0_xyy_0[i] * pb_x + g_0_xxxx_0_xyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xyz_0[i] = 4.0 * g_0_xxx_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyz_1[i] * fti_ab_0 + g_0_xxxx_0_yz_1[i] * fi_abcd_0 +
                               g_0_xxxx_0_xyz_0[i] * pb_x + g_0_xxxx_0_xyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xzz_0[i] = 4.0 * g_0_xxx_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxx_0_zz_1[i] * fi_abcd_0 +
                               g_0_xxxx_0_xzz_0[i] * pb_x + g_0_xxxx_0_xzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyy_0[i] =
            4.0 * g_0_xxx_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyy_0[i] * pb_x + g_0_xxxx_0_yyy_1[i] * wp_x[i];

        g_0_xxxxx_0_yyz_0[i] =
            4.0 * g_0_xxx_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyz_0[i] * pb_x + g_0_xxxx_0_yyz_1[i] * wp_x[i];

        g_0_xxxxx_0_yzz_0[i] =
            4.0 * g_0_xxx_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzz_0[i] * pb_x + g_0_xxxx_0_yzz_1[i] * wp_x[i];

        g_0_xxxxx_0_zzz_0[i] =
            4.0 * g_0_xxx_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzz_0[i] * pb_x + g_0_xxxx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : SHSF

    auto g_0_xxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 10);

    auto g_0_xxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 11);

    auto g_0_xxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 12);

    auto g_0_xxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 13);

    auto g_0_xxxxy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 14);

    auto g_0_xxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 15);

    auto g_0_xxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 16);

    auto g_0_xxxxy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 17);

    auto g_0_xxxxy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 18);

    auto g_0_xxxxy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 19);

#pragma omp simd aligned(g_0_xxxx_0_xx_1,       \
                             g_0_xxxx_0_xxx_0,  \
                             g_0_xxxx_0_xxx_1,  \
                             g_0_xxxx_0_xxy_0,  \
                             g_0_xxxx_0_xxy_1,  \
                             g_0_xxxx_0_xxz_0,  \
                             g_0_xxxx_0_xxz_1,  \
                             g_0_xxxx_0_xy_1,   \
                             g_0_xxxx_0_xyy_0,  \
                             g_0_xxxx_0_xyy_1,  \
                             g_0_xxxx_0_xyz_0,  \
                             g_0_xxxx_0_xyz_1,  \
                             g_0_xxxx_0_xz_1,   \
                             g_0_xxxx_0_xzz_0,  \
                             g_0_xxxx_0_xzz_1,  \
                             g_0_xxxx_0_yy_1,   \
                             g_0_xxxx_0_yyy_0,  \
                             g_0_xxxx_0_yyy_1,  \
                             g_0_xxxx_0_yyz_0,  \
                             g_0_xxxx_0_yyz_1,  \
                             g_0_xxxx_0_yz_1,   \
                             g_0_xxxx_0_yzz_0,  \
                             g_0_xxxx_0_yzz_1,  \
                             g_0_xxxx_0_zz_1,   \
                             g_0_xxxx_0_zzz_0,  \
                             g_0_xxxx_0_zzz_1,  \
                             g_0_xxxxy_0_xxx_0, \
                             g_0_xxxxy_0_xxy_0, \
                             g_0_xxxxy_0_xxz_0, \
                             g_0_xxxxy_0_xyy_0, \
                             g_0_xxxxy_0_xyz_0, \
                             g_0_xxxxy_0_xzz_0, \
                             g_0_xxxxy_0_yyy_0, \
                             g_0_xxxxy_0_yyz_0, \
                             g_0_xxxxy_0_yzz_0, \
                             g_0_xxxxy_0_zzz_0, \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xxx_0[i] = g_0_xxxx_0_xxx_0[i] * pb_y + g_0_xxxx_0_xxx_1[i] * wp_y[i];

        g_0_xxxxy_0_xxy_0[i] = g_0_xxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxy_0[i] * pb_y + g_0_xxxx_0_xxy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxz_0[i] = g_0_xxxx_0_xxz_0[i] * pb_y + g_0_xxxx_0_xxz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyy_0[i] = 2.0 * g_0_xxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyy_0[i] * pb_y + g_0_xxxx_0_xyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xyz_0[i] = g_0_xxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyz_0[i] * pb_y + g_0_xxxx_0_xyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xzz_0[i] = g_0_xxxx_0_xzz_0[i] * pb_y + g_0_xxxx_0_xzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyy_0[i] = 3.0 * g_0_xxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyy_0[i] * pb_y + g_0_xxxx_0_yyy_1[i] * wp_y[i];

        g_0_xxxxy_0_yyz_0[i] = 2.0 * g_0_xxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyz_0[i] * pb_y + g_0_xxxx_0_yyz_1[i] * wp_y[i];

        g_0_xxxxy_0_yzz_0[i] = g_0_xxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzz_0[i] * pb_y + g_0_xxxx_0_yzz_1[i] * wp_y[i];

        g_0_xxxxy_0_zzz_0[i] = g_0_xxxx_0_zzz_0[i] * pb_y + g_0_xxxx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : SHSF

    auto g_0_xxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 20);

    auto g_0_xxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 21);

    auto g_0_xxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 22);

    auto g_0_xxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 23);

    auto g_0_xxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 24);

    auto g_0_xxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 25);

    auto g_0_xxxxz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 26);

    auto g_0_xxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 27);

    auto g_0_xxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 28);

    auto g_0_xxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 29);

#pragma omp simd aligned(g_0_xxxx_0_xx_1,       \
                             g_0_xxxx_0_xxx_0,  \
                             g_0_xxxx_0_xxx_1,  \
                             g_0_xxxx_0_xxy_0,  \
                             g_0_xxxx_0_xxy_1,  \
                             g_0_xxxx_0_xxz_0,  \
                             g_0_xxxx_0_xxz_1,  \
                             g_0_xxxx_0_xy_1,   \
                             g_0_xxxx_0_xyy_0,  \
                             g_0_xxxx_0_xyy_1,  \
                             g_0_xxxx_0_xyz_0,  \
                             g_0_xxxx_0_xyz_1,  \
                             g_0_xxxx_0_xz_1,   \
                             g_0_xxxx_0_xzz_0,  \
                             g_0_xxxx_0_xzz_1,  \
                             g_0_xxxx_0_yy_1,   \
                             g_0_xxxx_0_yyy_0,  \
                             g_0_xxxx_0_yyy_1,  \
                             g_0_xxxx_0_yyz_0,  \
                             g_0_xxxx_0_yyz_1,  \
                             g_0_xxxx_0_yz_1,   \
                             g_0_xxxx_0_yzz_0,  \
                             g_0_xxxx_0_yzz_1,  \
                             g_0_xxxx_0_zz_1,   \
                             g_0_xxxx_0_zzz_0,  \
                             g_0_xxxx_0_zzz_1,  \
                             g_0_xxxxz_0_xxx_0, \
                             g_0_xxxxz_0_xxy_0, \
                             g_0_xxxxz_0_xxz_0, \
                             g_0_xxxxz_0_xyy_0, \
                             g_0_xxxxz_0_xyz_0, \
                             g_0_xxxxz_0_xzz_0, \
                             g_0_xxxxz_0_yyy_0, \
                             g_0_xxxxz_0_yyz_0, \
                             g_0_xxxxz_0_yzz_0, \
                             g_0_xxxxz_0_zzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xxx_0[i] = g_0_xxxx_0_xxx_0[i] * pb_z + g_0_xxxx_0_xxx_1[i] * wp_z[i];

        g_0_xxxxz_0_xxy_0[i] = g_0_xxxx_0_xxy_0[i] * pb_z + g_0_xxxx_0_xxy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxz_0[i] = g_0_xxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxz_0[i] * pb_z + g_0_xxxx_0_xxz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyy_0[i] = g_0_xxxx_0_xyy_0[i] * pb_z + g_0_xxxx_0_xyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xyz_0[i] = g_0_xxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyz_0[i] * pb_z + g_0_xxxx_0_xyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xzz_0[i] = 2.0 * g_0_xxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzz_0[i] * pb_z + g_0_xxxx_0_xzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyy_0[i] = g_0_xxxx_0_yyy_0[i] * pb_z + g_0_xxxx_0_yyy_1[i] * wp_z[i];

        g_0_xxxxz_0_yyz_0[i] = g_0_xxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyz_0[i] * pb_z + g_0_xxxx_0_yyz_1[i] * wp_z[i];

        g_0_xxxxz_0_yzz_0[i] = 2.0 * g_0_xxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzz_0[i] * pb_z + g_0_xxxx_0_yzz_1[i] * wp_z[i];

        g_0_xxxxz_0_zzz_0[i] = 3.0 * g_0_xxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxx_0_zzz_0[i] * pb_z + g_0_xxxx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : SHSF

    auto g_0_xxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 30);

    auto g_0_xxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 31);

    auto g_0_xxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 32);

    auto g_0_xxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 33);

    auto g_0_xxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 34);

    auto g_0_xxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 35);

    auto g_0_xxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 36);

    auto g_0_xxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 37);

    auto g_0_xxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 38);

    auto g_0_xxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 39);

#pragma omp simd aligned(g_0_xxx_0_xxx_0,       \
                             g_0_xxx_0_xxx_1,   \
                             g_0_xxx_0_xxz_0,   \
                             g_0_xxx_0_xxz_1,   \
                             g_0_xxx_0_xzz_0,   \
                             g_0_xxx_0_xzz_1,   \
                             g_0_xxxy_0_xxx_0,  \
                             g_0_xxxy_0_xxx_1,  \
                             g_0_xxxy_0_xxz_0,  \
                             g_0_xxxy_0_xxz_1,  \
                             g_0_xxxy_0_xzz_0,  \
                             g_0_xxxy_0_xzz_1,  \
                             g_0_xxxyy_0_xxx_0, \
                             g_0_xxxyy_0_xxy_0, \
                             g_0_xxxyy_0_xxz_0, \
                             g_0_xxxyy_0_xyy_0, \
                             g_0_xxxyy_0_xyz_0, \
                             g_0_xxxyy_0_xzz_0, \
                             g_0_xxxyy_0_yyy_0, \
                             g_0_xxxyy_0_yyz_0, \
                             g_0_xxxyy_0_yzz_0, \
                             g_0_xxxyy_0_zzz_0, \
                             g_0_xxyy_0_xxy_0,  \
                             g_0_xxyy_0_xxy_1,  \
                             g_0_xxyy_0_xy_1,   \
                             g_0_xxyy_0_xyy_0,  \
                             g_0_xxyy_0_xyy_1,  \
                             g_0_xxyy_0_xyz_0,  \
                             g_0_xxyy_0_xyz_1,  \
                             g_0_xxyy_0_yy_1,   \
                             g_0_xxyy_0_yyy_0,  \
                             g_0_xxyy_0_yyy_1,  \
                             g_0_xxyy_0_yyz_0,  \
                             g_0_xxyy_0_yyz_1,  \
                             g_0_xxyy_0_yz_1,   \
                             g_0_xxyy_0_yzz_0,  \
                             g_0_xxyy_0_yzz_1,  \
                             g_0_xxyy_0_zzz_0,  \
                             g_0_xxyy_0_zzz_1,  \
                             g_0_xyy_0_xxy_0,   \
                             g_0_xyy_0_xxy_1,   \
                             g_0_xyy_0_xyy_0,   \
                             g_0_xyy_0_xyy_1,   \
                             g_0_xyy_0_xyz_0,   \
                             g_0_xyy_0_xyz_1,   \
                             g_0_xyy_0_yyy_0,   \
                             g_0_xyy_0_yyy_1,   \
                             g_0_xyy_0_yyz_0,   \
                             g_0_xyy_0_yyz_1,   \
                             g_0_xyy_0_yzz_0,   \
                             g_0_xyy_0_yzz_1,   \
                             g_0_xyy_0_zzz_0,   \
                             g_0_xyy_0_zzz_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xxx_0[i] =
            g_0_xxx_0_xxx_0[i] * fi_ab_0 - g_0_xxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxy_0_xxx_0[i] * pb_y + g_0_xxxy_0_xxx_1[i] * wp_y[i];

        g_0_xxxyy_0_xxy_0[i] = 2.0 * g_0_xyy_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xy_1[i] * fi_abcd_0 +
                               g_0_xxyy_0_xxy_0[i] * pb_x + g_0_xxyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxz_0[i] =
            g_0_xxx_0_xxz_0[i] * fi_ab_0 - g_0_xxx_0_xxz_1[i] * fti_ab_0 + g_0_xxxy_0_xxz_0[i] * pb_y + g_0_xxxy_0_xxz_1[i] * wp_y[i];

        g_0_xxxyy_0_xyy_0[i] = 2.0 * g_0_xyy_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyy_0_yy_1[i] * fi_abcd_0 +
                               g_0_xxyy_0_xyy_0[i] * pb_x + g_0_xxyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xyz_0[i] = 2.0 * g_0_xyy_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyz_1[i] * fti_ab_0 + g_0_xxyy_0_yz_1[i] * fi_abcd_0 +
                               g_0_xxyy_0_xyz_0[i] * pb_x + g_0_xxyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xzz_0[i] =
            g_0_xxx_0_xzz_0[i] * fi_ab_0 - g_0_xxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxy_0_xzz_0[i] * pb_y + g_0_xxxy_0_xzz_1[i] * wp_y[i];

        g_0_xxxyy_0_yyy_0[i] =
            2.0 * g_0_xyy_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyy_0[i] * pb_x + g_0_xxyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxyy_0_yyz_0[i] =
            2.0 * g_0_xyy_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyz_0[i] * pb_x + g_0_xxyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxyy_0_yzz_0[i] =
            2.0 * g_0_xyy_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzz_0[i] * pb_x + g_0_xxyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxyy_0_zzz_0[i] =
            2.0 * g_0_xyy_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zzz_1[i] * fti_ab_0 + g_0_xxyy_0_zzz_0[i] * pb_x + g_0_xxyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : SHSF

    auto g_0_xxxyz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 40);

    auto g_0_xxxyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 41);

    auto g_0_xxxyz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 42);

    auto g_0_xxxyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 43);

    auto g_0_xxxyz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 44);

    auto g_0_xxxyz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 45);

    auto g_0_xxxyz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 46);

    auto g_0_xxxyz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 47);

    auto g_0_xxxyz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 48);

    auto g_0_xxxyz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 49);

#pragma omp simd aligned(g_0_xxxy_0_xxy_0,      \
                             g_0_xxxy_0_xxy_1,  \
                             g_0_xxxy_0_xyy_0,  \
                             g_0_xxxy_0_xyy_1,  \
                             g_0_xxxy_0_yyy_0,  \
                             g_0_xxxy_0_yyy_1,  \
                             g_0_xxxyz_0_xxx_0, \
                             g_0_xxxyz_0_xxy_0, \
                             g_0_xxxyz_0_xxz_0, \
                             g_0_xxxyz_0_xyy_0, \
                             g_0_xxxyz_0_xyz_0, \
                             g_0_xxxyz_0_xzz_0, \
                             g_0_xxxyz_0_yyy_0, \
                             g_0_xxxyz_0_yyz_0, \
                             g_0_xxxyz_0_yzz_0, \
                             g_0_xxxyz_0_zzz_0, \
                             g_0_xxxz_0_xxx_0,  \
                             g_0_xxxz_0_xxx_1,  \
                             g_0_xxxz_0_xxz_0,  \
                             g_0_xxxz_0_xxz_1,  \
                             g_0_xxxz_0_xyz_0,  \
                             g_0_xxxz_0_xyz_1,  \
                             g_0_xxxz_0_xz_1,   \
                             g_0_xxxz_0_xzz_0,  \
                             g_0_xxxz_0_xzz_1,  \
                             g_0_xxxz_0_yyz_0,  \
                             g_0_xxxz_0_yyz_1,  \
                             g_0_xxxz_0_yz_1,   \
                             g_0_xxxz_0_yzz_0,  \
                             g_0_xxxz_0_yzz_1,  \
                             g_0_xxxz_0_zz_1,   \
                             g_0_xxxz_0_zzz_0,  \
                             g_0_xxxz_0_zzz_1,  \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xxx_0[i] = g_0_xxxz_0_xxx_0[i] * pb_y + g_0_xxxz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyz_0_xxy_0[i] = g_0_xxxy_0_xxy_0[i] * pb_z + g_0_xxxy_0_xxy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxz_0[i] = g_0_xxxz_0_xxz_0[i] * pb_y + g_0_xxxz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyy_0[i] = g_0_xxxy_0_xyy_0[i] * pb_z + g_0_xxxy_0_xyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xyz_0[i] = g_0_xxxz_0_xz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyz_0[i] * pb_y + g_0_xxxz_0_xyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xzz_0[i] = g_0_xxxz_0_xzz_0[i] * pb_y + g_0_xxxz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyy_0[i] = g_0_xxxy_0_yyy_0[i] * pb_z + g_0_xxxy_0_yyy_1[i] * wp_z[i];

        g_0_xxxyz_0_yyz_0[i] = 2.0 * g_0_xxxz_0_yz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyz_0[i] * pb_y + g_0_xxxz_0_yyz_1[i] * wp_y[i];

        g_0_xxxyz_0_yzz_0[i] = g_0_xxxz_0_zz_1[i] * fi_abcd_0 + g_0_xxxz_0_yzz_0[i] * pb_y + g_0_xxxz_0_yzz_1[i] * wp_y[i];

        g_0_xxxyz_0_zzz_0[i] = g_0_xxxz_0_zzz_0[i] * pb_y + g_0_xxxz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 50-60 components of targeted buffer : SHSF

    auto g_0_xxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 50);

    auto g_0_xxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 51);

    auto g_0_xxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 52);

    auto g_0_xxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 53);

    auto g_0_xxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 54);

    auto g_0_xxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 55);

    auto g_0_xxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 56);

    auto g_0_xxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 57);

    auto g_0_xxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 58);

    auto g_0_xxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 59);

#pragma omp simd aligned(g_0_xxx_0_xxx_0,       \
                             g_0_xxx_0_xxx_1,   \
                             g_0_xxx_0_xxy_0,   \
                             g_0_xxx_0_xxy_1,   \
                             g_0_xxx_0_xyy_0,   \
                             g_0_xxx_0_xyy_1,   \
                             g_0_xxxz_0_xxx_0,  \
                             g_0_xxxz_0_xxx_1,  \
                             g_0_xxxz_0_xxy_0,  \
                             g_0_xxxz_0_xxy_1,  \
                             g_0_xxxz_0_xyy_0,  \
                             g_0_xxxz_0_xyy_1,  \
                             g_0_xxxzz_0_xxx_0, \
                             g_0_xxxzz_0_xxy_0, \
                             g_0_xxxzz_0_xxz_0, \
                             g_0_xxxzz_0_xyy_0, \
                             g_0_xxxzz_0_xyz_0, \
                             g_0_xxxzz_0_xzz_0, \
                             g_0_xxxzz_0_yyy_0, \
                             g_0_xxxzz_0_yyz_0, \
                             g_0_xxxzz_0_yzz_0, \
                             g_0_xxxzz_0_zzz_0, \
                             g_0_xxzz_0_xxz_0,  \
                             g_0_xxzz_0_xxz_1,  \
                             g_0_xxzz_0_xyz_0,  \
                             g_0_xxzz_0_xyz_1,  \
                             g_0_xxzz_0_xz_1,   \
                             g_0_xxzz_0_xzz_0,  \
                             g_0_xxzz_0_xzz_1,  \
                             g_0_xxzz_0_yyy_0,  \
                             g_0_xxzz_0_yyy_1,  \
                             g_0_xxzz_0_yyz_0,  \
                             g_0_xxzz_0_yyz_1,  \
                             g_0_xxzz_0_yz_1,   \
                             g_0_xxzz_0_yzz_0,  \
                             g_0_xxzz_0_yzz_1,  \
                             g_0_xxzz_0_zz_1,   \
                             g_0_xxzz_0_zzz_0,  \
                             g_0_xxzz_0_zzz_1,  \
                             g_0_xzz_0_xxz_0,   \
                             g_0_xzz_0_xxz_1,   \
                             g_0_xzz_0_xyz_0,   \
                             g_0_xzz_0_xyz_1,   \
                             g_0_xzz_0_xzz_0,   \
                             g_0_xzz_0_xzz_1,   \
                             g_0_xzz_0_yyy_0,   \
                             g_0_xzz_0_yyy_1,   \
                             g_0_xzz_0_yyz_0,   \
                             g_0_xzz_0_yyz_1,   \
                             g_0_xzz_0_yzz_0,   \
                             g_0_xzz_0_yzz_1,   \
                             g_0_xzz_0_zzz_0,   \
                             g_0_xzz_0_zzz_1,   \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xxx_0[i] =
            g_0_xxx_0_xxx_0[i] * fi_ab_0 - g_0_xxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxz_0_xxx_0[i] * pb_z + g_0_xxxz_0_xxx_1[i] * wp_z[i];

        g_0_xxxzz_0_xxy_0[i] =
            g_0_xxx_0_xxy_0[i] * fi_ab_0 - g_0_xxx_0_xxy_1[i] * fti_ab_0 + g_0_xxxz_0_xxy_0[i] * pb_z + g_0_xxxz_0_xxy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxz_0[i] = 2.0 * g_0_xzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xz_1[i] * fi_abcd_0 +
                               g_0_xxzz_0_xxz_0[i] * pb_x + g_0_xxzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyy_0[i] =
            g_0_xxx_0_xyy_0[i] * fi_ab_0 - g_0_xxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxz_0_xyy_0[i] * pb_z + g_0_xxxz_0_xyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xyz_0[i] = 2.0 * g_0_xzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyz_1[i] * fti_ab_0 + g_0_xxzz_0_yz_1[i] * fi_abcd_0 +
                               g_0_xxzz_0_xyz_0[i] * pb_x + g_0_xxzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xzz_0[i] = 2.0 * g_0_xzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xzz_1[i] * fti_ab_0 + g_0_xxzz_0_zz_1[i] * fi_abcd_0 +
                               g_0_xxzz_0_xzz_0[i] * pb_x + g_0_xxzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyy_0[i] =
            2.0 * g_0_xzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyy_1[i] * fti_ab_0 + g_0_xxzz_0_yyy_0[i] * pb_x + g_0_xxzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxzz_0_yyz_0[i] =
            2.0 * g_0_xzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyz_0[i] * pb_x + g_0_xxzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxzz_0_yzz_0[i] =
            2.0 * g_0_xzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzz_0[i] * pb_x + g_0_xxzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxzz_0_zzz_0[i] =
            2.0 * g_0_xzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzz_0[i] * pb_x + g_0_xxzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : SHSF

    auto g_0_xxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 60);

    auto g_0_xxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 61);

    auto g_0_xxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 62);

    auto g_0_xxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 63);

    auto g_0_xxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 64);

    auto g_0_xxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 65);

    auto g_0_xxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 66);

    auto g_0_xxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 67);

    auto g_0_xxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 68);

    auto g_0_xxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 69);

#pragma omp simd aligned(g_0_xxy_0_xxx_0,       \
                             g_0_xxy_0_xxx_1,   \
                             g_0_xxy_0_xxz_0,   \
                             g_0_xxy_0_xxz_1,   \
                             g_0_xxy_0_xzz_0,   \
                             g_0_xxy_0_xzz_1,   \
                             g_0_xxyy_0_xxx_0,  \
                             g_0_xxyy_0_xxx_1,  \
                             g_0_xxyy_0_xxz_0,  \
                             g_0_xxyy_0_xxz_1,  \
                             g_0_xxyy_0_xzz_0,  \
                             g_0_xxyy_0_xzz_1,  \
                             g_0_xxyyy_0_xxx_0, \
                             g_0_xxyyy_0_xxy_0, \
                             g_0_xxyyy_0_xxz_0, \
                             g_0_xxyyy_0_xyy_0, \
                             g_0_xxyyy_0_xyz_0, \
                             g_0_xxyyy_0_xzz_0, \
                             g_0_xxyyy_0_yyy_0, \
                             g_0_xxyyy_0_yyz_0, \
                             g_0_xxyyy_0_yzz_0, \
                             g_0_xxyyy_0_zzz_0, \
                             g_0_xyyy_0_xxy_0,  \
                             g_0_xyyy_0_xxy_1,  \
                             g_0_xyyy_0_xy_1,   \
                             g_0_xyyy_0_xyy_0,  \
                             g_0_xyyy_0_xyy_1,  \
                             g_0_xyyy_0_xyz_0,  \
                             g_0_xyyy_0_xyz_1,  \
                             g_0_xyyy_0_yy_1,   \
                             g_0_xyyy_0_yyy_0,  \
                             g_0_xyyy_0_yyy_1,  \
                             g_0_xyyy_0_yyz_0,  \
                             g_0_xyyy_0_yyz_1,  \
                             g_0_xyyy_0_yz_1,   \
                             g_0_xyyy_0_yzz_0,  \
                             g_0_xyyy_0_yzz_1,  \
                             g_0_xyyy_0_zzz_0,  \
                             g_0_xyyy_0_zzz_1,  \
                             g_0_yyy_0_xxy_0,   \
                             g_0_yyy_0_xxy_1,   \
                             g_0_yyy_0_xyy_0,   \
                             g_0_yyy_0_xyy_1,   \
                             g_0_yyy_0_xyz_0,   \
                             g_0_yyy_0_xyz_1,   \
                             g_0_yyy_0_yyy_0,   \
                             g_0_yyy_0_yyy_1,   \
                             g_0_yyy_0_yyz_0,   \
                             g_0_yyy_0_yyz_1,   \
                             g_0_yyy_0_yzz_0,   \
                             g_0_yyy_0_yzz_1,   \
                             g_0_yyy_0_zzz_0,   \
                             g_0_yyy_0_zzz_1,   \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xxx_0[i] =
            2.0 * g_0_xxy_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxx_1[i] * fti_ab_0 + g_0_xxyy_0_xxx_0[i] * pb_y + g_0_xxyy_0_xxx_1[i] * wp_y[i];

        g_0_xxyyy_0_xxy_0[i] = g_0_yyy_0_xxy_0[i] * fi_ab_0 - g_0_yyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xy_1[i] * fi_abcd_0 +
                               g_0_xyyy_0_xxy_0[i] * pb_x + g_0_xyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxz_0[i] =
            2.0 * g_0_xxy_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxz_1[i] * fti_ab_0 + g_0_xxyy_0_xxz_0[i] * pb_y + g_0_xxyy_0_xxz_1[i] * wp_y[i];

        g_0_xxyyy_0_xyy_0[i] = g_0_yyy_0_xyy_0[i] * fi_ab_0 - g_0_yyy_0_xyy_1[i] * fti_ab_0 + g_0_xyyy_0_yy_1[i] * fi_abcd_0 +
                               g_0_xyyy_0_xyy_0[i] * pb_x + g_0_xyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xyz_0[i] = g_0_yyy_0_xyz_0[i] * fi_ab_0 - g_0_yyy_0_xyz_1[i] * fti_ab_0 + g_0_xyyy_0_yz_1[i] * fi_abcd_0 +
                               g_0_xyyy_0_xyz_0[i] * pb_x + g_0_xyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xzz_0[i] =
            2.0 * g_0_xxy_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xzz_1[i] * fti_ab_0 + g_0_xxyy_0_xzz_0[i] * pb_y + g_0_xxyy_0_xzz_1[i] * wp_y[i];

        g_0_xxyyy_0_yyy_0[i] =
            g_0_yyy_0_yyy_0[i] * fi_ab_0 - g_0_yyy_0_yyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyy_0[i] * pb_x + g_0_xyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxyyy_0_yyz_0[i] =
            g_0_yyy_0_yyz_0[i] * fi_ab_0 - g_0_yyy_0_yyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyz_0[i] * pb_x + g_0_xyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxyyy_0_yzz_0[i] =
            g_0_yyy_0_yzz_0[i] * fi_ab_0 - g_0_yyy_0_yzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzz_0[i] * pb_x + g_0_xyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxyyy_0_zzz_0[i] =
            g_0_yyy_0_zzz_0[i] * fi_ab_0 - g_0_yyy_0_zzz_1[i] * fti_ab_0 + g_0_xyyy_0_zzz_0[i] * pb_x + g_0_xyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 70-80 components of targeted buffer : SHSF

    auto g_0_xxyyz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 70);

    auto g_0_xxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 71);

    auto g_0_xxyyz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 72);

    auto g_0_xxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 73);

    auto g_0_xxyyz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 74);

    auto g_0_xxyyz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 75);

    auto g_0_xxyyz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 76);

    auto g_0_xxyyz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 77);

    auto g_0_xxyyz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 78);

    auto g_0_xxyyz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 79);

#pragma omp simd aligned(g_0_xxyy_0_xx_1,       \
                             g_0_xxyy_0_xxx_0,  \
                             g_0_xxyy_0_xxx_1,  \
                             g_0_xxyy_0_xxy_0,  \
                             g_0_xxyy_0_xxy_1,  \
                             g_0_xxyy_0_xxz_0,  \
                             g_0_xxyy_0_xxz_1,  \
                             g_0_xxyy_0_xy_1,   \
                             g_0_xxyy_0_xyy_0,  \
                             g_0_xxyy_0_xyy_1,  \
                             g_0_xxyy_0_xyz_0,  \
                             g_0_xxyy_0_xyz_1,  \
                             g_0_xxyy_0_xz_1,   \
                             g_0_xxyy_0_xzz_0,  \
                             g_0_xxyy_0_xzz_1,  \
                             g_0_xxyy_0_yy_1,   \
                             g_0_xxyy_0_yyy_0,  \
                             g_0_xxyy_0_yyy_1,  \
                             g_0_xxyy_0_yyz_0,  \
                             g_0_xxyy_0_yyz_1,  \
                             g_0_xxyy_0_yz_1,   \
                             g_0_xxyy_0_yzz_0,  \
                             g_0_xxyy_0_yzz_1,  \
                             g_0_xxyy_0_zz_1,   \
                             g_0_xxyy_0_zzz_0,  \
                             g_0_xxyy_0_zzz_1,  \
                             g_0_xxyyz_0_xxx_0, \
                             g_0_xxyyz_0_xxy_0, \
                             g_0_xxyyz_0_xxz_0, \
                             g_0_xxyyz_0_xyy_0, \
                             g_0_xxyyz_0_xyz_0, \
                             g_0_xxyyz_0_xzz_0, \
                             g_0_xxyyz_0_yyy_0, \
                             g_0_xxyyz_0_yyz_0, \
                             g_0_xxyyz_0_yzz_0, \
                             g_0_xxyyz_0_zzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xxx_0[i] = g_0_xxyy_0_xxx_0[i] * pb_z + g_0_xxyy_0_xxx_1[i] * wp_z[i];

        g_0_xxyyz_0_xxy_0[i] = g_0_xxyy_0_xxy_0[i] * pb_z + g_0_xxyy_0_xxy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxz_0[i] = g_0_xxyy_0_xx_1[i] * fi_abcd_0 + g_0_xxyy_0_xxz_0[i] * pb_z + g_0_xxyy_0_xxz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyy_0[i] = g_0_xxyy_0_xyy_0[i] * pb_z + g_0_xxyy_0_xyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xyz_0[i] = g_0_xxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyz_0[i] * pb_z + g_0_xxyy_0_xyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xzz_0[i] = 2.0 * g_0_xxyy_0_xz_1[i] * fi_abcd_0 + g_0_xxyy_0_xzz_0[i] * pb_z + g_0_xxyy_0_xzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyy_0[i] = g_0_xxyy_0_yyy_0[i] * pb_z + g_0_xxyy_0_yyy_1[i] * wp_z[i];

        g_0_xxyyz_0_yyz_0[i] = g_0_xxyy_0_yy_1[i] * fi_abcd_0 + g_0_xxyy_0_yyz_0[i] * pb_z + g_0_xxyy_0_yyz_1[i] * wp_z[i];

        g_0_xxyyz_0_yzz_0[i] = 2.0 * g_0_xxyy_0_yz_1[i] * fi_abcd_0 + g_0_xxyy_0_yzz_0[i] * pb_z + g_0_xxyy_0_yzz_1[i] * wp_z[i];

        g_0_xxyyz_0_zzz_0[i] = 3.0 * g_0_xxyy_0_zz_1[i] * fi_abcd_0 + g_0_xxyy_0_zzz_0[i] * pb_z + g_0_xxyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 80-90 components of targeted buffer : SHSF

    auto g_0_xxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 80);

    auto g_0_xxyzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 81);

    auto g_0_xxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 82);

    auto g_0_xxyzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 83);

    auto g_0_xxyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 84);

    auto g_0_xxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 85);

    auto g_0_xxyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 86);

    auto g_0_xxyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 87);

    auto g_0_xxyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 88);

    auto g_0_xxyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 89);

#pragma omp simd aligned(g_0_xxyzz_0_xxx_0,     \
                             g_0_xxyzz_0_xxy_0, \
                             g_0_xxyzz_0_xxz_0, \
                             g_0_xxyzz_0_xyy_0, \
                             g_0_xxyzz_0_xyz_0, \
                             g_0_xxyzz_0_xzz_0, \
                             g_0_xxyzz_0_yyy_0, \
                             g_0_xxyzz_0_yyz_0, \
                             g_0_xxyzz_0_yzz_0, \
                             g_0_xxyzz_0_zzz_0, \
                             g_0_xxzz_0_xx_1,   \
                             g_0_xxzz_0_xxx_0,  \
                             g_0_xxzz_0_xxx_1,  \
                             g_0_xxzz_0_xxy_0,  \
                             g_0_xxzz_0_xxy_1,  \
                             g_0_xxzz_0_xxz_0,  \
                             g_0_xxzz_0_xxz_1,  \
                             g_0_xxzz_0_xy_1,   \
                             g_0_xxzz_0_xyy_0,  \
                             g_0_xxzz_0_xyy_1,  \
                             g_0_xxzz_0_xyz_0,  \
                             g_0_xxzz_0_xyz_1,  \
                             g_0_xxzz_0_xz_1,   \
                             g_0_xxzz_0_xzz_0,  \
                             g_0_xxzz_0_xzz_1,  \
                             g_0_xxzz_0_yy_1,   \
                             g_0_xxzz_0_yyy_0,  \
                             g_0_xxzz_0_yyy_1,  \
                             g_0_xxzz_0_yyz_0,  \
                             g_0_xxzz_0_yyz_1,  \
                             g_0_xxzz_0_yz_1,   \
                             g_0_xxzz_0_yzz_0,  \
                             g_0_xxzz_0_yzz_1,  \
                             g_0_xxzz_0_zz_1,   \
                             g_0_xxzz_0_zzz_0,  \
                             g_0_xxzz_0_zzz_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xxx_0[i] = g_0_xxzz_0_xxx_0[i] * pb_y + g_0_xxzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyzz_0_xxy_0[i] = g_0_xxzz_0_xx_1[i] * fi_abcd_0 + g_0_xxzz_0_xxy_0[i] * pb_y + g_0_xxzz_0_xxy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxz_0[i] = g_0_xxzz_0_xxz_0[i] * pb_y + g_0_xxzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyy_0[i] = 2.0 * g_0_xxzz_0_xy_1[i] * fi_abcd_0 + g_0_xxzz_0_xyy_0[i] * pb_y + g_0_xxzz_0_xyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xyz_0[i] = g_0_xxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyz_0[i] * pb_y + g_0_xxzz_0_xyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xzz_0[i] = g_0_xxzz_0_xzz_0[i] * pb_y + g_0_xxzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyy_0[i] = 3.0 * g_0_xxzz_0_yy_1[i] * fi_abcd_0 + g_0_xxzz_0_yyy_0[i] * pb_y + g_0_xxzz_0_yyy_1[i] * wp_y[i];

        g_0_xxyzz_0_yyz_0[i] = 2.0 * g_0_xxzz_0_yz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyz_0[i] * pb_y + g_0_xxzz_0_yyz_1[i] * wp_y[i];

        g_0_xxyzz_0_yzz_0[i] = g_0_xxzz_0_zz_1[i] * fi_abcd_0 + g_0_xxzz_0_yzz_0[i] * pb_y + g_0_xxzz_0_yzz_1[i] * wp_y[i];

        g_0_xxyzz_0_zzz_0[i] = g_0_xxzz_0_zzz_0[i] * pb_y + g_0_xxzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 90-100 components of targeted buffer : SHSF

    auto g_0_xxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 90);

    auto g_0_xxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 91);

    auto g_0_xxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 92);

    auto g_0_xxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 93);

    auto g_0_xxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 94);

    auto g_0_xxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 95);

    auto g_0_xxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 96);

    auto g_0_xxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 97);

    auto g_0_xxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 98);

    auto g_0_xxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 99);

#pragma omp simd aligned(g_0_xxz_0_xxx_0,       \
                             g_0_xxz_0_xxx_1,   \
                             g_0_xxz_0_xxy_0,   \
                             g_0_xxz_0_xxy_1,   \
                             g_0_xxz_0_xyy_0,   \
                             g_0_xxz_0_xyy_1,   \
                             g_0_xxzz_0_xxx_0,  \
                             g_0_xxzz_0_xxx_1,  \
                             g_0_xxzz_0_xxy_0,  \
                             g_0_xxzz_0_xxy_1,  \
                             g_0_xxzz_0_xyy_0,  \
                             g_0_xxzz_0_xyy_1,  \
                             g_0_xxzzz_0_xxx_0, \
                             g_0_xxzzz_0_xxy_0, \
                             g_0_xxzzz_0_xxz_0, \
                             g_0_xxzzz_0_xyy_0, \
                             g_0_xxzzz_0_xyz_0, \
                             g_0_xxzzz_0_xzz_0, \
                             g_0_xxzzz_0_yyy_0, \
                             g_0_xxzzz_0_yyz_0, \
                             g_0_xxzzz_0_yzz_0, \
                             g_0_xxzzz_0_zzz_0, \
                             g_0_xzzz_0_xxz_0,  \
                             g_0_xzzz_0_xxz_1,  \
                             g_0_xzzz_0_xyz_0,  \
                             g_0_xzzz_0_xyz_1,  \
                             g_0_xzzz_0_xz_1,   \
                             g_0_xzzz_0_xzz_0,  \
                             g_0_xzzz_0_xzz_1,  \
                             g_0_xzzz_0_yyy_0,  \
                             g_0_xzzz_0_yyy_1,  \
                             g_0_xzzz_0_yyz_0,  \
                             g_0_xzzz_0_yyz_1,  \
                             g_0_xzzz_0_yz_1,   \
                             g_0_xzzz_0_yzz_0,  \
                             g_0_xzzz_0_yzz_1,  \
                             g_0_xzzz_0_zz_1,   \
                             g_0_xzzz_0_zzz_0,  \
                             g_0_xzzz_0_zzz_1,  \
                             g_0_zzz_0_xxz_0,   \
                             g_0_zzz_0_xxz_1,   \
                             g_0_zzz_0_xyz_0,   \
                             g_0_zzz_0_xyz_1,   \
                             g_0_zzz_0_xzz_0,   \
                             g_0_zzz_0_xzz_1,   \
                             g_0_zzz_0_yyy_0,   \
                             g_0_zzz_0_yyy_1,   \
                             g_0_zzz_0_yyz_0,   \
                             g_0_zzz_0_yyz_1,   \
                             g_0_zzz_0_yzz_0,   \
                             g_0_zzz_0_yzz_1,   \
                             g_0_zzz_0_zzz_0,   \
                             g_0_zzz_0_zzz_1,   \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xxx_0[i] =
            2.0 * g_0_xxz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxx_1[i] * fti_ab_0 + g_0_xxzz_0_xxx_0[i] * pb_z + g_0_xxzz_0_xxx_1[i] * wp_z[i];

        g_0_xxzzz_0_xxy_0[i] =
            2.0 * g_0_xxz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxy_1[i] * fti_ab_0 + g_0_xxzz_0_xxy_0[i] * pb_z + g_0_xxzz_0_xxy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxz_0[i] = g_0_zzz_0_xxz_0[i] * fi_ab_0 - g_0_zzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xz_1[i] * fi_abcd_0 +
                               g_0_xzzz_0_xxz_0[i] * pb_x + g_0_xzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyy_0[i] =
            2.0 * g_0_xxz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xyy_1[i] * fti_ab_0 + g_0_xxzz_0_xyy_0[i] * pb_z + g_0_xxzz_0_xyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xyz_0[i] = g_0_zzz_0_xyz_0[i] * fi_ab_0 - g_0_zzz_0_xyz_1[i] * fti_ab_0 + g_0_xzzz_0_yz_1[i] * fi_abcd_0 +
                               g_0_xzzz_0_xyz_0[i] * pb_x + g_0_xzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xzz_0[i] = g_0_zzz_0_xzz_0[i] * fi_ab_0 - g_0_zzz_0_xzz_1[i] * fti_ab_0 + g_0_xzzz_0_zz_1[i] * fi_abcd_0 +
                               g_0_xzzz_0_xzz_0[i] * pb_x + g_0_xzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyy_0[i] =
            g_0_zzz_0_yyy_0[i] * fi_ab_0 - g_0_zzz_0_yyy_1[i] * fti_ab_0 + g_0_xzzz_0_yyy_0[i] * pb_x + g_0_xzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxzzz_0_yyz_0[i] =
            g_0_zzz_0_yyz_0[i] * fi_ab_0 - g_0_zzz_0_yyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyz_0[i] * pb_x + g_0_xzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxzzz_0_yzz_0[i] =
            g_0_zzz_0_yzz_0[i] * fi_ab_0 - g_0_zzz_0_yzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzz_0[i] * pb_x + g_0_xzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxzzz_0_zzz_0[i] =
            g_0_zzz_0_zzz_0[i] * fi_ab_0 - g_0_zzz_0_zzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzz_0[i] * pb_x + g_0_xzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 100-110 components of targeted buffer : SHSF

    auto g_0_xyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 100);

    auto g_0_xyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 101);

    auto g_0_xyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 102);

    auto g_0_xyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 103);

    auto g_0_xyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 104);

    auto g_0_xyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 105);

    auto g_0_xyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 106);

    auto g_0_xyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 107);

    auto g_0_xyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 108);

    auto g_0_xyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 109);

#pragma omp simd aligned(g_0_xyyyy_0_xxx_0,     \
                             g_0_xyyyy_0_xxy_0, \
                             g_0_xyyyy_0_xxz_0, \
                             g_0_xyyyy_0_xyy_0, \
                             g_0_xyyyy_0_xyz_0, \
                             g_0_xyyyy_0_xzz_0, \
                             g_0_xyyyy_0_yyy_0, \
                             g_0_xyyyy_0_yyz_0, \
                             g_0_xyyyy_0_yzz_0, \
                             g_0_xyyyy_0_zzz_0, \
                             g_0_yyyy_0_xx_1,   \
                             g_0_yyyy_0_xxx_0,  \
                             g_0_yyyy_0_xxx_1,  \
                             g_0_yyyy_0_xxy_0,  \
                             g_0_yyyy_0_xxy_1,  \
                             g_0_yyyy_0_xxz_0,  \
                             g_0_yyyy_0_xxz_1,  \
                             g_0_yyyy_0_xy_1,   \
                             g_0_yyyy_0_xyy_0,  \
                             g_0_yyyy_0_xyy_1,  \
                             g_0_yyyy_0_xyz_0,  \
                             g_0_yyyy_0_xyz_1,  \
                             g_0_yyyy_0_xz_1,   \
                             g_0_yyyy_0_xzz_0,  \
                             g_0_yyyy_0_xzz_1,  \
                             g_0_yyyy_0_yy_1,   \
                             g_0_yyyy_0_yyy_0,  \
                             g_0_yyyy_0_yyy_1,  \
                             g_0_yyyy_0_yyz_0,  \
                             g_0_yyyy_0_yyz_1,  \
                             g_0_yyyy_0_yz_1,   \
                             g_0_yyyy_0_yzz_0,  \
                             g_0_yyyy_0_yzz_1,  \
                             g_0_yyyy_0_zz_1,   \
                             g_0_yyyy_0_zzz_0,  \
                             g_0_yyyy_0_zzz_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xxx_0[i] = 3.0 * g_0_yyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxx_0[i] * pb_x + g_0_yyyy_0_xxx_1[i] * wp_x[i];

        g_0_xyyyy_0_xxy_0[i] = 2.0 * g_0_yyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxy_0[i] * pb_x + g_0_yyyy_0_xxy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxz_0[i] = 2.0 * g_0_yyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxz_0[i] * pb_x + g_0_yyyy_0_xxz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyy_0[i] = g_0_yyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyy_0[i] * pb_x + g_0_yyyy_0_xyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xyz_0[i] = g_0_yyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyz_0[i] * pb_x + g_0_yyyy_0_xyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xzz_0[i] = g_0_yyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzz_0[i] * pb_x + g_0_yyyy_0_xzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyy_0[i] = g_0_yyyy_0_yyy_0[i] * pb_x + g_0_yyyy_0_yyy_1[i] * wp_x[i];

        g_0_xyyyy_0_yyz_0[i] = g_0_yyyy_0_yyz_0[i] * pb_x + g_0_yyyy_0_yyz_1[i] * wp_x[i];

        g_0_xyyyy_0_yzz_0[i] = g_0_yyyy_0_yzz_0[i] * pb_x + g_0_yyyy_0_yzz_1[i] * wp_x[i];

        g_0_xyyyy_0_zzz_0[i] = g_0_yyyy_0_zzz_0[i] * pb_x + g_0_yyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 110-120 components of targeted buffer : SHSF

    auto g_0_xyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 110);

    auto g_0_xyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 111);

    auto g_0_xyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 112);

    auto g_0_xyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 113);

    auto g_0_xyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 114);

    auto g_0_xyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 115);

    auto g_0_xyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 116);

    auto g_0_xyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 117);

    auto g_0_xyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 118);

    auto g_0_xyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 119);

#pragma omp simd aligned(g_0_xyyy_0_xxx_0,      \
                             g_0_xyyy_0_xxx_1,  \
                             g_0_xyyy_0_xxy_0,  \
                             g_0_xyyy_0_xxy_1,  \
                             g_0_xyyy_0_xyy_0,  \
                             g_0_xyyy_0_xyy_1,  \
                             g_0_xyyyz_0_xxx_0, \
                             g_0_xyyyz_0_xxy_0, \
                             g_0_xyyyz_0_xxz_0, \
                             g_0_xyyyz_0_xyy_0, \
                             g_0_xyyyz_0_xyz_0, \
                             g_0_xyyyz_0_xzz_0, \
                             g_0_xyyyz_0_yyy_0, \
                             g_0_xyyyz_0_yyz_0, \
                             g_0_xyyyz_0_yzz_0, \
                             g_0_xyyyz_0_zzz_0, \
                             g_0_yyyz_0_xxz_0,  \
                             g_0_yyyz_0_xxz_1,  \
                             g_0_yyyz_0_xyz_0,  \
                             g_0_yyyz_0_xyz_1,  \
                             g_0_yyyz_0_xz_1,   \
                             g_0_yyyz_0_xzz_0,  \
                             g_0_yyyz_0_xzz_1,  \
                             g_0_yyyz_0_yyy_0,  \
                             g_0_yyyz_0_yyy_1,  \
                             g_0_yyyz_0_yyz_0,  \
                             g_0_yyyz_0_yyz_1,  \
                             g_0_yyyz_0_yz_1,   \
                             g_0_yyyz_0_yzz_0,  \
                             g_0_yyyz_0_yzz_1,  \
                             g_0_yyyz_0_zz_1,   \
                             g_0_yyyz_0_zzz_0,  \
                             g_0_yyyz_0_zzz_1,  \
                             wp_x,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xxx_0[i] = g_0_xyyy_0_xxx_0[i] * pb_z + g_0_xyyy_0_xxx_1[i] * wp_z[i];

        g_0_xyyyz_0_xxy_0[i] = g_0_xyyy_0_xxy_0[i] * pb_z + g_0_xyyy_0_xxy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxz_0[i] = 2.0 * g_0_yyyz_0_xz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxz_0[i] * pb_x + g_0_yyyz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyy_0[i] = g_0_xyyy_0_xyy_0[i] * pb_z + g_0_xyyy_0_xyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xyz_0[i] = g_0_yyyz_0_yz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyz_0[i] * pb_x + g_0_yyyz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xzz_0[i] = g_0_yyyz_0_zz_1[i] * fi_abcd_0 + g_0_yyyz_0_xzz_0[i] * pb_x + g_0_yyyz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyy_0[i] = g_0_yyyz_0_yyy_0[i] * pb_x + g_0_yyyz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyz_0_yyz_0[i] = g_0_yyyz_0_yyz_0[i] * pb_x + g_0_yyyz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyz_0_yzz_0[i] = g_0_yyyz_0_yzz_0[i] * pb_x + g_0_yyyz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyz_0_zzz_0[i] = g_0_yyyz_0_zzz_0[i] * pb_x + g_0_yyyz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 120-130 components of targeted buffer : SHSF

    auto g_0_xyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 120);

    auto g_0_xyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 121);

    auto g_0_xyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 122);

    auto g_0_xyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 123);

    auto g_0_xyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 124);

    auto g_0_xyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 125);

    auto g_0_xyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 126);

    auto g_0_xyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 127);

    auto g_0_xyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 128);

    auto g_0_xyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 129);

#pragma omp simd aligned(g_0_xyyzz_0_xxx_0,     \
                             g_0_xyyzz_0_xxy_0, \
                             g_0_xyyzz_0_xxz_0, \
                             g_0_xyyzz_0_xyy_0, \
                             g_0_xyyzz_0_xyz_0, \
                             g_0_xyyzz_0_xzz_0, \
                             g_0_xyyzz_0_yyy_0, \
                             g_0_xyyzz_0_yyz_0, \
                             g_0_xyyzz_0_yzz_0, \
                             g_0_xyyzz_0_zzz_0, \
                             g_0_yyzz_0_xx_1,   \
                             g_0_yyzz_0_xxx_0,  \
                             g_0_yyzz_0_xxx_1,  \
                             g_0_yyzz_0_xxy_0,  \
                             g_0_yyzz_0_xxy_1,  \
                             g_0_yyzz_0_xxz_0,  \
                             g_0_yyzz_0_xxz_1,  \
                             g_0_yyzz_0_xy_1,   \
                             g_0_yyzz_0_xyy_0,  \
                             g_0_yyzz_0_xyy_1,  \
                             g_0_yyzz_0_xyz_0,  \
                             g_0_yyzz_0_xyz_1,  \
                             g_0_yyzz_0_xz_1,   \
                             g_0_yyzz_0_xzz_0,  \
                             g_0_yyzz_0_xzz_1,  \
                             g_0_yyzz_0_yy_1,   \
                             g_0_yyzz_0_yyy_0,  \
                             g_0_yyzz_0_yyy_1,  \
                             g_0_yyzz_0_yyz_0,  \
                             g_0_yyzz_0_yyz_1,  \
                             g_0_yyzz_0_yz_1,   \
                             g_0_yyzz_0_yzz_0,  \
                             g_0_yyzz_0_yzz_1,  \
                             g_0_yyzz_0_zz_1,   \
                             g_0_yyzz_0_zzz_0,  \
                             g_0_yyzz_0_zzz_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xxx_0[i] = 3.0 * g_0_yyzz_0_xx_1[i] * fi_abcd_0 + g_0_yyzz_0_xxx_0[i] * pb_x + g_0_yyzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyzz_0_xxy_0[i] = 2.0 * g_0_yyzz_0_xy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxy_0[i] * pb_x + g_0_yyzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxz_0[i] = 2.0 * g_0_yyzz_0_xz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxz_0[i] * pb_x + g_0_yyzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyy_0[i] = g_0_yyzz_0_yy_1[i] * fi_abcd_0 + g_0_yyzz_0_xyy_0[i] * pb_x + g_0_yyzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xyz_0[i] = g_0_yyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyz_0[i] * pb_x + g_0_yyzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xzz_0[i] = g_0_yyzz_0_zz_1[i] * fi_abcd_0 + g_0_yyzz_0_xzz_0[i] * pb_x + g_0_yyzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyy_0[i] = g_0_yyzz_0_yyy_0[i] * pb_x + g_0_yyzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyzz_0_yyz_0[i] = g_0_yyzz_0_yyz_0[i] * pb_x + g_0_yyzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyzz_0_yzz_0[i] = g_0_yyzz_0_yzz_0[i] * pb_x + g_0_yyzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyzz_0_zzz_0[i] = g_0_yyzz_0_zzz_0[i] * pb_x + g_0_yyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 130-140 components of targeted buffer : SHSF

    auto g_0_xyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 130);

    auto g_0_xyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 131);

    auto g_0_xyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 132);

    auto g_0_xyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 133);

    auto g_0_xyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 134);

    auto g_0_xyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 135);

    auto g_0_xyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 136);

    auto g_0_xyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 137);

    auto g_0_xyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 138);

    auto g_0_xyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 139);

#pragma omp simd aligned(g_0_xyzzz_0_xxx_0,     \
                             g_0_xyzzz_0_xxy_0, \
                             g_0_xyzzz_0_xxz_0, \
                             g_0_xyzzz_0_xyy_0, \
                             g_0_xyzzz_0_xyz_0, \
                             g_0_xyzzz_0_xzz_0, \
                             g_0_xyzzz_0_yyy_0, \
                             g_0_xyzzz_0_yyz_0, \
                             g_0_xyzzz_0_yzz_0, \
                             g_0_xyzzz_0_zzz_0, \
                             g_0_xzzz_0_xxx_0,  \
                             g_0_xzzz_0_xxx_1,  \
                             g_0_xzzz_0_xxz_0,  \
                             g_0_xzzz_0_xxz_1,  \
                             g_0_xzzz_0_xzz_0,  \
                             g_0_xzzz_0_xzz_1,  \
                             g_0_yzzz_0_xxy_0,  \
                             g_0_yzzz_0_xxy_1,  \
                             g_0_yzzz_0_xy_1,   \
                             g_0_yzzz_0_xyy_0,  \
                             g_0_yzzz_0_xyy_1,  \
                             g_0_yzzz_0_xyz_0,  \
                             g_0_yzzz_0_xyz_1,  \
                             g_0_yzzz_0_yy_1,   \
                             g_0_yzzz_0_yyy_0,  \
                             g_0_yzzz_0_yyy_1,  \
                             g_0_yzzz_0_yyz_0,  \
                             g_0_yzzz_0_yyz_1,  \
                             g_0_yzzz_0_yz_1,   \
                             g_0_yzzz_0_yzz_0,  \
                             g_0_yzzz_0_yzz_1,  \
                             g_0_yzzz_0_zzz_0,  \
                             g_0_yzzz_0_zzz_1,  \
                             wp_x,              \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xxx_0[i] = g_0_xzzz_0_xxx_0[i] * pb_y + g_0_xzzz_0_xxx_1[i] * wp_y[i];

        g_0_xyzzz_0_xxy_0[i] = 2.0 * g_0_yzzz_0_xy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxy_0[i] * pb_x + g_0_yzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxz_0[i] = g_0_xzzz_0_xxz_0[i] * pb_y + g_0_xzzz_0_xxz_1[i] * wp_y[i];

        g_0_xyzzz_0_xyy_0[i] = g_0_yzzz_0_yy_1[i] * fi_abcd_0 + g_0_yzzz_0_xyy_0[i] * pb_x + g_0_yzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xyz_0[i] = g_0_yzzz_0_yz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyz_0[i] * pb_x + g_0_yzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xzz_0[i] = g_0_xzzz_0_xzz_0[i] * pb_y + g_0_xzzz_0_xzz_1[i] * wp_y[i];

        g_0_xyzzz_0_yyy_0[i] = g_0_yzzz_0_yyy_0[i] * pb_x + g_0_yzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyzzz_0_yyz_0[i] = g_0_yzzz_0_yyz_0[i] * pb_x + g_0_yzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyzzz_0_yzz_0[i] = g_0_yzzz_0_yzz_0[i] * pb_x + g_0_yzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyzzz_0_zzz_0[i] = g_0_yzzz_0_zzz_0[i] * pb_x + g_0_yzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 140-150 components of targeted buffer : SHSF

    auto g_0_xzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 140);

    auto g_0_xzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 141);

    auto g_0_xzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 142);

    auto g_0_xzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 143);

    auto g_0_xzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 144);

    auto g_0_xzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 145);

    auto g_0_xzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 146);

    auto g_0_xzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 147);

    auto g_0_xzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 148);

    auto g_0_xzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 149);

#pragma omp simd aligned(g_0_xzzzz_0_xxx_0,     \
                             g_0_xzzzz_0_xxy_0, \
                             g_0_xzzzz_0_xxz_0, \
                             g_0_xzzzz_0_xyy_0, \
                             g_0_xzzzz_0_xyz_0, \
                             g_0_xzzzz_0_xzz_0, \
                             g_0_xzzzz_0_yyy_0, \
                             g_0_xzzzz_0_yyz_0, \
                             g_0_xzzzz_0_yzz_0, \
                             g_0_xzzzz_0_zzz_0, \
                             g_0_zzzz_0_xx_1,   \
                             g_0_zzzz_0_xxx_0,  \
                             g_0_zzzz_0_xxx_1,  \
                             g_0_zzzz_0_xxy_0,  \
                             g_0_zzzz_0_xxy_1,  \
                             g_0_zzzz_0_xxz_0,  \
                             g_0_zzzz_0_xxz_1,  \
                             g_0_zzzz_0_xy_1,   \
                             g_0_zzzz_0_xyy_0,  \
                             g_0_zzzz_0_xyy_1,  \
                             g_0_zzzz_0_xyz_0,  \
                             g_0_zzzz_0_xyz_1,  \
                             g_0_zzzz_0_xz_1,   \
                             g_0_zzzz_0_xzz_0,  \
                             g_0_zzzz_0_xzz_1,  \
                             g_0_zzzz_0_yy_1,   \
                             g_0_zzzz_0_yyy_0,  \
                             g_0_zzzz_0_yyy_1,  \
                             g_0_zzzz_0_yyz_0,  \
                             g_0_zzzz_0_yyz_1,  \
                             g_0_zzzz_0_yz_1,   \
                             g_0_zzzz_0_yzz_0,  \
                             g_0_zzzz_0_yzz_1,  \
                             g_0_zzzz_0_zz_1,   \
                             g_0_zzzz_0_zzz_0,  \
                             g_0_zzzz_0_zzz_1,  \
                             wp_x,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xxx_0[i] = 3.0 * g_0_zzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxx_0[i] * pb_x + g_0_zzzz_0_xxx_1[i] * wp_x[i];

        g_0_xzzzz_0_xxy_0[i] = 2.0 * g_0_zzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxy_0[i] * pb_x + g_0_zzzz_0_xxy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxz_0[i] = 2.0 * g_0_zzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxz_0[i] * pb_x + g_0_zzzz_0_xxz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyy_0[i] = g_0_zzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyy_0[i] * pb_x + g_0_zzzz_0_xyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xyz_0[i] = g_0_zzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyz_0[i] * pb_x + g_0_zzzz_0_xyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xzz_0[i] = g_0_zzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzz_0[i] * pb_x + g_0_zzzz_0_xzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyy_0[i] = g_0_zzzz_0_yyy_0[i] * pb_x + g_0_zzzz_0_yyy_1[i] * wp_x[i];

        g_0_xzzzz_0_yyz_0[i] = g_0_zzzz_0_yyz_0[i] * pb_x + g_0_zzzz_0_yyz_1[i] * wp_x[i];

        g_0_xzzzz_0_yzz_0[i] = g_0_zzzz_0_yzz_0[i] * pb_x + g_0_zzzz_0_yzz_1[i] * wp_x[i];

        g_0_xzzzz_0_zzz_0[i] = g_0_zzzz_0_zzz_0[i] * pb_x + g_0_zzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 150-160 components of targeted buffer : SHSF

    auto g_0_yyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 150);

    auto g_0_yyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 151);

    auto g_0_yyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 152);

    auto g_0_yyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 153);

    auto g_0_yyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 154);

    auto g_0_yyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 155);

    auto g_0_yyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 156);

    auto g_0_yyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 157);

    auto g_0_yyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 158);

    auto g_0_yyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 159);

#pragma omp simd aligned(g_0_yyy_0_xxx_0,       \
                             g_0_yyy_0_xxx_1,   \
                             g_0_yyy_0_xxy_0,   \
                             g_0_yyy_0_xxy_1,   \
                             g_0_yyy_0_xxz_0,   \
                             g_0_yyy_0_xxz_1,   \
                             g_0_yyy_0_xyy_0,   \
                             g_0_yyy_0_xyy_1,   \
                             g_0_yyy_0_xyz_0,   \
                             g_0_yyy_0_xyz_1,   \
                             g_0_yyy_0_xzz_0,   \
                             g_0_yyy_0_xzz_1,   \
                             g_0_yyy_0_yyy_0,   \
                             g_0_yyy_0_yyy_1,   \
                             g_0_yyy_0_yyz_0,   \
                             g_0_yyy_0_yyz_1,   \
                             g_0_yyy_0_yzz_0,   \
                             g_0_yyy_0_yzz_1,   \
                             g_0_yyy_0_zzz_0,   \
                             g_0_yyy_0_zzz_1,   \
                             g_0_yyyy_0_xx_1,   \
                             g_0_yyyy_0_xxx_0,  \
                             g_0_yyyy_0_xxx_1,  \
                             g_0_yyyy_0_xxy_0,  \
                             g_0_yyyy_0_xxy_1,  \
                             g_0_yyyy_0_xxz_0,  \
                             g_0_yyyy_0_xxz_1,  \
                             g_0_yyyy_0_xy_1,   \
                             g_0_yyyy_0_xyy_0,  \
                             g_0_yyyy_0_xyy_1,  \
                             g_0_yyyy_0_xyz_0,  \
                             g_0_yyyy_0_xyz_1,  \
                             g_0_yyyy_0_xz_1,   \
                             g_0_yyyy_0_xzz_0,  \
                             g_0_yyyy_0_xzz_1,  \
                             g_0_yyyy_0_yy_1,   \
                             g_0_yyyy_0_yyy_0,  \
                             g_0_yyyy_0_yyy_1,  \
                             g_0_yyyy_0_yyz_0,  \
                             g_0_yyyy_0_yyz_1,  \
                             g_0_yyyy_0_yz_1,   \
                             g_0_yyyy_0_yzz_0,  \
                             g_0_yyyy_0_yzz_1,  \
                             g_0_yyyy_0_zz_1,   \
                             g_0_yyyy_0_zzz_0,  \
                             g_0_yyyy_0_zzz_1,  \
                             g_0_yyyyy_0_xxx_0, \
                             g_0_yyyyy_0_xxy_0, \
                             g_0_yyyyy_0_xxz_0, \
                             g_0_yyyyy_0_xyy_0, \
                             g_0_yyyyy_0_xyz_0, \
                             g_0_yyyyy_0_xzz_0, \
                             g_0_yyyyy_0_yyy_0, \
                             g_0_yyyyy_0_yyz_0, \
                             g_0_yyyyy_0_yzz_0, \
                             g_0_yyyyy_0_zzz_0, \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xxx_0[i] =
            4.0 * g_0_yyy_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxx_1[i] * fti_ab_0 + g_0_yyyy_0_xxx_0[i] * pb_y + g_0_yyyy_0_xxx_1[i] * wp_y[i];

        g_0_yyyyy_0_xxy_0[i] = 4.0 * g_0_yyy_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyy_0_xx_1[i] * fi_abcd_0 +
                               g_0_yyyy_0_xxy_0[i] * pb_y + g_0_yyyy_0_xxy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxz_0[i] =
            4.0 * g_0_yyy_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxz_1[i] * fti_ab_0 + g_0_yyyy_0_xxz_0[i] * pb_y + g_0_yyyy_0_xxz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyy_0[i] = 4.0 * g_0_yyy_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xy_1[i] * fi_abcd_0 +
                               g_0_yyyy_0_xyy_0[i] * pb_y + g_0_yyyy_0_xyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xyz_0[i] = 4.0 * g_0_yyy_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyz_1[i] * fti_ab_0 + g_0_yyyy_0_xz_1[i] * fi_abcd_0 +
                               g_0_yyyy_0_xyz_0[i] * pb_y + g_0_yyyy_0_xyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xzz_0[i] =
            4.0 * g_0_yyy_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzz_0[i] * pb_y + g_0_yyyy_0_xzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyy_0[i] = 4.0 * g_0_yyy_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_yy_1[i] * fi_abcd_0 +
                               g_0_yyyy_0_yyy_0[i] * pb_y + g_0_yyyy_0_yyy_1[i] * wp_y[i];

        g_0_yyyyy_0_yyz_0[i] = 4.0 * g_0_yyy_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_yz_1[i] * fi_abcd_0 +
                               g_0_yyyy_0_yyz_0[i] * pb_y + g_0_yyyy_0_yyz_1[i] * wp_y[i];

        g_0_yyyyy_0_yzz_0[i] = 4.0 * g_0_yyy_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yzz_1[i] * fti_ab_0 + g_0_yyyy_0_zz_1[i] * fi_abcd_0 +
                               g_0_yyyy_0_yzz_0[i] * pb_y + g_0_yyyy_0_yzz_1[i] * wp_y[i];

        g_0_yyyyy_0_zzz_0[i] =
            4.0 * g_0_yyy_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzz_0[i] * pb_y + g_0_yyyy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 160-170 components of targeted buffer : SHSF

    auto g_0_yyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 160);

    auto g_0_yyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 161);

    auto g_0_yyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 162);

    auto g_0_yyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 163);

    auto g_0_yyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 164);

    auto g_0_yyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 165);

    auto g_0_yyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 166);

    auto g_0_yyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 167);

    auto g_0_yyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 168);

    auto g_0_yyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 169);

#pragma omp simd aligned(g_0_yyyy_0_xx_1,       \
                             g_0_yyyy_0_xxx_0,  \
                             g_0_yyyy_0_xxx_1,  \
                             g_0_yyyy_0_xxy_0,  \
                             g_0_yyyy_0_xxy_1,  \
                             g_0_yyyy_0_xxz_0,  \
                             g_0_yyyy_0_xxz_1,  \
                             g_0_yyyy_0_xy_1,   \
                             g_0_yyyy_0_xyy_0,  \
                             g_0_yyyy_0_xyy_1,  \
                             g_0_yyyy_0_xyz_0,  \
                             g_0_yyyy_0_xyz_1,  \
                             g_0_yyyy_0_xz_1,   \
                             g_0_yyyy_0_xzz_0,  \
                             g_0_yyyy_0_xzz_1,  \
                             g_0_yyyy_0_yy_1,   \
                             g_0_yyyy_0_yyy_0,  \
                             g_0_yyyy_0_yyy_1,  \
                             g_0_yyyy_0_yyz_0,  \
                             g_0_yyyy_0_yyz_1,  \
                             g_0_yyyy_0_yz_1,   \
                             g_0_yyyy_0_yzz_0,  \
                             g_0_yyyy_0_yzz_1,  \
                             g_0_yyyy_0_zz_1,   \
                             g_0_yyyy_0_zzz_0,  \
                             g_0_yyyy_0_zzz_1,  \
                             g_0_yyyyz_0_xxx_0, \
                             g_0_yyyyz_0_xxy_0, \
                             g_0_yyyyz_0_xxz_0, \
                             g_0_yyyyz_0_xyy_0, \
                             g_0_yyyyz_0_xyz_0, \
                             g_0_yyyyz_0_xzz_0, \
                             g_0_yyyyz_0_yyy_0, \
                             g_0_yyyyz_0_yyz_0, \
                             g_0_yyyyz_0_yzz_0, \
                             g_0_yyyyz_0_zzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xxx_0[i] = g_0_yyyy_0_xxx_0[i] * pb_z + g_0_yyyy_0_xxx_1[i] * wp_z[i];

        g_0_yyyyz_0_xxy_0[i] = g_0_yyyy_0_xxy_0[i] * pb_z + g_0_yyyy_0_xxy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxz_0[i] = g_0_yyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxz_0[i] * pb_z + g_0_yyyy_0_xxz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyy_0[i] = g_0_yyyy_0_xyy_0[i] * pb_z + g_0_yyyy_0_xyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xyz_0[i] = g_0_yyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyz_0[i] * pb_z + g_0_yyyy_0_xyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xzz_0[i] = 2.0 * g_0_yyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzz_0[i] * pb_z + g_0_yyyy_0_xzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyy_0[i] = g_0_yyyy_0_yyy_0[i] * pb_z + g_0_yyyy_0_yyy_1[i] * wp_z[i];

        g_0_yyyyz_0_yyz_0[i] = g_0_yyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyz_0[i] * pb_z + g_0_yyyy_0_yyz_1[i] * wp_z[i];

        g_0_yyyyz_0_yzz_0[i] = 2.0 * g_0_yyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzz_0[i] * pb_z + g_0_yyyy_0_yzz_1[i] * wp_z[i];

        g_0_yyyyz_0_zzz_0[i] = 3.0 * g_0_yyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyy_0_zzz_0[i] * pb_z + g_0_yyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 170-180 components of targeted buffer : SHSF

    auto g_0_yyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 170);

    auto g_0_yyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 171);

    auto g_0_yyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 172);

    auto g_0_yyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 173);

    auto g_0_yyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 174);

    auto g_0_yyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 175);

    auto g_0_yyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 176);

    auto g_0_yyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 177);

    auto g_0_yyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 178);

    auto g_0_yyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 179);

#pragma omp simd aligned(g_0_yyy_0_xxy_0,       \
                             g_0_yyy_0_xxy_1,   \
                             g_0_yyy_0_xyy_0,   \
                             g_0_yyy_0_xyy_1,   \
                             g_0_yyy_0_yyy_0,   \
                             g_0_yyy_0_yyy_1,   \
                             g_0_yyyz_0_xxy_0,  \
                             g_0_yyyz_0_xxy_1,  \
                             g_0_yyyz_0_xyy_0,  \
                             g_0_yyyz_0_xyy_1,  \
                             g_0_yyyz_0_yyy_0,  \
                             g_0_yyyz_0_yyy_1,  \
                             g_0_yyyzz_0_xxx_0, \
                             g_0_yyyzz_0_xxy_0, \
                             g_0_yyyzz_0_xxz_0, \
                             g_0_yyyzz_0_xyy_0, \
                             g_0_yyyzz_0_xyz_0, \
                             g_0_yyyzz_0_xzz_0, \
                             g_0_yyyzz_0_yyy_0, \
                             g_0_yyyzz_0_yyz_0, \
                             g_0_yyyzz_0_yzz_0, \
                             g_0_yyyzz_0_zzz_0, \
                             g_0_yyzz_0_xxx_0,  \
                             g_0_yyzz_0_xxx_1,  \
                             g_0_yyzz_0_xxz_0,  \
                             g_0_yyzz_0_xxz_1,  \
                             g_0_yyzz_0_xyz_0,  \
                             g_0_yyzz_0_xyz_1,  \
                             g_0_yyzz_0_xz_1,   \
                             g_0_yyzz_0_xzz_0,  \
                             g_0_yyzz_0_xzz_1,  \
                             g_0_yyzz_0_yyz_0,  \
                             g_0_yyzz_0_yyz_1,  \
                             g_0_yyzz_0_yz_1,   \
                             g_0_yyzz_0_yzz_0,  \
                             g_0_yyzz_0_yzz_1,  \
                             g_0_yyzz_0_zz_1,   \
                             g_0_yyzz_0_zzz_0,  \
                             g_0_yyzz_0_zzz_1,  \
                             g_0_yzz_0_xxx_0,   \
                             g_0_yzz_0_xxx_1,   \
                             g_0_yzz_0_xxz_0,   \
                             g_0_yzz_0_xxz_1,   \
                             g_0_yzz_0_xyz_0,   \
                             g_0_yzz_0_xyz_1,   \
                             g_0_yzz_0_xzz_0,   \
                             g_0_yzz_0_xzz_1,   \
                             g_0_yzz_0_yyz_0,   \
                             g_0_yzz_0_yyz_1,   \
                             g_0_yzz_0_yzz_0,   \
                             g_0_yzz_0_yzz_1,   \
                             g_0_yzz_0_zzz_0,   \
                             g_0_yzz_0_zzz_1,   \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xxx_0[i] =
            2.0 * g_0_yzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxx_1[i] * fti_ab_0 + g_0_yyzz_0_xxx_0[i] * pb_y + g_0_yyzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyzz_0_xxy_0[i] =
            g_0_yyy_0_xxy_0[i] * fi_ab_0 - g_0_yyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyz_0_xxy_0[i] * pb_z + g_0_yyyz_0_xxy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxz_0[i] =
            2.0 * g_0_yzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxz_1[i] * fti_ab_0 + g_0_yyzz_0_xxz_0[i] * pb_y + g_0_yyzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyy_0[i] =
            g_0_yyy_0_xyy_0[i] * fi_ab_0 - g_0_yyy_0_xyy_1[i] * fti_ab_0 + g_0_yyyz_0_xyy_0[i] * pb_z + g_0_yyyz_0_xyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xyz_0[i] = 2.0 * g_0_yzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyz_1[i] * fti_ab_0 + g_0_yyzz_0_xz_1[i] * fi_abcd_0 +
                               g_0_yyzz_0_xyz_0[i] * pb_y + g_0_yyzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xzz_0[i] =
            2.0 * g_0_yzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzz_0[i] * pb_y + g_0_yyzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyy_0[i] =
            g_0_yyy_0_yyy_0[i] * fi_ab_0 - g_0_yyy_0_yyy_1[i] * fti_ab_0 + g_0_yyyz_0_yyy_0[i] * pb_z + g_0_yyyz_0_yyy_1[i] * wp_z[i];

        g_0_yyyzz_0_yyz_0[i] = 2.0 * g_0_yzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_yz_1[i] * fi_abcd_0 +
                               g_0_yyzz_0_yyz_0[i] * pb_y + g_0_yyzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyzz_0_yzz_0[i] = 2.0 * g_0_yzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yzz_1[i] * fti_ab_0 + g_0_yyzz_0_zz_1[i] * fi_abcd_0 +
                               g_0_yyzz_0_yzz_0[i] * pb_y + g_0_yyzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyzz_0_zzz_0[i] =
            2.0 * g_0_yzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzz_0[i] * pb_y + g_0_yyzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 180-190 components of targeted buffer : SHSF

    auto g_0_yyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 180);

    auto g_0_yyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 181);

    auto g_0_yyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 182);

    auto g_0_yyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 183);

    auto g_0_yyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 184);

    auto g_0_yyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 185);

    auto g_0_yyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 186);

    auto g_0_yyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 187);

    auto g_0_yyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 188);

    auto g_0_yyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 189);

#pragma omp simd aligned(g_0_yyz_0_xxy_0,       \
                             g_0_yyz_0_xxy_1,   \
                             g_0_yyz_0_xyy_0,   \
                             g_0_yyz_0_xyy_1,   \
                             g_0_yyz_0_yyy_0,   \
                             g_0_yyz_0_yyy_1,   \
                             g_0_yyzz_0_xxy_0,  \
                             g_0_yyzz_0_xxy_1,  \
                             g_0_yyzz_0_xyy_0,  \
                             g_0_yyzz_0_xyy_1,  \
                             g_0_yyzz_0_yyy_0,  \
                             g_0_yyzz_0_yyy_1,  \
                             g_0_yyzzz_0_xxx_0, \
                             g_0_yyzzz_0_xxy_0, \
                             g_0_yyzzz_0_xxz_0, \
                             g_0_yyzzz_0_xyy_0, \
                             g_0_yyzzz_0_xyz_0, \
                             g_0_yyzzz_0_xzz_0, \
                             g_0_yyzzz_0_yyy_0, \
                             g_0_yyzzz_0_yyz_0, \
                             g_0_yyzzz_0_yzz_0, \
                             g_0_yyzzz_0_zzz_0, \
                             g_0_yzzz_0_xxx_0,  \
                             g_0_yzzz_0_xxx_1,  \
                             g_0_yzzz_0_xxz_0,  \
                             g_0_yzzz_0_xxz_1,  \
                             g_0_yzzz_0_xyz_0,  \
                             g_0_yzzz_0_xyz_1,  \
                             g_0_yzzz_0_xz_1,   \
                             g_0_yzzz_0_xzz_0,  \
                             g_0_yzzz_0_xzz_1,  \
                             g_0_yzzz_0_yyz_0,  \
                             g_0_yzzz_0_yyz_1,  \
                             g_0_yzzz_0_yz_1,   \
                             g_0_yzzz_0_yzz_0,  \
                             g_0_yzzz_0_yzz_1,  \
                             g_0_yzzz_0_zz_1,   \
                             g_0_yzzz_0_zzz_0,  \
                             g_0_yzzz_0_zzz_1,  \
                             g_0_zzz_0_xxx_0,   \
                             g_0_zzz_0_xxx_1,   \
                             g_0_zzz_0_xxz_0,   \
                             g_0_zzz_0_xxz_1,   \
                             g_0_zzz_0_xyz_0,   \
                             g_0_zzz_0_xyz_1,   \
                             g_0_zzz_0_xzz_0,   \
                             g_0_zzz_0_xzz_1,   \
                             g_0_zzz_0_yyz_0,   \
                             g_0_zzz_0_yyz_1,   \
                             g_0_zzz_0_yzz_0,   \
                             g_0_zzz_0_yzz_1,   \
                             g_0_zzz_0_zzz_0,   \
                             g_0_zzz_0_zzz_1,   \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xxx_0[i] =
            g_0_zzz_0_xxx_0[i] * fi_ab_0 - g_0_zzz_0_xxx_1[i] * fti_ab_0 + g_0_yzzz_0_xxx_0[i] * pb_y + g_0_yzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyzzz_0_xxy_0[i] =
            2.0 * g_0_yyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxy_1[i] * fti_ab_0 + g_0_yyzz_0_xxy_0[i] * pb_z + g_0_yyzz_0_xxy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxz_0[i] =
            g_0_zzz_0_xxz_0[i] * fi_ab_0 - g_0_zzz_0_xxz_1[i] * fti_ab_0 + g_0_yzzz_0_xxz_0[i] * pb_y + g_0_yzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyy_0[i] =
            2.0 * g_0_yyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xyy_1[i] * fti_ab_0 + g_0_yyzz_0_xyy_0[i] * pb_z + g_0_yyzz_0_xyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xyz_0[i] = g_0_zzz_0_xyz_0[i] * fi_ab_0 - g_0_zzz_0_xyz_1[i] * fti_ab_0 + g_0_yzzz_0_xz_1[i] * fi_abcd_0 +
                               g_0_yzzz_0_xyz_0[i] * pb_y + g_0_yzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xzz_0[i] =
            g_0_zzz_0_xzz_0[i] * fi_ab_0 - g_0_zzz_0_xzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzz_0[i] * pb_y + g_0_yzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyy_0[i] =
            2.0 * g_0_yyz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yyy_1[i] * fti_ab_0 + g_0_yyzz_0_yyy_0[i] * pb_z + g_0_yyzz_0_yyy_1[i] * wp_z[i];

        g_0_yyzzz_0_yyz_0[i] = g_0_zzz_0_yyz_0[i] * fi_ab_0 - g_0_zzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_yz_1[i] * fi_abcd_0 +
                               g_0_yzzz_0_yyz_0[i] * pb_y + g_0_yzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyzzz_0_yzz_0[i] = g_0_zzz_0_yzz_0[i] * fi_ab_0 - g_0_zzz_0_yzz_1[i] * fti_ab_0 + g_0_yzzz_0_zz_1[i] * fi_abcd_0 +
                               g_0_yzzz_0_yzz_0[i] * pb_y + g_0_yzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyzzz_0_zzz_0[i] =
            g_0_zzz_0_zzz_0[i] * fi_ab_0 - g_0_zzz_0_zzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzz_0[i] * pb_y + g_0_yzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 190-200 components of targeted buffer : SHSF

    auto g_0_yzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 190);

    auto g_0_yzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 191);

    auto g_0_yzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 192);

    auto g_0_yzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 193);

    auto g_0_yzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 194);

    auto g_0_yzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 195);

    auto g_0_yzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 196);

    auto g_0_yzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 197);

    auto g_0_yzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 198);

    auto g_0_yzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 199);

#pragma omp simd aligned(g_0_yzzzz_0_xxx_0,     \
                             g_0_yzzzz_0_xxy_0, \
                             g_0_yzzzz_0_xxz_0, \
                             g_0_yzzzz_0_xyy_0, \
                             g_0_yzzzz_0_xyz_0, \
                             g_0_yzzzz_0_xzz_0, \
                             g_0_yzzzz_0_yyy_0, \
                             g_0_yzzzz_0_yyz_0, \
                             g_0_yzzzz_0_yzz_0, \
                             g_0_yzzzz_0_zzz_0, \
                             g_0_zzzz_0_xx_1,   \
                             g_0_zzzz_0_xxx_0,  \
                             g_0_zzzz_0_xxx_1,  \
                             g_0_zzzz_0_xxy_0,  \
                             g_0_zzzz_0_xxy_1,  \
                             g_0_zzzz_0_xxz_0,  \
                             g_0_zzzz_0_xxz_1,  \
                             g_0_zzzz_0_xy_1,   \
                             g_0_zzzz_0_xyy_0,  \
                             g_0_zzzz_0_xyy_1,  \
                             g_0_zzzz_0_xyz_0,  \
                             g_0_zzzz_0_xyz_1,  \
                             g_0_zzzz_0_xz_1,   \
                             g_0_zzzz_0_xzz_0,  \
                             g_0_zzzz_0_xzz_1,  \
                             g_0_zzzz_0_yy_1,   \
                             g_0_zzzz_0_yyy_0,  \
                             g_0_zzzz_0_yyy_1,  \
                             g_0_zzzz_0_yyz_0,  \
                             g_0_zzzz_0_yyz_1,  \
                             g_0_zzzz_0_yz_1,   \
                             g_0_zzzz_0_yzz_0,  \
                             g_0_zzzz_0_yzz_1,  \
                             g_0_zzzz_0_zz_1,   \
                             g_0_zzzz_0_zzz_0,  \
                             g_0_zzzz_0_zzz_1,  \
                             wp_y,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xxx_0[i] = g_0_zzzz_0_xxx_0[i] * pb_y + g_0_zzzz_0_xxx_1[i] * wp_y[i];

        g_0_yzzzz_0_xxy_0[i] = g_0_zzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxy_0[i] * pb_y + g_0_zzzz_0_xxy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxz_0[i] = g_0_zzzz_0_xxz_0[i] * pb_y + g_0_zzzz_0_xxz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyy_0[i] = 2.0 * g_0_zzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyy_0[i] * pb_y + g_0_zzzz_0_xyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xyz_0[i] = g_0_zzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyz_0[i] * pb_y + g_0_zzzz_0_xyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xzz_0[i] = g_0_zzzz_0_xzz_0[i] * pb_y + g_0_zzzz_0_xzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyy_0[i] = 3.0 * g_0_zzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyy_0[i] * pb_y + g_0_zzzz_0_yyy_1[i] * wp_y[i];

        g_0_yzzzz_0_yyz_0[i] = 2.0 * g_0_zzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyz_0[i] * pb_y + g_0_zzzz_0_yyz_1[i] * wp_y[i];

        g_0_yzzzz_0_yzz_0[i] = g_0_zzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzz_0[i] * pb_y + g_0_zzzz_0_yzz_1[i] * wp_y[i];

        g_0_yzzzz_0_zzz_0[i] = g_0_zzzz_0_zzz_0[i] * pb_y + g_0_zzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 200-210 components of targeted buffer : SHSF

    auto g_0_zzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 200);

    auto g_0_zzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 201);

    auto g_0_zzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 202);

    auto g_0_zzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 203);

    auto g_0_zzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 204);

    auto g_0_zzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 205);

    auto g_0_zzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 206);

    auto g_0_zzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 207);

    auto g_0_zzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 208);

    auto g_0_zzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 209);

#pragma omp simd aligned(g_0_zzz_0_xxx_0,       \
                             g_0_zzz_0_xxx_1,   \
                             g_0_zzz_0_xxy_0,   \
                             g_0_zzz_0_xxy_1,   \
                             g_0_zzz_0_xxz_0,   \
                             g_0_zzz_0_xxz_1,   \
                             g_0_zzz_0_xyy_0,   \
                             g_0_zzz_0_xyy_1,   \
                             g_0_zzz_0_xyz_0,   \
                             g_0_zzz_0_xyz_1,   \
                             g_0_zzz_0_xzz_0,   \
                             g_0_zzz_0_xzz_1,   \
                             g_0_zzz_0_yyy_0,   \
                             g_0_zzz_0_yyy_1,   \
                             g_0_zzz_0_yyz_0,   \
                             g_0_zzz_0_yyz_1,   \
                             g_0_zzz_0_yzz_0,   \
                             g_0_zzz_0_yzz_1,   \
                             g_0_zzz_0_zzz_0,   \
                             g_0_zzz_0_zzz_1,   \
                             g_0_zzzz_0_xx_1,   \
                             g_0_zzzz_0_xxx_0,  \
                             g_0_zzzz_0_xxx_1,  \
                             g_0_zzzz_0_xxy_0,  \
                             g_0_zzzz_0_xxy_1,  \
                             g_0_zzzz_0_xxz_0,  \
                             g_0_zzzz_0_xxz_1,  \
                             g_0_zzzz_0_xy_1,   \
                             g_0_zzzz_0_xyy_0,  \
                             g_0_zzzz_0_xyy_1,  \
                             g_0_zzzz_0_xyz_0,  \
                             g_0_zzzz_0_xyz_1,  \
                             g_0_zzzz_0_xz_1,   \
                             g_0_zzzz_0_xzz_0,  \
                             g_0_zzzz_0_xzz_1,  \
                             g_0_zzzz_0_yy_1,   \
                             g_0_zzzz_0_yyy_0,  \
                             g_0_zzzz_0_yyy_1,  \
                             g_0_zzzz_0_yyz_0,  \
                             g_0_zzzz_0_yyz_1,  \
                             g_0_zzzz_0_yz_1,   \
                             g_0_zzzz_0_yzz_0,  \
                             g_0_zzzz_0_yzz_1,  \
                             g_0_zzzz_0_zz_1,   \
                             g_0_zzzz_0_zzz_0,  \
                             g_0_zzzz_0_zzz_1,  \
                             g_0_zzzzz_0_xxx_0, \
                             g_0_zzzzz_0_xxy_0, \
                             g_0_zzzzz_0_xxz_0, \
                             g_0_zzzzz_0_xyy_0, \
                             g_0_zzzzz_0_xyz_0, \
                             g_0_zzzzz_0_xzz_0, \
                             g_0_zzzzz_0_yyy_0, \
                             g_0_zzzzz_0_yyz_0, \
                             g_0_zzzzz_0_yzz_0, \
                             g_0_zzzzz_0_zzz_0, \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xxx_0[i] =
            4.0 * g_0_zzz_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxx_1[i] * fti_ab_0 + g_0_zzzz_0_xxx_0[i] * pb_z + g_0_zzzz_0_xxx_1[i] * wp_z[i];

        g_0_zzzzz_0_xxy_0[i] =
            4.0 * g_0_zzz_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxy_1[i] * fti_ab_0 + g_0_zzzz_0_xxy_0[i] * pb_z + g_0_zzzz_0_xxy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxz_0[i] = 4.0 * g_0_zzz_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxz_1[i] * fti_ab_0 + g_0_zzzz_0_xx_1[i] * fi_abcd_0 +
                               g_0_zzzz_0_xxz_0[i] * pb_z + g_0_zzzz_0_xxz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyy_0[i] =
            4.0 * g_0_zzz_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyy_1[i] * fti_ab_0 + g_0_zzzz_0_xyy_0[i] * pb_z + g_0_zzzz_0_xyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xyz_0[i] = 4.0 * g_0_zzz_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyz_1[i] * fti_ab_0 + g_0_zzzz_0_xy_1[i] * fi_abcd_0 +
                               g_0_zzzz_0_xyz_0[i] * pb_z + g_0_zzzz_0_xyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xzz_0[i] = 4.0 * g_0_zzz_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xz_1[i] * fi_abcd_0 +
                               g_0_zzzz_0_xzz_0[i] * pb_z + g_0_zzzz_0_xzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyy_0[i] =
            4.0 * g_0_zzz_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyy_1[i] * fti_ab_0 + g_0_zzzz_0_yyy_0[i] * pb_z + g_0_zzzz_0_yyy_1[i] * wp_z[i];

        g_0_zzzzz_0_yyz_0[i] = 4.0 * g_0_zzz_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyz_1[i] * fti_ab_0 + g_0_zzzz_0_yy_1[i] * fi_abcd_0 +
                               g_0_zzzz_0_yyz_0[i] * pb_z + g_0_zzzz_0_yyz_1[i] * wp_z[i];

        g_0_zzzzz_0_yzz_0[i] = 4.0 * g_0_zzz_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_yz_1[i] * fi_abcd_0 +
                               g_0_zzzz_0_yzz_0[i] * pb_z + g_0_zzzz_0_yzz_1[i] * wp_z[i];

        g_0_zzzzz_0_zzz_0[i] = 4.0 * g_0_zzz_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_zz_1[i] * fi_abcd_0 +
                               g_0_zzzz_0_zzz_0[i] * pb_z + g_0_zzzz_0_zzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
