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

#include "KineticEnergyPrimRecGF.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gf(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gf,
                            const size_t              idx_ovl_df,
                            const size_t              idx_kin_df,
                            const size_t              idx_kin_fd,
                            const size_t              idx_kin_ff,
                            const size_t              idx_ovl_gf,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_ovl_df);

    auto ts_xx_xxy = pbuffer.data(idx_ovl_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_ovl_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_ovl_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_ovl_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_ovl_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_ovl_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_ovl_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_ovl_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_ovl_df + 9);

    auto ts_yy_xxx = pbuffer.data(idx_ovl_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_ovl_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_ovl_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_ovl_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_ovl_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_ovl_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_ovl_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_ovl_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_ovl_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_ovl_df + 39);

    auto ts_zz_xxx = pbuffer.data(idx_ovl_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_ovl_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_ovl_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_ovl_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_ovl_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_ovl_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_ovl_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_ovl_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_ovl_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_ovl_df + 59);

    // Set up components of auxiliary buffer : DF

    auto tk_xx_xxx = pbuffer.data(idx_kin_df);

    auto tk_xx_xxy = pbuffer.data(idx_kin_df + 1);

    auto tk_xx_xxz = pbuffer.data(idx_kin_df + 2);

    auto tk_xx_xyy = pbuffer.data(idx_kin_df + 3);

    auto tk_xx_xyz = pbuffer.data(idx_kin_df + 4);

    auto tk_xx_xzz = pbuffer.data(idx_kin_df + 5);

    auto tk_xx_yyy = pbuffer.data(idx_kin_df + 6);

    auto tk_xx_yyz = pbuffer.data(idx_kin_df + 7);

    auto tk_xx_yzz = pbuffer.data(idx_kin_df + 8);

    auto tk_xx_zzz = pbuffer.data(idx_kin_df + 9);

    auto tk_yy_xxx = pbuffer.data(idx_kin_df + 30);

    auto tk_yy_xxy = pbuffer.data(idx_kin_df + 31);

    auto tk_yy_xxz = pbuffer.data(idx_kin_df + 32);

    auto tk_yy_xyy = pbuffer.data(idx_kin_df + 33);

    auto tk_yy_xyz = pbuffer.data(idx_kin_df + 34);

    auto tk_yy_xzz = pbuffer.data(idx_kin_df + 35);

    auto tk_yy_yyy = pbuffer.data(idx_kin_df + 36);

    auto tk_yy_yyz = pbuffer.data(idx_kin_df + 37);

    auto tk_yy_yzz = pbuffer.data(idx_kin_df + 38);

    auto tk_yy_zzz = pbuffer.data(idx_kin_df + 39);

    auto tk_zz_xxx = pbuffer.data(idx_kin_df + 50);

    auto tk_zz_xxy = pbuffer.data(idx_kin_df + 51);

    auto tk_zz_xxz = pbuffer.data(idx_kin_df + 52);

    auto tk_zz_xyy = pbuffer.data(idx_kin_df + 53);

    auto tk_zz_xyz = pbuffer.data(idx_kin_df + 54);

    auto tk_zz_xzz = pbuffer.data(idx_kin_df + 55);

    auto tk_zz_yyy = pbuffer.data(idx_kin_df + 56);

    auto tk_zz_yyz = pbuffer.data(idx_kin_df + 57);

    auto tk_zz_yzz = pbuffer.data(idx_kin_df + 58);

    auto tk_zz_zzz = pbuffer.data(idx_kin_df + 59);

    // Set up components of auxiliary buffer : FD

    auto tk_xxx_xx = pbuffer.data(idx_kin_fd);

    auto tk_xxx_xy = pbuffer.data(idx_kin_fd + 1);

    auto tk_xxx_xz = pbuffer.data(idx_kin_fd + 2);

    auto tk_xxx_yy = pbuffer.data(idx_kin_fd + 3);

    auto tk_xxx_yz = pbuffer.data(idx_kin_fd + 4);

    auto tk_xxx_zz = pbuffer.data(idx_kin_fd + 5);

    auto tk_xxz_xz = pbuffer.data(idx_kin_fd + 14);

    auto tk_xxz_yz = pbuffer.data(idx_kin_fd + 16);

    auto tk_xxz_zz = pbuffer.data(idx_kin_fd + 17);

    auto tk_xyy_xy = pbuffer.data(idx_kin_fd + 19);

    auto tk_xyy_yy = pbuffer.data(idx_kin_fd + 21);

    auto tk_xyy_yz = pbuffer.data(idx_kin_fd + 22);

    auto tk_xzz_xz = pbuffer.data(idx_kin_fd + 32);

    auto tk_xzz_yz = pbuffer.data(idx_kin_fd + 34);

    auto tk_xzz_zz = pbuffer.data(idx_kin_fd + 35);

    auto tk_yyy_xx = pbuffer.data(idx_kin_fd + 36);

    auto tk_yyy_xy = pbuffer.data(idx_kin_fd + 37);

    auto tk_yyy_xz = pbuffer.data(idx_kin_fd + 38);

    auto tk_yyy_yy = pbuffer.data(idx_kin_fd + 39);

    auto tk_yyy_yz = pbuffer.data(idx_kin_fd + 40);

    auto tk_yyy_zz = pbuffer.data(idx_kin_fd + 41);

    auto tk_yyz_xz = pbuffer.data(idx_kin_fd + 44);

    auto tk_yyz_yz = pbuffer.data(idx_kin_fd + 46);

    auto tk_yyz_zz = pbuffer.data(idx_kin_fd + 47);

    auto tk_yzz_xy = pbuffer.data(idx_kin_fd + 49);

    auto tk_yzz_xz = pbuffer.data(idx_kin_fd + 50);

    auto tk_yzz_yy = pbuffer.data(idx_kin_fd + 51);

    auto tk_yzz_yz = pbuffer.data(idx_kin_fd + 52);

    auto tk_yzz_zz = pbuffer.data(idx_kin_fd + 53);

    auto tk_zzz_xx = pbuffer.data(idx_kin_fd + 54);

    auto tk_zzz_xy = pbuffer.data(idx_kin_fd + 55);

    auto tk_zzz_xz = pbuffer.data(idx_kin_fd + 56);

    auto tk_zzz_yy = pbuffer.data(idx_kin_fd + 57);

    auto tk_zzz_yz = pbuffer.data(idx_kin_fd + 58);

    auto tk_zzz_zz = pbuffer.data(idx_kin_fd + 59);

    // Set up components of auxiliary buffer : FF

    auto tk_xxx_xxx = pbuffer.data(idx_kin_ff);

    auto tk_xxx_xxy = pbuffer.data(idx_kin_ff + 1);

    auto tk_xxx_xxz = pbuffer.data(idx_kin_ff + 2);

    auto tk_xxx_xyy = pbuffer.data(idx_kin_ff + 3);

    auto tk_xxx_xyz = pbuffer.data(idx_kin_ff + 4);

    auto tk_xxx_xzz = pbuffer.data(idx_kin_ff + 5);

    auto tk_xxx_yyy = pbuffer.data(idx_kin_ff + 6);

    auto tk_xxx_yyz = pbuffer.data(idx_kin_ff + 7);

    auto tk_xxx_yzz = pbuffer.data(idx_kin_ff + 8);

    auto tk_xxx_zzz = pbuffer.data(idx_kin_ff + 9);

    auto tk_xxy_xxx = pbuffer.data(idx_kin_ff + 10);

    auto tk_xxy_xxy = pbuffer.data(idx_kin_ff + 11);

    auto tk_xxy_xxz = pbuffer.data(idx_kin_ff + 12);

    auto tk_xxy_xyy = pbuffer.data(idx_kin_ff + 13);

    auto tk_xxy_xzz = pbuffer.data(idx_kin_ff + 15);

    auto tk_xxy_yyy = pbuffer.data(idx_kin_ff + 16);

    auto tk_xxz_xxx = pbuffer.data(idx_kin_ff + 20);

    auto tk_xxz_xxy = pbuffer.data(idx_kin_ff + 21);

    auto tk_xxz_xxz = pbuffer.data(idx_kin_ff + 22);

    auto tk_xxz_xyy = pbuffer.data(idx_kin_ff + 23);

    auto tk_xxz_xyz = pbuffer.data(idx_kin_ff + 24);

    auto tk_xxz_xzz = pbuffer.data(idx_kin_ff + 25);

    auto tk_xxz_yyz = pbuffer.data(idx_kin_ff + 27);

    auto tk_xxz_yzz = pbuffer.data(idx_kin_ff + 28);

    auto tk_xxz_zzz = pbuffer.data(idx_kin_ff + 29);

    auto tk_xyy_xxx = pbuffer.data(idx_kin_ff + 30);

    auto tk_xyy_xxy = pbuffer.data(idx_kin_ff + 31);

    auto tk_xyy_xyy = pbuffer.data(idx_kin_ff + 33);

    auto tk_xyy_xyz = pbuffer.data(idx_kin_ff + 34);

    auto tk_xyy_yyy = pbuffer.data(idx_kin_ff + 36);

    auto tk_xyy_yyz = pbuffer.data(idx_kin_ff + 37);

    auto tk_xyy_yzz = pbuffer.data(idx_kin_ff + 38);

    auto tk_xyy_zzz = pbuffer.data(idx_kin_ff + 39);

    auto tk_xzz_xxx = pbuffer.data(idx_kin_ff + 50);

    auto tk_xzz_xxz = pbuffer.data(idx_kin_ff + 52);

    auto tk_xzz_xyz = pbuffer.data(idx_kin_ff + 54);

    auto tk_xzz_xzz = pbuffer.data(idx_kin_ff + 55);

    auto tk_xzz_yyy = pbuffer.data(idx_kin_ff + 56);

    auto tk_xzz_yyz = pbuffer.data(idx_kin_ff + 57);

    auto tk_xzz_yzz = pbuffer.data(idx_kin_ff + 58);

    auto tk_xzz_zzz = pbuffer.data(idx_kin_ff + 59);

    auto tk_yyy_xxx = pbuffer.data(idx_kin_ff + 60);

    auto tk_yyy_xxy = pbuffer.data(idx_kin_ff + 61);

    auto tk_yyy_xxz = pbuffer.data(idx_kin_ff + 62);

    auto tk_yyy_xyy = pbuffer.data(idx_kin_ff + 63);

    auto tk_yyy_xyz = pbuffer.data(idx_kin_ff + 64);

    auto tk_yyy_xzz = pbuffer.data(idx_kin_ff + 65);

    auto tk_yyy_yyy = pbuffer.data(idx_kin_ff + 66);

    auto tk_yyy_yyz = pbuffer.data(idx_kin_ff + 67);

    auto tk_yyy_yzz = pbuffer.data(idx_kin_ff + 68);

    auto tk_yyy_zzz = pbuffer.data(idx_kin_ff + 69);

    auto tk_yyz_xxy = pbuffer.data(idx_kin_ff + 71);

    auto tk_yyz_xxz = pbuffer.data(idx_kin_ff + 72);

    auto tk_yyz_xyy = pbuffer.data(idx_kin_ff + 73);

    auto tk_yyz_xyz = pbuffer.data(idx_kin_ff + 74);

    auto tk_yyz_xzz = pbuffer.data(idx_kin_ff + 75);

    auto tk_yyz_yyy = pbuffer.data(idx_kin_ff + 76);

    auto tk_yyz_yyz = pbuffer.data(idx_kin_ff + 77);

    auto tk_yyz_yzz = pbuffer.data(idx_kin_ff + 78);

    auto tk_yyz_zzz = pbuffer.data(idx_kin_ff + 79);

    auto tk_yzz_xxx = pbuffer.data(idx_kin_ff + 80);

    auto tk_yzz_xxy = pbuffer.data(idx_kin_ff + 81);

    auto tk_yzz_xxz = pbuffer.data(idx_kin_ff + 82);

    auto tk_yzz_xyy = pbuffer.data(idx_kin_ff + 83);

    auto tk_yzz_xyz = pbuffer.data(idx_kin_ff + 84);

    auto tk_yzz_xzz = pbuffer.data(idx_kin_ff + 85);

    auto tk_yzz_yyy = pbuffer.data(idx_kin_ff + 86);

    auto tk_yzz_yyz = pbuffer.data(idx_kin_ff + 87);

    auto tk_yzz_yzz = pbuffer.data(idx_kin_ff + 88);

    auto tk_yzz_zzz = pbuffer.data(idx_kin_ff + 89);

    auto tk_zzz_xxx = pbuffer.data(idx_kin_ff + 90);

    auto tk_zzz_xxy = pbuffer.data(idx_kin_ff + 91);

    auto tk_zzz_xxz = pbuffer.data(idx_kin_ff + 92);

    auto tk_zzz_xyy = pbuffer.data(idx_kin_ff + 93);

    auto tk_zzz_xyz = pbuffer.data(idx_kin_ff + 94);

    auto tk_zzz_xzz = pbuffer.data(idx_kin_ff + 95);

    auto tk_zzz_yyy = pbuffer.data(idx_kin_ff + 96);

    auto tk_zzz_yyz = pbuffer.data(idx_kin_ff + 97);

    auto tk_zzz_yzz = pbuffer.data(idx_kin_ff + 98);

    auto tk_zzz_zzz = pbuffer.data(idx_kin_ff + 99);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_ovl_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_ovl_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_ovl_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_ovl_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_ovl_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_ovl_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_ovl_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_ovl_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_ovl_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_ovl_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_ovl_gf + 10);

    auto ts_xxxy_xxy = pbuffer.data(idx_ovl_gf + 11);

    auto ts_xxxy_xxz = pbuffer.data(idx_ovl_gf + 12);

    auto ts_xxxy_xyy = pbuffer.data(idx_ovl_gf + 13);

    auto ts_xxxy_xyz = pbuffer.data(idx_ovl_gf + 14);

    auto ts_xxxy_xzz = pbuffer.data(idx_ovl_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_ovl_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_ovl_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_ovl_gf + 18);

    auto ts_xxxy_zzz = pbuffer.data(idx_ovl_gf + 19);

    auto ts_xxxz_xxx = pbuffer.data(idx_ovl_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_ovl_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_ovl_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_ovl_gf + 23);

    auto ts_xxxz_xyz = pbuffer.data(idx_ovl_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_ovl_gf + 25);

    auto ts_xxxz_yyy = pbuffer.data(idx_ovl_gf + 26);

    auto ts_xxxz_yyz = pbuffer.data(idx_ovl_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_ovl_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_ovl_gf + 29);

    auto ts_xxyy_xxx = pbuffer.data(idx_ovl_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_ovl_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_ovl_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_ovl_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_ovl_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_ovl_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_ovl_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_ovl_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_ovl_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_ovl_gf + 39);

    auto ts_xxyz_xxx = pbuffer.data(idx_ovl_gf + 40);

    auto ts_xxyz_xxy = pbuffer.data(idx_ovl_gf + 41);

    auto ts_xxyz_xxz = pbuffer.data(idx_ovl_gf + 42);

    auto ts_xxyz_xyy = pbuffer.data(idx_ovl_gf + 43);

    auto ts_xxyz_xyz = pbuffer.data(idx_ovl_gf + 44);

    auto ts_xxyz_xzz = pbuffer.data(idx_ovl_gf + 45);

    auto ts_xxyz_yyy = pbuffer.data(idx_ovl_gf + 46);

    auto ts_xxyz_yyz = pbuffer.data(idx_ovl_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_ovl_gf + 48);

    auto ts_xxyz_zzz = pbuffer.data(idx_ovl_gf + 49);

    auto ts_xxzz_xxx = pbuffer.data(idx_ovl_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_ovl_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_ovl_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_ovl_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_ovl_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_ovl_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_ovl_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_ovl_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_ovl_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_ovl_gf + 59);

    auto ts_xyyy_xxx = pbuffer.data(idx_ovl_gf + 60);

    auto ts_xyyy_xxy = pbuffer.data(idx_ovl_gf + 61);

    auto ts_xyyy_xxz = pbuffer.data(idx_ovl_gf + 62);

    auto ts_xyyy_xyy = pbuffer.data(idx_ovl_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_ovl_gf + 64);

    auto ts_xyyy_xzz = pbuffer.data(idx_ovl_gf + 65);

    auto ts_xyyy_yyy = pbuffer.data(idx_ovl_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_ovl_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_ovl_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_ovl_gf + 69);

    auto ts_xyyz_xxx = pbuffer.data(idx_ovl_gf + 70);

    auto ts_xyyz_xxy = pbuffer.data(idx_ovl_gf + 71);

    auto ts_xyyz_xxz = pbuffer.data(idx_ovl_gf + 72);

    auto ts_xyyz_xyy = pbuffer.data(idx_ovl_gf + 73);

    auto ts_xyyz_xyz = pbuffer.data(idx_ovl_gf + 74);

    auto ts_xyyz_xzz = pbuffer.data(idx_ovl_gf + 75);

    auto ts_xyyz_yyy = pbuffer.data(idx_ovl_gf + 76);

    auto ts_xyyz_yyz = pbuffer.data(idx_ovl_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_ovl_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_ovl_gf + 79);

    auto ts_xyzz_xxx = pbuffer.data(idx_ovl_gf + 80);

    auto ts_xyzz_xxy = pbuffer.data(idx_ovl_gf + 81);

    auto ts_xyzz_xxz = pbuffer.data(idx_ovl_gf + 82);

    auto ts_xyzz_xyy = pbuffer.data(idx_ovl_gf + 83);

    auto ts_xyzz_xyz = pbuffer.data(idx_ovl_gf + 84);

    auto ts_xyzz_xzz = pbuffer.data(idx_ovl_gf + 85);

    auto ts_xyzz_yyy = pbuffer.data(idx_ovl_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_ovl_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_ovl_gf + 88);

    auto ts_xyzz_zzz = pbuffer.data(idx_ovl_gf + 89);

    auto ts_xzzz_xxx = pbuffer.data(idx_ovl_gf + 90);

    auto ts_xzzz_xxy = pbuffer.data(idx_ovl_gf + 91);

    auto ts_xzzz_xxz = pbuffer.data(idx_ovl_gf + 92);

    auto ts_xzzz_xyy = pbuffer.data(idx_ovl_gf + 93);

    auto ts_xzzz_xyz = pbuffer.data(idx_ovl_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_ovl_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_ovl_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_ovl_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_ovl_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_ovl_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_ovl_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_ovl_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_ovl_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_ovl_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_ovl_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_ovl_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_ovl_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_ovl_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_ovl_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_ovl_gf + 109);

    auto ts_yyyz_xxx = pbuffer.data(idx_ovl_gf + 110);

    auto ts_yyyz_xxy = pbuffer.data(idx_ovl_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_ovl_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_ovl_gf + 113);

    auto ts_yyyz_xyz = pbuffer.data(idx_ovl_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_ovl_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_ovl_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_ovl_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_ovl_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_ovl_gf + 119);

    auto ts_yyzz_xxx = pbuffer.data(idx_ovl_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_ovl_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_ovl_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_ovl_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_ovl_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_ovl_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_ovl_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_ovl_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_ovl_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_ovl_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_ovl_gf + 130);

    auto ts_yzzz_xxy = pbuffer.data(idx_ovl_gf + 131);

    auto ts_yzzz_xxz = pbuffer.data(idx_ovl_gf + 132);

    auto ts_yzzz_xyy = pbuffer.data(idx_ovl_gf + 133);

    auto ts_yzzz_xyz = pbuffer.data(idx_ovl_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_ovl_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_ovl_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_ovl_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_ovl_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_ovl_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_ovl_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_ovl_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_ovl_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_ovl_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_ovl_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_ovl_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_ovl_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_ovl_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_ovl_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_ovl_gf + 149);

    // Set up 0-10 components of targeted buffer : GF

    auto tk_xxxx_xxx = pbuffer.data(idx_kin_gf);

    auto tk_xxxx_xxy = pbuffer.data(idx_kin_gf + 1);

    auto tk_xxxx_xxz = pbuffer.data(idx_kin_gf + 2);

    auto tk_xxxx_xyy = pbuffer.data(idx_kin_gf + 3);

    auto tk_xxxx_xyz = pbuffer.data(idx_kin_gf + 4);

    auto tk_xxxx_xzz = pbuffer.data(idx_kin_gf + 5);

    auto tk_xxxx_yyy = pbuffer.data(idx_kin_gf + 6);

    auto tk_xxxx_yyz = pbuffer.data(idx_kin_gf + 7);

    auto tk_xxxx_yzz = pbuffer.data(idx_kin_gf + 8);

    auto tk_xxxx_zzz = pbuffer.data(idx_kin_gf + 9);

#pragma omp simd aligned(pa_x,            \
                             tk_xx_xxx,   \
                             tk_xx_xxy,   \
                             tk_xx_xxz,   \
                             tk_xx_xyy,   \
                             tk_xx_xyz,   \
                             tk_xx_xzz,   \
                             tk_xx_yyy,   \
                             tk_xx_yyz,   \
                             tk_xx_yzz,   \
                             tk_xx_zzz,   \
                             tk_xxx_xx,   \
                             tk_xxx_xxx,  \
                             tk_xxx_xxy,  \
                             tk_xxx_xxz,  \
                             tk_xxx_xy,   \
                             tk_xxx_xyy,  \
                             tk_xxx_xyz,  \
                             tk_xxx_xz,   \
                             tk_xxx_xzz,  \
                             tk_xxx_yy,   \
                             tk_xxx_yyy,  \
                             tk_xxx_yyz,  \
                             tk_xxx_yz,   \
                             tk_xxx_yzz,  \
                             tk_xxx_zz,   \
                             tk_xxx_zzz,  \
                             tk_xxxx_xxx, \
                             tk_xxxx_xxy, \
                             tk_xxxx_xxz, \
                             tk_xxxx_xyy, \
                             tk_xxxx_xyz, \
                             tk_xxxx_xzz, \
                             tk_xxxx_yyy, \
                             tk_xxxx_yyz, \
                             tk_xxxx_yzz, \
                             tk_xxxx_zzz, \
                             ts_xx_xxx,   \
                             ts_xx_xxy,   \
                             ts_xx_xxz,   \
                             ts_xx_xyy,   \
                             ts_xx_xyz,   \
                             ts_xx_xzz,   \
                             ts_xx_yyy,   \
                             ts_xx_yyz,   \
                             ts_xx_yzz,   \
                             ts_xx_zzz,   \
                             ts_xxxx_xxx, \
                             ts_xxxx_xxy, \
                             ts_xxxx_xxz, \
                             ts_xxxx_xyy, \
                             ts_xxxx_xyz, \
                             ts_xxxx_xzz, \
                             ts_xxxx_yyy, \
                             ts_xxxx_yyz, \
                             ts_xxxx_yzz, \
                             ts_xxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_xxx[i] = -6.0 * ts_xx_xxx[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxx[i] * fe_0 + 3.0 * tk_xxx_xx[i] * fe_0 + tk_xxx_xxx[i] * pa_x[i] +
                         2.0 * ts_xxxx_xxx[i] * fz_0;

        tk_xxxx_xxy[i] = -6.0 * ts_xx_xxy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxy[i] * fe_0 + 2.0 * tk_xxx_xy[i] * fe_0 + tk_xxx_xxy[i] * pa_x[i] +
                         2.0 * ts_xxxx_xxy[i] * fz_0;

        tk_xxxx_xxz[i] = -6.0 * ts_xx_xxz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxz[i] * fe_0 + 2.0 * tk_xxx_xz[i] * fe_0 + tk_xxx_xxz[i] * pa_x[i] +
                         2.0 * ts_xxxx_xxz[i] * fz_0;

        tk_xxxx_xyy[i] = -6.0 * ts_xx_xyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyy[i] * fe_0 + tk_xxx_yy[i] * fe_0 + tk_xxx_xyy[i] * pa_x[i] +
                         2.0 * ts_xxxx_xyy[i] * fz_0;

        tk_xxxx_xyz[i] = -6.0 * ts_xx_xyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyz[i] * fe_0 + tk_xxx_yz[i] * fe_0 + tk_xxx_xyz[i] * pa_x[i] +
                         2.0 * ts_xxxx_xyz[i] * fz_0;

        tk_xxxx_xzz[i] = -6.0 * ts_xx_xzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xzz[i] * fe_0 + tk_xxx_zz[i] * fe_0 + tk_xxx_xzz[i] * pa_x[i] +
                         2.0 * ts_xxxx_xzz[i] * fz_0;

        tk_xxxx_yyy[i] = -6.0 * ts_xx_yyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyy[i] * fe_0 + tk_xxx_yyy[i] * pa_x[i] + 2.0 * ts_xxxx_yyy[i] * fz_0;

        tk_xxxx_yyz[i] = -6.0 * ts_xx_yyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyz[i] * fe_0 + tk_xxx_yyz[i] * pa_x[i] + 2.0 * ts_xxxx_yyz[i] * fz_0;

        tk_xxxx_yzz[i] = -6.0 * ts_xx_yzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yzz[i] * fe_0 + tk_xxx_yzz[i] * pa_x[i] + 2.0 * ts_xxxx_yzz[i] * fz_0;

        tk_xxxx_zzz[i] = -6.0 * ts_xx_zzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_zzz[i] * fe_0 + tk_xxx_zzz[i] * pa_x[i] + 2.0 * ts_xxxx_zzz[i] * fz_0;
    }

    // Set up 10-20 components of targeted buffer : GF

    auto tk_xxxy_xxx = pbuffer.data(idx_kin_gf + 10);

    auto tk_xxxy_xxy = pbuffer.data(idx_kin_gf + 11);

    auto tk_xxxy_xxz = pbuffer.data(idx_kin_gf + 12);

    auto tk_xxxy_xyy = pbuffer.data(idx_kin_gf + 13);

    auto tk_xxxy_xyz = pbuffer.data(idx_kin_gf + 14);

    auto tk_xxxy_xzz = pbuffer.data(idx_kin_gf + 15);

    auto tk_xxxy_yyy = pbuffer.data(idx_kin_gf + 16);

    auto tk_xxxy_yyz = pbuffer.data(idx_kin_gf + 17);

    auto tk_xxxy_yzz = pbuffer.data(idx_kin_gf + 18);

    auto tk_xxxy_zzz = pbuffer.data(idx_kin_gf + 19);

#pragma omp simd aligned(pa_y,            \
                             tk_xxx_xx,   \
                             tk_xxx_xxx,  \
                             tk_xxx_xxy,  \
                             tk_xxx_xxz,  \
                             tk_xxx_xy,   \
                             tk_xxx_xyy,  \
                             tk_xxx_xyz,  \
                             tk_xxx_xz,   \
                             tk_xxx_xzz,  \
                             tk_xxx_yy,   \
                             tk_xxx_yyy,  \
                             tk_xxx_yyz,  \
                             tk_xxx_yz,   \
                             tk_xxx_yzz,  \
                             tk_xxx_zz,   \
                             tk_xxx_zzz,  \
                             tk_xxxy_xxx, \
                             tk_xxxy_xxy, \
                             tk_xxxy_xxz, \
                             tk_xxxy_xyy, \
                             tk_xxxy_xyz, \
                             tk_xxxy_xzz, \
                             tk_xxxy_yyy, \
                             tk_xxxy_yyz, \
                             tk_xxxy_yzz, \
                             tk_xxxy_zzz, \
                             ts_xxxy_xxx, \
                             ts_xxxy_xxy, \
                             ts_xxxy_xxz, \
                             ts_xxxy_xyy, \
                             ts_xxxy_xyz, \
                             ts_xxxy_xzz, \
                             ts_xxxy_yyy, \
                             ts_xxxy_yyz, \
                             ts_xxxy_yzz, \
                             ts_xxxy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxy_xxx[i] = tk_xxx_xxx[i] * pa_y[i] + 2.0 * ts_xxxy_xxx[i] * fz_0;

        tk_xxxy_xxy[i] = tk_xxx_xx[i] * fe_0 + tk_xxx_xxy[i] * pa_y[i] + 2.0 * ts_xxxy_xxy[i] * fz_0;

        tk_xxxy_xxz[i] = tk_xxx_xxz[i] * pa_y[i] + 2.0 * ts_xxxy_xxz[i] * fz_0;

        tk_xxxy_xyy[i] = 2.0 * tk_xxx_xy[i] * fe_0 + tk_xxx_xyy[i] * pa_y[i] + 2.0 * ts_xxxy_xyy[i] * fz_0;

        tk_xxxy_xyz[i] = tk_xxx_xz[i] * fe_0 + tk_xxx_xyz[i] * pa_y[i] + 2.0 * ts_xxxy_xyz[i] * fz_0;

        tk_xxxy_xzz[i] = tk_xxx_xzz[i] * pa_y[i] + 2.0 * ts_xxxy_xzz[i] * fz_0;

        tk_xxxy_yyy[i] = 3.0 * tk_xxx_yy[i] * fe_0 + tk_xxx_yyy[i] * pa_y[i] + 2.0 * ts_xxxy_yyy[i] * fz_0;

        tk_xxxy_yyz[i] = 2.0 * tk_xxx_yz[i] * fe_0 + tk_xxx_yyz[i] * pa_y[i] + 2.0 * ts_xxxy_yyz[i] * fz_0;

        tk_xxxy_yzz[i] = tk_xxx_zz[i] * fe_0 + tk_xxx_yzz[i] * pa_y[i] + 2.0 * ts_xxxy_yzz[i] * fz_0;

        tk_xxxy_zzz[i] = tk_xxx_zzz[i] * pa_y[i] + 2.0 * ts_xxxy_zzz[i] * fz_0;
    }

    // Set up 20-30 components of targeted buffer : GF

    auto tk_xxxz_xxx = pbuffer.data(idx_kin_gf + 20);

    auto tk_xxxz_xxy = pbuffer.data(idx_kin_gf + 21);

    auto tk_xxxz_xxz = pbuffer.data(idx_kin_gf + 22);

    auto tk_xxxz_xyy = pbuffer.data(idx_kin_gf + 23);

    auto tk_xxxz_xyz = pbuffer.data(idx_kin_gf + 24);

    auto tk_xxxz_xzz = pbuffer.data(idx_kin_gf + 25);

    auto tk_xxxz_yyy = pbuffer.data(idx_kin_gf + 26);

    auto tk_xxxz_yyz = pbuffer.data(idx_kin_gf + 27);

    auto tk_xxxz_yzz = pbuffer.data(idx_kin_gf + 28);

    auto tk_xxxz_zzz = pbuffer.data(idx_kin_gf + 29);

#pragma omp simd aligned(pa_z,            \
                             tk_xxx_xx,   \
                             tk_xxx_xxx,  \
                             tk_xxx_xxy,  \
                             tk_xxx_xxz,  \
                             tk_xxx_xy,   \
                             tk_xxx_xyy,  \
                             tk_xxx_xyz,  \
                             tk_xxx_xz,   \
                             tk_xxx_xzz,  \
                             tk_xxx_yy,   \
                             tk_xxx_yyy,  \
                             tk_xxx_yyz,  \
                             tk_xxx_yz,   \
                             tk_xxx_yzz,  \
                             tk_xxx_zz,   \
                             tk_xxx_zzz,  \
                             tk_xxxz_xxx, \
                             tk_xxxz_xxy, \
                             tk_xxxz_xxz, \
                             tk_xxxz_xyy, \
                             tk_xxxz_xyz, \
                             tk_xxxz_xzz, \
                             tk_xxxz_yyy, \
                             tk_xxxz_yyz, \
                             tk_xxxz_yzz, \
                             tk_xxxz_zzz, \
                             ts_xxxz_xxx, \
                             ts_xxxz_xxy, \
                             ts_xxxz_xxz, \
                             ts_xxxz_xyy, \
                             ts_xxxz_xyz, \
                             ts_xxxz_xzz, \
                             ts_xxxz_yyy, \
                             ts_xxxz_yyz, \
                             ts_xxxz_yzz, \
                             ts_xxxz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxz_xxx[i] = tk_xxx_xxx[i] * pa_z[i] + 2.0 * ts_xxxz_xxx[i] * fz_0;

        tk_xxxz_xxy[i] = tk_xxx_xxy[i] * pa_z[i] + 2.0 * ts_xxxz_xxy[i] * fz_0;

        tk_xxxz_xxz[i] = tk_xxx_xx[i] * fe_0 + tk_xxx_xxz[i] * pa_z[i] + 2.0 * ts_xxxz_xxz[i] * fz_0;

        tk_xxxz_xyy[i] = tk_xxx_xyy[i] * pa_z[i] + 2.0 * ts_xxxz_xyy[i] * fz_0;

        tk_xxxz_xyz[i] = tk_xxx_xy[i] * fe_0 + tk_xxx_xyz[i] * pa_z[i] + 2.0 * ts_xxxz_xyz[i] * fz_0;

        tk_xxxz_xzz[i] = 2.0 * tk_xxx_xz[i] * fe_0 + tk_xxx_xzz[i] * pa_z[i] + 2.0 * ts_xxxz_xzz[i] * fz_0;

        tk_xxxz_yyy[i] = tk_xxx_yyy[i] * pa_z[i] + 2.0 * ts_xxxz_yyy[i] * fz_0;

        tk_xxxz_yyz[i] = tk_xxx_yy[i] * fe_0 + tk_xxx_yyz[i] * pa_z[i] + 2.0 * ts_xxxz_yyz[i] * fz_0;

        tk_xxxz_yzz[i] = 2.0 * tk_xxx_yz[i] * fe_0 + tk_xxx_yzz[i] * pa_z[i] + 2.0 * ts_xxxz_yzz[i] * fz_0;

        tk_xxxz_zzz[i] = 3.0 * tk_xxx_zz[i] * fe_0 + tk_xxx_zzz[i] * pa_z[i] + 2.0 * ts_xxxz_zzz[i] * fz_0;
    }

    // Set up 30-40 components of targeted buffer : GF

    auto tk_xxyy_xxx = pbuffer.data(idx_kin_gf + 30);

    auto tk_xxyy_xxy = pbuffer.data(idx_kin_gf + 31);

    auto tk_xxyy_xxz = pbuffer.data(idx_kin_gf + 32);

    auto tk_xxyy_xyy = pbuffer.data(idx_kin_gf + 33);

    auto tk_xxyy_xyz = pbuffer.data(idx_kin_gf + 34);

    auto tk_xxyy_xzz = pbuffer.data(idx_kin_gf + 35);

    auto tk_xxyy_yyy = pbuffer.data(idx_kin_gf + 36);

    auto tk_xxyy_yyz = pbuffer.data(idx_kin_gf + 37);

    auto tk_xxyy_yzz = pbuffer.data(idx_kin_gf + 38);

    auto tk_xxyy_zzz = pbuffer.data(idx_kin_gf + 39);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tk_xx_xxx,   \
                             tk_xx_xxz,   \
                             tk_xx_xzz,   \
                             tk_xxy_xxx,  \
                             tk_xxy_xxz,  \
                             tk_xxy_xzz,  \
                             tk_xxyy_xxx, \
                             tk_xxyy_xxy, \
                             tk_xxyy_xxz, \
                             tk_xxyy_xyy, \
                             tk_xxyy_xyz, \
                             tk_xxyy_xzz, \
                             tk_xxyy_yyy, \
                             tk_xxyy_yyz, \
                             tk_xxyy_yzz, \
                             tk_xxyy_zzz, \
                             tk_xyy_xxy,  \
                             tk_xyy_xy,   \
                             tk_xyy_xyy,  \
                             tk_xyy_xyz,  \
                             tk_xyy_yy,   \
                             tk_xyy_yyy,  \
                             tk_xyy_yyz,  \
                             tk_xyy_yz,   \
                             tk_xyy_yzz,  \
                             tk_xyy_zzz,  \
                             tk_yy_xxy,   \
                             tk_yy_xyy,   \
                             tk_yy_xyz,   \
                             tk_yy_yyy,   \
                             tk_yy_yyz,   \
                             tk_yy_yzz,   \
                             tk_yy_zzz,   \
                             ts_xx_xxx,   \
                             ts_xx_xxz,   \
                             ts_xx_xzz,   \
                             ts_xxyy_xxx, \
                             ts_xxyy_xxy, \
                             ts_xxyy_xxz, \
                             ts_xxyy_xyy, \
                             ts_xxyy_xyz, \
                             ts_xxyy_xzz, \
                             ts_xxyy_yyy, \
                             ts_xxyy_yyz, \
                             ts_xxyy_yzz, \
                             ts_xxyy_zzz, \
                             ts_yy_xxy,   \
                             ts_yy_xyy,   \
                             ts_yy_xyz,   \
                             ts_yy_yyy,   \
                             ts_yy_yyz,   \
                             ts_yy_yzz,   \
                             ts_yy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyy_xxx[i] = -2.0 * ts_xx_xxx[i] * fbe_0 * fz_0 + tk_xx_xxx[i] * fe_0 + tk_xxy_xxx[i] * pa_y[i] + 2.0 * ts_xxyy_xxx[i] * fz_0;

        tk_xxyy_xxy[i] = -2.0 * ts_yy_xxy[i] * fbe_0 * fz_0 + tk_yy_xxy[i] * fe_0 + 2.0 * tk_xyy_xy[i] * fe_0 + tk_xyy_xxy[i] * pa_x[i] +
                         2.0 * ts_xxyy_xxy[i] * fz_0;

        tk_xxyy_xxz[i] = -2.0 * ts_xx_xxz[i] * fbe_0 * fz_0 + tk_xx_xxz[i] * fe_0 + tk_xxy_xxz[i] * pa_y[i] + 2.0 * ts_xxyy_xxz[i] * fz_0;

        tk_xxyy_xyy[i] =
            -2.0 * ts_yy_xyy[i] * fbe_0 * fz_0 + tk_yy_xyy[i] * fe_0 + tk_xyy_yy[i] * fe_0 + tk_xyy_xyy[i] * pa_x[i] + 2.0 * ts_xxyy_xyy[i] * fz_0;

        tk_xxyy_xyz[i] =
            -2.0 * ts_yy_xyz[i] * fbe_0 * fz_0 + tk_yy_xyz[i] * fe_0 + tk_xyy_yz[i] * fe_0 + tk_xyy_xyz[i] * pa_x[i] + 2.0 * ts_xxyy_xyz[i] * fz_0;

        tk_xxyy_xzz[i] = -2.0 * ts_xx_xzz[i] * fbe_0 * fz_0 + tk_xx_xzz[i] * fe_0 + tk_xxy_xzz[i] * pa_y[i] + 2.0 * ts_xxyy_xzz[i] * fz_0;

        tk_xxyy_yyy[i] = -2.0 * ts_yy_yyy[i] * fbe_0 * fz_0 + tk_yy_yyy[i] * fe_0 + tk_xyy_yyy[i] * pa_x[i] + 2.0 * ts_xxyy_yyy[i] * fz_0;

        tk_xxyy_yyz[i] = -2.0 * ts_yy_yyz[i] * fbe_0 * fz_0 + tk_yy_yyz[i] * fe_0 + tk_xyy_yyz[i] * pa_x[i] + 2.0 * ts_xxyy_yyz[i] * fz_0;

        tk_xxyy_yzz[i] = -2.0 * ts_yy_yzz[i] * fbe_0 * fz_0 + tk_yy_yzz[i] * fe_0 + tk_xyy_yzz[i] * pa_x[i] + 2.0 * ts_xxyy_yzz[i] * fz_0;

        tk_xxyy_zzz[i] = -2.0 * ts_yy_zzz[i] * fbe_0 * fz_0 + tk_yy_zzz[i] * fe_0 + tk_xyy_zzz[i] * pa_x[i] + 2.0 * ts_xxyy_zzz[i] * fz_0;
    }

    // Set up 40-50 components of targeted buffer : GF

    auto tk_xxyz_xxx = pbuffer.data(idx_kin_gf + 40);

    auto tk_xxyz_xxy = pbuffer.data(idx_kin_gf + 41);

    auto tk_xxyz_xxz = pbuffer.data(idx_kin_gf + 42);

    auto tk_xxyz_xyy = pbuffer.data(idx_kin_gf + 43);

    auto tk_xxyz_xyz = pbuffer.data(idx_kin_gf + 44);

    auto tk_xxyz_xzz = pbuffer.data(idx_kin_gf + 45);

    auto tk_xxyz_yyy = pbuffer.data(idx_kin_gf + 46);

    auto tk_xxyz_yyz = pbuffer.data(idx_kin_gf + 47);

    auto tk_xxyz_yzz = pbuffer.data(idx_kin_gf + 48);

    auto tk_xxyz_zzz = pbuffer.data(idx_kin_gf + 49);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tk_xxy_xxy,  \
                             tk_xxy_xyy,  \
                             tk_xxy_yyy,  \
                             tk_xxyz_xxx, \
                             tk_xxyz_xxy, \
                             tk_xxyz_xxz, \
                             tk_xxyz_xyy, \
                             tk_xxyz_xyz, \
                             tk_xxyz_xzz, \
                             tk_xxyz_yyy, \
                             tk_xxyz_yyz, \
                             tk_xxyz_yzz, \
                             tk_xxyz_zzz, \
                             tk_xxz_xxx,  \
                             tk_xxz_xxz,  \
                             tk_xxz_xyz,  \
                             tk_xxz_xz,   \
                             tk_xxz_xzz,  \
                             tk_xxz_yyz,  \
                             tk_xxz_yz,   \
                             tk_xxz_yzz,  \
                             tk_xxz_zz,   \
                             tk_xxz_zzz,  \
                             ts_xxyz_xxx, \
                             ts_xxyz_xxy, \
                             ts_xxyz_xxz, \
                             ts_xxyz_xyy, \
                             ts_xxyz_xyz, \
                             ts_xxyz_xzz, \
                             ts_xxyz_yyy, \
                             ts_xxyz_yyz, \
                             ts_xxyz_yzz, \
                             ts_xxyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyz_xxx[i] = tk_xxz_xxx[i] * pa_y[i] + 2.0 * ts_xxyz_xxx[i] * fz_0;

        tk_xxyz_xxy[i] = tk_xxy_xxy[i] * pa_z[i] + 2.0 * ts_xxyz_xxy[i] * fz_0;

        tk_xxyz_xxz[i] = tk_xxz_xxz[i] * pa_y[i] + 2.0 * ts_xxyz_xxz[i] * fz_0;

        tk_xxyz_xyy[i] = tk_xxy_xyy[i] * pa_z[i] + 2.0 * ts_xxyz_xyy[i] * fz_0;

        tk_xxyz_xyz[i] = tk_xxz_xz[i] * fe_0 + tk_xxz_xyz[i] * pa_y[i] + 2.0 * ts_xxyz_xyz[i] * fz_0;

        tk_xxyz_xzz[i] = tk_xxz_xzz[i] * pa_y[i] + 2.0 * ts_xxyz_xzz[i] * fz_0;

        tk_xxyz_yyy[i] = tk_xxy_yyy[i] * pa_z[i] + 2.0 * ts_xxyz_yyy[i] * fz_0;

        tk_xxyz_yyz[i] = 2.0 * tk_xxz_yz[i] * fe_0 + tk_xxz_yyz[i] * pa_y[i] + 2.0 * ts_xxyz_yyz[i] * fz_0;

        tk_xxyz_yzz[i] = tk_xxz_zz[i] * fe_0 + tk_xxz_yzz[i] * pa_y[i] + 2.0 * ts_xxyz_yzz[i] * fz_0;

        tk_xxyz_zzz[i] = tk_xxz_zzz[i] * pa_y[i] + 2.0 * ts_xxyz_zzz[i] * fz_0;
    }

    // Set up 50-60 components of targeted buffer : GF

    auto tk_xxzz_xxx = pbuffer.data(idx_kin_gf + 50);

    auto tk_xxzz_xxy = pbuffer.data(idx_kin_gf + 51);

    auto tk_xxzz_xxz = pbuffer.data(idx_kin_gf + 52);

    auto tk_xxzz_xyy = pbuffer.data(idx_kin_gf + 53);

    auto tk_xxzz_xyz = pbuffer.data(idx_kin_gf + 54);

    auto tk_xxzz_xzz = pbuffer.data(idx_kin_gf + 55);

    auto tk_xxzz_yyy = pbuffer.data(idx_kin_gf + 56);

    auto tk_xxzz_yyz = pbuffer.data(idx_kin_gf + 57);

    auto tk_xxzz_yzz = pbuffer.data(idx_kin_gf + 58);

    auto tk_xxzz_zzz = pbuffer.data(idx_kin_gf + 59);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tk_xx_xxx,   \
                             tk_xx_xxy,   \
                             tk_xx_xyy,   \
                             tk_xxz_xxx,  \
                             tk_xxz_xxy,  \
                             tk_xxz_xyy,  \
                             tk_xxzz_xxx, \
                             tk_xxzz_xxy, \
                             tk_xxzz_xxz, \
                             tk_xxzz_xyy, \
                             tk_xxzz_xyz, \
                             tk_xxzz_xzz, \
                             tk_xxzz_yyy, \
                             tk_xxzz_yyz, \
                             tk_xxzz_yzz, \
                             tk_xxzz_zzz, \
                             tk_xzz_xxz,  \
                             tk_xzz_xyz,  \
                             tk_xzz_xz,   \
                             tk_xzz_xzz,  \
                             tk_xzz_yyy,  \
                             tk_xzz_yyz,  \
                             tk_xzz_yz,   \
                             tk_xzz_yzz,  \
                             tk_xzz_zz,   \
                             tk_xzz_zzz,  \
                             tk_zz_xxz,   \
                             tk_zz_xyz,   \
                             tk_zz_xzz,   \
                             tk_zz_yyy,   \
                             tk_zz_yyz,   \
                             tk_zz_yzz,   \
                             tk_zz_zzz,   \
                             ts_xx_xxx,   \
                             ts_xx_xxy,   \
                             ts_xx_xyy,   \
                             ts_xxzz_xxx, \
                             ts_xxzz_xxy, \
                             ts_xxzz_xxz, \
                             ts_xxzz_xyy, \
                             ts_xxzz_xyz, \
                             ts_xxzz_xzz, \
                             ts_xxzz_yyy, \
                             ts_xxzz_yyz, \
                             ts_xxzz_yzz, \
                             ts_xxzz_zzz, \
                             ts_zz_xxz,   \
                             ts_zz_xyz,   \
                             ts_zz_xzz,   \
                             ts_zz_yyy,   \
                             ts_zz_yyz,   \
                             ts_zz_yzz,   \
                             ts_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzz_xxx[i] = -2.0 * ts_xx_xxx[i] * fbe_0 * fz_0 + tk_xx_xxx[i] * fe_0 + tk_xxz_xxx[i] * pa_z[i] + 2.0 * ts_xxzz_xxx[i] * fz_0;

        tk_xxzz_xxy[i] = -2.0 * ts_xx_xxy[i] * fbe_0 * fz_0 + tk_xx_xxy[i] * fe_0 + tk_xxz_xxy[i] * pa_z[i] + 2.0 * ts_xxzz_xxy[i] * fz_0;

        tk_xxzz_xxz[i] = -2.0 * ts_zz_xxz[i] * fbe_0 * fz_0 + tk_zz_xxz[i] * fe_0 + 2.0 * tk_xzz_xz[i] * fe_0 + tk_xzz_xxz[i] * pa_x[i] +
                         2.0 * ts_xxzz_xxz[i] * fz_0;

        tk_xxzz_xyy[i] = -2.0 * ts_xx_xyy[i] * fbe_0 * fz_0 + tk_xx_xyy[i] * fe_0 + tk_xxz_xyy[i] * pa_z[i] + 2.0 * ts_xxzz_xyy[i] * fz_0;

        tk_xxzz_xyz[i] =
            -2.0 * ts_zz_xyz[i] * fbe_0 * fz_0 + tk_zz_xyz[i] * fe_0 + tk_xzz_yz[i] * fe_0 + tk_xzz_xyz[i] * pa_x[i] + 2.0 * ts_xxzz_xyz[i] * fz_0;

        tk_xxzz_xzz[i] =
            -2.0 * ts_zz_xzz[i] * fbe_0 * fz_0 + tk_zz_xzz[i] * fe_0 + tk_xzz_zz[i] * fe_0 + tk_xzz_xzz[i] * pa_x[i] + 2.0 * ts_xxzz_xzz[i] * fz_0;

        tk_xxzz_yyy[i] = -2.0 * ts_zz_yyy[i] * fbe_0 * fz_0 + tk_zz_yyy[i] * fe_0 + tk_xzz_yyy[i] * pa_x[i] + 2.0 * ts_xxzz_yyy[i] * fz_0;

        tk_xxzz_yyz[i] = -2.0 * ts_zz_yyz[i] * fbe_0 * fz_0 + tk_zz_yyz[i] * fe_0 + tk_xzz_yyz[i] * pa_x[i] + 2.0 * ts_xxzz_yyz[i] * fz_0;

        tk_xxzz_yzz[i] = -2.0 * ts_zz_yzz[i] * fbe_0 * fz_0 + tk_zz_yzz[i] * fe_0 + tk_xzz_yzz[i] * pa_x[i] + 2.0 * ts_xxzz_yzz[i] * fz_0;

        tk_xxzz_zzz[i] = -2.0 * ts_zz_zzz[i] * fbe_0 * fz_0 + tk_zz_zzz[i] * fe_0 + tk_xzz_zzz[i] * pa_x[i] + 2.0 * ts_xxzz_zzz[i] * fz_0;
    }

    // Set up 60-70 components of targeted buffer : GF

    auto tk_xyyy_xxx = pbuffer.data(idx_kin_gf + 60);

    auto tk_xyyy_xxy = pbuffer.data(idx_kin_gf + 61);

    auto tk_xyyy_xxz = pbuffer.data(idx_kin_gf + 62);

    auto tk_xyyy_xyy = pbuffer.data(idx_kin_gf + 63);

    auto tk_xyyy_xyz = pbuffer.data(idx_kin_gf + 64);

    auto tk_xyyy_xzz = pbuffer.data(idx_kin_gf + 65);

    auto tk_xyyy_yyy = pbuffer.data(idx_kin_gf + 66);

    auto tk_xyyy_yyz = pbuffer.data(idx_kin_gf + 67);

    auto tk_xyyy_yzz = pbuffer.data(idx_kin_gf + 68);

    auto tk_xyyy_zzz = pbuffer.data(idx_kin_gf + 69);

#pragma omp simd aligned(pa_x,            \
                             tk_xyyy_xxx, \
                             tk_xyyy_xxy, \
                             tk_xyyy_xxz, \
                             tk_xyyy_xyy, \
                             tk_xyyy_xyz, \
                             tk_xyyy_xzz, \
                             tk_xyyy_yyy, \
                             tk_xyyy_yyz, \
                             tk_xyyy_yzz, \
                             tk_xyyy_zzz, \
                             tk_yyy_xx,   \
                             tk_yyy_xxx,  \
                             tk_yyy_xxy,  \
                             tk_yyy_xxz,  \
                             tk_yyy_xy,   \
                             tk_yyy_xyy,  \
                             tk_yyy_xyz,  \
                             tk_yyy_xz,   \
                             tk_yyy_xzz,  \
                             tk_yyy_yy,   \
                             tk_yyy_yyy,  \
                             tk_yyy_yyz,  \
                             tk_yyy_yz,   \
                             tk_yyy_yzz,  \
                             tk_yyy_zz,   \
                             tk_yyy_zzz,  \
                             ts_xyyy_xxx, \
                             ts_xyyy_xxy, \
                             ts_xyyy_xxz, \
                             ts_xyyy_xyy, \
                             ts_xyyy_xyz, \
                             ts_xyyy_xzz, \
                             ts_xyyy_yyy, \
                             ts_xyyy_yyz, \
                             ts_xyyy_yzz, \
                             ts_xyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyy_xxx[i] = 3.0 * tk_yyy_xx[i] * fe_0 + tk_yyy_xxx[i] * pa_x[i] + 2.0 * ts_xyyy_xxx[i] * fz_0;

        tk_xyyy_xxy[i] = 2.0 * tk_yyy_xy[i] * fe_0 + tk_yyy_xxy[i] * pa_x[i] + 2.0 * ts_xyyy_xxy[i] * fz_0;

        tk_xyyy_xxz[i] = 2.0 * tk_yyy_xz[i] * fe_0 + tk_yyy_xxz[i] * pa_x[i] + 2.0 * ts_xyyy_xxz[i] * fz_0;

        tk_xyyy_xyy[i] = tk_yyy_yy[i] * fe_0 + tk_yyy_xyy[i] * pa_x[i] + 2.0 * ts_xyyy_xyy[i] * fz_0;

        tk_xyyy_xyz[i] = tk_yyy_yz[i] * fe_0 + tk_yyy_xyz[i] * pa_x[i] + 2.0 * ts_xyyy_xyz[i] * fz_0;

        tk_xyyy_xzz[i] = tk_yyy_zz[i] * fe_0 + tk_yyy_xzz[i] * pa_x[i] + 2.0 * ts_xyyy_xzz[i] * fz_0;

        tk_xyyy_yyy[i] = tk_yyy_yyy[i] * pa_x[i] + 2.0 * ts_xyyy_yyy[i] * fz_0;

        tk_xyyy_yyz[i] = tk_yyy_yyz[i] * pa_x[i] + 2.0 * ts_xyyy_yyz[i] * fz_0;

        tk_xyyy_yzz[i] = tk_yyy_yzz[i] * pa_x[i] + 2.0 * ts_xyyy_yzz[i] * fz_0;

        tk_xyyy_zzz[i] = tk_yyy_zzz[i] * pa_x[i] + 2.0 * ts_xyyy_zzz[i] * fz_0;
    }

    // Set up 70-80 components of targeted buffer : GF

    auto tk_xyyz_xxx = pbuffer.data(idx_kin_gf + 70);

    auto tk_xyyz_xxy = pbuffer.data(idx_kin_gf + 71);

    auto tk_xyyz_xxz = pbuffer.data(idx_kin_gf + 72);

    auto tk_xyyz_xyy = pbuffer.data(idx_kin_gf + 73);

    auto tk_xyyz_xyz = pbuffer.data(idx_kin_gf + 74);

    auto tk_xyyz_xzz = pbuffer.data(idx_kin_gf + 75);

    auto tk_xyyz_yyy = pbuffer.data(idx_kin_gf + 76);

    auto tk_xyyz_yyz = pbuffer.data(idx_kin_gf + 77);

    auto tk_xyyz_yzz = pbuffer.data(idx_kin_gf + 78);

    auto tk_xyyz_zzz = pbuffer.data(idx_kin_gf + 79);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tk_xyy_xxx,  \
                             tk_xyy_xxy,  \
                             tk_xyy_xyy,  \
                             tk_xyyz_xxx, \
                             tk_xyyz_xxy, \
                             tk_xyyz_xxz, \
                             tk_xyyz_xyy, \
                             tk_xyyz_xyz, \
                             tk_xyyz_xzz, \
                             tk_xyyz_yyy, \
                             tk_xyyz_yyz, \
                             tk_xyyz_yzz, \
                             tk_xyyz_zzz, \
                             tk_yyz_xxz,  \
                             tk_yyz_xyz,  \
                             tk_yyz_xz,   \
                             tk_yyz_xzz,  \
                             tk_yyz_yyy,  \
                             tk_yyz_yyz,  \
                             tk_yyz_yz,   \
                             tk_yyz_yzz,  \
                             tk_yyz_zz,   \
                             tk_yyz_zzz,  \
                             ts_xyyz_xxx, \
                             ts_xyyz_xxy, \
                             ts_xyyz_xxz, \
                             ts_xyyz_xyy, \
                             ts_xyyz_xyz, \
                             ts_xyyz_xzz, \
                             ts_xyyz_yyy, \
                             ts_xyyz_yyz, \
                             ts_xyyz_yzz, \
                             ts_xyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyz_xxx[i] = tk_xyy_xxx[i] * pa_z[i] + 2.0 * ts_xyyz_xxx[i] * fz_0;

        tk_xyyz_xxy[i] = tk_xyy_xxy[i] * pa_z[i] + 2.0 * ts_xyyz_xxy[i] * fz_0;

        tk_xyyz_xxz[i] = 2.0 * tk_yyz_xz[i] * fe_0 + tk_yyz_xxz[i] * pa_x[i] + 2.0 * ts_xyyz_xxz[i] * fz_0;

        tk_xyyz_xyy[i] = tk_xyy_xyy[i] * pa_z[i] + 2.0 * ts_xyyz_xyy[i] * fz_0;

        tk_xyyz_xyz[i] = tk_yyz_yz[i] * fe_0 + tk_yyz_xyz[i] * pa_x[i] + 2.0 * ts_xyyz_xyz[i] * fz_0;

        tk_xyyz_xzz[i] = tk_yyz_zz[i] * fe_0 + tk_yyz_xzz[i] * pa_x[i] + 2.0 * ts_xyyz_xzz[i] * fz_0;

        tk_xyyz_yyy[i] = tk_yyz_yyy[i] * pa_x[i] + 2.0 * ts_xyyz_yyy[i] * fz_0;

        tk_xyyz_yyz[i] = tk_yyz_yyz[i] * pa_x[i] + 2.0 * ts_xyyz_yyz[i] * fz_0;

        tk_xyyz_yzz[i] = tk_yyz_yzz[i] * pa_x[i] + 2.0 * ts_xyyz_yzz[i] * fz_0;

        tk_xyyz_zzz[i] = tk_yyz_zzz[i] * pa_x[i] + 2.0 * ts_xyyz_zzz[i] * fz_0;
    }

    // Set up 80-90 components of targeted buffer : GF

    auto tk_xyzz_xxx = pbuffer.data(idx_kin_gf + 80);

    auto tk_xyzz_xxy = pbuffer.data(idx_kin_gf + 81);

    auto tk_xyzz_xxz = pbuffer.data(idx_kin_gf + 82);

    auto tk_xyzz_xyy = pbuffer.data(idx_kin_gf + 83);

    auto tk_xyzz_xyz = pbuffer.data(idx_kin_gf + 84);

    auto tk_xyzz_xzz = pbuffer.data(idx_kin_gf + 85);

    auto tk_xyzz_yyy = pbuffer.data(idx_kin_gf + 86);

    auto tk_xyzz_yyz = pbuffer.data(idx_kin_gf + 87);

    auto tk_xyzz_yzz = pbuffer.data(idx_kin_gf + 88);

    auto tk_xyzz_zzz = pbuffer.data(idx_kin_gf + 89);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tk_xyzz_xxx, \
                             tk_xyzz_xxy, \
                             tk_xyzz_xxz, \
                             tk_xyzz_xyy, \
                             tk_xyzz_xyz, \
                             tk_xyzz_xzz, \
                             tk_xyzz_yyy, \
                             tk_xyzz_yyz, \
                             tk_xyzz_yzz, \
                             tk_xyzz_zzz, \
                             tk_xzz_xxx,  \
                             tk_xzz_xxz,  \
                             tk_xzz_xzz,  \
                             tk_yzz_xxy,  \
                             tk_yzz_xy,   \
                             tk_yzz_xyy,  \
                             tk_yzz_xyz,  \
                             tk_yzz_yy,   \
                             tk_yzz_yyy,  \
                             tk_yzz_yyz,  \
                             tk_yzz_yz,   \
                             tk_yzz_yzz,  \
                             tk_yzz_zzz,  \
                             ts_xyzz_xxx, \
                             ts_xyzz_xxy, \
                             ts_xyzz_xxz, \
                             ts_xyzz_xyy, \
                             ts_xyzz_xyz, \
                             ts_xyzz_xzz, \
                             ts_xyzz_yyy, \
                             ts_xyzz_yyz, \
                             ts_xyzz_yzz, \
                             ts_xyzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzz_xxx[i] = tk_xzz_xxx[i] * pa_y[i] + 2.0 * ts_xyzz_xxx[i] * fz_0;

        tk_xyzz_xxy[i] = 2.0 * tk_yzz_xy[i] * fe_0 + tk_yzz_xxy[i] * pa_x[i] + 2.0 * ts_xyzz_xxy[i] * fz_0;

        tk_xyzz_xxz[i] = tk_xzz_xxz[i] * pa_y[i] + 2.0 * ts_xyzz_xxz[i] * fz_0;

        tk_xyzz_xyy[i] = tk_yzz_yy[i] * fe_0 + tk_yzz_xyy[i] * pa_x[i] + 2.0 * ts_xyzz_xyy[i] * fz_0;

        tk_xyzz_xyz[i] = tk_yzz_yz[i] * fe_0 + tk_yzz_xyz[i] * pa_x[i] + 2.0 * ts_xyzz_xyz[i] * fz_0;

        tk_xyzz_xzz[i] = tk_xzz_xzz[i] * pa_y[i] + 2.0 * ts_xyzz_xzz[i] * fz_0;

        tk_xyzz_yyy[i] = tk_yzz_yyy[i] * pa_x[i] + 2.0 * ts_xyzz_yyy[i] * fz_0;

        tk_xyzz_yyz[i] = tk_yzz_yyz[i] * pa_x[i] + 2.0 * ts_xyzz_yyz[i] * fz_0;

        tk_xyzz_yzz[i] = tk_yzz_yzz[i] * pa_x[i] + 2.0 * ts_xyzz_yzz[i] * fz_0;

        tk_xyzz_zzz[i] = tk_yzz_zzz[i] * pa_x[i] + 2.0 * ts_xyzz_zzz[i] * fz_0;
    }

    // Set up 90-100 components of targeted buffer : GF

    auto tk_xzzz_xxx = pbuffer.data(idx_kin_gf + 90);

    auto tk_xzzz_xxy = pbuffer.data(idx_kin_gf + 91);

    auto tk_xzzz_xxz = pbuffer.data(idx_kin_gf + 92);

    auto tk_xzzz_xyy = pbuffer.data(idx_kin_gf + 93);

    auto tk_xzzz_xyz = pbuffer.data(idx_kin_gf + 94);

    auto tk_xzzz_xzz = pbuffer.data(idx_kin_gf + 95);

    auto tk_xzzz_yyy = pbuffer.data(idx_kin_gf + 96);

    auto tk_xzzz_yyz = pbuffer.data(idx_kin_gf + 97);

    auto tk_xzzz_yzz = pbuffer.data(idx_kin_gf + 98);

    auto tk_xzzz_zzz = pbuffer.data(idx_kin_gf + 99);

#pragma omp simd aligned(pa_x,            \
                             tk_xzzz_xxx, \
                             tk_xzzz_xxy, \
                             tk_xzzz_xxz, \
                             tk_xzzz_xyy, \
                             tk_xzzz_xyz, \
                             tk_xzzz_xzz, \
                             tk_xzzz_yyy, \
                             tk_xzzz_yyz, \
                             tk_xzzz_yzz, \
                             tk_xzzz_zzz, \
                             tk_zzz_xx,   \
                             tk_zzz_xxx,  \
                             tk_zzz_xxy,  \
                             tk_zzz_xxz,  \
                             tk_zzz_xy,   \
                             tk_zzz_xyy,  \
                             tk_zzz_xyz,  \
                             tk_zzz_xz,   \
                             tk_zzz_xzz,  \
                             tk_zzz_yy,   \
                             tk_zzz_yyy,  \
                             tk_zzz_yyz,  \
                             tk_zzz_yz,   \
                             tk_zzz_yzz,  \
                             tk_zzz_zz,   \
                             tk_zzz_zzz,  \
                             ts_xzzz_xxx, \
                             ts_xzzz_xxy, \
                             ts_xzzz_xxz, \
                             ts_xzzz_xyy, \
                             ts_xzzz_xyz, \
                             ts_xzzz_xzz, \
                             ts_xzzz_yyy, \
                             ts_xzzz_yyz, \
                             ts_xzzz_yzz, \
                             ts_xzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzz_xxx[i] = 3.0 * tk_zzz_xx[i] * fe_0 + tk_zzz_xxx[i] * pa_x[i] + 2.0 * ts_xzzz_xxx[i] * fz_0;

        tk_xzzz_xxy[i] = 2.0 * tk_zzz_xy[i] * fe_0 + tk_zzz_xxy[i] * pa_x[i] + 2.0 * ts_xzzz_xxy[i] * fz_0;

        tk_xzzz_xxz[i] = 2.0 * tk_zzz_xz[i] * fe_0 + tk_zzz_xxz[i] * pa_x[i] + 2.0 * ts_xzzz_xxz[i] * fz_0;

        tk_xzzz_xyy[i] = tk_zzz_yy[i] * fe_0 + tk_zzz_xyy[i] * pa_x[i] + 2.0 * ts_xzzz_xyy[i] * fz_0;

        tk_xzzz_xyz[i] = tk_zzz_yz[i] * fe_0 + tk_zzz_xyz[i] * pa_x[i] + 2.0 * ts_xzzz_xyz[i] * fz_0;

        tk_xzzz_xzz[i] = tk_zzz_zz[i] * fe_0 + tk_zzz_xzz[i] * pa_x[i] + 2.0 * ts_xzzz_xzz[i] * fz_0;

        tk_xzzz_yyy[i] = tk_zzz_yyy[i] * pa_x[i] + 2.0 * ts_xzzz_yyy[i] * fz_0;

        tk_xzzz_yyz[i] = tk_zzz_yyz[i] * pa_x[i] + 2.0 * ts_xzzz_yyz[i] * fz_0;

        tk_xzzz_yzz[i] = tk_zzz_yzz[i] * pa_x[i] + 2.0 * ts_xzzz_yzz[i] * fz_0;

        tk_xzzz_zzz[i] = tk_zzz_zzz[i] * pa_x[i] + 2.0 * ts_xzzz_zzz[i] * fz_0;
    }

    // Set up 100-110 components of targeted buffer : GF

    auto tk_yyyy_xxx = pbuffer.data(idx_kin_gf + 100);

    auto tk_yyyy_xxy = pbuffer.data(idx_kin_gf + 101);

    auto tk_yyyy_xxz = pbuffer.data(idx_kin_gf + 102);

    auto tk_yyyy_xyy = pbuffer.data(idx_kin_gf + 103);

    auto tk_yyyy_xyz = pbuffer.data(idx_kin_gf + 104);

    auto tk_yyyy_xzz = pbuffer.data(idx_kin_gf + 105);

    auto tk_yyyy_yyy = pbuffer.data(idx_kin_gf + 106);

    auto tk_yyyy_yyz = pbuffer.data(idx_kin_gf + 107);

    auto tk_yyyy_yzz = pbuffer.data(idx_kin_gf + 108);

    auto tk_yyyy_zzz = pbuffer.data(idx_kin_gf + 109);

#pragma omp simd aligned(pa_y,            \
                             tk_yy_xxx,   \
                             tk_yy_xxy,   \
                             tk_yy_xxz,   \
                             tk_yy_xyy,   \
                             tk_yy_xyz,   \
                             tk_yy_xzz,   \
                             tk_yy_yyy,   \
                             tk_yy_yyz,   \
                             tk_yy_yzz,   \
                             tk_yy_zzz,   \
                             tk_yyy_xx,   \
                             tk_yyy_xxx,  \
                             tk_yyy_xxy,  \
                             tk_yyy_xxz,  \
                             tk_yyy_xy,   \
                             tk_yyy_xyy,  \
                             tk_yyy_xyz,  \
                             tk_yyy_xz,   \
                             tk_yyy_xzz,  \
                             tk_yyy_yy,   \
                             tk_yyy_yyy,  \
                             tk_yyy_yyz,  \
                             tk_yyy_yz,   \
                             tk_yyy_yzz,  \
                             tk_yyy_zz,   \
                             tk_yyy_zzz,  \
                             tk_yyyy_xxx, \
                             tk_yyyy_xxy, \
                             tk_yyyy_xxz, \
                             tk_yyyy_xyy, \
                             tk_yyyy_xyz, \
                             tk_yyyy_xzz, \
                             tk_yyyy_yyy, \
                             tk_yyyy_yyz, \
                             tk_yyyy_yzz, \
                             tk_yyyy_zzz, \
                             ts_yy_xxx,   \
                             ts_yy_xxy,   \
                             ts_yy_xxz,   \
                             ts_yy_xyy,   \
                             ts_yy_xyz,   \
                             ts_yy_xzz,   \
                             ts_yy_yyy,   \
                             ts_yy_yyz,   \
                             ts_yy_yzz,   \
                             ts_yy_zzz,   \
                             ts_yyyy_xxx, \
                             ts_yyyy_xxy, \
                             ts_yyyy_xxz, \
                             ts_yyyy_xyy, \
                             ts_yyyy_xyz, \
                             ts_yyyy_xzz, \
                             ts_yyyy_yyy, \
                             ts_yyyy_yyz, \
                             ts_yyyy_yzz, \
                             ts_yyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyy_xxx[i] = -6.0 * ts_yy_xxx[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxx[i] * fe_0 + tk_yyy_xxx[i] * pa_y[i] + 2.0 * ts_yyyy_xxx[i] * fz_0;

        tk_yyyy_xxy[i] = -6.0 * ts_yy_xxy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxy[i] * fe_0 + tk_yyy_xx[i] * fe_0 + tk_yyy_xxy[i] * pa_y[i] +
                         2.0 * ts_yyyy_xxy[i] * fz_0;

        tk_yyyy_xxz[i] = -6.0 * ts_yy_xxz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxz[i] * fe_0 + tk_yyy_xxz[i] * pa_y[i] + 2.0 * ts_yyyy_xxz[i] * fz_0;

        tk_yyyy_xyy[i] = -6.0 * ts_yy_xyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyy[i] * fe_0 + 2.0 * tk_yyy_xy[i] * fe_0 + tk_yyy_xyy[i] * pa_y[i] +
                         2.0 * ts_yyyy_xyy[i] * fz_0;

        tk_yyyy_xyz[i] = -6.0 * ts_yy_xyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyz[i] * fe_0 + tk_yyy_xz[i] * fe_0 + tk_yyy_xyz[i] * pa_y[i] +
                         2.0 * ts_yyyy_xyz[i] * fz_0;

        tk_yyyy_xzz[i] = -6.0 * ts_yy_xzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xzz[i] * fe_0 + tk_yyy_xzz[i] * pa_y[i] + 2.0 * ts_yyyy_xzz[i] * fz_0;

        tk_yyyy_yyy[i] = -6.0 * ts_yy_yyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyy[i] * fe_0 + 3.0 * tk_yyy_yy[i] * fe_0 + tk_yyy_yyy[i] * pa_y[i] +
                         2.0 * ts_yyyy_yyy[i] * fz_0;

        tk_yyyy_yyz[i] = -6.0 * ts_yy_yyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyz[i] * fe_0 + 2.0 * tk_yyy_yz[i] * fe_0 + tk_yyy_yyz[i] * pa_y[i] +
                         2.0 * ts_yyyy_yyz[i] * fz_0;

        tk_yyyy_yzz[i] = -6.0 * ts_yy_yzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yzz[i] * fe_0 + tk_yyy_zz[i] * fe_0 + tk_yyy_yzz[i] * pa_y[i] +
                         2.0 * ts_yyyy_yzz[i] * fz_0;

        tk_yyyy_zzz[i] = -6.0 * ts_yy_zzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_zzz[i] * fe_0 + tk_yyy_zzz[i] * pa_y[i] + 2.0 * ts_yyyy_zzz[i] * fz_0;
    }

    // Set up 110-120 components of targeted buffer : GF

    auto tk_yyyz_xxx = pbuffer.data(idx_kin_gf + 110);

    auto tk_yyyz_xxy = pbuffer.data(idx_kin_gf + 111);

    auto tk_yyyz_xxz = pbuffer.data(idx_kin_gf + 112);

    auto tk_yyyz_xyy = pbuffer.data(idx_kin_gf + 113);

    auto tk_yyyz_xyz = pbuffer.data(idx_kin_gf + 114);

    auto tk_yyyz_xzz = pbuffer.data(idx_kin_gf + 115);

    auto tk_yyyz_yyy = pbuffer.data(idx_kin_gf + 116);

    auto tk_yyyz_yyz = pbuffer.data(idx_kin_gf + 117);

    auto tk_yyyz_yzz = pbuffer.data(idx_kin_gf + 118);

    auto tk_yyyz_zzz = pbuffer.data(idx_kin_gf + 119);

#pragma omp simd aligned(pa_z,            \
                             tk_yyy_xx,   \
                             tk_yyy_xxx,  \
                             tk_yyy_xxy,  \
                             tk_yyy_xxz,  \
                             tk_yyy_xy,   \
                             tk_yyy_xyy,  \
                             tk_yyy_xyz,  \
                             tk_yyy_xz,   \
                             tk_yyy_xzz,  \
                             tk_yyy_yy,   \
                             tk_yyy_yyy,  \
                             tk_yyy_yyz,  \
                             tk_yyy_yz,   \
                             tk_yyy_yzz,  \
                             tk_yyy_zz,   \
                             tk_yyy_zzz,  \
                             tk_yyyz_xxx, \
                             tk_yyyz_xxy, \
                             tk_yyyz_xxz, \
                             tk_yyyz_xyy, \
                             tk_yyyz_xyz, \
                             tk_yyyz_xzz, \
                             tk_yyyz_yyy, \
                             tk_yyyz_yyz, \
                             tk_yyyz_yzz, \
                             tk_yyyz_zzz, \
                             ts_yyyz_xxx, \
                             ts_yyyz_xxy, \
                             ts_yyyz_xxz, \
                             ts_yyyz_xyy, \
                             ts_yyyz_xyz, \
                             ts_yyyz_xzz, \
                             ts_yyyz_yyy, \
                             ts_yyyz_yyz, \
                             ts_yyyz_yzz, \
                             ts_yyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyz_xxx[i] = tk_yyy_xxx[i] * pa_z[i] + 2.0 * ts_yyyz_xxx[i] * fz_0;

        tk_yyyz_xxy[i] = tk_yyy_xxy[i] * pa_z[i] + 2.0 * ts_yyyz_xxy[i] * fz_0;

        tk_yyyz_xxz[i] = tk_yyy_xx[i] * fe_0 + tk_yyy_xxz[i] * pa_z[i] + 2.0 * ts_yyyz_xxz[i] * fz_0;

        tk_yyyz_xyy[i] = tk_yyy_xyy[i] * pa_z[i] + 2.0 * ts_yyyz_xyy[i] * fz_0;

        tk_yyyz_xyz[i] = tk_yyy_xy[i] * fe_0 + tk_yyy_xyz[i] * pa_z[i] + 2.0 * ts_yyyz_xyz[i] * fz_0;

        tk_yyyz_xzz[i] = 2.0 * tk_yyy_xz[i] * fe_0 + tk_yyy_xzz[i] * pa_z[i] + 2.0 * ts_yyyz_xzz[i] * fz_0;

        tk_yyyz_yyy[i] = tk_yyy_yyy[i] * pa_z[i] + 2.0 * ts_yyyz_yyy[i] * fz_0;

        tk_yyyz_yyz[i] = tk_yyy_yy[i] * fe_0 + tk_yyy_yyz[i] * pa_z[i] + 2.0 * ts_yyyz_yyz[i] * fz_0;

        tk_yyyz_yzz[i] = 2.0 * tk_yyy_yz[i] * fe_0 + tk_yyy_yzz[i] * pa_z[i] + 2.0 * ts_yyyz_yzz[i] * fz_0;

        tk_yyyz_zzz[i] = 3.0 * tk_yyy_zz[i] * fe_0 + tk_yyy_zzz[i] * pa_z[i] + 2.0 * ts_yyyz_zzz[i] * fz_0;
    }

    // Set up 120-130 components of targeted buffer : GF

    auto tk_yyzz_xxx = pbuffer.data(idx_kin_gf + 120);

    auto tk_yyzz_xxy = pbuffer.data(idx_kin_gf + 121);

    auto tk_yyzz_xxz = pbuffer.data(idx_kin_gf + 122);

    auto tk_yyzz_xyy = pbuffer.data(idx_kin_gf + 123);

    auto tk_yyzz_xyz = pbuffer.data(idx_kin_gf + 124);

    auto tk_yyzz_xzz = pbuffer.data(idx_kin_gf + 125);

    auto tk_yyzz_yyy = pbuffer.data(idx_kin_gf + 126);

    auto tk_yyzz_yyz = pbuffer.data(idx_kin_gf + 127);

    auto tk_yyzz_yzz = pbuffer.data(idx_kin_gf + 128);

    auto tk_yyzz_zzz = pbuffer.data(idx_kin_gf + 129);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tk_yy_xxy,   \
                             tk_yy_xyy,   \
                             tk_yy_yyy,   \
                             tk_yyz_xxy,  \
                             tk_yyz_xyy,  \
                             tk_yyz_yyy,  \
                             tk_yyzz_xxx, \
                             tk_yyzz_xxy, \
                             tk_yyzz_xxz, \
                             tk_yyzz_xyy, \
                             tk_yyzz_xyz, \
                             tk_yyzz_xzz, \
                             tk_yyzz_yyy, \
                             tk_yyzz_yyz, \
                             tk_yyzz_yzz, \
                             tk_yyzz_zzz, \
                             tk_yzz_xxx,  \
                             tk_yzz_xxz,  \
                             tk_yzz_xyz,  \
                             tk_yzz_xz,   \
                             tk_yzz_xzz,  \
                             tk_yzz_yyz,  \
                             tk_yzz_yz,   \
                             tk_yzz_yzz,  \
                             tk_yzz_zz,   \
                             tk_yzz_zzz,  \
                             tk_zz_xxx,   \
                             tk_zz_xxz,   \
                             tk_zz_xyz,   \
                             tk_zz_xzz,   \
                             tk_zz_yyz,   \
                             tk_zz_yzz,   \
                             tk_zz_zzz,   \
                             ts_yy_xxy,   \
                             ts_yy_xyy,   \
                             ts_yy_yyy,   \
                             ts_yyzz_xxx, \
                             ts_yyzz_xxy, \
                             ts_yyzz_xxz, \
                             ts_yyzz_xyy, \
                             ts_yyzz_xyz, \
                             ts_yyzz_xzz, \
                             ts_yyzz_yyy, \
                             ts_yyzz_yyz, \
                             ts_yyzz_yzz, \
                             ts_yyzz_zzz, \
                             ts_zz_xxx,   \
                             ts_zz_xxz,   \
                             ts_zz_xyz,   \
                             ts_zz_xzz,   \
                             ts_zz_yyz,   \
                             ts_zz_yzz,   \
                             ts_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzz_xxx[i] = -2.0 * ts_zz_xxx[i] * fbe_0 * fz_0 + tk_zz_xxx[i] * fe_0 + tk_yzz_xxx[i] * pa_y[i] + 2.0 * ts_yyzz_xxx[i] * fz_0;

        tk_yyzz_xxy[i] = -2.0 * ts_yy_xxy[i] * fbe_0 * fz_0 + tk_yy_xxy[i] * fe_0 + tk_yyz_xxy[i] * pa_z[i] + 2.0 * ts_yyzz_xxy[i] * fz_0;

        tk_yyzz_xxz[i] = -2.0 * ts_zz_xxz[i] * fbe_0 * fz_0 + tk_zz_xxz[i] * fe_0 + tk_yzz_xxz[i] * pa_y[i] + 2.0 * ts_yyzz_xxz[i] * fz_0;

        tk_yyzz_xyy[i] = -2.0 * ts_yy_xyy[i] * fbe_0 * fz_0 + tk_yy_xyy[i] * fe_0 + tk_yyz_xyy[i] * pa_z[i] + 2.0 * ts_yyzz_xyy[i] * fz_0;

        tk_yyzz_xyz[i] =
            -2.0 * ts_zz_xyz[i] * fbe_0 * fz_0 + tk_zz_xyz[i] * fe_0 + tk_yzz_xz[i] * fe_0 + tk_yzz_xyz[i] * pa_y[i] + 2.0 * ts_yyzz_xyz[i] * fz_0;

        tk_yyzz_xzz[i] = -2.0 * ts_zz_xzz[i] * fbe_0 * fz_0 + tk_zz_xzz[i] * fe_0 + tk_yzz_xzz[i] * pa_y[i] + 2.0 * ts_yyzz_xzz[i] * fz_0;

        tk_yyzz_yyy[i] = -2.0 * ts_yy_yyy[i] * fbe_0 * fz_0 + tk_yy_yyy[i] * fe_0 + tk_yyz_yyy[i] * pa_z[i] + 2.0 * ts_yyzz_yyy[i] * fz_0;

        tk_yyzz_yyz[i] = -2.0 * ts_zz_yyz[i] * fbe_0 * fz_0 + tk_zz_yyz[i] * fe_0 + 2.0 * tk_yzz_yz[i] * fe_0 + tk_yzz_yyz[i] * pa_y[i] +
                         2.0 * ts_yyzz_yyz[i] * fz_0;

        tk_yyzz_yzz[i] =
            -2.0 * ts_zz_yzz[i] * fbe_0 * fz_0 + tk_zz_yzz[i] * fe_0 + tk_yzz_zz[i] * fe_0 + tk_yzz_yzz[i] * pa_y[i] + 2.0 * ts_yyzz_yzz[i] * fz_0;

        tk_yyzz_zzz[i] = -2.0 * ts_zz_zzz[i] * fbe_0 * fz_0 + tk_zz_zzz[i] * fe_0 + tk_yzz_zzz[i] * pa_y[i] + 2.0 * ts_yyzz_zzz[i] * fz_0;
    }

    // Set up 130-140 components of targeted buffer : GF

    auto tk_yzzz_xxx = pbuffer.data(idx_kin_gf + 130);

    auto tk_yzzz_xxy = pbuffer.data(idx_kin_gf + 131);

    auto tk_yzzz_xxz = pbuffer.data(idx_kin_gf + 132);

    auto tk_yzzz_xyy = pbuffer.data(idx_kin_gf + 133);

    auto tk_yzzz_xyz = pbuffer.data(idx_kin_gf + 134);

    auto tk_yzzz_xzz = pbuffer.data(idx_kin_gf + 135);

    auto tk_yzzz_yyy = pbuffer.data(idx_kin_gf + 136);

    auto tk_yzzz_yyz = pbuffer.data(idx_kin_gf + 137);

    auto tk_yzzz_yzz = pbuffer.data(idx_kin_gf + 138);

    auto tk_yzzz_zzz = pbuffer.data(idx_kin_gf + 139);

#pragma omp simd aligned(pa_y,            \
                             tk_yzzz_xxx, \
                             tk_yzzz_xxy, \
                             tk_yzzz_xxz, \
                             tk_yzzz_xyy, \
                             tk_yzzz_xyz, \
                             tk_yzzz_xzz, \
                             tk_yzzz_yyy, \
                             tk_yzzz_yyz, \
                             tk_yzzz_yzz, \
                             tk_yzzz_zzz, \
                             tk_zzz_xx,   \
                             tk_zzz_xxx,  \
                             tk_zzz_xxy,  \
                             tk_zzz_xxz,  \
                             tk_zzz_xy,   \
                             tk_zzz_xyy,  \
                             tk_zzz_xyz,  \
                             tk_zzz_xz,   \
                             tk_zzz_xzz,  \
                             tk_zzz_yy,   \
                             tk_zzz_yyy,  \
                             tk_zzz_yyz,  \
                             tk_zzz_yz,   \
                             tk_zzz_yzz,  \
                             tk_zzz_zz,   \
                             tk_zzz_zzz,  \
                             ts_yzzz_xxx, \
                             ts_yzzz_xxy, \
                             ts_yzzz_xxz, \
                             ts_yzzz_xyy, \
                             ts_yzzz_xyz, \
                             ts_yzzz_xzz, \
                             ts_yzzz_yyy, \
                             ts_yzzz_yyz, \
                             ts_yzzz_yzz, \
                             ts_yzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzz_xxx[i] = tk_zzz_xxx[i] * pa_y[i] + 2.0 * ts_yzzz_xxx[i] * fz_0;

        tk_yzzz_xxy[i] = tk_zzz_xx[i] * fe_0 + tk_zzz_xxy[i] * pa_y[i] + 2.0 * ts_yzzz_xxy[i] * fz_0;

        tk_yzzz_xxz[i] = tk_zzz_xxz[i] * pa_y[i] + 2.0 * ts_yzzz_xxz[i] * fz_0;

        tk_yzzz_xyy[i] = 2.0 * tk_zzz_xy[i] * fe_0 + tk_zzz_xyy[i] * pa_y[i] + 2.0 * ts_yzzz_xyy[i] * fz_0;

        tk_yzzz_xyz[i] = tk_zzz_xz[i] * fe_0 + tk_zzz_xyz[i] * pa_y[i] + 2.0 * ts_yzzz_xyz[i] * fz_0;

        tk_yzzz_xzz[i] = tk_zzz_xzz[i] * pa_y[i] + 2.0 * ts_yzzz_xzz[i] * fz_0;

        tk_yzzz_yyy[i] = 3.0 * tk_zzz_yy[i] * fe_0 + tk_zzz_yyy[i] * pa_y[i] + 2.0 * ts_yzzz_yyy[i] * fz_0;

        tk_yzzz_yyz[i] = 2.0 * tk_zzz_yz[i] * fe_0 + tk_zzz_yyz[i] * pa_y[i] + 2.0 * ts_yzzz_yyz[i] * fz_0;

        tk_yzzz_yzz[i] = tk_zzz_zz[i] * fe_0 + tk_zzz_yzz[i] * pa_y[i] + 2.0 * ts_yzzz_yzz[i] * fz_0;

        tk_yzzz_zzz[i] = tk_zzz_zzz[i] * pa_y[i] + 2.0 * ts_yzzz_zzz[i] * fz_0;
    }

    // Set up 140-150 components of targeted buffer : GF

    auto tk_zzzz_xxx = pbuffer.data(idx_kin_gf + 140);

    auto tk_zzzz_xxy = pbuffer.data(idx_kin_gf + 141);

    auto tk_zzzz_xxz = pbuffer.data(idx_kin_gf + 142);

    auto tk_zzzz_xyy = pbuffer.data(idx_kin_gf + 143);

    auto tk_zzzz_xyz = pbuffer.data(idx_kin_gf + 144);

    auto tk_zzzz_xzz = pbuffer.data(idx_kin_gf + 145);

    auto tk_zzzz_yyy = pbuffer.data(idx_kin_gf + 146);

    auto tk_zzzz_yyz = pbuffer.data(idx_kin_gf + 147);

    auto tk_zzzz_yzz = pbuffer.data(idx_kin_gf + 148);

    auto tk_zzzz_zzz = pbuffer.data(idx_kin_gf + 149);

#pragma omp simd aligned(pa_z,            \
                             tk_zz_xxx,   \
                             tk_zz_xxy,   \
                             tk_zz_xxz,   \
                             tk_zz_xyy,   \
                             tk_zz_xyz,   \
                             tk_zz_xzz,   \
                             tk_zz_yyy,   \
                             tk_zz_yyz,   \
                             tk_zz_yzz,   \
                             tk_zz_zzz,   \
                             tk_zzz_xx,   \
                             tk_zzz_xxx,  \
                             tk_zzz_xxy,  \
                             tk_zzz_xxz,  \
                             tk_zzz_xy,   \
                             tk_zzz_xyy,  \
                             tk_zzz_xyz,  \
                             tk_zzz_xz,   \
                             tk_zzz_xzz,  \
                             tk_zzz_yy,   \
                             tk_zzz_yyy,  \
                             tk_zzz_yyz,  \
                             tk_zzz_yz,   \
                             tk_zzz_yzz,  \
                             tk_zzz_zz,   \
                             tk_zzz_zzz,  \
                             tk_zzzz_xxx, \
                             tk_zzzz_xxy, \
                             tk_zzzz_xxz, \
                             tk_zzzz_xyy, \
                             tk_zzzz_xyz, \
                             tk_zzzz_xzz, \
                             tk_zzzz_yyy, \
                             tk_zzzz_yyz, \
                             tk_zzzz_yzz, \
                             tk_zzzz_zzz, \
                             ts_zz_xxx,   \
                             ts_zz_xxy,   \
                             ts_zz_xxz,   \
                             ts_zz_xyy,   \
                             ts_zz_xyz,   \
                             ts_zz_xzz,   \
                             ts_zz_yyy,   \
                             ts_zz_yyz,   \
                             ts_zz_yzz,   \
                             ts_zz_zzz,   \
                             ts_zzzz_xxx, \
                             ts_zzzz_xxy, \
                             ts_zzzz_xxz, \
                             ts_zzzz_xyy, \
                             ts_zzzz_xyz, \
                             ts_zzzz_xzz, \
                             ts_zzzz_yyy, \
                             ts_zzzz_yyz, \
                             ts_zzzz_yzz, \
                             ts_zzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzz_xxx[i] = -6.0 * ts_zz_xxx[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxx[i] * fe_0 + tk_zzz_xxx[i] * pa_z[i] + 2.0 * ts_zzzz_xxx[i] * fz_0;

        tk_zzzz_xxy[i] = -6.0 * ts_zz_xxy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxy[i] * fe_0 + tk_zzz_xxy[i] * pa_z[i] + 2.0 * ts_zzzz_xxy[i] * fz_0;

        tk_zzzz_xxz[i] = -6.0 * ts_zz_xxz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxz[i] * fe_0 + tk_zzz_xx[i] * fe_0 + tk_zzz_xxz[i] * pa_z[i] +
                         2.0 * ts_zzzz_xxz[i] * fz_0;

        tk_zzzz_xyy[i] = -6.0 * ts_zz_xyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyy[i] * fe_0 + tk_zzz_xyy[i] * pa_z[i] + 2.0 * ts_zzzz_xyy[i] * fz_0;

        tk_zzzz_xyz[i] = -6.0 * ts_zz_xyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyz[i] * fe_0 + tk_zzz_xy[i] * fe_0 + tk_zzz_xyz[i] * pa_z[i] +
                         2.0 * ts_zzzz_xyz[i] * fz_0;

        tk_zzzz_xzz[i] = -6.0 * ts_zz_xzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xzz[i] * fe_0 + 2.0 * tk_zzz_xz[i] * fe_0 + tk_zzz_xzz[i] * pa_z[i] +
                         2.0 * ts_zzzz_xzz[i] * fz_0;

        tk_zzzz_yyy[i] = -6.0 * ts_zz_yyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyy[i] * fe_0 + tk_zzz_yyy[i] * pa_z[i] + 2.0 * ts_zzzz_yyy[i] * fz_0;

        tk_zzzz_yyz[i] = -6.0 * ts_zz_yyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyz[i] * fe_0 + tk_zzz_yy[i] * fe_0 + tk_zzz_yyz[i] * pa_z[i] +
                         2.0 * ts_zzzz_yyz[i] * fz_0;

        tk_zzzz_yzz[i] = -6.0 * ts_zz_yzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yzz[i] * fe_0 + 2.0 * tk_zzz_yz[i] * fe_0 + tk_zzz_yzz[i] * pa_z[i] +
                         2.0 * ts_zzzz_yzz[i] * fz_0;

        tk_zzzz_zzz[i] = -6.0 * ts_zz_zzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_zzz[i] * fe_0 + 3.0 * tk_zzz_zz[i] * fe_0 + tk_zzz_zzz[i] * pa_z[i] +
                         2.0 * ts_zzzz_zzz[i] * fz_0;
    }
}

}  // namespace kinrec
