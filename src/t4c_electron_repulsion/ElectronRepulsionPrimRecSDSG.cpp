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

#include "ElectronRepulsionPrimRecSDSG.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sdsg(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sdsg,
                                  size_t                idx_eri_0_sssg,
                                  size_t                idx_eri_1_sssg,
                                  size_t                idx_eri_1_spsf,
                                  size_t                idx_eri_0_spsg,
                                  size_t                idx_eri_1_spsg,
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

    /// Set up components of auxilary buffer : SSSG

    auto g_0_0_0_xxxx_0 = pbuffer.data(idx_eri_0_sssg);

    auto g_0_0_0_xxxy_0 = pbuffer.data(idx_eri_0_sssg + 1);

    auto g_0_0_0_xxxz_0 = pbuffer.data(idx_eri_0_sssg + 2);

    auto g_0_0_0_xxyy_0 = pbuffer.data(idx_eri_0_sssg + 3);

    auto g_0_0_0_xxyz_0 = pbuffer.data(idx_eri_0_sssg + 4);

    auto g_0_0_0_xxzz_0 = pbuffer.data(idx_eri_0_sssg + 5);

    auto g_0_0_0_xyyy_0 = pbuffer.data(idx_eri_0_sssg + 6);

    auto g_0_0_0_xyyz_0 = pbuffer.data(idx_eri_0_sssg + 7);

    auto g_0_0_0_xyzz_0 = pbuffer.data(idx_eri_0_sssg + 8);

    auto g_0_0_0_xzzz_0 = pbuffer.data(idx_eri_0_sssg + 9);

    auto g_0_0_0_yyyy_0 = pbuffer.data(idx_eri_0_sssg + 10);

    auto g_0_0_0_yyyz_0 = pbuffer.data(idx_eri_0_sssg + 11);

    auto g_0_0_0_yyzz_0 = pbuffer.data(idx_eri_0_sssg + 12);

    auto g_0_0_0_yzzz_0 = pbuffer.data(idx_eri_0_sssg + 13);

    auto g_0_0_0_zzzz_0 = pbuffer.data(idx_eri_0_sssg + 14);

    /// Set up components of auxilary buffer : SSSG

    auto g_0_0_0_xxxx_1 = pbuffer.data(idx_eri_1_sssg);

    auto g_0_0_0_xxxy_1 = pbuffer.data(idx_eri_1_sssg + 1);

    auto g_0_0_0_xxxz_1 = pbuffer.data(idx_eri_1_sssg + 2);

    auto g_0_0_0_xxyy_1 = pbuffer.data(idx_eri_1_sssg + 3);

    auto g_0_0_0_xxyz_1 = pbuffer.data(idx_eri_1_sssg + 4);

    auto g_0_0_0_xxzz_1 = pbuffer.data(idx_eri_1_sssg + 5);

    auto g_0_0_0_xyyy_1 = pbuffer.data(idx_eri_1_sssg + 6);

    auto g_0_0_0_xyyz_1 = pbuffer.data(idx_eri_1_sssg + 7);

    auto g_0_0_0_xyzz_1 = pbuffer.data(idx_eri_1_sssg + 8);

    auto g_0_0_0_xzzz_1 = pbuffer.data(idx_eri_1_sssg + 9);

    auto g_0_0_0_yyyy_1 = pbuffer.data(idx_eri_1_sssg + 10);

    auto g_0_0_0_yyyz_1 = pbuffer.data(idx_eri_1_sssg + 11);

    auto g_0_0_0_yyzz_1 = pbuffer.data(idx_eri_1_sssg + 12);

    auto g_0_0_0_yzzz_1 = pbuffer.data(idx_eri_1_sssg + 13);

    auto g_0_0_0_zzzz_1 = pbuffer.data(idx_eri_1_sssg + 14);

    /// Set up components of auxilary buffer : SPSF

    auto g_0_x_0_xxx_1 = pbuffer.data(idx_eri_1_spsf);

    auto g_0_x_0_xxy_1 = pbuffer.data(idx_eri_1_spsf + 1);

    auto g_0_x_0_xxz_1 = pbuffer.data(idx_eri_1_spsf + 2);

    auto g_0_x_0_xyy_1 = pbuffer.data(idx_eri_1_spsf + 3);

    auto g_0_x_0_xyz_1 = pbuffer.data(idx_eri_1_spsf + 4);

    auto g_0_x_0_xzz_1 = pbuffer.data(idx_eri_1_spsf + 5);

    auto g_0_x_0_yyy_1 = pbuffer.data(idx_eri_1_spsf + 6);

    auto g_0_x_0_yyz_1 = pbuffer.data(idx_eri_1_spsf + 7);

    auto g_0_x_0_yzz_1 = pbuffer.data(idx_eri_1_spsf + 8);

    auto g_0_x_0_zzz_1 = pbuffer.data(idx_eri_1_spsf + 9);

    auto g_0_y_0_xxx_1 = pbuffer.data(idx_eri_1_spsf + 10);

    auto g_0_y_0_xxy_1 = pbuffer.data(idx_eri_1_spsf + 11);

    auto g_0_y_0_xxz_1 = pbuffer.data(idx_eri_1_spsf + 12);

    auto g_0_y_0_xyy_1 = pbuffer.data(idx_eri_1_spsf + 13);

    auto g_0_y_0_xyz_1 = pbuffer.data(idx_eri_1_spsf + 14);

    auto g_0_y_0_xzz_1 = pbuffer.data(idx_eri_1_spsf + 15);

    auto g_0_y_0_yyy_1 = pbuffer.data(idx_eri_1_spsf + 16);

    auto g_0_y_0_yyz_1 = pbuffer.data(idx_eri_1_spsf + 17);

    auto g_0_y_0_yzz_1 = pbuffer.data(idx_eri_1_spsf + 18);

    auto g_0_y_0_zzz_1 = pbuffer.data(idx_eri_1_spsf + 19);

    auto g_0_z_0_xxx_1 = pbuffer.data(idx_eri_1_spsf + 20);

    auto g_0_z_0_xxy_1 = pbuffer.data(idx_eri_1_spsf + 21);

    auto g_0_z_0_xxz_1 = pbuffer.data(idx_eri_1_spsf + 22);

    auto g_0_z_0_xyy_1 = pbuffer.data(idx_eri_1_spsf + 23);

    auto g_0_z_0_xyz_1 = pbuffer.data(idx_eri_1_spsf + 24);

    auto g_0_z_0_xzz_1 = pbuffer.data(idx_eri_1_spsf + 25);

    auto g_0_z_0_yyy_1 = pbuffer.data(idx_eri_1_spsf + 26);

    auto g_0_z_0_yyz_1 = pbuffer.data(idx_eri_1_spsf + 27);

    auto g_0_z_0_yzz_1 = pbuffer.data(idx_eri_1_spsf + 28);

    auto g_0_z_0_zzz_1 = pbuffer.data(idx_eri_1_spsf + 29);

    /// Set up components of auxilary buffer : SPSG

    auto g_0_x_0_xxxx_0 = pbuffer.data(idx_eri_0_spsg);

    auto g_0_x_0_xxxy_0 = pbuffer.data(idx_eri_0_spsg + 1);

    auto g_0_x_0_xxxz_0 = pbuffer.data(idx_eri_0_spsg + 2);

    auto g_0_x_0_xxyy_0 = pbuffer.data(idx_eri_0_spsg + 3);

    auto g_0_x_0_xxyz_0 = pbuffer.data(idx_eri_0_spsg + 4);

    auto g_0_x_0_xxzz_0 = pbuffer.data(idx_eri_0_spsg + 5);

    auto g_0_x_0_xyyy_0 = pbuffer.data(idx_eri_0_spsg + 6);

    auto g_0_x_0_xyyz_0 = pbuffer.data(idx_eri_0_spsg + 7);

    auto g_0_x_0_xyzz_0 = pbuffer.data(idx_eri_0_spsg + 8);

    auto g_0_x_0_xzzz_0 = pbuffer.data(idx_eri_0_spsg + 9);

    auto g_0_x_0_yyyy_0 = pbuffer.data(idx_eri_0_spsg + 10);

    auto g_0_x_0_yyyz_0 = pbuffer.data(idx_eri_0_spsg + 11);

    auto g_0_x_0_yyzz_0 = pbuffer.data(idx_eri_0_spsg + 12);

    auto g_0_x_0_yzzz_0 = pbuffer.data(idx_eri_0_spsg + 13);

    auto g_0_x_0_zzzz_0 = pbuffer.data(idx_eri_0_spsg + 14);

    auto g_0_y_0_xxxx_0 = pbuffer.data(idx_eri_0_spsg + 15);

    auto g_0_y_0_xxxy_0 = pbuffer.data(idx_eri_0_spsg + 16);

    auto g_0_y_0_xxxz_0 = pbuffer.data(idx_eri_0_spsg + 17);

    auto g_0_y_0_xxyy_0 = pbuffer.data(idx_eri_0_spsg + 18);

    auto g_0_y_0_xxyz_0 = pbuffer.data(idx_eri_0_spsg + 19);

    auto g_0_y_0_xxzz_0 = pbuffer.data(idx_eri_0_spsg + 20);

    auto g_0_y_0_xyyy_0 = pbuffer.data(idx_eri_0_spsg + 21);

    auto g_0_y_0_xyyz_0 = pbuffer.data(idx_eri_0_spsg + 22);

    auto g_0_y_0_xyzz_0 = pbuffer.data(idx_eri_0_spsg + 23);

    auto g_0_y_0_xzzz_0 = pbuffer.data(idx_eri_0_spsg + 24);

    auto g_0_y_0_yyyy_0 = pbuffer.data(idx_eri_0_spsg + 25);

    auto g_0_y_0_yyyz_0 = pbuffer.data(idx_eri_0_spsg + 26);

    auto g_0_y_0_yyzz_0 = pbuffer.data(idx_eri_0_spsg + 27);

    auto g_0_y_0_yzzz_0 = pbuffer.data(idx_eri_0_spsg + 28);

    auto g_0_y_0_zzzz_0 = pbuffer.data(idx_eri_0_spsg + 29);

    auto g_0_z_0_xxxx_0 = pbuffer.data(idx_eri_0_spsg + 30);

    auto g_0_z_0_xxxy_0 = pbuffer.data(idx_eri_0_spsg + 31);

    auto g_0_z_0_xxxz_0 = pbuffer.data(idx_eri_0_spsg + 32);

    auto g_0_z_0_xxyy_0 = pbuffer.data(idx_eri_0_spsg + 33);

    auto g_0_z_0_xxyz_0 = pbuffer.data(idx_eri_0_spsg + 34);

    auto g_0_z_0_xxzz_0 = pbuffer.data(idx_eri_0_spsg + 35);

    auto g_0_z_0_xyyy_0 = pbuffer.data(idx_eri_0_spsg + 36);

    auto g_0_z_0_xyyz_0 = pbuffer.data(idx_eri_0_spsg + 37);

    auto g_0_z_0_xyzz_0 = pbuffer.data(idx_eri_0_spsg + 38);

    auto g_0_z_0_xzzz_0 = pbuffer.data(idx_eri_0_spsg + 39);

    auto g_0_z_0_yyyy_0 = pbuffer.data(idx_eri_0_spsg + 40);

    auto g_0_z_0_yyyz_0 = pbuffer.data(idx_eri_0_spsg + 41);

    auto g_0_z_0_yyzz_0 = pbuffer.data(idx_eri_0_spsg + 42);

    auto g_0_z_0_yzzz_0 = pbuffer.data(idx_eri_0_spsg + 43);

    auto g_0_z_0_zzzz_0 = pbuffer.data(idx_eri_0_spsg + 44);

    /// Set up components of auxilary buffer : SPSG

    auto g_0_x_0_xxxx_1 = pbuffer.data(idx_eri_1_spsg);

    auto g_0_x_0_xxxy_1 = pbuffer.data(idx_eri_1_spsg + 1);

    auto g_0_x_0_xxxz_1 = pbuffer.data(idx_eri_1_spsg + 2);

    auto g_0_x_0_xxyy_1 = pbuffer.data(idx_eri_1_spsg + 3);

    auto g_0_x_0_xxyz_1 = pbuffer.data(idx_eri_1_spsg + 4);

    auto g_0_x_0_xxzz_1 = pbuffer.data(idx_eri_1_spsg + 5);

    auto g_0_x_0_xyyy_1 = pbuffer.data(idx_eri_1_spsg + 6);

    auto g_0_x_0_xyyz_1 = pbuffer.data(idx_eri_1_spsg + 7);

    auto g_0_x_0_xyzz_1 = pbuffer.data(idx_eri_1_spsg + 8);

    auto g_0_x_0_xzzz_1 = pbuffer.data(idx_eri_1_spsg + 9);

    auto g_0_x_0_yyyy_1 = pbuffer.data(idx_eri_1_spsg + 10);

    auto g_0_x_0_yyyz_1 = pbuffer.data(idx_eri_1_spsg + 11);

    auto g_0_x_0_yyzz_1 = pbuffer.data(idx_eri_1_spsg + 12);

    auto g_0_x_0_yzzz_1 = pbuffer.data(idx_eri_1_spsg + 13);

    auto g_0_x_0_zzzz_1 = pbuffer.data(idx_eri_1_spsg + 14);

    auto g_0_y_0_xxxx_1 = pbuffer.data(idx_eri_1_spsg + 15);

    auto g_0_y_0_xxxy_1 = pbuffer.data(idx_eri_1_spsg + 16);

    auto g_0_y_0_xxxz_1 = pbuffer.data(idx_eri_1_spsg + 17);

    auto g_0_y_0_xxyy_1 = pbuffer.data(idx_eri_1_spsg + 18);

    auto g_0_y_0_xxyz_1 = pbuffer.data(idx_eri_1_spsg + 19);

    auto g_0_y_0_xxzz_1 = pbuffer.data(idx_eri_1_spsg + 20);

    auto g_0_y_0_xyyy_1 = pbuffer.data(idx_eri_1_spsg + 21);

    auto g_0_y_0_xyyz_1 = pbuffer.data(idx_eri_1_spsg + 22);

    auto g_0_y_0_xyzz_1 = pbuffer.data(idx_eri_1_spsg + 23);

    auto g_0_y_0_xzzz_1 = pbuffer.data(idx_eri_1_spsg + 24);

    auto g_0_y_0_yyyy_1 = pbuffer.data(idx_eri_1_spsg + 25);

    auto g_0_y_0_yyyz_1 = pbuffer.data(idx_eri_1_spsg + 26);

    auto g_0_y_0_yyzz_1 = pbuffer.data(idx_eri_1_spsg + 27);

    auto g_0_y_0_yzzz_1 = pbuffer.data(idx_eri_1_spsg + 28);

    auto g_0_y_0_zzzz_1 = pbuffer.data(idx_eri_1_spsg + 29);

    auto g_0_z_0_xxxx_1 = pbuffer.data(idx_eri_1_spsg + 30);

    auto g_0_z_0_xxxy_1 = pbuffer.data(idx_eri_1_spsg + 31);

    auto g_0_z_0_xxxz_1 = pbuffer.data(idx_eri_1_spsg + 32);

    auto g_0_z_0_xxyy_1 = pbuffer.data(idx_eri_1_spsg + 33);

    auto g_0_z_0_xxyz_1 = pbuffer.data(idx_eri_1_spsg + 34);

    auto g_0_z_0_xxzz_1 = pbuffer.data(idx_eri_1_spsg + 35);

    auto g_0_z_0_xyyy_1 = pbuffer.data(idx_eri_1_spsg + 36);

    auto g_0_z_0_xyyz_1 = pbuffer.data(idx_eri_1_spsg + 37);

    auto g_0_z_0_xyzz_1 = pbuffer.data(idx_eri_1_spsg + 38);

    auto g_0_z_0_xzzz_1 = pbuffer.data(idx_eri_1_spsg + 39);

    auto g_0_z_0_yyyy_1 = pbuffer.data(idx_eri_1_spsg + 40);

    auto g_0_z_0_yyyz_1 = pbuffer.data(idx_eri_1_spsg + 41);

    auto g_0_z_0_yyzz_1 = pbuffer.data(idx_eri_1_spsg + 42);

    auto g_0_z_0_yzzz_1 = pbuffer.data(idx_eri_1_spsg + 43);

    auto g_0_z_0_zzzz_1 = pbuffer.data(idx_eri_1_spsg + 44);

    /// Set up 0-15 components of targeted buffer : SDSG

    auto g_0_xx_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg);

    auto g_0_xx_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 1);

    auto g_0_xx_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 2);

    auto g_0_xx_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 3);

    auto g_0_xx_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 4);

    auto g_0_xx_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 5);

    auto g_0_xx_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 6);

    auto g_0_xx_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 7);

    auto g_0_xx_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 8);

    auto g_0_xx_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 9);

    auto g_0_xx_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 10);

    auto g_0_xx_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 11);

    auto g_0_xx_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 12);

    auto g_0_xx_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 13);

    auto g_0_xx_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 14);

#pragma omp simd aligned(g_0_0_0_xxxx_0,      \
                             g_0_0_0_xxxx_1,  \
                             g_0_0_0_xxxy_0,  \
                             g_0_0_0_xxxy_1,  \
                             g_0_0_0_xxxz_0,  \
                             g_0_0_0_xxxz_1,  \
                             g_0_0_0_xxyy_0,  \
                             g_0_0_0_xxyy_1,  \
                             g_0_0_0_xxyz_0,  \
                             g_0_0_0_xxyz_1,  \
                             g_0_0_0_xxzz_0,  \
                             g_0_0_0_xxzz_1,  \
                             g_0_0_0_xyyy_0,  \
                             g_0_0_0_xyyy_1,  \
                             g_0_0_0_xyyz_0,  \
                             g_0_0_0_xyyz_1,  \
                             g_0_0_0_xyzz_0,  \
                             g_0_0_0_xyzz_1,  \
                             g_0_0_0_xzzz_0,  \
                             g_0_0_0_xzzz_1,  \
                             g_0_0_0_yyyy_0,  \
                             g_0_0_0_yyyy_1,  \
                             g_0_0_0_yyyz_0,  \
                             g_0_0_0_yyyz_1,  \
                             g_0_0_0_yyzz_0,  \
                             g_0_0_0_yyzz_1,  \
                             g_0_0_0_yzzz_0,  \
                             g_0_0_0_yzzz_1,  \
                             g_0_0_0_zzzz_0,  \
                             g_0_0_0_zzzz_1,  \
                             g_0_x_0_xxx_1,   \
                             g_0_x_0_xxxx_0,  \
                             g_0_x_0_xxxx_1,  \
                             g_0_x_0_xxxy_0,  \
                             g_0_x_0_xxxy_1,  \
                             g_0_x_0_xxxz_0,  \
                             g_0_x_0_xxxz_1,  \
                             g_0_x_0_xxy_1,   \
                             g_0_x_0_xxyy_0,  \
                             g_0_x_0_xxyy_1,  \
                             g_0_x_0_xxyz_0,  \
                             g_0_x_0_xxyz_1,  \
                             g_0_x_0_xxz_1,   \
                             g_0_x_0_xxzz_0,  \
                             g_0_x_0_xxzz_1,  \
                             g_0_x_0_xyy_1,   \
                             g_0_x_0_xyyy_0,  \
                             g_0_x_0_xyyy_1,  \
                             g_0_x_0_xyyz_0,  \
                             g_0_x_0_xyyz_1,  \
                             g_0_x_0_xyz_1,   \
                             g_0_x_0_xyzz_0,  \
                             g_0_x_0_xyzz_1,  \
                             g_0_x_0_xzz_1,   \
                             g_0_x_0_xzzz_0,  \
                             g_0_x_0_xzzz_1,  \
                             g_0_x_0_yyy_1,   \
                             g_0_x_0_yyyy_0,  \
                             g_0_x_0_yyyy_1,  \
                             g_0_x_0_yyyz_0,  \
                             g_0_x_0_yyyz_1,  \
                             g_0_x_0_yyz_1,   \
                             g_0_x_0_yyzz_0,  \
                             g_0_x_0_yyzz_1,  \
                             g_0_x_0_yzz_1,   \
                             g_0_x_0_yzzz_0,  \
                             g_0_x_0_yzzz_1,  \
                             g_0_x_0_zzz_1,   \
                             g_0_x_0_zzzz_0,  \
                             g_0_x_0_zzzz_1,  \
                             g_0_xx_0_xxxx_0, \
                             g_0_xx_0_xxxy_0, \
                             g_0_xx_0_xxxz_0, \
                             g_0_xx_0_xxyy_0, \
                             g_0_xx_0_xxyz_0, \
                             g_0_xx_0_xxzz_0, \
                             g_0_xx_0_xyyy_0, \
                             g_0_xx_0_xyyz_0, \
                             g_0_xx_0_xyzz_0, \
                             g_0_xx_0_xzzz_0, \
                             g_0_xx_0_yyyy_0, \
                             g_0_xx_0_yyyz_0, \
                             g_0_xx_0_yyzz_0, \
                             g_0_xx_0_yzzz_0, \
                             g_0_xx_0_zzzz_0, \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xxxx_0[i] = g_0_0_0_xxxx_0[i] * fi_ab_0 - g_0_0_0_xxxx_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxx_1[i] * fi_abcd_0 +
                             g_0_x_0_xxxx_0[i] * pb_x + g_0_x_0_xxxx_1[i] * wp_x[i];

        g_0_xx_0_xxxy_0[i] = g_0_0_0_xxxy_0[i] * fi_ab_0 - g_0_0_0_xxxy_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxy_1[i] * fi_abcd_0 +
                             g_0_x_0_xxxy_0[i] * pb_x + g_0_x_0_xxxy_1[i] * wp_x[i];

        g_0_xx_0_xxxz_0[i] = g_0_0_0_xxxz_0[i] * fi_ab_0 - g_0_0_0_xxxz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxz_1[i] * fi_abcd_0 +
                             g_0_x_0_xxxz_0[i] * pb_x + g_0_x_0_xxxz_1[i] * wp_x[i];

        g_0_xx_0_xxyy_0[i] = g_0_0_0_xxyy_0[i] * fi_ab_0 - g_0_0_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyy_1[i] * fi_abcd_0 +
                             g_0_x_0_xxyy_0[i] * pb_x + g_0_x_0_xxyy_1[i] * wp_x[i];

        g_0_xx_0_xxyz_0[i] = g_0_0_0_xxyz_0[i] * fi_ab_0 - g_0_0_0_xxyz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyz_1[i] * fi_abcd_0 +
                             g_0_x_0_xxyz_0[i] * pb_x + g_0_x_0_xxyz_1[i] * wp_x[i];

        g_0_xx_0_xxzz_0[i] = g_0_0_0_xxzz_0[i] * fi_ab_0 - g_0_0_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xzz_1[i] * fi_abcd_0 +
                             g_0_x_0_xxzz_0[i] * pb_x + g_0_x_0_xxzz_1[i] * wp_x[i];

        g_0_xx_0_xyyy_0[i] = g_0_0_0_xyyy_0[i] * fi_ab_0 - g_0_0_0_xyyy_1[i] * fti_ab_0 + g_0_x_0_yyy_1[i] * fi_abcd_0 + g_0_x_0_xyyy_0[i] * pb_x +
                             g_0_x_0_xyyy_1[i] * wp_x[i];

        g_0_xx_0_xyyz_0[i] = g_0_0_0_xyyz_0[i] * fi_ab_0 - g_0_0_0_xyyz_1[i] * fti_ab_0 + g_0_x_0_yyz_1[i] * fi_abcd_0 + g_0_x_0_xyyz_0[i] * pb_x +
                             g_0_x_0_xyyz_1[i] * wp_x[i];

        g_0_xx_0_xyzz_0[i] = g_0_0_0_xyzz_0[i] * fi_ab_0 - g_0_0_0_xyzz_1[i] * fti_ab_0 + g_0_x_0_yzz_1[i] * fi_abcd_0 + g_0_x_0_xyzz_0[i] * pb_x +
                             g_0_x_0_xyzz_1[i] * wp_x[i];

        g_0_xx_0_xzzz_0[i] = g_0_0_0_xzzz_0[i] * fi_ab_0 - g_0_0_0_xzzz_1[i] * fti_ab_0 + g_0_x_0_zzz_1[i] * fi_abcd_0 + g_0_x_0_xzzz_0[i] * pb_x +
                             g_0_x_0_xzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyy_0[i] = g_0_0_0_yyyy_0[i] * fi_ab_0 - g_0_0_0_yyyy_1[i] * fti_ab_0 + g_0_x_0_yyyy_0[i] * pb_x + g_0_x_0_yyyy_1[i] * wp_x[i];

        g_0_xx_0_yyyz_0[i] = g_0_0_0_yyyz_0[i] * fi_ab_0 - g_0_0_0_yyyz_1[i] * fti_ab_0 + g_0_x_0_yyyz_0[i] * pb_x + g_0_x_0_yyyz_1[i] * wp_x[i];

        g_0_xx_0_yyzz_0[i] = g_0_0_0_yyzz_0[i] * fi_ab_0 - g_0_0_0_yyzz_1[i] * fti_ab_0 + g_0_x_0_yyzz_0[i] * pb_x + g_0_x_0_yyzz_1[i] * wp_x[i];

        g_0_xx_0_yzzz_0[i] = g_0_0_0_yzzz_0[i] * fi_ab_0 - g_0_0_0_yzzz_1[i] * fti_ab_0 + g_0_x_0_yzzz_0[i] * pb_x + g_0_x_0_yzzz_1[i] * wp_x[i];

        g_0_xx_0_zzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_ab_0 - g_0_0_0_zzzz_1[i] * fti_ab_0 + g_0_x_0_zzzz_0[i] * pb_x + g_0_x_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 15-30 components of targeted buffer : SDSG

    auto g_0_xy_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 15);

    auto g_0_xy_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 16);

    auto g_0_xy_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 17);

    auto g_0_xy_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 18);

    auto g_0_xy_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 19);

    auto g_0_xy_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 20);

    auto g_0_xy_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 21);

    auto g_0_xy_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 22);

    auto g_0_xy_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 23);

    auto g_0_xy_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 24);

    auto g_0_xy_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 25);

    auto g_0_xy_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 26);

    auto g_0_xy_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 27);

    auto g_0_xy_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 28);

    auto g_0_xy_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 29);

#pragma omp simd aligned(g_0_x_0_xxxx_0,      \
                             g_0_x_0_xxxx_1,  \
                             g_0_x_0_xxxz_0,  \
                             g_0_x_0_xxxz_1,  \
                             g_0_x_0_xxzz_0,  \
                             g_0_x_0_xxzz_1,  \
                             g_0_x_0_xzzz_0,  \
                             g_0_x_0_xzzz_1,  \
                             g_0_xy_0_xxxx_0, \
                             g_0_xy_0_xxxy_0, \
                             g_0_xy_0_xxxz_0, \
                             g_0_xy_0_xxyy_0, \
                             g_0_xy_0_xxyz_0, \
                             g_0_xy_0_xxzz_0, \
                             g_0_xy_0_xyyy_0, \
                             g_0_xy_0_xyyz_0, \
                             g_0_xy_0_xyzz_0, \
                             g_0_xy_0_xzzz_0, \
                             g_0_xy_0_yyyy_0, \
                             g_0_xy_0_yyyz_0, \
                             g_0_xy_0_yyzz_0, \
                             g_0_xy_0_yzzz_0, \
                             g_0_xy_0_zzzz_0, \
                             g_0_y_0_xxxy_0,  \
                             g_0_y_0_xxxy_1,  \
                             g_0_y_0_xxy_1,   \
                             g_0_y_0_xxyy_0,  \
                             g_0_y_0_xxyy_1,  \
                             g_0_y_0_xxyz_0,  \
                             g_0_y_0_xxyz_1,  \
                             g_0_y_0_xyy_1,   \
                             g_0_y_0_xyyy_0,  \
                             g_0_y_0_xyyy_1,  \
                             g_0_y_0_xyyz_0,  \
                             g_0_y_0_xyyz_1,  \
                             g_0_y_0_xyz_1,   \
                             g_0_y_0_xyzz_0,  \
                             g_0_y_0_xyzz_1,  \
                             g_0_y_0_yyy_1,   \
                             g_0_y_0_yyyy_0,  \
                             g_0_y_0_yyyy_1,  \
                             g_0_y_0_yyyz_0,  \
                             g_0_y_0_yyyz_1,  \
                             g_0_y_0_yyz_1,   \
                             g_0_y_0_yyzz_0,  \
                             g_0_y_0_yyzz_1,  \
                             g_0_y_0_yzz_1,   \
                             g_0_y_0_yzzz_0,  \
                             g_0_y_0_yzzz_1,  \
                             g_0_y_0_zzzz_0,  \
                             g_0_y_0_zzzz_1,  \
                             wp_x,            \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xxxx_0[i] = g_0_x_0_xxxx_0[i] * pb_y + g_0_x_0_xxxx_1[i] * wp_y[i];

        g_0_xy_0_xxxy_0[i] = 3.0 * g_0_y_0_xxy_1[i] * fi_abcd_0 + g_0_y_0_xxxy_0[i] * pb_x + g_0_y_0_xxxy_1[i] * wp_x[i];

        g_0_xy_0_xxxz_0[i] = g_0_x_0_xxxz_0[i] * pb_y + g_0_x_0_xxxz_1[i] * wp_y[i];

        g_0_xy_0_xxyy_0[i] = 2.0 * g_0_y_0_xyy_1[i] * fi_abcd_0 + g_0_y_0_xxyy_0[i] * pb_x + g_0_y_0_xxyy_1[i] * wp_x[i];

        g_0_xy_0_xxyz_0[i] = 2.0 * g_0_y_0_xyz_1[i] * fi_abcd_0 + g_0_y_0_xxyz_0[i] * pb_x + g_0_y_0_xxyz_1[i] * wp_x[i];

        g_0_xy_0_xxzz_0[i] = g_0_x_0_xxzz_0[i] * pb_y + g_0_x_0_xxzz_1[i] * wp_y[i];

        g_0_xy_0_xyyy_0[i] = g_0_y_0_yyy_1[i] * fi_abcd_0 + g_0_y_0_xyyy_0[i] * pb_x + g_0_y_0_xyyy_1[i] * wp_x[i];

        g_0_xy_0_xyyz_0[i] = g_0_y_0_yyz_1[i] * fi_abcd_0 + g_0_y_0_xyyz_0[i] * pb_x + g_0_y_0_xyyz_1[i] * wp_x[i];

        g_0_xy_0_xyzz_0[i] = g_0_y_0_yzz_1[i] * fi_abcd_0 + g_0_y_0_xyzz_0[i] * pb_x + g_0_y_0_xyzz_1[i] * wp_x[i];

        g_0_xy_0_xzzz_0[i] = g_0_x_0_xzzz_0[i] * pb_y + g_0_x_0_xzzz_1[i] * wp_y[i];

        g_0_xy_0_yyyy_0[i] = g_0_y_0_yyyy_0[i] * pb_x + g_0_y_0_yyyy_1[i] * wp_x[i];

        g_0_xy_0_yyyz_0[i] = g_0_y_0_yyyz_0[i] * pb_x + g_0_y_0_yyyz_1[i] * wp_x[i];

        g_0_xy_0_yyzz_0[i] = g_0_y_0_yyzz_0[i] * pb_x + g_0_y_0_yyzz_1[i] * wp_x[i];

        g_0_xy_0_yzzz_0[i] = g_0_y_0_yzzz_0[i] * pb_x + g_0_y_0_yzzz_1[i] * wp_x[i];

        g_0_xy_0_zzzz_0[i] = g_0_y_0_zzzz_0[i] * pb_x + g_0_y_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 30-45 components of targeted buffer : SDSG

    auto g_0_xz_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 30);

    auto g_0_xz_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 31);

    auto g_0_xz_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 32);

    auto g_0_xz_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 33);

    auto g_0_xz_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 34);

    auto g_0_xz_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 35);

    auto g_0_xz_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 36);

    auto g_0_xz_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 37);

    auto g_0_xz_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 38);

    auto g_0_xz_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 39);

    auto g_0_xz_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 40);

    auto g_0_xz_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 41);

    auto g_0_xz_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 42);

    auto g_0_xz_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 43);

    auto g_0_xz_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 44);

#pragma omp simd aligned(g_0_x_0_xxxx_0,      \
                             g_0_x_0_xxxx_1,  \
                             g_0_x_0_xxxy_0,  \
                             g_0_x_0_xxxy_1,  \
                             g_0_x_0_xxyy_0,  \
                             g_0_x_0_xxyy_1,  \
                             g_0_x_0_xyyy_0,  \
                             g_0_x_0_xyyy_1,  \
                             g_0_xz_0_xxxx_0, \
                             g_0_xz_0_xxxy_0, \
                             g_0_xz_0_xxxz_0, \
                             g_0_xz_0_xxyy_0, \
                             g_0_xz_0_xxyz_0, \
                             g_0_xz_0_xxzz_0, \
                             g_0_xz_0_xyyy_0, \
                             g_0_xz_0_xyyz_0, \
                             g_0_xz_0_xyzz_0, \
                             g_0_xz_0_xzzz_0, \
                             g_0_xz_0_yyyy_0, \
                             g_0_xz_0_yyyz_0, \
                             g_0_xz_0_yyzz_0, \
                             g_0_xz_0_yzzz_0, \
                             g_0_xz_0_zzzz_0, \
                             g_0_z_0_xxxz_0,  \
                             g_0_z_0_xxxz_1,  \
                             g_0_z_0_xxyz_0,  \
                             g_0_z_0_xxyz_1,  \
                             g_0_z_0_xxz_1,   \
                             g_0_z_0_xxzz_0,  \
                             g_0_z_0_xxzz_1,  \
                             g_0_z_0_xyyz_0,  \
                             g_0_z_0_xyyz_1,  \
                             g_0_z_0_xyz_1,   \
                             g_0_z_0_xyzz_0,  \
                             g_0_z_0_xyzz_1,  \
                             g_0_z_0_xzz_1,   \
                             g_0_z_0_xzzz_0,  \
                             g_0_z_0_xzzz_1,  \
                             g_0_z_0_yyyy_0,  \
                             g_0_z_0_yyyy_1,  \
                             g_0_z_0_yyyz_0,  \
                             g_0_z_0_yyyz_1,  \
                             g_0_z_0_yyz_1,   \
                             g_0_z_0_yyzz_0,  \
                             g_0_z_0_yyzz_1,  \
                             g_0_z_0_yzz_1,   \
                             g_0_z_0_yzzz_0,  \
                             g_0_z_0_yzzz_1,  \
                             g_0_z_0_zzz_1,   \
                             g_0_z_0_zzzz_0,  \
                             g_0_z_0_zzzz_1,  \
                             wp_x,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xxxx_0[i] = g_0_x_0_xxxx_0[i] * pb_z + g_0_x_0_xxxx_1[i] * wp_z[i];

        g_0_xz_0_xxxy_0[i] = g_0_x_0_xxxy_0[i] * pb_z + g_0_x_0_xxxy_1[i] * wp_z[i];

        g_0_xz_0_xxxz_0[i] = 3.0 * g_0_z_0_xxz_1[i] * fi_abcd_0 + g_0_z_0_xxxz_0[i] * pb_x + g_0_z_0_xxxz_1[i] * wp_x[i];

        g_0_xz_0_xxyy_0[i] = g_0_x_0_xxyy_0[i] * pb_z + g_0_x_0_xxyy_1[i] * wp_z[i];

        g_0_xz_0_xxyz_0[i] = 2.0 * g_0_z_0_xyz_1[i] * fi_abcd_0 + g_0_z_0_xxyz_0[i] * pb_x + g_0_z_0_xxyz_1[i] * wp_x[i];

        g_0_xz_0_xxzz_0[i] = 2.0 * g_0_z_0_xzz_1[i] * fi_abcd_0 + g_0_z_0_xxzz_0[i] * pb_x + g_0_z_0_xxzz_1[i] * wp_x[i];

        g_0_xz_0_xyyy_0[i] = g_0_x_0_xyyy_0[i] * pb_z + g_0_x_0_xyyy_1[i] * wp_z[i];

        g_0_xz_0_xyyz_0[i] = g_0_z_0_yyz_1[i] * fi_abcd_0 + g_0_z_0_xyyz_0[i] * pb_x + g_0_z_0_xyyz_1[i] * wp_x[i];

        g_0_xz_0_xyzz_0[i] = g_0_z_0_yzz_1[i] * fi_abcd_0 + g_0_z_0_xyzz_0[i] * pb_x + g_0_z_0_xyzz_1[i] * wp_x[i];

        g_0_xz_0_xzzz_0[i] = g_0_z_0_zzz_1[i] * fi_abcd_0 + g_0_z_0_xzzz_0[i] * pb_x + g_0_z_0_xzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyy_0[i] = g_0_z_0_yyyy_0[i] * pb_x + g_0_z_0_yyyy_1[i] * wp_x[i];

        g_0_xz_0_yyyz_0[i] = g_0_z_0_yyyz_0[i] * pb_x + g_0_z_0_yyyz_1[i] * wp_x[i];

        g_0_xz_0_yyzz_0[i] = g_0_z_0_yyzz_0[i] * pb_x + g_0_z_0_yyzz_1[i] * wp_x[i];

        g_0_xz_0_yzzz_0[i] = g_0_z_0_yzzz_0[i] * pb_x + g_0_z_0_yzzz_1[i] * wp_x[i];

        g_0_xz_0_zzzz_0[i] = g_0_z_0_zzzz_0[i] * pb_x + g_0_z_0_zzzz_1[i] * wp_x[i];
    }

    /// Set up 45-60 components of targeted buffer : SDSG

    auto g_0_yy_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 45);

    auto g_0_yy_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 46);

    auto g_0_yy_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 47);

    auto g_0_yy_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 48);

    auto g_0_yy_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 49);

    auto g_0_yy_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 50);

    auto g_0_yy_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 51);

    auto g_0_yy_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 52);

    auto g_0_yy_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 53);

    auto g_0_yy_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 54);

    auto g_0_yy_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 55);

    auto g_0_yy_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 56);

    auto g_0_yy_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 57);

    auto g_0_yy_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 58);

    auto g_0_yy_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 59);

#pragma omp simd aligned(g_0_0_0_xxxx_0,      \
                             g_0_0_0_xxxx_1,  \
                             g_0_0_0_xxxy_0,  \
                             g_0_0_0_xxxy_1,  \
                             g_0_0_0_xxxz_0,  \
                             g_0_0_0_xxxz_1,  \
                             g_0_0_0_xxyy_0,  \
                             g_0_0_0_xxyy_1,  \
                             g_0_0_0_xxyz_0,  \
                             g_0_0_0_xxyz_1,  \
                             g_0_0_0_xxzz_0,  \
                             g_0_0_0_xxzz_1,  \
                             g_0_0_0_xyyy_0,  \
                             g_0_0_0_xyyy_1,  \
                             g_0_0_0_xyyz_0,  \
                             g_0_0_0_xyyz_1,  \
                             g_0_0_0_xyzz_0,  \
                             g_0_0_0_xyzz_1,  \
                             g_0_0_0_xzzz_0,  \
                             g_0_0_0_xzzz_1,  \
                             g_0_0_0_yyyy_0,  \
                             g_0_0_0_yyyy_1,  \
                             g_0_0_0_yyyz_0,  \
                             g_0_0_0_yyyz_1,  \
                             g_0_0_0_yyzz_0,  \
                             g_0_0_0_yyzz_1,  \
                             g_0_0_0_yzzz_0,  \
                             g_0_0_0_yzzz_1,  \
                             g_0_0_0_zzzz_0,  \
                             g_0_0_0_zzzz_1,  \
                             g_0_y_0_xxx_1,   \
                             g_0_y_0_xxxx_0,  \
                             g_0_y_0_xxxx_1,  \
                             g_0_y_0_xxxy_0,  \
                             g_0_y_0_xxxy_1,  \
                             g_0_y_0_xxxz_0,  \
                             g_0_y_0_xxxz_1,  \
                             g_0_y_0_xxy_1,   \
                             g_0_y_0_xxyy_0,  \
                             g_0_y_0_xxyy_1,  \
                             g_0_y_0_xxyz_0,  \
                             g_0_y_0_xxyz_1,  \
                             g_0_y_0_xxz_1,   \
                             g_0_y_0_xxzz_0,  \
                             g_0_y_0_xxzz_1,  \
                             g_0_y_0_xyy_1,   \
                             g_0_y_0_xyyy_0,  \
                             g_0_y_0_xyyy_1,  \
                             g_0_y_0_xyyz_0,  \
                             g_0_y_0_xyyz_1,  \
                             g_0_y_0_xyz_1,   \
                             g_0_y_0_xyzz_0,  \
                             g_0_y_0_xyzz_1,  \
                             g_0_y_0_xzz_1,   \
                             g_0_y_0_xzzz_0,  \
                             g_0_y_0_xzzz_1,  \
                             g_0_y_0_yyy_1,   \
                             g_0_y_0_yyyy_0,  \
                             g_0_y_0_yyyy_1,  \
                             g_0_y_0_yyyz_0,  \
                             g_0_y_0_yyyz_1,  \
                             g_0_y_0_yyz_1,   \
                             g_0_y_0_yyzz_0,  \
                             g_0_y_0_yyzz_1,  \
                             g_0_y_0_yzz_1,   \
                             g_0_y_0_yzzz_0,  \
                             g_0_y_0_yzzz_1,  \
                             g_0_y_0_zzz_1,   \
                             g_0_y_0_zzzz_0,  \
                             g_0_y_0_zzzz_1,  \
                             g_0_yy_0_xxxx_0, \
                             g_0_yy_0_xxxy_0, \
                             g_0_yy_0_xxxz_0, \
                             g_0_yy_0_xxyy_0, \
                             g_0_yy_0_xxyz_0, \
                             g_0_yy_0_xxzz_0, \
                             g_0_yy_0_xyyy_0, \
                             g_0_yy_0_xyyz_0, \
                             g_0_yy_0_xyzz_0, \
                             g_0_yy_0_xzzz_0, \
                             g_0_yy_0_yyyy_0, \
                             g_0_yy_0_yyyz_0, \
                             g_0_yy_0_yyzz_0, \
                             g_0_yy_0_yzzz_0, \
                             g_0_yy_0_zzzz_0, \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xxxx_0[i] = g_0_0_0_xxxx_0[i] * fi_ab_0 - g_0_0_0_xxxx_1[i] * fti_ab_0 + g_0_y_0_xxxx_0[i] * pb_y + g_0_y_0_xxxx_1[i] * wp_y[i];

        g_0_yy_0_xxxy_0[i] = g_0_0_0_xxxy_0[i] * fi_ab_0 - g_0_0_0_xxxy_1[i] * fti_ab_0 + g_0_y_0_xxx_1[i] * fi_abcd_0 + g_0_y_0_xxxy_0[i] * pb_y +
                             g_0_y_0_xxxy_1[i] * wp_y[i];

        g_0_yy_0_xxxz_0[i] = g_0_0_0_xxxz_0[i] * fi_ab_0 - g_0_0_0_xxxz_1[i] * fti_ab_0 + g_0_y_0_xxxz_0[i] * pb_y + g_0_y_0_xxxz_1[i] * wp_y[i];

        g_0_yy_0_xxyy_0[i] = g_0_0_0_xxyy_0[i] * fi_ab_0 - g_0_0_0_xxyy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxy_1[i] * fi_abcd_0 +
                             g_0_y_0_xxyy_0[i] * pb_y + g_0_y_0_xxyy_1[i] * wp_y[i];

        g_0_yy_0_xxyz_0[i] = g_0_0_0_xxyz_0[i] * fi_ab_0 - g_0_0_0_xxyz_1[i] * fti_ab_0 + g_0_y_0_xxz_1[i] * fi_abcd_0 + g_0_y_0_xxyz_0[i] * pb_y +
                             g_0_y_0_xxyz_1[i] * wp_y[i];

        g_0_yy_0_xxzz_0[i] = g_0_0_0_xxzz_0[i] * fi_ab_0 - g_0_0_0_xxzz_1[i] * fti_ab_0 + g_0_y_0_xxzz_0[i] * pb_y + g_0_y_0_xxzz_1[i] * wp_y[i];

        g_0_yy_0_xyyy_0[i] = g_0_0_0_xyyy_0[i] * fi_ab_0 - g_0_0_0_xyyy_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xyy_1[i] * fi_abcd_0 +
                             g_0_y_0_xyyy_0[i] * pb_y + g_0_y_0_xyyy_1[i] * wp_y[i];

        g_0_yy_0_xyyz_0[i] = g_0_0_0_xyyz_0[i] * fi_ab_0 - g_0_0_0_xyyz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xyz_1[i] * fi_abcd_0 +
                             g_0_y_0_xyyz_0[i] * pb_y + g_0_y_0_xyyz_1[i] * wp_y[i];

        g_0_yy_0_xyzz_0[i] = g_0_0_0_xyzz_0[i] * fi_ab_0 - g_0_0_0_xyzz_1[i] * fti_ab_0 + g_0_y_0_xzz_1[i] * fi_abcd_0 + g_0_y_0_xyzz_0[i] * pb_y +
                             g_0_y_0_xyzz_1[i] * wp_y[i];

        g_0_yy_0_xzzz_0[i] = g_0_0_0_xzzz_0[i] * fi_ab_0 - g_0_0_0_xzzz_1[i] * fti_ab_0 + g_0_y_0_xzzz_0[i] * pb_y + g_0_y_0_xzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyy_0[i] = g_0_0_0_yyyy_0[i] * fi_ab_0 - g_0_0_0_yyyy_1[i] * fti_ab_0 + 4.0 * g_0_y_0_yyy_1[i] * fi_abcd_0 +
                             g_0_y_0_yyyy_0[i] * pb_y + g_0_y_0_yyyy_1[i] * wp_y[i];

        g_0_yy_0_yyyz_0[i] = g_0_0_0_yyyz_0[i] * fi_ab_0 - g_0_0_0_yyyz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_yyz_1[i] * fi_abcd_0 +
                             g_0_y_0_yyyz_0[i] * pb_y + g_0_y_0_yyyz_1[i] * wp_y[i];

        g_0_yy_0_yyzz_0[i] = g_0_0_0_yyzz_0[i] * fi_ab_0 - g_0_0_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_yzz_1[i] * fi_abcd_0 +
                             g_0_y_0_yyzz_0[i] * pb_y + g_0_y_0_yyzz_1[i] * wp_y[i];

        g_0_yy_0_yzzz_0[i] = g_0_0_0_yzzz_0[i] * fi_ab_0 - g_0_0_0_yzzz_1[i] * fti_ab_0 + g_0_y_0_zzz_1[i] * fi_abcd_0 + g_0_y_0_yzzz_0[i] * pb_y +
                             g_0_y_0_yzzz_1[i] * wp_y[i];

        g_0_yy_0_zzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_ab_0 - g_0_0_0_zzzz_1[i] * fti_ab_0 + g_0_y_0_zzzz_0[i] * pb_y + g_0_y_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 60-75 components of targeted buffer : SDSG

    auto g_0_yz_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 60);

    auto g_0_yz_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 61);

    auto g_0_yz_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 62);

    auto g_0_yz_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 63);

    auto g_0_yz_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 64);

    auto g_0_yz_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 65);

    auto g_0_yz_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 66);

    auto g_0_yz_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 67);

    auto g_0_yz_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 68);

    auto g_0_yz_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 69);

    auto g_0_yz_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 70);

    auto g_0_yz_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 71);

    auto g_0_yz_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 72);

    auto g_0_yz_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 73);

    auto g_0_yz_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 74);

#pragma omp simd aligned(g_0_y_0_xxxy_0,      \
                             g_0_y_0_xxxy_1,  \
                             g_0_y_0_xxyy_0,  \
                             g_0_y_0_xxyy_1,  \
                             g_0_y_0_xyyy_0,  \
                             g_0_y_0_xyyy_1,  \
                             g_0_y_0_yyyy_0,  \
                             g_0_y_0_yyyy_1,  \
                             g_0_yz_0_xxxx_0, \
                             g_0_yz_0_xxxy_0, \
                             g_0_yz_0_xxxz_0, \
                             g_0_yz_0_xxyy_0, \
                             g_0_yz_0_xxyz_0, \
                             g_0_yz_0_xxzz_0, \
                             g_0_yz_0_xyyy_0, \
                             g_0_yz_0_xyyz_0, \
                             g_0_yz_0_xyzz_0, \
                             g_0_yz_0_xzzz_0, \
                             g_0_yz_0_yyyy_0, \
                             g_0_yz_0_yyyz_0, \
                             g_0_yz_0_yyzz_0, \
                             g_0_yz_0_yzzz_0, \
                             g_0_yz_0_zzzz_0, \
                             g_0_z_0_xxxx_0,  \
                             g_0_z_0_xxxx_1,  \
                             g_0_z_0_xxxz_0,  \
                             g_0_z_0_xxxz_1,  \
                             g_0_z_0_xxyz_0,  \
                             g_0_z_0_xxyz_1,  \
                             g_0_z_0_xxz_1,   \
                             g_0_z_0_xxzz_0,  \
                             g_0_z_0_xxzz_1,  \
                             g_0_z_0_xyyz_0,  \
                             g_0_z_0_xyyz_1,  \
                             g_0_z_0_xyz_1,   \
                             g_0_z_0_xyzz_0,  \
                             g_0_z_0_xyzz_1,  \
                             g_0_z_0_xzz_1,   \
                             g_0_z_0_xzzz_0,  \
                             g_0_z_0_xzzz_1,  \
                             g_0_z_0_yyyz_0,  \
                             g_0_z_0_yyyz_1,  \
                             g_0_z_0_yyz_1,   \
                             g_0_z_0_yyzz_0,  \
                             g_0_z_0_yyzz_1,  \
                             g_0_z_0_yzz_1,   \
                             g_0_z_0_yzzz_0,  \
                             g_0_z_0_yzzz_1,  \
                             g_0_z_0_zzz_1,   \
                             g_0_z_0_zzzz_0,  \
                             g_0_z_0_zzzz_1,  \
                             wp_y,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xxxx_0[i] = g_0_z_0_xxxx_0[i] * pb_y + g_0_z_0_xxxx_1[i] * wp_y[i];

        g_0_yz_0_xxxy_0[i] = g_0_y_0_xxxy_0[i] * pb_z + g_0_y_0_xxxy_1[i] * wp_z[i];

        g_0_yz_0_xxxz_0[i] = g_0_z_0_xxxz_0[i] * pb_y + g_0_z_0_xxxz_1[i] * wp_y[i];

        g_0_yz_0_xxyy_0[i] = g_0_y_0_xxyy_0[i] * pb_z + g_0_y_0_xxyy_1[i] * wp_z[i];

        g_0_yz_0_xxyz_0[i] = g_0_z_0_xxz_1[i] * fi_abcd_0 + g_0_z_0_xxyz_0[i] * pb_y + g_0_z_0_xxyz_1[i] * wp_y[i];

        g_0_yz_0_xxzz_0[i] = g_0_z_0_xxzz_0[i] * pb_y + g_0_z_0_xxzz_1[i] * wp_y[i];

        g_0_yz_0_xyyy_0[i] = g_0_y_0_xyyy_0[i] * pb_z + g_0_y_0_xyyy_1[i] * wp_z[i];

        g_0_yz_0_xyyz_0[i] = 2.0 * g_0_z_0_xyz_1[i] * fi_abcd_0 + g_0_z_0_xyyz_0[i] * pb_y + g_0_z_0_xyyz_1[i] * wp_y[i];

        g_0_yz_0_xyzz_0[i] = g_0_z_0_xzz_1[i] * fi_abcd_0 + g_0_z_0_xyzz_0[i] * pb_y + g_0_z_0_xyzz_1[i] * wp_y[i];

        g_0_yz_0_xzzz_0[i] = g_0_z_0_xzzz_0[i] * pb_y + g_0_z_0_xzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyy_0[i] = g_0_y_0_yyyy_0[i] * pb_z + g_0_y_0_yyyy_1[i] * wp_z[i];

        g_0_yz_0_yyyz_0[i] = 3.0 * g_0_z_0_yyz_1[i] * fi_abcd_0 + g_0_z_0_yyyz_0[i] * pb_y + g_0_z_0_yyyz_1[i] * wp_y[i];

        g_0_yz_0_yyzz_0[i] = 2.0 * g_0_z_0_yzz_1[i] * fi_abcd_0 + g_0_z_0_yyzz_0[i] * pb_y + g_0_z_0_yyzz_1[i] * wp_y[i];

        g_0_yz_0_yzzz_0[i] = g_0_z_0_zzz_1[i] * fi_abcd_0 + g_0_z_0_yzzz_0[i] * pb_y + g_0_z_0_yzzz_1[i] * wp_y[i];

        g_0_yz_0_zzzz_0[i] = g_0_z_0_zzzz_0[i] * pb_y + g_0_z_0_zzzz_1[i] * wp_y[i];
    }

    /// Set up 75-90 components of targeted buffer : SDSG

    auto g_0_zz_0_xxxx_0 = pbuffer.data(idx_eri_0_sdsg + 75);

    auto g_0_zz_0_xxxy_0 = pbuffer.data(idx_eri_0_sdsg + 76);

    auto g_0_zz_0_xxxz_0 = pbuffer.data(idx_eri_0_sdsg + 77);

    auto g_0_zz_0_xxyy_0 = pbuffer.data(idx_eri_0_sdsg + 78);

    auto g_0_zz_0_xxyz_0 = pbuffer.data(idx_eri_0_sdsg + 79);

    auto g_0_zz_0_xxzz_0 = pbuffer.data(idx_eri_0_sdsg + 80);

    auto g_0_zz_0_xyyy_0 = pbuffer.data(idx_eri_0_sdsg + 81);

    auto g_0_zz_0_xyyz_0 = pbuffer.data(idx_eri_0_sdsg + 82);

    auto g_0_zz_0_xyzz_0 = pbuffer.data(idx_eri_0_sdsg + 83);

    auto g_0_zz_0_xzzz_0 = pbuffer.data(idx_eri_0_sdsg + 84);

    auto g_0_zz_0_yyyy_0 = pbuffer.data(idx_eri_0_sdsg + 85);

    auto g_0_zz_0_yyyz_0 = pbuffer.data(idx_eri_0_sdsg + 86);

    auto g_0_zz_0_yyzz_0 = pbuffer.data(idx_eri_0_sdsg + 87);

    auto g_0_zz_0_yzzz_0 = pbuffer.data(idx_eri_0_sdsg + 88);

    auto g_0_zz_0_zzzz_0 = pbuffer.data(idx_eri_0_sdsg + 89);

#pragma omp simd aligned(g_0_0_0_xxxx_0,      \
                             g_0_0_0_xxxx_1,  \
                             g_0_0_0_xxxy_0,  \
                             g_0_0_0_xxxy_1,  \
                             g_0_0_0_xxxz_0,  \
                             g_0_0_0_xxxz_1,  \
                             g_0_0_0_xxyy_0,  \
                             g_0_0_0_xxyy_1,  \
                             g_0_0_0_xxyz_0,  \
                             g_0_0_0_xxyz_1,  \
                             g_0_0_0_xxzz_0,  \
                             g_0_0_0_xxzz_1,  \
                             g_0_0_0_xyyy_0,  \
                             g_0_0_0_xyyy_1,  \
                             g_0_0_0_xyyz_0,  \
                             g_0_0_0_xyyz_1,  \
                             g_0_0_0_xyzz_0,  \
                             g_0_0_0_xyzz_1,  \
                             g_0_0_0_xzzz_0,  \
                             g_0_0_0_xzzz_1,  \
                             g_0_0_0_yyyy_0,  \
                             g_0_0_0_yyyy_1,  \
                             g_0_0_0_yyyz_0,  \
                             g_0_0_0_yyyz_1,  \
                             g_0_0_0_yyzz_0,  \
                             g_0_0_0_yyzz_1,  \
                             g_0_0_0_yzzz_0,  \
                             g_0_0_0_yzzz_1,  \
                             g_0_0_0_zzzz_0,  \
                             g_0_0_0_zzzz_1,  \
                             g_0_z_0_xxx_1,   \
                             g_0_z_0_xxxx_0,  \
                             g_0_z_0_xxxx_1,  \
                             g_0_z_0_xxxy_0,  \
                             g_0_z_0_xxxy_1,  \
                             g_0_z_0_xxxz_0,  \
                             g_0_z_0_xxxz_1,  \
                             g_0_z_0_xxy_1,   \
                             g_0_z_0_xxyy_0,  \
                             g_0_z_0_xxyy_1,  \
                             g_0_z_0_xxyz_0,  \
                             g_0_z_0_xxyz_1,  \
                             g_0_z_0_xxz_1,   \
                             g_0_z_0_xxzz_0,  \
                             g_0_z_0_xxzz_1,  \
                             g_0_z_0_xyy_1,   \
                             g_0_z_0_xyyy_0,  \
                             g_0_z_0_xyyy_1,  \
                             g_0_z_0_xyyz_0,  \
                             g_0_z_0_xyyz_1,  \
                             g_0_z_0_xyz_1,   \
                             g_0_z_0_xyzz_0,  \
                             g_0_z_0_xyzz_1,  \
                             g_0_z_0_xzz_1,   \
                             g_0_z_0_xzzz_0,  \
                             g_0_z_0_xzzz_1,  \
                             g_0_z_0_yyy_1,   \
                             g_0_z_0_yyyy_0,  \
                             g_0_z_0_yyyy_1,  \
                             g_0_z_0_yyyz_0,  \
                             g_0_z_0_yyyz_1,  \
                             g_0_z_0_yyz_1,   \
                             g_0_z_0_yyzz_0,  \
                             g_0_z_0_yyzz_1,  \
                             g_0_z_0_yzz_1,   \
                             g_0_z_0_yzzz_0,  \
                             g_0_z_0_yzzz_1,  \
                             g_0_z_0_zzz_1,   \
                             g_0_z_0_zzzz_0,  \
                             g_0_z_0_zzzz_1,  \
                             g_0_zz_0_xxxx_0, \
                             g_0_zz_0_xxxy_0, \
                             g_0_zz_0_xxxz_0, \
                             g_0_zz_0_xxyy_0, \
                             g_0_zz_0_xxyz_0, \
                             g_0_zz_0_xxzz_0, \
                             g_0_zz_0_xyyy_0, \
                             g_0_zz_0_xyyz_0, \
                             g_0_zz_0_xyzz_0, \
                             g_0_zz_0_xzzz_0, \
                             g_0_zz_0_yyyy_0, \
                             g_0_zz_0_yyyz_0, \
                             g_0_zz_0_yyzz_0, \
                             g_0_zz_0_yzzz_0, \
                             g_0_zz_0_zzzz_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xxxx_0[i] = g_0_0_0_xxxx_0[i] * fi_ab_0 - g_0_0_0_xxxx_1[i] * fti_ab_0 + g_0_z_0_xxxx_0[i] * pb_z + g_0_z_0_xxxx_1[i] * wp_z[i];

        g_0_zz_0_xxxy_0[i] = g_0_0_0_xxxy_0[i] * fi_ab_0 - g_0_0_0_xxxy_1[i] * fti_ab_0 + g_0_z_0_xxxy_0[i] * pb_z + g_0_z_0_xxxy_1[i] * wp_z[i];

        g_0_zz_0_xxxz_0[i] = g_0_0_0_xxxz_0[i] * fi_ab_0 - g_0_0_0_xxxz_1[i] * fti_ab_0 + g_0_z_0_xxx_1[i] * fi_abcd_0 + g_0_z_0_xxxz_0[i] * pb_z +
                             g_0_z_0_xxxz_1[i] * wp_z[i];

        g_0_zz_0_xxyy_0[i] = g_0_0_0_xxyy_0[i] * fi_ab_0 - g_0_0_0_xxyy_1[i] * fti_ab_0 + g_0_z_0_xxyy_0[i] * pb_z + g_0_z_0_xxyy_1[i] * wp_z[i];

        g_0_zz_0_xxyz_0[i] = g_0_0_0_xxyz_0[i] * fi_ab_0 - g_0_0_0_xxyz_1[i] * fti_ab_0 + g_0_z_0_xxy_1[i] * fi_abcd_0 + g_0_z_0_xxyz_0[i] * pb_z +
                             g_0_z_0_xxyz_1[i] * wp_z[i];

        g_0_zz_0_xxzz_0[i] = g_0_0_0_xxzz_0[i] * fi_ab_0 - g_0_0_0_xxzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxz_1[i] * fi_abcd_0 +
                             g_0_z_0_xxzz_0[i] * pb_z + g_0_z_0_xxzz_1[i] * wp_z[i];

        g_0_zz_0_xyyy_0[i] = g_0_0_0_xyyy_0[i] * fi_ab_0 - g_0_0_0_xyyy_1[i] * fti_ab_0 + g_0_z_0_xyyy_0[i] * pb_z + g_0_z_0_xyyy_1[i] * wp_z[i];

        g_0_zz_0_xyyz_0[i] = g_0_0_0_xyyz_0[i] * fi_ab_0 - g_0_0_0_xyyz_1[i] * fti_ab_0 + g_0_z_0_xyy_1[i] * fi_abcd_0 + g_0_z_0_xyyz_0[i] * pb_z +
                             g_0_z_0_xyyz_1[i] * wp_z[i];

        g_0_zz_0_xyzz_0[i] = g_0_0_0_xyzz_0[i] * fi_ab_0 - g_0_0_0_xyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xyz_1[i] * fi_abcd_0 +
                             g_0_z_0_xyzz_0[i] * pb_z + g_0_z_0_xyzz_1[i] * wp_z[i];

        g_0_zz_0_xzzz_0[i] = g_0_0_0_xzzz_0[i] * fi_ab_0 - g_0_0_0_xzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xzz_1[i] * fi_abcd_0 +
                             g_0_z_0_xzzz_0[i] * pb_z + g_0_z_0_xzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyy_0[i] = g_0_0_0_yyyy_0[i] * fi_ab_0 - g_0_0_0_yyyy_1[i] * fti_ab_0 + g_0_z_0_yyyy_0[i] * pb_z + g_0_z_0_yyyy_1[i] * wp_z[i];

        g_0_zz_0_yyyz_0[i] = g_0_0_0_yyyz_0[i] * fi_ab_0 - g_0_0_0_yyyz_1[i] * fti_ab_0 + g_0_z_0_yyy_1[i] * fi_abcd_0 + g_0_z_0_yyyz_0[i] * pb_z +
                             g_0_z_0_yyyz_1[i] * wp_z[i];

        g_0_zz_0_yyzz_0[i] = g_0_0_0_yyzz_0[i] * fi_ab_0 - g_0_0_0_yyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_yyz_1[i] * fi_abcd_0 +
                             g_0_z_0_yyzz_0[i] * pb_z + g_0_z_0_yyzz_1[i] * wp_z[i];

        g_0_zz_0_yzzz_0[i] = g_0_0_0_yzzz_0[i] * fi_ab_0 - g_0_0_0_yzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_yzz_1[i] * fi_abcd_0 +
                             g_0_z_0_yzzz_0[i] * pb_z + g_0_z_0_yzzz_1[i] * wp_z[i];

        g_0_zz_0_zzzz_0[i] = g_0_0_0_zzzz_0[i] * fi_ab_0 - g_0_0_0_zzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_zzz_1[i] * fi_abcd_0 +
                             g_0_z_0_zzzz_0[i] * pb_z + g_0_z_0_zzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
