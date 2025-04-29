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

#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsf,
                                 size_t idx_eri_0_dsf,
                                 size_t idx_eri_1_dsf,
                                 size_t idx_eri_1_fsd,
                                 size_t idx_eri_1_fsf,
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

    /// Set up components of auxilary buffer : DSF

    auto g_xx_0_xxx_0 = pbuffer.data(idx_eri_0_dsf);

    auto g_xx_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 1);

    auto g_xx_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 2);

    auto g_xx_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 3);

    auto g_xx_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 4);

    auto g_xx_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 5);

    auto g_xx_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 6);

    auto g_xx_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 7);

    auto g_xx_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 8);

    auto g_xx_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 9);

    auto g_yy_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 30);

    auto g_yy_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 31);

    auto g_yy_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 32);

    auto g_yy_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 33);

    auto g_yy_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 34);

    auto g_yy_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 35);

    auto g_yy_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 36);

    auto g_yy_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 37);

    auto g_yy_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 38);

    auto g_yy_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 39);

    auto g_zz_0_xxx_0 = pbuffer.data(idx_eri_0_dsf + 50);

    auto g_zz_0_xxy_0 = pbuffer.data(idx_eri_0_dsf + 51);

    auto g_zz_0_xxz_0 = pbuffer.data(idx_eri_0_dsf + 52);

    auto g_zz_0_xyy_0 = pbuffer.data(idx_eri_0_dsf + 53);

    auto g_zz_0_xyz_0 = pbuffer.data(idx_eri_0_dsf + 54);

    auto g_zz_0_xzz_0 = pbuffer.data(idx_eri_0_dsf + 55);

    auto g_zz_0_yyy_0 = pbuffer.data(idx_eri_0_dsf + 56);

    auto g_zz_0_yyz_0 = pbuffer.data(idx_eri_0_dsf + 57);

    auto g_zz_0_yzz_0 = pbuffer.data(idx_eri_0_dsf + 58);

    auto g_zz_0_zzz_0 = pbuffer.data(idx_eri_0_dsf + 59);

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

    /// Set up components of auxilary buffer : FSD

    auto g_xxx_0_xx_1 = pbuffer.data(idx_eri_1_fsd);

    auto g_xxx_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 1);

    auto g_xxx_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 2);

    auto g_xxx_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 3);

    auto g_xxx_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 4);

    auto g_xxx_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 5);

    auto g_xxz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 14);

    auto g_xxz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 16);

    auto g_xxz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 17);

    auto g_xyy_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 19);

    auto g_xyy_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 21);

    auto g_xyy_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 22);

    auto g_xzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 32);

    auto g_xzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 34);

    auto g_xzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 35);

    auto g_yyy_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 36);

    auto g_yyy_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 37);

    auto g_yyy_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 38);

    auto g_yyy_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 39);

    auto g_yyy_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 40);

    auto g_yyy_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 41);

    auto g_yyz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 44);

    auto g_yyz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 46);

    auto g_yyz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 47);

    auto g_yzz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 49);

    auto g_yzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 50);

    auto g_yzz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 51);

    auto g_yzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 52);

    auto g_yzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 53);

    auto g_zzz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 54);

    auto g_zzz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 55);

    auto g_zzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 56);

    auto g_zzz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 57);

    auto g_zzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 58);

    auto g_zzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 59);

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

    auto g_xxy_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 11);

    auto g_xxy_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 12);

    auto g_xxy_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 13);

    auto g_xxy_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 15);

    auto g_xxy_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 16);

    auto g_xxz_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 20);

    auto g_xxz_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 21);

    auto g_xxz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 22);

    auto g_xxz_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 23);

    auto g_xxz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 24);

    auto g_xxz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 25);

    auto g_xxz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 27);

    auto g_xxz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 28);

    auto g_xxz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 29);

    auto g_xyy_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 30);

    auto g_xyy_0_xxy_1 = pbuffer.data(idx_eri_1_fsf + 31);

    auto g_xyy_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 33);

    auto g_xyy_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 34);

    auto g_xyy_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 36);

    auto g_xyy_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 37);

    auto g_xyy_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 38);

    auto g_xyy_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 39);

    auto g_xzz_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 50);

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

    auto g_yyz_0_xxz_1 = pbuffer.data(idx_eri_1_fsf + 72);

    auto g_yyz_0_xyy_1 = pbuffer.data(idx_eri_1_fsf + 73);

    auto g_yyz_0_xyz_1 = pbuffer.data(idx_eri_1_fsf + 74);

    auto g_yyz_0_xzz_1 = pbuffer.data(idx_eri_1_fsf + 75);

    auto g_yyz_0_yyy_1 = pbuffer.data(idx_eri_1_fsf + 76);

    auto g_yyz_0_yyz_1 = pbuffer.data(idx_eri_1_fsf + 77);

    auto g_yyz_0_yzz_1 = pbuffer.data(idx_eri_1_fsf + 78);

    auto g_yyz_0_zzz_1 = pbuffer.data(idx_eri_1_fsf + 79);

    auto g_yzz_0_xxx_1 = pbuffer.data(idx_eri_1_fsf + 80);

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

    /// Set up 0-10 components of targeted buffer : GSF

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

    #pragma omp simd aligned(g_xx_0_xxx_0, g_xx_0_xxx_1, g_xx_0_xxy_0, g_xx_0_xxy_1, g_xx_0_xxz_0, g_xx_0_xxz_1, g_xx_0_xyy_0, g_xx_0_xyy_1, g_xx_0_xyz_0, g_xx_0_xyz_1, g_xx_0_xzz_0, g_xx_0_xzz_1, g_xx_0_yyy_0, g_xx_0_yyy_1, g_xx_0_yyz_0, g_xx_0_yyz_1, g_xx_0_yzz_0, g_xx_0_yzz_1, g_xx_0_zzz_0, g_xx_0_zzz_1, g_xxx_0_xx_1, g_xxx_0_xxx_1, g_xxx_0_xxy_1, g_xxx_0_xxz_1, g_xxx_0_xy_1, g_xxx_0_xyy_1, g_xxx_0_xyz_1, g_xxx_0_xz_1, g_xxx_0_xzz_1, g_xxx_0_yy_1, g_xxx_0_yyy_1, g_xxx_0_yyz_1, g_xxx_0_yz_1, g_xxx_0_yzz_1, g_xxx_0_zz_1, g_xxx_0_zzz_1, g_xxxx_0_xxx_0, g_xxxx_0_xxy_0, g_xxxx_0_xxz_0, g_xxxx_0_xyy_0, g_xxxx_0_xyz_0, g_xxxx_0_xzz_0, g_xxxx_0_yyy_0, g_xxxx_0_yyz_0, g_xxxx_0_yzz_0, g_xxxx_0_zzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxx_0[i] = 3.0 * g_xx_0_xxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxx_1[i] * fz_be_0 + 3.0 * g_xxx_0_xx_1[i] * fi_acd_0 + g_xxx_0_xxx_1[i] * wa_x[i];

        g_xxxx_0_xxy_0[i] = 3.0 * g_xx_0_xxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xy_1[i] * fi_acd_0 + g_xxx_0_xxy_1[i] * wa_x[i];

        g_xxxx_0_xxz_0[i] = 3.0 * g_xx_0_xxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xz_1[i] * fi_acd_0 + g_xxx_0_xxz_1[i] * wa_x[i];

        g_xxxx_0_xyy_0[i] = 3.0 * g_xx_0_xyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyy_1[i] * fz_be_0 + g_xxx_0_yy_1[i] * fi_acd_0 + g_xxx_0_xyy_1[i] * wa_x[i];

        g_xxxx_0_xyz_0[i] = 3.0 * g_xx_0_xyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyz_1[i] * fz_be_0 + g_xxx_0_yz_1[i] * fi_acd_0 + g_xxx_0_xyz_1[i] * wa_x[i];

        g_xxxx_0_xzz_0[i] = 3.0 * g_xx_0_xzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzz_1[i] * fz_be_0 + g_xxx_0_zz_1[i] * fi_acd_0 + g_xxx_0_xzz_1[i] * wa_x[i];

        g_xxxx_0_yyy_0[i] = 3.0 * g_xx_0_yyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyy_1[i] * fz_be_0 + g_xxx_0_yyy_1[i] * wa_x[i];

        g_xxxx_0_yyz_0[i] = 3.0 * g_xx_0_yyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyz_1[i] * fz_be_0 + g_xxx_0_yyz_1[i] * wa_x[i];

        g_xxxx_0_yzz_0[i] = 3.0 * g_xx_0_yzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzz_1[i] * fz_be_0 + g_xxx_0_yzz_1[i] * wa_x[i];

        g_xxxx_0_zzz_0[i] = 3.0 * g_xx_0_zzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzz_1[i] * fz_be_0 + g_xxx_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : GSF

    auto g_xxxy_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 10);

    auto g_xxxy_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 11);

    auto g_xxxy_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 12);

    auto g_xxxy_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 13);

    auto g_xxxy_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 14);

    auto g_xxxy_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 15);

    auto g_xxxy_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 16);

    auto g_xxxy_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 17);

    auto g_xxxy_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 18);

    auto g_xxxy_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 19);

    #pragma omp simd aligned(g_xxx_0_xx_1, g_xxx_0_xxx_1, g_xxx_0_xxy_1, g_xxx_0_xxz_1, g_xxx_0_xy_1, g_xxx_0_xyy_1, g_xxx_0_xyz_1, g_xxx_0_xz_1, g_xxx_0_xzz_1, g_xxx_0_yy_1, g_xxx_0_yyy_1, g_xxx_0_yyz_1, g_xxx_0_yz_1, g_xxx_0_yzz_1, g_xxx_0_zz_1, g_xxx_0_zzz_1, g_xxxy_0_xxx_0, g_xxxy_0_xxy_0, g_xxxy_0_xxz_0, g_xxxy_0_xyy_0, g_xxxy_0_xyz_0, g_xxxy_0_xzz_0, g_xxxy_0_yyy_0, g_xxxy_0_yyz_0, g_xxxy_0_yzz_0, g_xxxy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxx_0[i] = g_xxx_0_xxx_1[i] * wa_y[i];

        g_xxxy_0_xxy_0[i] = g_xxx_0_xx_1[i] * fi_acd_0 + g_xxx_0_xxy_1[i] * wa_y[i];

        g_xxxy_0_xxz_0[i] = g_xxx_0_xxz_1[i] * wa_y[i];

        g_xxxy_0_xyy_0[i] = 2.0 * g_xxx_0_xy_1[i] * fi_acd_0 + g_xxx_0_xyy_1[i] * wa_y[i];

        g_xxxy_0_xyz_0[i] = g_xxx_0_xz_1[i] * fi_acd_0 + g_xxx_0_xyz_1[i] * wa_y[i];

        g_xxxy_0_xzz_0[i] = g_xxx_0_xzz_1[i] * wa_y[i];

        g_xxxy_0_yyy_0[i] = 3.0 * g_xxx_0_yy_1[i] * fi_acd_0 + g_xxx_0_yyy_1[i] * wa_y[i];

        g_xxxy_0_yyz_0[i] = 2.0 * g_xxx_0_yz_1[i] * fi_acd_0 + g_xxx_0_yyz_1[i] * wa_y[i];

        g_xxxy_0_yzz_0[i] = g_xxx_0_zz_1[i] * fi_acd_0 + g_xxx_0_yzz_1[i] * wa_y[i];

        g_xxxy_0_zzz_0[i] = g_xxx_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 20-30 components of targeted buffer : GSF

    auto g_xxxz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 20);

    auto g_xxxz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 21);

    auto g_xxxz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 22);

    auto g_xxxz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 23);

    auto g_xxxz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 24);

    auto g_xxxz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 25);

    auto g_xxxz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 26);

    auto g_xxxz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 27);

    auto g_xxxz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 28);

    auto g_xxxz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 29);

    #pragma omp simd aligned(g_xxx_0_xx_1, g_xxx_0_xxx_1, g_xxx_0_xxy_1, g_xxx_0_xxz_1, g_xxx_0_xy_1, g_xxx_0_xyy_1, g_xxx_0_xyz_1, g_xxx_0_xz_1, g_xxx_0_xzz_1, g_xxx_0_yy_1, g_xxx_0_yyy_1, g_xxx_0_yyz_1, g_xxx_0_yz_1, g_xxx_0_yzz_1, g_xxx_0_zz_1, g_xxx_0_zzz_1, g_xxxz_0_xxx_0, g_xxxz_0_xxy_0, g_xxxz_0_xxz_0, g_xxxz_0_xyy_0, g_xxxz_0_xyz_0, g_xxxz_0_xzz_0, g_xxxz_0_yyy_0, g_xxxz_0_yyz_0, g_xxxz_0_yzz_0, g_xxxz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxx_0[i] = g_xxx_0_xxx_1[i] * wa_z[i];

        g_xxxz_0_xxy_0[i] = g_xxx_0_xxy_1[i] * wa_z[i];

        g_xxxz_0_xxz_0[i] = g_xxx_0_xx_1[i] * fi_acd_0 + g_xxx_0_xxz_1[i] * wa_z[i];

        g_xxxz_0_xyy_0[i] = g_xxx_0_xyy_1[i] * wa_z[i];

        g_xxxz_0_xyz_0[i] = g_xxx_0_xy_1[i] * fi_acd_0 + g_xxx_0_xyz_1[i] * wa_z[i];

        g_xxxz_0_xzz_0[i] = 2.0 * g_xxx_0_xz_1[i] * fi_acd_0 + g_xxx_0_xzz_1[i] * wa_z[i];

        g_xxxz_0_yyy_0[i] = g_xxx_0_yyy_1[i] * wa_z[i];

        g_xxxz_0_yyz_0[i] = g_xxx_0_yy_1[i] * fi_acd_0 + g_xxx_0_yyz_1[i] * wa_z[i];

        g_xxxz_0_yzz_0[i] = 2.0 * g_xxx_0_yz_1[i] * fi_acd_0 + g_xxx_0_yzz_1[i] * wa_z[i];

        g_xxxz_0_zzz_0[i] = 3.0 * g_xxx_0_zz_1[i] * fi_acd_0 + g_xxx_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 30-40 components of targeted buffer : GSF

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

    #pragma omp simd aligned(g_xx_0_xxx_0, g_xx_0_xxx_1, g_xx_0_xxz_0, g_xx_0_xxz_1, g_xx_0_xzz_0, g_xx_0_xzz_1, g_xxy_0_xxx_1, g_xxy_0_xxz_1, g_xxy_0_xzz_1, g_xxyy_0_xxx_0, g_xxyy_0_xxy_0, g_xxyy_0_xxz_0, g_xxyy_0_xyy_0, g_xxyy_0_xyz_0, g_xxyy_0_xzz_0, g_xxyy_0_yyy_0, g_xxyy_0_yyz_0, g_xxyy_0_yzz_0, g_xxyy_0_zzz_0, g_xyy_0_xxy_1, g_xyy_0_xy_1, g_xyy_0_xyy_1, g_xyy_0_xyz_1, g_xyy_0_yy_1, g_xyy_0_yyy_1, g_xyy_0_yyz_1, g_xyy_0_yz_1, g_xyy_0_yzz_1, g_xyy_0_zzz_1, g_yy_0_xxy_0, g_yy_0_xxy_1, g_yy_0_xyy_0, g_yy_0_xyy_1, g_yy_0_xyz_0, g_yy_0_xyz_1, g_yy_0_yyy_0, g_yy_0_yyy_1, g_yy_0_yyz_0, g_yy_0_yyz_1, g_yy_0_yzz_0, g_yy_0_yzz_1, g_yy_0_zzz_0, g_yy_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxx_0[i] = g_xx_0_xxx_0[i] * fbe_0 - g_xx_0_xxx_1[i] * fz_be_0 + g_xxy_0_xxx_1[i] * wa_y[i];

        g_xxyy_0_xxy_0[i] = g_yy_0_xxy_0[i] * fbe_0 - g_yy_0_xxy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xy_1[i] * fi_acd_0 + g_xyy_0_xxy_1[i] * wa_x[i];

        g_xxyy_0_xxz_0[i] = g_xx_0_xxz_0[i] * fbe_0 - g_xx_0_xxz_1[i] * fz_be_0 + g_xxy_0_xxz_1[i] * wa_y[i];

        g_xxyy_0_xyy_0[i] = g_yy_0_xyy_0[i] * fbe_0 - g_yy_0_xyy_1[i] * fz_be_0 + g_xyy_0_yy_1[i] * fi_acd_0 + g_xyy_0_xyy_1[i] * wa_x[i];

        g_xxyy_0_xyz_0[i] = g_yy_0_xyz_0[i] * fbe_0 - g_yy_0_xyz_1[i] * fz_be_0 + g_xyy_0_yz_1[i] * fi_acd_0 + g_xyy_0_xyz_1[i] * wa_x[i];

        g_xxyy_0_xzz_0[i] = g_xx_0_xzz_0[i] * fbe_0 - g_xx_0_xzz_1[i] * fz_be_0 + g_xxy_0_xzz_1[i] * wa_y[i];

        g_xxyy_0_yyy_0[i] = g_yy_0_yyy_0[i] * fbe_0 - g_yy_0_yyy_1[i] * fz_be_0 + g_xyy_0_yyy_1[i] * wa_x[i];

        g_xxyy_0_yyz_0[i] = g_yy_0_yyz_0[i] * fbe_0 - g_yy_0_yyz_1[i] * fz_be_0 + g_xyy_0_yyz_1[i] * wa_x[i];

        g_xxyy_0_yzz_0[i] = g_yy_0_yzz_0[i] * fbe_0 - g_yy_0_yzz_1[i] * fz_be_0 + g_xyy_0_yzz_1[i] * wa_x[i];

        g_xxyy_0_zzz_0[i] = g_yy_0_zzz_0[i] * fbe_0 - g_yy_0_zzz_1[i] * fz_be_0 + g_xyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 40-50 components of targeted buffer : GSF

    auto g_xxyz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 40);

    auto g_xxyz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 41);

    auto g_xxyz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 42);

    auto g_xxyz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 43);

    auto g_xxyz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 44);

    auto g_xxyz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 45);

    auto g_xxyz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 46);

    auto g_xxyz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 47);

    auto g_xxyz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 48);

    auto g_xxyz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 49);

    #pragma omp simd aligned(g_xxy_0_xxy_1, g_xxy_0_xyy_1, g_xxy_0_yyy_1, g_xxyz_0_xxx_0, g_xxyz_0_xxy_0, g_xxyz_0_xxz_0, g_xxyz_0_xyy_0, g_xxyz_0_xyz_0, g_xxyz_0_xzz_0, g_xxyz_0_yyy_0, g_xxyz_0_yyz_0, g_xxyz_0_yzz_0, g_xxyz_0_zzz_0, g_xxz_0_xxx_1, g_xxz_0_xxz_1, g_xxz_0_xyz_1, g_xxz_0_xz_1, g_xxz_0_xzz_1, g_xxz_0_yyz_1, g_xxz_0_yz_1, g_xxz_0_yzz_1, g_xxz_0_zz_1, g_xxz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxx_0[i] = g_xxz_0_xxx_1[i] * wa_y[i];

        g_xxyz_0_xxy_0[i] = g_xxy_0_xxy_1[i] * wa_z[i];

        g_xxyz_0_xxz_0[i] = g_xxz_0_xxz_1[i] * wa_y[i];

        g_xxyz_0_xyy_0[i] = g_xxy_0_xyy_1[i] * wa_z[i];

        g_xxyz_0_xyz_0[i] = g_xxz_0_xz_1[i] * fi_acd_0 + g_xxz_0_xyz_1[i] * wa_y[i];

        g_xxyz_0_xzz_0[i] = g_xxz_0_xzz_1[i] * wa_y[i];

        g_xxyz_0_yyy_0[i] = g_xxy_0_yyy_1[i] * wa_z[i];

        g_xxyz_0_yyz_0[i] = 2.0 * g_xxz_0_yz_1[i] * fi_acd_0 + g_xxz_0_yyz_1[i] * wa_y[i];

        g_xxyz_0_yzz_0[i] = g_xxz_0_zz_1[i] * fi_acd_0 + g_xxz_0_yzz_1[i] * wa_y[i];

        g_xxyz_0_zzz_0[i] = g_xxz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 50-60 components of targeted buffer : GSF

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

    #pragma omp simd aligned(g_xx_0_xxx_0, g_xx_0_xxx_1, g_xx_0_xxy_0, g_xx_0_xxy_1, g_xx_0_xyy_0, g_xx_0_xyy_1, g_xxz_0_xxx_1, g_xxz_0_xxy_1, g_xxz_0_xyy_1, g_xxzz_0_xxx_0, g_xxzz_0_xxy_0, g_xxzz_0_xxz_0, g_xxzz_0_xyy_0, g_xxzz_0_xyz_0, g_xxzz_0_xzz_0, g_xxzz_0_yyy_0, g_xxzz_0_yyz_0, g_xxzz_0_yzz_0, g_xxzz_0_zzz_0, g_xzz_0_xxz_1, g_xzz_0_xyz_1, g_xzz_0_xz_1, g_xzz_0_xzz_1, g_xzz_0_yyy_1, g_xzz_0_yyz_1, g_xzz_0_yz_1, g_xzz_0_yzz_1, g_xzz_0_zz_1, g_xzz_0_zzz_1, g_zz_0_xxz_0, g_zz_0_xxz_1, g_zz_0_xyz_0, g_zz_0_xyz_1, g_zz_0_xzz_0, g_zz_0_xzz_1, g_zz_0_yyy_0, g_zz_0_yyy_1, g_zz_0_yyz_0, g_zz_0_yyz_1, g_zz_0_yzz_0, g_zz_0_yzz_1, g_zz_0_zzz_0, g_zz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxx_0[i] = g_xx_0_xxx_0[i] * fbe_0 - g_xx_0_xxx_1[i] * fz_be_0 + g_xxz_0_xxx_1[i] * wa_z[i];

        g_xxzz_0_xxy_0[i] = g_xx_0_xxy_0[i] * fbe_0 - g_xx_0_xxy_1[i] * fz_be_0 + g_xxz_0_xxy_1[i] * wa_z[i];

        g_xxzz_0_xxz_0[i] = g_zz_0_xxz_0[i] * fbe_0 - g_zz_0_xxz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xz_1[i] * fi_acd_0 + g_xzz_0_xxz_1[i] * wa_x[i];

        g_xxzz_0_xyy_0[i] = g_xx_0_xyy_0[i] * fbe_0 - g_xx_0_xyy_1[i] * fz_be_0 + g_xxz_0_xyy_1[i] * wa_z[i];

        g_xxzz_0_xyz_0[i] = g_zz_0_xyz_0[i] * fbe_0 - g_zz_0_xyz_1[i] * fz_be_0 + g_xzz_0_yz_1[i] * fi_acd_0 + g_xzz_0_xyz_1[i] * wa_x[i];

        g_xxzz_0_xzz_0[i] = g_zz_0_xzz_0[i] * fbe_0 - g_zz_0_xzz_1[i] * fz_be_0 + g_xzz_0_zz_1[i] * fi_acd_0 + g_xzz_0_xzz_1[i] * wa_x[i];

        g_xxzz_0_yyy_0[i] = g_zz_0_yyy_0[i] * fbe_0 - g_zz_0_yyy_1[i] * fz_be_0 + g_xzz_0_yyy_1[i] * wa_x[i];

        g_xxzz_0_yyz_0[i] = g_zz_0_yyz_0[i] * fbe_0 - g_zz_0_yyz_1[i] * fz_be_0 + g_xzz_0_yyz_1[i] * wa_x[i];

        g_xxzz_0_yzz_0[i] = g_zz_0_yzz_0[i] * fbe_0 - g_zz_0_yzz_1[i] * fz_be_0 + g_xzz_0_yzz_1[i] * wa_x[i];

        g_xxzz_0_zzz_0[i] = g_zz_0_zzz_0[i] * fbe_0 - g_zz_0_zzz_1[i] * fz_be_0 + g_xzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 60-70 components of targeted buffer : GSF

    auto g_xyyy_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 60);

    auto g_xyyy_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 61);

    auto g_xyyy_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 62);

    auto g_xyyy_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 63);

    auto g_xyyy_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 64);

    auto g_xyyy_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 65);

    auto g_xyyy_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 66);

    auto g_xyyy_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 67);

    auto g_xyyy_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 68);

    auto g_xyyy_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 69);

    #pragma omp simd aligned(g_xyyy_0_xxx_0, g_xyyy_0_xxy_0, g_xyyy_0_xxz_0, g_xyyy_0_xyy_0, g_xyyy_0_xyz_0, g_xyyy_0_xzz_0, g_xyyy_0_yyy_0, g_xyyy_0_yyz_0, g_xyyy_0_yzz_0, g_xyyy_0_zzz_0, g_yyy_0_xx_1, g_yyy_0_xxx_1, g_yyy_0_xxy_1, g_yyy_0_xxz_1, g_yyy_0_xy_1, g_yyy_0_xyy_1, g_yyy_0_xyz_1, g_yyy_0_xz_1, g_yyy_0_xzz_1, g_yyy_0_yy_1, g_yyy_0_yyy_1, g_yyy_0_yyz_1, g_yyy_0_yz_1, g_yyy_0_yzz_1, g_yyy_0_zz_1, g_yyy_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxx_0[i] = 3.0 * g_yyy_0_xx_1[i] * fi_acd_0 + g_yyy_0_xxx_1[i] * wa_x[i];

        g_xyyy_0_xxy_0[i] = 2.0 * g_yyy_0_xy_1[i] * fi_acd_0 + g_yyy_0_xxy_1[i] * wa_x[i];

        g_xyyy_0_xxz_0[i] = 2.0 * g_yyy_0_xz_1[i] * fi_acd_0 + g_yyy_0_xxz_1[i] * wa_x[i];

        g_xyyy_0_xyy_0[i] = g_yyy_0_yy_1[i] * fi_acd_0 + g_yyy_0_xyy_1[i] * wa_x[i];

        g_xyyy_0_xyz_0[i] = g_yyy_0_yz_1[i] * fi_acd_0 + g_yyy_0_xyz_1[i] * wa_x[i];

        g_xyyy_0_xzz_0[i] = g_yyy_0_zz_1[i] * fi_acd_0 + g_yyy_0_xzz_1[i] * wa_x[i];

        g_xyyy_0_yyy_0[i] = g_yyy_0_yyy_1[i] * wa_x[i];

        g_xyyy_0_yyz_0[i] = g_yyy_0_yyz_1[i] * wa_x[i];

        g_xyyy_0_yzz_0[i] = g_yyy_0_yzz_1[i] * wa_x[i];

        g_xyyy_0_zzz_0[i] = g_yyy_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 70-80 components of targeted buffer : GSF

    auto g_xyyz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 70);

    auto g_xyyz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 71);

    auto g_xyyz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 72);

    auto g_xyyz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 73);

    auto g_xyyz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 74);

    auto g_xyyz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 75);

    auto g_xyyz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 76);

    auto g_xyyz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 77);

    auto g_xyyz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 78);

    auto g_xyyz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 79);

    #pragma omp simd aligned(g_xyy_0_xxx_1, g_xyy_0_xxy_1, g_xyy_0_xyy_1, g_xyyz_0_xxx_0, g_xyyz_0_xxy_0, g_xyyz_0_xxz_0, g_xyyz_0_xyy_0, g_xyyz_0_xyz_0, g_xyyz_0_xzz_0, g_xyyz_0_yyy_0, g_xyyz_0_yyz_0, g_xyyz_0_yzz_0, g_xyyz_0_zzz_0, g_yyz_0_xxz_1, g_yyz_0_xyz_1, g_yyz_0_xz_1, g_yyz_0_xzz_1, g_yyz_0_yyy_1, g_yyz_0_yyz_1, g_yyz_0_yz_1, g_yyz_0_yzz_1, g_yyz_0_zz_1, g_yyz_0_zzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxx_0[i] = g_xyy_0_xxx_1[i] * wa_z[i];

        g_xyyz_0_xxy_0[i] = g_xyy_0_xxy_1[i] * wa_z[i];

        g_xyyz_0_xxz_0[i] = 2.0 * g_yyz_0_xz_1[i] * fi_acd_0 + g_yyz_0_xxz_1[i] * wa_x[i];

        g_xyyz_0_xyy_0[i] = g_xyy_0_xyy_1[i] * wa_z[i];

        g_xyyz_0_xyz_0[i] = g_yyz_0_yz_1[i] * fi_acd_0 + g_yyz_0_xyz_1[i] * wa_x[i];

        g_xyyz_0_xzz_0[i] = g_yyz_0_zz_1[i] * fi_acd_0 + g_yyz_0_xzz_1[i] * wa_x[i];

        g_xyyz_0_yyy_0[i] = g_yyz_0_yyy_1[i] * wa_x[i];

        g_xyyz_0_yyz_0[i] = g_yyz_0_yyz_1[i] * wa_x[i];

        g_xyyz_0_yzz_0[i] = g_yyz_0_yzz_1[i] * wa_x[i];

        g_xyyz_0_zzz_0[i] = g_yyz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 80-90 components of targeted buffer : GSF

    auto g_xyzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 80);

    auto g_xyzz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 81);

    auto g_xyzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 82);

    auto g_xyzz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 83);

    auto g_xyzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 84);

    auto g_xyzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 85);

    auto g_xyzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 86);

    auto g_xyzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 87);

    auto g_xyzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 88);

    auto g_xyzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 89);

    #pragma omp simd aligned(g_xyzz_0_xxx_0, g_xyzz_0_xxy_0, g_xyzz_0_xxz_0, g_xyzz_0_xyy_0, g_xyzz_0_xyz_0, g_xyzz_0_xzz_0, g_xyzz_0_yyy_0, g_xyzz_0_yyz_0, g_xyzz_0_yzz_0, g_xyzz_0_zzz_0, g_xzz_0_xxx_1, g_xzz_0_xxz_1, g_xzz_0_xzz_1, g_yzz_0_xxy_1, g_yzz_0_xy_1, g_yzz_0_xyy_1, g_yzz_0_xyz_1, g_yzz_0_yy_1, g_yzz_0_yyy_1, g_yzz_0_yyz_1, g_yzz_0_yz_1, g_yzz_0_yzz_1, g_yzz_0_zzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxx_0[i] = g_xzz_0_xxx_1[i] * wa_y[i];

        g_xyzz_0_xxy_0[i] = 2.0 * g_yzz_0_xy_1[i] * fi_acd_0 + g_yzz_0_xxy_1[i] * wa_x[i];

        g_xyzz_0_xxz_0[i] = g_xzz_0_xxz_1[i] * wa_y[i];

        g_xyzz_0_xyy_0[i] = g_yzz_0_yy_1[i] * fi_acd_0 + g_yzz_0_xyy_1[i] * wa_x[i];

        g_xyzz_0_xyz_0[i] = g_yzz_0_yz_1[i] * fi_acd_0 + g_yzz_0_xyz_1[i] * wa_x[i];

        g_xyzz_0_xzz_0[i] = g_xzz_0_xzz_1[i] * wa_y[i];

        g_xyzz_0_yyy_0[i] = g_yzz_0_yyy_1[i] * wa_x[i];

        g_xyzz_0_yyz_0[i] = g_yzz_0_yyz_1[i] * wa_x[i];

        g_xyzz_0_yzz_0[i] = g_yzz_0_yzz_1[i] * wa_x[i];

        g_xyzz_0_zzz_0[i] = g_yzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 90-100 components of targeted buffer : GSF

    auto g_xzzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 90);

    auto g_xzzz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 91);

    auto g_xzzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 92);

    auto g_xzzz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 93);

    auto g_xzzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 94);

    auto g_xzzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 95);

    auto g_xzzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 96);

    auto g_xzzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 97);

    auto g_xzzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 98);

    auto g_xzzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 99);

    #pragma omp simd aligned(g_xzzz_0_xxx_0, g_xzzz_0_xxy_0, g_xzzz_0_xxz_0, g_xzzz_0_xyy_0, g_xzzz_0_xyz_0, g_xzzz_0_xzz_0, g_xzzz_0_yyy_0, g_xzzz_0_yyz_0, g_xzzz_0_yzz_0, g_xzzz_0_zzz_0, g_zzz_0_xx_1, g_zzz_0_xxx_1, g_zzz_0_xxy_1, g_zzz_0_xxz_1, g_zzz_0_xy_1, g_zzz_0_xyy_1, g_zzz_0_xyz_1, g_zzz_0_xz_1, g_zzz_0_xzz_1, g_zzz_0_yy_1, g_zzz_0_yyy_1, g_zzz_0_yyz_1, g_zzz_0_yz_1, g_zzz_0_yzz_1, g_zzz_0_zz_1, g_zzz_0_zzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxx_0[i] = 3.0 * g_zzz_0_xx_1[i] * fi_acd_0 + g_zzz_0_xxx_1[i] * wa_x[i];

        g_xzzz_0_xxy_0[i] = 2.0 * g_zzz_0_xy_1[i] * fi_acd_0 + g_zzz_0_xxy_1[i] * wa_x[i];

        g_xzzz_0_xxz_0[i] = 2.0 * g_zzz_0_xz_1[i] * fi_acd_0 + g_zzz_0_xxz_1[i] * wa_x[i];

        g_xzzz_0_xyy_0[i] = g_zzz_0_yy_1[i] * fi_acd_0 + g_zzz_0_xyy_1[i] * wa_x[i];

        g_xzzz_0_xyz_0[i] = g_zzz_0_yz_1[i] * fi_acd_0 + g_zzz_0_xyz_1[i] * wa_x[i];

        g_xzzz_0_xzz_0[i] = g_zzz_0_zz_1[i] * fi_acd_0 + g_zzz_0_xzz_1[i] * wa_x[i];

        g_xzzz_0_yyy_0[i] = g_zzz_0_yyy_1[i] * wa_x[i];

        g_xzzz_0_yyz_0[i] = g_zzz_0_yyz_1[i] * wa_x[i];

        g_xzzz_0_yzz_0[i] = g_zzz_0_yzz_1[i] * wa_x[i];

        g_xzzz_0_zzz_0[i] = g_zzz_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 100-110 components of targeted buffer : GSF

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

    #pragma omp simd aligned(g_yy_0_xxx_0, g_yy_0_xxx_1, g_yy_0_xxy_0, g_yy_0_xxy_1, g_yy_0_xxz_0, g_yy_0_xxz_1, g_yy_0_xyy_0, g_yy_0_xyy_1, g_yy_0_xyz_0, g_yy_0_xyz_1, g_yy_0_xzz_0, g_yy_0_xzz_1, g_yy_0_yyy_0, g_yy_0_yyy_1, g_yy_0_yyz_0, g_yy_0_yyz_1, g_yy_0_yzz_0, g_yy_0_yzz_1, g_yy_0_zzz_0, g_yy_0_zzz_1, g_yyy_0_xx_1, g_yyy_0_xxx_1, g_yyy_0_xxy_1, g_yyy_0_xxz_1, g_yyy_0_xy_1, g_yyy_0_xyy_1, g_yyy_0_xyz_1, g_yyy_0_xz_1, g_yyy_0_xzz_1, g_yyy_0_yy_1, g_yyy_0_yyy_1, g_yyy_0_yyz_1, g_yyy_0_yz_1, g_yyy_0_yzz_1, g_yyy_0_zz_1, g_yyy_0_zzz_1, g_yyyy_0_xxx_0, g_yyyy_0_xxy_0, g_yyyy_0_xxz_0, g_yyyy_0_xyy_0, g_yyyy_0_xyz_0, g_yyyy_0_xzz_0, g_yyyy_0_yyy_0, g_yyyy_0_yyz_0, g_yyyy_0_yzz_0, g_yyyy_0_zzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxx_0[i] = 3.0 * g_yy_0_xxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxx_1[i] * fz_be_0 + g_yyy_0_xxx_1[i] * wa_y[i];

        g_yyyy_0_xxy_0[i] = 3.0 * g_yy_0_xxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxy_1[i] * fz_be_0 + g_yyy_0_xx_1[i] * fi_acd_0 + g_yyy_0_xxy_1[i] * wa_y[i];

        g_yyyy_0_xxz_0[i] = 3.0 * g_yy_0_xxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxz_1[i] * fz_be_0 + g_yyy_0_xxz_1[i] * wa_y[i];

        g_yyyy_0_xyy_0[i] = 3.0 * g_yy_0_xyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xy_1[i] * fi_acd_0 + g_yyy_0_xyy_1[i] * wa_y[i];

        g_yyyy_0_xyz_0[i] = 3.0 * g_yy_0_xyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyz_1[i] * fz_be_0 + g_yyy_0_xz_1[i] * fi_acd_0 + g_yyy_0_xyz_1[i] * wa_y[i];

        g_yyyy_0_xzz_0[i] = 3.0 * g_yy_0_xzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzz_1[i] * fz_be_0 + g_yyy_0_xzz_1[i] * wa_y[i];

        g_yyyy_0_yyy_0[i] = 3.0 * g_yy_0_yyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_yy_1[i] * fi_acd_0 + g_yyy_0_yyy_1[i] * wa_y[i];

        g_yyyy_0_yyz_0[i] = 3.0 * g_yy_0_yyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yz_1[i] * fi_acd_0 + g_yyy_0_yyz_1[i] * wa_y[i];

        g_yyyy_0_yzz_0[i] = 3.0 * g_yy_0_yzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzz_1[i] * fz_be_0 + g_yyy_0_zz_1[i] * fi_acd_0 + g_yyy_0_yzz_1[i] * wa_y[i];

        g_yyyy_0_zzz_0[i] = 3.0 * g_yy_0_zzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzz_1[i] * fz_be_0 + g_yyy_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 110-120 components of targeted buffer : GSF

    auto g_yyyz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 110);

    auto g_yyyz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 111);

    auto g_yyyz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 112);

    auto g_yyyz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 113);

    auto g_yyyz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 114);

    auto g_yyyz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 115);

    auto g_yyyz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 116);

    auto g_yyyz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 117);

    auto g_yyyz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 118);

    auto g_yyyz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 119);

    #pragma omp simd aligned(g_yyy_0_xx_1, g_yyy_0_xxx_1, g_yyy_0_xxy_1, g_yyy_0_xxz_1, g_yyy_0_xy_1, g_yyy_0_xyy_1, g_yyy_0_xyz_1, g_yyy_0_xz_1, g_yyy_0_xzz_1, g_yyy_0_yy_1, g_yyy_0_yyy_1, g_yyy_0_yyz_1, g_yyy_0_yz_1, g_yyy_0_yzz_1, g_yyy_0_zz_1, g_yyy_0_zzz_1, g_yyyz_0_xxx_0, g_yyyz_0_xxy_0, g_yyyz_0_xxz_0, g_yyyz_0_xyy_0, g_yyyz_0_xyz_0, g_yyyz_0_xzz_0, g_yyyz_0_yyy_0, g_yyyz_0_yyz_0, g_yyyz_0_yzz_0, g_yyyz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxx_0[i] = g_yyy_0_xxx_1[i] * wa_z[i];

        g_yyyz_0_xxy_0[i] = g_yyy_0_xxy_1[i] * wa_z[i];

        g_yyyz_0_xxz_0[i] = g_yyy_0_xx_1[i] * fi_acd_0 + g_yyy_0_xxz_1[i] * wa_z[i];

        g_yyyz_0_xyy_0[i] = g_yyy_0_xyy_1[i] * wa_z[i];

        g_yyyz_0_xyz_0[i] = g_yyy_0_xy_1[i] * fi_acd_0 + g_yyy_0_xyz_1[i] * wa_z[i];

        g_yyyz_0_xzz_0[i] = 2.0 * g_yyy_0_xz_1[i] * fi_acd_0 + g_yyy_0_xzz_1[i] * wa_z[i];

        g_yyyz_0_yyy_0[i] = g_yyy_0_yyy_1[i] * wa_z[i];

        g_yyyz_0_yyz_0[i] = g_yyy_0_yy_1[i] * fi_acd_0 + g_yyy_0_yyz_1[i] * wa_z[i];

        g_yyyz_0_yzz_0[i] = 2.0 * g_yyy_0_yz_1[i] * fi_acd_0 + g_yyy_0_yzz_1[i] * wa_z[i];

        g_yyyz_0_zzz_0[i] = 3.0 * g_yyy_0_zz_1[i] * fi_acd_0 + g_yyy_0_zzz_1[i] * wa_z[i];
    }

    /// Set up 120-130 components of targeted buffer : GSF

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

    #pragma omp simd aligned(g_yy_0_xxy_0, g_yy_0_xxy_1, g_yy_0_xyy_0, g_yy_0_xyy_1, g_yy_0_yyy_0, g_yy_0_yyy_1, g_yyz_0_xxy_1, g_yyz_0_xyy_1, g_yyz_0_yyy_1, g_yyzz_0_xxx_0, g_yyzz_0_xxy_0, g_yyzz_0_xxz_0, g_yyzz_0_xyy_0, g_yyzz_0_xyz_0, g_yyzz_0_xzz_0, g_yyzz_0_yyy_0, g_yyzz_0_yyz_0, g_yyzz_0_yzz_0, g_yyzz_0_zzz_0, g_yzz_0_xxx_1, g_yzz_0_xxz_1, g_yzz_0_xyz_1, g_yzz_0_xz_1, g_yzz_0_xzz_1, g_yzz_0_yyz_1, g_yzz_0_yz_1, g_yzz_0_yzz_1, g_yzz_0_zz_1, g_yzz_0_zzz_1, g_zz_0_xxx_0, g_zz_0_xxx_1, g_zz_0_xxz_0, g_zz_0_xxz_1, g_zz_0_xyz_0, g_zz_0_xyz_1, g_zz_0_xzz_0, g_zz_0_xzz_1, g_zz_0_yyz_0, g_zz_0_yyz_1, g_zz_0_yzz_0, g_zz_0_yzz_1, g_zz_0_zzz_0, g_zz_0_zzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxx_0[i] = g_zz_0_xxx_0[i] * fbe_0 - g_zz_0_xxx_1[i] * fz_be_0 + g_yzz_0_xxx_1[i] * wa_y[i];

        g_yyzz_0_xxy_0[i] = g_yy_0_xxy_0[i] * fbe_0 - g_yy_0_xxy_1[i] * fz_be_0 + g_yyz_0_xxy_1[i] * wa_z[i];

        g_yyzz_0_xxz_0[i] = g_zz_0_xxz_0[i] * fbe_0 - g_zz_0_xxz_1[i] * fz_be_0 + g_yzz_0_xxz_1[i] * wa_y[i];

        g_yyzz_0_xyy_0[i] = g_yy_0_xyy_0[i] * fbe_0 - g_yy_0_xyy_1[i] * fz_be_0 + g_yyz_0_xyy_1[i] * wa_z[i];

        g_yyzz_0_xyz_0[i] = g_zz_0_xyz_0[i] * fbe_0 - g_zz_0_xyz_1[i] * fz_be_0 + g_yzz_0_xz_1[i] * fi_acd_0 + g_yzz_0_xyz_1[i] * wa_y[i];

        g_yyzz_0_xzz_0[i] = g_zz_0_xzz_0[i] * fbe_0 - g_zz_0_xzz_1[i] * fz_be_0 + g_yzz_0_xzz_1[i] * wa_y[i];

        g_yyzz_0_yyy_0[i] = g_yy_0_yyy_0[i] * fbe_0 - g_yy_0_yyy_1[i] * fz_be_0 + g_yyz_0_yyy_1[i] * wa_z[i];

        g_yyzz_0_yyz_0[i] = g_zz_0_yyz_0[i] * fbe_0 - g_zz_0_yyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yz_1[i] * fi_acd_0 + g_yzz_0_yyz_1[i] * wa_y[i];

        g_yyzz_0_yzz_0[i] = g_zz_0_yzz_0[i] * fbe_0 - g_zz_0_yzz_1[i] * fz_be_0 + g_yzz_0_zz_1[i] * fi_acd_0 + g_yzz_0_yzz_1[i] * wa_y[i];

        g_yyzz_0_zzz_0[i] = g_zz_0_zzz_0[i] * fbe_0 - g_zz_0_zzz_1[i] * fz_be_0 + g_yzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 130-140 components of targeted buffer : GSF

    auto g_yzzz_0_xxx_0 = pbuffer.data(idx_eri_0_gsf + 130);

    auto g_yzzz_0_xxy_0 = pbuffer.data(idx_eri_0_gsf + 131);

    auto g_yzzz_0_xxz_0 = pbuffer.data(idx_eri_0_gsf + 132);

    auto g_yzzz_0_xyy_0 = pbuffer.data(idx_eri_0_gsf + 133);

    auto g_yzzz_0_xyz_0 = pbuffer.data(idx_eri_0_gsf + 134);

    auto g_yzzz_0_xzz_0 = pbuffer.data(idx_eri_0_gsf + 135);

    auto g_yzzz_0_yyy_0 = pbuffer.data(idx_eri_0_gsf + 136);

    auto g_yzzz_0_yyz_0 = pbuffer.data(idx_eri_0_gsf + 137);

    auto g_yzzz_0_yzz_0 = pbuffer.data(idx_eri_0_gsf + 138);

    auto g_yzzz_0_zzz_0 = pbuffer.data(idx_eri_0_gsf + 139);

    #pragma omp simd aligned(g_yzzz_0_xxx_0, g_yzzz_0_xxy_0, g_yzzz_0_xxz_0, g_yzzz_0_xyy_0, g_yzzz_0_xyz_0, g_yzzz_0_xzz_0, g_yzzz_0_yyy_0, g_yzzz_0_yyz_0, g_yzzz_0_yzz_0, g_yzzz_0_zzz_0, g_zzz_0_xx_1, g_zzz_0_xxx_1, g_zzz_0_xxy_1, g_zzz_0_xxz_1, g_zzz_0_xy_1, g_zzz_0_xyy_1, g_zzz_0_xyz_1, g_zzz_0_xz_1, g_zzz_0_xzz_1, g_zzz_0_yy_1, g_zzz_0_yyy_1, g_zzz_0_yyz_1, g_zzz_0_yz_1, g_zzz_0_yzz_1, g_zzz_0_zz_1, g_zzz_0_zzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxx_0[i] = g_zzz_0_xxx_1[i] * wa_y[i];

        g_yzzz_0_xxy_0[i] = g_zzz_0_xx_1[i] * fi_acd_0 + g_zzz_0_xxy_1[i] * wa_y[i];

        g_yzzz_0_xxz_0[i] = g_zzz_0_xxz_1[i] * wa_y[i];

        g_yzzz_0_xyy_0[i] = 2.0 * g_zzz_0_xy_1[i] * fi_acd_0 + g_zzz_0_xyy_1[i] * wa_y[i];

        g_yzzz_0_xyz_0[i] = g_zzz_0_xz_1[i] * fi_acd_0 + g_zzz_0_xyz_1[i] * wa_y[i];

        g_yzzz_0_xzz_0[i] = g_zzz_0_xzz_1[i] * wa_y[i];

        g_yzzz_0_yyy_0[i] = 3.0 * g_zzz_0_yy_1[i] * fi_acd_0 + g_zzz_0_yyy_1[i] * wa_y[i];

        g_yzzz_0_yyz_0[i] = 2.0 * g_zzz_0_yz_1[i] * fi_acd_0 + g_zzz_0_yyz_1[i] * wa_y[i];

        g_yzzz_0_yzz_0[i] = g_zzz_0_zz_1[i] * fi_acd_0 + g_zzz_0_yzz_1[i] * wa_y[i];

        g_yzzz_0_zzz_0[i] = g_zzz_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 140-150 components of targeted buffer : GSF

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

    #pragma omp simd aligned(g_zz_0_xxx_0, g_zz_0_xxx_1, g_zz_0_xxy_0, g_zz_0_xxy_1, g_zz_0_xxz_0, g_zz_0_xxz_1, g_zz_0_xyy_0, g_zz_0_xyy_1, g_zz_0_xyz_0, g_zz_0_xyz_1, g_zz_0_xzz_0, g_zz_0_xzz_1, g_zz_0_yyy_0, g_zz_0_yyy_1, g_zz_0_yyz_0, g_zz_0_yyz_1, g_zz_0_yzz_0, g_zz_0_yzz_1, g_zz_0_zzz_0, g_zz_0_zzz_1, g_zzz_0_xx_1, g_zzz_0_xxx_1, g_zzz_0_xxy_1, g_zzz_0_xxz_1, g_zzz_0_xy_1, g_zzz_0_xyy_1, g_zzz_0_xyz_1, g_zzz_0_xz_1, g_zzz_0_xzz_1, g_zzz_0_yy_1, g_zzz_0_yyy_1, g_zzz_0_yyz_1, g_zzz_0_yz_1, g_zzz_0_yzz_1, g_zzz_0_zz_1, g_zzz_0_zzz_1, g_zzzz_0_xxx_0, g_zzzz_0_xxy_0, g_zzzz_0_xxz_0, g_zzzz_0_xyy_0, g_zzzz_0_xyz_0, g_zzzz_0_xzz_0, g_zzzz_0_yyy_0, g_zzzz_0_yyz_0, g_zzzz_0_yzz_0, g_zzzz_0_zzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxx_0[i] = 3.0 * g_zz_0_xxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxx_1[i] * fz_be_0 + g_zzz_0_xxx_1[i] * wa_z[i];

        g_zzzz_0_xxy_0[i] = 3.0 * g_zz_0_xxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxy_1[i] * fz_be_0 + g_zzz_0_xxy_1[i] * wa_z[i];

        g_zzzz_0_xxz_0[i] = 3.0 * g_zz_0_xxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxz_1[i] * fz_be_0 + g_zzz_0_xx_1[i] * fi_acd_0 + g_zzz_0_xxz_1[i] * wa_z[i];

        g_zzzz_0_xyy_0[i] = 3.0 * g_zz_0_xyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyy_1[i] * fz_be_0 + g_zzz_0_xyy_1[i] * wa_z[i];

        g_zzzz_0_xyz_0[i] = 3.0 * g_zz_0_xyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyz_1[i] * fz_be_0 + g_zzz_0_xy_1[i] * fi_acd_0 + g_zzz_0_xyz_1[i] * wa_z[i];

        g_zzzz_0_xzz_0[i] = 3.0 * g_zz_0_xzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xz_1[i] * fi_acd_0 + g_zzz_0_xzz_1[i] * wa_z[i];

        g_zzzz_0_yyy_0[i] = 3.0 * g_zz_0_yyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyy_1[i] * fz_be_0 + g_zzz_0_yyy_1[i] * wa_z[i];

        g_zzzz_0_yyz_0[i] = 3.0 * g_zz_0_yyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyz_1[i] * fz_be_0 + g_zzz_0_yy_1[i] * fi_acd_0 + g_zzz_0_yyz_1[i] * wa_z[i];

        g_zzzz_0_yzz_0[i] = 3.0 * g_zz_0_yzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yz_1[i] * fi_acd_0 + g_zzz_0_yzz_1[i] * wa_z[i];

        g_zzzz_0_zzz_0[i] = 3.0 * g_zz_0_zzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_zz_1[i] * fi_acd_0 + g_zzz_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

