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

#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psf,
                                 size_t idx_eri_1_ssd,
                                 size_t idx_eri_1_ssf,
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

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_1 = pbuffer.data(idx_eri_1_ssd);

    auto g_0_0_xy_1 = pbuffer.data(idx_eri_1_ssd + 1);

    auto g_0_0_xz_1 = pbuffer.data(idx_eri_1_ssd + 2);

    auto g_0_0_yy_1 = pbuffer.data(idx_eri_1_ssd + 3);

    auto g_0_0_yz_1 = pbuffer.data(idx_eri_1_ssd + 4);

    auto g_0_0_zz_1 = pbuffer.data(idx_eri_1_ssd + 5);

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_1 = pbuffer.data(idx_eri_1_ssf);

    auto g_0_0_xxy_1 = pbuffer.data(idx_eri_1_ssf + 1);

    auto g_0_0_xxz_1 = pbuffer.data(idx_eri_1_ssf + 2);

    auto g_0_0_xyy_1 = pbuffer.data(idx_eri_1_ssf + 3);

    auto g_0_0_xyz_1 = pbuffer.data(idx_eri_1_ssf + 4);

    auto g_0_0_xzz_1 = pbuffer.data(idx_eri_1_ssf + 5);

    auto g_0_0_yyy_1 = pbuffer.data(idx_eri_1_ssf + 6);

    auto g_0_0_yyz_1 = pbuffer.data(idx_eri_1_ssf + 7);

    auto g_0_0_yzz_1 = pbuffer.data(idx_eri_1_ssf + 8);

    auto g_0_0_zzz_1 = pbuffer.data(idx_eri_1_ssf + 9);

    /// Set up 0-10 components of targeted buffer : PSF

    auto g_x_0_xxx_0 = pbuffer.data(idx_eri_0_psf);

    auto g_x_0_xxy_0 = pbuffer.data(idx_eri_0_psf + 1);

    auto g_x_0_xxz_0 = pbuffer.data(idx_eri_0_psf + 2);

    auto g_x_0_xyy_0 = pbuffer.data(idx_eri_0_psf + 3);

    auto g_x_0_xyz_0 = pbuffer.data(idx_eri_0_psf + 4);

    auto g_x_0_xzz_0 = pbuffer.data(idx_eri_0_psf + 5);

    auto g_x_0_yyy_0 = pbuffer.data(idx_eri_0_psf + 6);

    auto g_x_0_yyz_0 = pbuffer.data(idx_eri_0_psf + 7);

    auto g_x_0_yzz_0 = pbuffer.data(idx_eri_0_psf + 8);

    auto g_x_0_zzz_0 = pbuffer.data(idx_eri_0_psf + 9);

    #pragma omp simd aligned(g_0_0_xx_1, g_0_0_xxx_1, g_0_0_xxy_1, g_0_0_xxz_1, g_0_0_xy_1, g_0_0_xyy_1, g_0_0_xyz_1, g_0_0_xz_1, g_0_0_xzz_1, g_0_0_yy_1, g_0_0_yyy_1, g_0_0_yyz_1, g_0_0_yz_1, g_0_0_yzz_1, g_0_0_zz_1, g_0_0_zzz_1, g_x_0_xxx_0, g_x_0_xxy_0, g_x_0_xxz_0, g_x_0_xyy_0, g_x_0_xyz_0, g_x_0_xzz_0, g_x_0_yyy_0, g_x_0_yyz_0, g_x_0_yzz_0, g_x_0_zzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxx_0[i] = 3.0 * g_0_0_xx_1[i] * fi_acd_0 + g_0_0_xxx_1[i] * wa_x[i];

        g_x_0_xxy_0[i] = 2.0 * g_0_0_xy_1[i] * fi_acd_0 + g_0_0_xxy_1[i] * wa_x[i];

        g_x_0_xxz_0[i] = 2.0 * g_0_0_xz_1[i] * fi_acd_0 + g_0_0_xxz_1[i] * wa_x[i];

        g_x_0_xyy_0[i] = g_0_0_yy_1[i] * fi_acd_0 + g_0_0_xyy_1[i] * wa_x[i];

        g_x_0_xyz_0[i] = g_0_0_yz_1[i] * fi_acd_0 + g_0_0_xyz_1[i] * wa_x[i];

        g_x_0_xzz_0[i] = g_0_0_zz_1[i] * fi_acd_0 + g_0_0_xzz_1[i] * wa_x[i];

        g_x_0_yyy_0[i] = g_0_0_yyy_1[i] * wa_x[i];

        g_x_0_yyz_0[i] = g_0_0_yyz_1[i] * wa_x[i];

        g_x_0_yzz_0[i] = g_0_0_yzz_1[i] * wa_x[i];

        g_x_0_zzz_0[i] = g_0_0_zzz_1[i] * wa_x[i];
    }

    /// Set up 10-20 components of targeted buffer : PSF

    auto g_y_0_xxx_0 = pbuffer.data(idx_eri_0_psf + 10);

    auto g_y_0_xxy_0 = pbuffer.data(idx_eri_0_psf + 11);

    auto g_y_0_xxz_0 = pbuffer.data(idx_eri_0_psf + 12);

    auto g_y_0_xyy_0 = pbuffer.data(idx_eri_0_psf + 13);

    auto g_y_0_xyz_0 = pbuffer.data(idx_eri_0_psf + 14);

    auto g_y_0_xzz_0 = pbuffer.data(idx_eri_0_psf + 15);

    auto g_y_0_yyy_0 = pbuffer.data(idx_eri_0_psf + 16);

    auto g_y_0_yyz_0 = pbuffer.data(idx_eri_0_psf + 17);

    auto g_y_0_yzz_0 = pbuffer.data(idx_eri_0_psf + 18);

    auto g_y_0_zzz_0 = pbuffer.data(idx_eri_0_psf + 19);

    #pragma omp simd aligned(g_0_0_xx_1, g_0_0_xxx_1, g_0_0_xxy_1, g_0_0_xxz_1, g_0_0_xy_1, g_0_0_xyy_1, g_0_0_xyz_1, g_0_0_xz_1, g_0_0_xzz_1, g_0_0_yy_1, g_0_0_yyy_1, g_0_0_yyz_1, g_0_0_yz_1, g_0_0_yzz_1, g_0_0_zz_1, g_0_0_zzz_1, g_y_0_xxx_0, g_y_0_xxy_0, g_y_0_xxz_0, g_y_0_xyy_0, g_y_0_xyz_0, g_y_0_xzz_0, g_y_0_yyy_0, g_y_0_yyz_0, g_y_0_yzz_0, g_y_0_zzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxx_0[i] = g_0_0_xxx_1[i] * wa_y[i];

        g_y_0_xxy_0[i] = g_0_0_xx_1[i] * fi_acd_0 + g_0_0_xxy_1[i] * wa_y[i];

        g_y_0_xxz_0[i] = g_0_0_xxz_1[i] * wa_y[i];

        g_y_0_xyy_0[i] = 2.0 * g_0_0_xy_1[i] * fi_acd_0 + g_0_0_xyy_1[i] * wa_y[i];

        g_y_0_xyz_0[i] = g_0_0_xz_1[i] * fi_acd_0 + g_0_0_xyz_1[i] * wa_y[i];

        g_y_0_xzz_0[i] = g_0_0_xzz_1[i] * wa_y[i];

        g_y_0_yyy_0[i] = 3.0 * g_0_0_yy_1[i] * fi_acd_0 + g_0_0_yyy_1[i] * wa_y[i];

        g_y_0_yyz_0[i] = 2.0 * g_0_0_yz_1[i] * fi_acd_0 + g_0_0_yyz_1[i] * wa_y[i];

        g_y_0_yzz_0[i] = g_0_0_zz_1[i] * fi_acd_0 + g_0_0_yzz_1[i] * wa_y[i];

        g_y_0_zzz_0[i] = g_0_0_zzz_1[i] * wa_y[i];
    }

    /// Set up 20-30 components of targeted buffer : PSF

    auto g_z_0_xxx_0 = pbuffer.data(idx_eri_0_psf + 20);

    auto g_z_0_xxy_0 = pbuffer.data(idx_eri_0_psf + 21);

    auto g_z_0_xxz_0 = pbuffer.data(idx_eri_0_psf + 22);

    auto g_z_0_xyy_0 = pbuffer.data(idx_eri_0_psf + 23);

    auto g_z_0_xyz_0 = pbuffer.data(idx_eri_0_psf + 24);

    auto g_z_0_xzz_0 = pbuffer.data(idx_eri_0_psf + 25);

    auto g_z_0_yyy_0 = pbuffer.data(idx_eri_0_psf + 26);

    auto g_z_0_yyz_0 = pbuffer.data(idx_eri_0_psf + 27);

    auto g_z_0_yzz_0 = pbuffer.data(idx_eri_0_psf + 28);

    auto g_z_0_zzz_0 = pbuffer.data(idx_eri_0_psf + 29);

    #pragma omp simd aligned(g_0_0_xx_1, g_0_0_xxx_1, g_0_0_xxy_1, g_0_0_xxz_1, g_0_0_xy_1, g_0_0_xyy_1, g_0_0_xyz_1, g_0_0_xz_1, g_0_0_xzz_1, g_0_0_yy_1, g_0_0_yyy_1, g_0_0_yyz_1, g_0_0_yz_1, g_0_0_yzz_1, g_0_0_zz_1, g_0_0_zzz_1, g_z_0_xxx_0, g_z_0_xxy_0, g_z_0_xxz_0, g_z_0_xyy_0, g_z_0_xyz_0, g_z_0_xzz_0, g_z_0_yyy_0, g_z_0_yyz_0, g_z_0_yzz_0, g_z_0_zzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxx_0[i] = g_0_0_xxx_1[i] * wa_z[i];

        g_z_0_xxy_0[i] = g_0_0_xxy_1[i] * wa_z[i];

        g_z_0_xxz_0[i] = g_0_0_xx_1[i] * fi_acd_0 + g_0_0_xxz_1[i] * wa_z[i];

        g_z_0_xyy_0[i] = g_0_0_xyy_1[i] * wa_z[i];

        g_z_0_xyz_0[i] = g_0_0_xy_1[i] * fi_acd_0 + g_0_0_xyz_1[i] * wa_z[i];

        g_z_0_xzz_0[i] = 2.0 * g_0_0_xz_1[i] * fi_acd_0 + g_0_0_xzz_1[i] * wa_z[i];

        g_z_0_yyy_0[i] = g_0_0_yyy_1[i] * wa_z[i];

        g_z_0_yyz_0[i] = g_0_0_yy_1[i] * fi_acd_0 + g_0_0_yyz_1[i] * wa_z[i];

        g_z_0_yzz_0[i] = 2.0 * g_0_0_yz_1[i] * fi_acd_0 + g_0_0_yzz_1[i] * wa_z[i];

        g_z_0_zzz_0[i] = 3.0 * g_0_0_zz_1[i] * fi_acd_0 + g_0_0_zzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

