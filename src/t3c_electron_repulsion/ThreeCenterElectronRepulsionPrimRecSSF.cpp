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

#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssf(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssf,
                                 size_t idx_eri_0_ssp,
                                 size_t idx_eri_1_ssp,
                                 size_t idx_eri_0_ssd,
                                 size_t idx_eri_1_ssd,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSP

    auto g_0_0_x_0 = pbuffer.data(idx_eri_0_ssp);

    auto g_0_0_y_0 = pbuffer.data(idx_eri_0_ssp + 1);

    auto g_0_0_z_0 = pbuffer.data(idx_eri_0_ssp + 2);

    /// Set up components of auxilary buffer : SSP

    auto g_0_0_x_1 = pbuffer.data(idx_eri_1_ssp);

    auto g_0_0_y_1 = pbuffer.data(idx_eri_1_ssp + 1);

    auto g_0_0_z_1 = pbuffer.data(idx_eri_1_ssp + 2);

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_0 = pbuffer.data(idx_eri_0_ssd);

    auto g_0_0_yy_0 = pbuffer.data(idx_eri_0_ssd + 3);

    auto g_0_0_yz_0 = pbuffer.data(idx_eri_0_ssd + 4);

    auto g_0_0_zz_0 = pbuffer.data(idx_eri_0_ssd + 5);

    /// Set up components of auxilary buffer : SSD

    auto g_0_0_xx_1 = pbuffer.data(idx_eri_1_ssd);

    auto g_0_0_yy_1 = pbuffer.data(idx_eri_1_ssd + 3);

    auto g_0_0_yz_1 = pbuffer.data(idx_eri_1_ssd + 4);

    auto g_0_0_zz_1 = pbuffer.data(idx_eri_1_ssd + 5);

    /// Set up components of targeted buffer : SSF

    auto g_0_0_xxx_0 = pbuffer.data(idx_eri_0_ssf);

    auto g_0_0_xxy_0 = pbuffer.data(idx_eri_0_ssf + 1);

    auto g_0_0_xxz_0 = pbuffer.data(idx_eri_0_ssf + 2);

    auto g_0_0_xyy_0 = pbuffer.data(idx_eri_0_ssf + 3);

    auto g_0_0_xyz_0 = pbuffer.data(idx_eri_0_ssf + 4);

    auto g_0_0_xzz_0 = pbuffer.data(idx_eri_0_ssf + 5);

    auto g_0_0_yyy_0 = pbuffer.data(idx_eri_0_ssf + 6);

    auto g_0_0_yyz_0 = pbuffer.data(idx_eri_0_ssf + 7);

    auto g_0_0_yzz_0 = pbuffer.data(idx_eri_0_ssf + 8);

    auto g_0_0_zzz_0 = pbuffer.data(idx_eri_0_ssf + 9);

    #pragma omp simd aligned(g_0_0_x_0, g_0_0_x_1, g_0_0_xx_0, g_0_0_xx_1, g_0_0_xxx_0, g_0_0_xxy_0, g_0_0_xxz_0, g_0_0_xyy_0, g_0_0_xyz_0, g_0_0_xzz_0, g_0_0_y_0, g_0_0_y_1, g_0_0_yy_0, g_0_0_yy_1, g_0_0_yyy_0, g_0_0_yyz_0, g_0_0_yz_0, g_0_0_yz_1, g_0_0_yzz_0, g_0_0_z_0, g_0_0_z_1, g_0_0_zz_0, g_0_0_zz_1, g_0_0_zzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxx_0[i] = 2.0 * g_0_0_x_0[i] * fi_cd_0 - 2.0 * g_0_0_x_1[i] * fzi_cd_0 + g_0_0_xx_0[i] * qd_x[i] + g_0_0_xx_1[i] * wq_x[i];

        g_0_0_xxy_0[i] = g_0_0_xx_0[i] * qd_y[i] + g_0_0_xx_1[i] * wq_y[i];

        g_0_0_xxz_0[i] = g_0_0_xx_0[i] * qd_z[i] + g_0_0_xx_1[i] * wq_z[i];

        g_0_0_xyy_0[i] = g_0_0_yy_0[i] * qd_x[i] + g_0_0_yy_1[i] * wq_x[i];

        g_0_0_xyz_0[i] = g_0_0_yz_0[i] * qd_x[i] + g_0_0_yz_1[i] * wq_x[i];

        g_0_0_xzz_0[i] = g_0_0_zz_0[i] * qd_x[i] + g_0_0_zz_1[i] * wq_x[i];

        g_0_0_yyy_0[i] = 2.0 * g_0_0_y_0[i] * fi_cd_0 - 2.0 * g_0_0_y_1[i] * fzi_cd_0 + g_0_0_yy_0[i] * qd_y[i] + g_0_0_yy_1[i] * wq_y[i];

        g_0_0_yyz_0[i] = g_0_0_yy_0[i] * qd_z[i] + g_0_0_yy_1[i] * wq_z[i];

        g_0_0_yzz_0[i] = g_0_0_zz_0[i] * qd_y[i] + g_0_0_zz_1[i] * wq_y[i];

        g_0_0_zzz_0[i] = 2.0 * g_0_0_z_0[i] * fi_cd_0 - 2.0 * g_0_0_z_1[i] * fzi_cd_0 + g_0_0_zz_0[i] * qd_z[i] + g_0_0_zz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

