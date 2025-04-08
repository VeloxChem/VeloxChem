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

#include "ElectronRepulsionPrimRecSDSS.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sdss(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sdss,
                                  size_t                idx_eri_0_ssss,
                                  size_t                idx_eri_1_ssss,
                                  size_t                idx_eri_0_spss,
                                  size_t                idx_eri_1_spss,
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

    /// Set up components of auxilary buffer : SSSS

    auto g_0_0_0_0_0 = pbuffer.data(idx_eri_0_ssss);

    /// Set up components of auxilary buffer : SSSS

    auto g_0_0_0_0_1 = pbuffer.data(idx_eri_1_ssss);

    /// Set up components of auxilary buffer : SPSS

    auto g_0_x_0_0_0 = pbuffer.data(idx_eri_0_spss);

    auto g_0_y_0_0_0 = pbuffer.data(idx_eri_0_spss + 1);

    auto g_0_z_0_0_0 = pbuffer.data(idx_eri_0_spss + 2);

    /// Set up components of auxilary buffer : SPSS

    auto g_0_x_0_0_1 = pbuffer.data(idx_eri_1_spss);

    auto g_0_y_0_0_1 = pbuffer.data(idx_eri_1_spss + 1);

    auto g_0_z_0_0_1 = pbuffer.data(idx_eri_1_spss + 2);

    /// Set up components of targeted buffer : SDSS

    auto g_0_xx_0_0_0 = pbuffer.data(idx_eri_0_sdss);

    auto g_0_xy_0_0_0 = pbuffer.data(idx_eri_0_sdss + 1);

    auto g_0_xz_0_0_0 = pbuffer.data(idx_eri_0_sdss + 2);

    auto g_0_yy_0_0_0 = pbuffer.data(idx_eri_0_sdss + 3);

    auto g_0_yz_0_0_0 = pbuffer.data(idx_eri_0_sdss + 4);

    auto g_0_zz_0_0_0 = pbuffer.data(idx_eri_0_sdss + 5);

#pragma omp simd aligned(g_0_0_0_0_0,      \
                             g_0_0_0_0_1,  \
                             g_0_x_0_0_0,  \
                             g_0_x_0_0_1,  \
                             g_0_xx_0_0_0, \
                             g_0_xy_0_0_0, \
                             g_0_xz_0_0_0, \
                             g_0_y_0_0_0,  \
                             g_0_y_0_0_1,  \
                             g_0_yy_0_0_0, \
                             g_0_yz_0_0_0, \
                             g_0_z_0_0_0,  \
                             g_0_z_0_0_1,  \
                             g_0_zz_0_0_0, \
                             wp_x,         \
                             wp_y,         \
                             wp_z,         \
                             c_exps,       \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xx_0_0_0[i] = g_0_0_0_0_0[i] * fi_ab_0 - g_0_0_0_0_1[i] * fti_ab_0 + g_0_x_0_0_0[i] * pb_x + g_0_x_0_0_1[i] * wp_x[i];

        g_0_xy_0_0_0[i] = g_0_y_0_0_0[i] * pb_x + g_0_y_0_0_1[i] * wp_x[i];

        g_0_xz_0_0_0[i] = g_0_z_0_0_0[i] * pb_x + g_0_z_0_0_1[i] * wp_x[i];

        g_0_yy_0_0_0[i] = g_0_0_0_0_0[i] * fi_ab_0 - g_0_0_0_0_1[i] * fti_ab_0 + g_0_y_0_0_0[i] * pb_y + g_0_y_0_0_1[i] * wp_y[i];

        g_0_yz_0_0_0[i] = g_0_z_0_0_0[i] * pb_y + g_0_z_0_0_1[i] * wp_y[i];

        g_0_zz_0_0_0[i] = g_0_0_0_0_0[i] * fi_ab_0 - g_0_0_0_0_1[i] * fti_ab_0 + g_0_z_0_0_0[i] * pb_z + g_0_z_0_0_1[i] * wp_z[i];
    }
}

}  // namespace erirec
