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

#include "ElectronRepulsionContrRecDPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_dpxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_dpxx,
                                     const size_t idx_ppxx,
                                     const size_t idx_pdxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : PPSS

            const auto pp_off = idx_ppxx + i * dcomps + j;

            auto g_x_x = cbuffer.data(pp_off + 0 * ccomps * dcomps);

            auto g_x_y = cbuffer.data(pp_off + 1 * ccomps * dcomps);

            auto g_x_z = cbuffer.data(pp_off + 2 * ccomps * dcomps);

            auto g_y_x = cbuffer.data(pp_off + 3 * ccomps * dcomps);

            auto g_y_y = cbuffer.data(pp_off + 4 * ccomps * dcomps);

            auto g_y_z = cbuffer.data(pp_off + 5 * ccomps * dcomps);

            auto g_z_x = cbuffer.data(pp_off + 6 * ccomps * dcomps);

            auto g_z_y = cbuffer.data(pp_off + 7 * ccomps * dcomps);

            auto g_z_z = cbuffer.data(pp_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PDSS

            const auto pd_off = idx_pdxx + i * dcomps + j;

            auto g_x_xx = cbuffer.data(pd_off + 0 * ccomps * dcomps);

            auto g_x_xy = cbuffer.data(pd_off + 1 * ccomps * dcomps);

            auto g_x_xz = cbuffer.data(pd_off + 2 * ccomps * dcomps);

            auto g_y_xx = cbuffer.data(pd_off + 6 * ccomps * dcomps);

            auto g_y_xy = cbuffer.data(pd_off + 7 * ccomps * dcomps);

            auto g_y_xz = cbuffer.data(pd_off + 8 * ccomps * dcomps);

            auto g_y_yy = cbuffer.data(pd_off + 9 * ccomps * dcomps);

            auto g_y_yz = cbuffer.data(pd_off + 10 * ccomps * dcomps);

            auto g_z_xx = cbuffer.data(pd_off + 12 * ccomps * dcomps);

            auto g_z_xy = cbuffer.data(pd_off + 13 * ccomps * dcomps);

            auto g_z_xz = cbuffer.data(pd_off + 14 * ccomps * dcomps);

            auto g_z_yy = cbuffer.data(pd_off + 15 * ccomps * dcomps);

            auto g_z_yz = cbuffer.data(pd_off + 16 * ccomps * dcomps);

            auto g_z_zz = cbuffer.data(pd_off + 17 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dpxx

            const auto dp_off = idx_dpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_xx_x = cbuffer.data(dp_off + 0 * ccomps * dcomps);

            auto g_xx_y = cbuffer.data(dp_off + 1 * ccomps * dcomps);

            auto g_xx_z = cbuffer.data(dp_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x, g_x_xx, g_x_xy, g_x_xz, g_x_y, g_x_z, g_xx_x, g_xx_y, g_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_x[k] = -g_x_x[k] * ab_x + g_x_xx[k];

                g_xx_y[k] = -g_x_y[k] * ab_x + g_x_xy[k];

                g_xx_z[k] = -g_x_z[k] * ab_x + g_x_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_xy_x = cbuffer.data(dp_off + 3 * ccomps * dcomps);

            auto g_xy_y = cbuffer.data(dp_off + 4 * ccomps * dcomps);

            auto g_xy_z = cbuffer.data(dp_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_x, g_xy_y, g_xy_z, g_y_x, g_y_xx, g_y_xy, g_y_xz, g_y_y, g_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_x[k] = -g_y_x[k] * ab_x + g_y_xx[k];

                g_xy_y[k] = -g_y_y[k] * ab_x + g_y_xy[k];

                g_xy_z[k] = -g_y_z[k] * ab_x + g_y_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_xz_x = cbuffer.data(dp_off + 6 * ccomps * dcomps);

            auto g_xz_y = cbuffer.data(dp_off + 7 * ccomps * dcomps);

            auto g_xz_z = cbuffer.data(dp_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_x, g_xz_y, g_xz_z, g_z_x, g_z_xx, g_z_xy, g_z_xz, g_z_y, g_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_x[k] = -g_z_x[k] * ab_x + g_z_xx[k];

                g_xz_y[k] = -g_z_y[k] * ab_x + g_z_xy[k];

                g_xz_z[k] = -g_z_z[k] * ab_x + g_z_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_yy_x = cbuffer.data(dp_off + 9 * ccomps * dcomps);

            auto g_yy_y = cbuffer.data(dp_off + 10 * ccomps * dcomps);

            auto g_yy_z = cbuffer.data(dp_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x, g_y_xy, g_y_y, g_y_yy, g_y_yz, g_y_z, g_yy_x, g_yy_y, g_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_x[k] = -g_y_x[k] * ab_y + g_y_xy[k];

                g_yy_y[k] = -g_y_y[k] * ab_y + g_y_yy[k];

                g_yy_z[k] = -g_y_z[k] * ab_y + g_y_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_yz_x = cbuffer.data(dp_off + 12 * ccomps * dcomps);

            auto g_yz_y = cbuffer.data(dp_off + 13 * ccomps * dcomps);

            auto g_yz_z = cbuffer.data(dp_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_x, g_yz_y, g_yz_z, g_z_x, g_z_xy, g_z_y, g_z_yy, g_z_yz, g_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_x[k] = -g_z_x[k] * ab_y + g_z_xy[k];

                g_yz_y[k] = -g_z_y[k] * ab_y + g_z_yy[k];

                g_yz_z[k] = -g_z_z[k] * ab_y + g_z_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_zz_x = cbuffer.data(dp_off + 15 * ccomps * dcomps);

            auto g_zz_y = cbuffer.data(dp_off + 16 * ccomps * dcomps);

            auto g_zz_z = cbuffer.data(dp_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x, g_z_xz, g_z_y, g_z_yz, g_z_z, g_z_zz, g_zz_x, g_zz_y, g_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_x[k] = -g_z_x[k] * ab_z + g_z_xz[k];

                g_zz_y[k] = -g_z_y[k] * ab_z + g_z_yz[k];

                g_zz_z[k] = -g_z_z[k] * ab_z + g_z_zz[k];
            }
        }
    }
}

} // erirec namespace

