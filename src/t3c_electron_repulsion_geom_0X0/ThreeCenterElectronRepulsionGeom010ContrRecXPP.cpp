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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xpp(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xpp,
                                        const size_t idx_xsp,
                                        const size_t idx_geom_10_xsp,
                                        const size_t idx_geom_10_xsd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SSP

        const auto sp_off = idx_xsp + i * 3;

        auto g_0_x = cbuffer.data(sp_off + 0);

        auto g_0_y = cbuffer.data(sp_off + 1);

        auto g_0_z = cbuffer.data(sp_off + 2);

        /// Set up components of auxilary buffer : SSP

        const auto sp_geom_10_off = idx_geom_10_xsp + i * 3;

        auto g_x_0_0_x = cbuffer.data(sp_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_y = cbuffer.data(sp_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_z = cbuffer.data(sp_geom_10_off + 0 * acomps + 2);

        auto g_y_0_0_x = cbuffer.data(sp_geom_10_off + 3 * acomps + 0);

        auto g_y_0_0_y = cbuffer.data(sp_geom_10_off + 3 * acomps + 1);

        auto g_y_0_0_z = cbuffer.data(sp_geom_10_off + 3 * acomps + 2);

        auto g_z_0_0_x = cbuffer.data(sp_geom_10_off + 6 * acomps + 0);

        auto g_z_0_0_y = cbuffer.data(sp_geom_10_off + 6 * acomps + 1);

        auto g_z_0_0_z = cbuffer.data(sp_geom_10_off + 6 * acomps + 2);

        /// Set up components of auxilary buffer : SSD

        const auto sd_geom_10_off = idx_geom_10_xsd + i * 6;

        auto g_x_0_0_xx = cbuffer.data(sd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xy = cbuffer.data(sd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xz = cbuffer.data(sd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_yy = cbuffer.data(sd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_yz = cbuffer.data(sd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_zz = cbuffer.data(sd_geom_10_off + 0 * acomps + 5);

        auto g_y_0_0_xx = cbuffer.data(sd_geom_10_off + 6 * acomps + 0);

        auto g_y_0_0_xy = cbuffer.data(sd_geom_10_off + 6 * acomps + 1);

        auto g_y_0_0_xz = cbuffer.data(sd_geom_10_off + 6 * acomps + 2);

        auto g_y_0_0_yy = cbuffer.data(sd_geom_10_off + 6 * acomps + 3);

        auto g_y_0_0_yz = cbuffer.data(sd_geom_10_off + 6 * acomps + 4);

        auto g_y_0_0_zz = cbuffer.data(sd_geom_10_off + 6 * acomps + 5);

        auto g_z_0_0_xx = cbuffer.data(sd_geom_10_off + 12 * acomps + 0);

        auto g_z_0_0_xy = cbuffer.data(sd_geom_10_off + 12 * acomps + 1);

        auto g_z_0_0_xz = cbuffer.data(sd_geom_10_off + 12 * acomps + 2);

        auto g_z_0_0_yy = cbuffer.data(sd_geom_10_off + 12 * acomps + 3);

        auto g_z_0_0_yz = cbuffer.data(sd_geom_10_off + 12 * acomps + 4);

        auto g_z_0_0_zz = cbuffer.data(sd_geom_10_off + 12 * acomps + 5);

        /// set up bra offset for contr_buffer_xxpp

        const auto pp_geom_10_off = idx_geom_10_xpp + i * 9;

        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_x_0_x_x = cbuffer.data(pp_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_x_y = cbuffer.data(pp_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_x_z = cbuffer.data(pp_geom_10_off + 0 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_0_x, g_0_y, g_0_z, g_x_0_0_x, g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_y, g_x_0_0_z, g_x_0_x_x, g_x_0_x_y, g_x_0_x_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_x_x[k] = -g_0_x[k] - g_x_0_0_x[k] * cd_x[k] + g_x_0_0_xx[k];

            g_x_0_x_y[k] = -g_0_y[k] - g_x_0_0_y[k] * cd_x[k] + g_x_0_0_xy[k];

            g_x_0_x_z[k] = -g_0_z[k] - g_x_0_0_z[k] * cd_x[k] + g_x_0_0_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_y_x = cbuffer.data(pp_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_y_y = cbuffer.data(pp_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_y_z = cbuffer.data(pp_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_y, g_x_0_0_x, g_x_0_0_xy, g_x_0_0_y, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_z, g_x_0_y_x, g_x_0_y_y, g_x_0_y_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_y_x[k] = -g_x_0_0_x[k] * cd_y[k] + g_x_0_0_xy[k];

            g_x_0_y_y[k] = -g_x_0_0_y[k] * cd_y[k] + g_x_0_0_yy[k];

            g_x_0_y_z[k] = -g_x_0_0_z[k] * cd_y[k] + g_x_0_0_yz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_x_0_z_x = cbuffer.data(pp_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_z_y = cbuffer.data(pp_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_z_z = cbuffer.data(pp_geom_10_off + 0 * acomps  + 8);

        #pragma omp simd aligned(cd_z, g_x_0_0_x, g_x_0_0_xz, g_x_0_0_y, g_x_0_0_yz, g_x_0_0_z, g_x_0_0_zz, g_x_0_z_x, g_x_0_z_y, g_x_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_z_x[k] = -g_x_0_0_x[k] * cd_z[k] + g_x_0_0_xz[k];

            g_x_0_z_y[k] = -g_x_0_0_y[k] * cd_z[k] + g_x_0_0_yz[k];

            g_x_0_z_z[k] = -g_x_0_0_z[k] * cd_z[k] + g_x_0_0_zz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_y_0_x_x = cbuffer.data(pp_geom_10_off + 9 * acomps  + 0);

        auto g_y_0_x_y = cbuffer.data(pp_geom_10_off + 9 * acomps  + 1);

        auto g_y_0_x_z = cbuffer.data(pp_geom_10_off + 9 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_y_0_0_x, g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_y, g_y_0_0_z, g_y_0_x_x, g_y_0_x_y, g_y_0_x_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_x_x[k] = -g_y_0_0_x[k] * cd_x[k] + g_y_0_0_xx[k];

            g_y_0_x_y[k] = -g_y_0_0_y[k] * cd_x[k] + g_y_0_0_xy[k];

            g_y_0_x_z[k] = -g_y_0_0_z[k] * cd_x[k] + g_y_0_0_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_y_x = cbuffer.data(pp_geom_10_off + 9 * acomps  + 3);

        auto g_y_0_y_y = cbuffer.data(pp_geom_10_off + 9 * acomps  + 4);

        auto g_y_0_y_z = cbuffer.data(pp_geom_10_off + 9 * acomps  + 5);

        #pragma omp simd aligned(cd_y, g_0_x, g_0_y, g_0_z, g_y_0_0_x, g_y_0_0_xy, g_y_0_0_y, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_z, g_y_0_y_x, g_y_0_y_y, g_y_0_y_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_y_x[k] = -g_0_x[k] - g_y_0_0_x[k] * cd_y[k] + g_y_0_0_xy[k];

            g_y_0_y_y[k] = -g_0_y[k] - g_y_0_0_y[k] * cd_y[k] + g_y_0_0_yy[k];

            g_y_0_y_z[k] = -g_0_z[k] - g_y_0_0_z[k] * cd_y[k] + g_y_0_0_yz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_y_0_z_x = cbuffer.data(pp_geom_10_off + 9 * acomps  + 6);

        auto g_y_0_z_y = cbuffer.data(pp_geom_10_off + 9 * acomps  + 7);

        auto g_y_0_z_z = cbuffer.data(pp_geom_10_off + 9 * acomps  + 8);

        #pragma omp simd aligned(cd_z, g_y_0_0_x, g_y_0_0_xz, g_y_0_0_y, g_y_0_0_yz, g_y_0_0_z, g_y_0_0_zz, g_y_0_z_x, g_y_0_z_y, g_y_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_z_x[k] = -g_y_0_0_x[k] * cd_z[k] + g_y_0_0_xz[k];

            g_y_0_z_y[k] = -g_y_0_0_y[k] * cd_z[k] + g_y_0_0_yz[k];

            g_y_0_z_z[k] = -g_y_0_0_z[k] * cd_z[k] + g_y_0_0_zz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_z_0_x_x = cbuffer.data(pp_geom_10_off + 18 * acomps  + 0);

        auto g_z_0_x_y = cbuffer.data(pp_geom_10_off + 18 * acomps  + 1);

        auto g_z_0_x_z = cbuffer.data(pp_geom_10_off + 18 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_z_0_0_x, g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_y, g_z_0_0_z, g_z_0_x_x, g_z_0_x_y, g_z_0_x_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_x_x[k] = -g_z_0_0_x[k] * cd_x[k] + g_z_0_0_xx[k];

            g_z_0_x_y[k] = -g_z_0_0_y[k] * cd_x[k] + g_z_0_0_xy[k];

            g_z_0_x_z[k] = -g_z_0_0_z[k] * cd_x[k] + g_z_0_0_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_y_x = cbuffer.data(pp_geom_10_off + 18 * acomps  + 3);

        auto g_z_0_y_y = cbuffer.data(pp_geom_10_off + 18 * acomps  + 4);

        auto g_z_0_y_z = cbuffer.data(pp_geom_10_off + 18 * acomps  + 5);

        #pragma omp simd aligned(cd_y, g_z_0_0_x, g_z_0_0_xy, g_z_0_0_y, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_z, g_z_0_y_x, g_z_0_y_y, g_z_0_y_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_y_x[k] = -g_z_0_0_x[k] * cd_y[k] + g_z_0_0_xy[k];

            g_z_0_y_y[k] = -g_z_0_0_y[k] * cd_y[k] + g_z_0_0_yy[k];

            g_z_0_y_z[k] = -g_z_0_0_z[k] * cd_y[k] + g_z_0_0_yz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_z_0_z_x = cbuffer.data(pp_geom_10_off + 18 * acomps  + 6);

        auto g_z_0_z_y = cbuffer.data(pp_geom_10_off + 18 * acomps  + 7);

        auto g_z_0_z_z = cbuffer.data(pp_geom_10_off + 18 * acomps  + 8);

        #pragma omp simd aligned(cd_z, g_0_x, g_0_y, g_0_z, g_z_0_0_x, g_z_0_0_xz, g_z_0_0_y, g_z_0_0_yz, g_z_0_0_z, g_z_0_0_zz, g_z_0_z_x, g_z_0_z_y, g_z_0_z_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_z_x[k] = -g_0_x[k] - g_z_0_0_x[k] * cd_z[k] + g_z_0_0_xz[k];

            g_z_0_z_y[k] = -g_0_y[k] - g_z_0_0_y[k] * cd_z[k] + g_z_0_0_yz[k];

            g_z_0_z_z[k] = -g_0_z[k] - g_z_0_0_z[k] * cd_z[k] + g_z_0_0_zz[k];
        }
    }
}

} // t3ceri namespace

