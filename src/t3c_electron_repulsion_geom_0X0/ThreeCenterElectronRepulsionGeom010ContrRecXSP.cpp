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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xsp(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xsp,
                                        const size_t idx_xsp,
                                        const size_t idx_xsd,
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

        /// Set up components of auxilary buffer : SSD

        const auto sd_off = idx_xsd + i * 6;

        auto g_0_xx = cbuffer.data(sd_off + 0);

        auto g_0_xy = cbuffer.data(sd_off + 1);

        auto g_0_xz = cbuffer.data(sd_off + 2);

        auto g_0_yy = cbuffer.data(sd_off + 3);

        auto g_0_yz = cbuffer.data(sd_off + 4);

        auto g_0_zz = cbuffer.data(sd_off + 5);

        /// set up bra offset for contr_buffer_xxsp

        const auto sp_geom_10_off = idx_geom_10_xsp + i * 3;

        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_x_0_0_x = cbuffer.data(sp_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_0_y = cbuffer.data(sp_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_0_z = cbuffer.data(sp_geom_10_off + 0 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_0_x, g_0_xx, g_0_xy, g_0_xz, g_0_y, g_0_z, g_x_0_0_x, g_x_0_0_y, g_x_0_0_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_0_x[k] = -g_0_x[k] * cd_x[k] + g_0_xx[k];

            g_x_0_0_y[k] = -g_0_y[k] * cd_x[k] + g_0_xy[k];

            g_x_0_0_z[k] = -g_0_z[k] * cd_x[k] + g_0_xz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_y_0_0_x = cbuffer.data(sp_geom_10_off + 3 * acomps  + 0);

        auto g_y_0_0_y = cbuffer.data(sp_geom_10_off + 3 * acomps  + 1);

        auto g_y_0_0_z = cbuffer.data(sp_geom_10_off + 3 * acomps  + 2);

        #pragma omp simd aligned(cd_y, g_0_x, g_0_xy, g_0_y, g_0_yy, g_0_yz, g_0_z, g_y_0_0_x, g_y_0_0_y, g_y_0_0_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_0_x[k] = -g_0_x[k] * cd_y[k] + g_0_xy[k];

            g_y_0_0_y[k] = -g_0_y[k] * cd_y[k] + g_0_yy[k];

            g_y_0_0_z[k] = -g_0_z[k] * cd_y[k] + g_0_yz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_z_0_0_x = cbuffer.data(sp_geom_10_off + 6 * acomps  + 0);

        auto g_z_0_0_y = cbuffer.data(sp_geom_10_off + 6 * acomps  + 1);

        auto g_z_0_0_z = cbuffer.data(sp_geom_10_off + 6 * acomps  + 2);

        #pragma omp simd aligned(cd_z, g_0_x, g_0_xz, g_0_y, g_0_yz, g_0_z, g_0_zz, g_z_0_0_x, g_z_0_0_y, g_z_0_0_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_0_x[k] = -g_0_x[k] * cd_z[k] + g_0_xz[k];

            g_z_0_0_y[k] = -g_0_y[k] * cd_z[k] + g_0_yz[k];

            g_z_0_0_z[k] = -g_0_z[k] * cd_z[k] + g_0_zz[k];
        }
    }
}

} // t3ceri namespace

