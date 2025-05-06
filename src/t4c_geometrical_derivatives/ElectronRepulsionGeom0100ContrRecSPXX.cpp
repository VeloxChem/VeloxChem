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

#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_spxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_spxx,
                                            const size_t idx_ssxx,
                                            const size_t idx_sdxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSS

            const auto ss_off = idx_ssxx + i * dcomps + j;

            auto g_0_0 = cbuffer.data(ss_off + 0 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SDSS

            const auto sd_off = idx_sdxx + i * dcomps + j;

            auto g_0_xx = cbuffer.data(sd_off + 0 * ccomps * dcomps);

            auto g_0_xy = cbuffer.data(sd_off + 1 * ccomps * dcomps);

            auto g_0_xz = cbuffer.data(sd_off + 2 * ccomps * dcomps);

            auto g_0_yy = cbuffer.data(sd_off + 3 * ccomps * dcomps);

            auto g_0_yz = cbuffer.data(sd_off + 4 * ccomps * dcomps);

            auto g_0_zz = cbuffer.data(sd_off + 5 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_spxx

            const auto sp_geom_01_off = idx_geom_01_spxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_0_x = cbuffer.data(sp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_y = cbuffer.data(sp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_z = cbuffer.data(sp_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_0_xx, g_0_xy, g_0_xz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_x[k] = -g_0_0[k] + g_0_xx[k];

                g_0_x_0_y[k] = g_0_xy[k];

                g_0_x_0_z[k] = g_0_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_0_y_0_x = cbuffer.data(sp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_y_0_y = cbuffer.data(sp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_y_0_z = cbuffer.data(sp_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_xy, g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_0_yy, g_0_yz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_x[k] = g_0_xy[k];

                g_0_y_0_y[k] = -g_0_0[k] + g_0_yy[k];

                g_0_y_0_z[k] = g_0_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_0_z_0_x = cbuffer.data(sp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_z_0_y = cbuffer.data(sp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_z_0_z = cbuffer.data(sp_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_xz, g_0_yz, g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_x[k] = g_0_xz[k];

                g_0_z_0_y[k] = g_0_yz[k];

                g_0_z_0_z[k] = -g_0_0[k] + g_0_zz[k];
            }
        }
    }
}

} // erirec namespace

