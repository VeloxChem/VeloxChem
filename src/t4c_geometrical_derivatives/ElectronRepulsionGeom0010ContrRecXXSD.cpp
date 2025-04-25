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

#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxsd(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsd,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsd,
                                            const size_t idx_xxsf,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSD

            const auto sd_off = idx_xxsd + (i * bcomps + j) * 6;

            auto g_0_xx = pbuffer.data(sd_off + 0);

            auto g_0_xy = pbuffer.data(sd_off + 1);

            auto g_0_xz = pbuffer.data(sd_off + 2);

            auto g_0_yy = pbuffer.data(sd_off + 3);

            auto g_0_yz = pbuffer.data(sd_off + 4);

            auto g_0_zz = pbuffer.data(sd_off + 5);

            /// Set up components of auxilary buffer : SSSF

            const auto sf_off = idx_xxsf + (i * bcomps + j) * 10;

            auto g_0_xxx = pbuffer.data(sf_off + 0);

            auto g_0_xxy = pbuffer.data(sf_off + 1);

            auto g_0_xxz = pbuffer.data(sf_off + 2);

            auto g_0_xyy = pbuffer.data(sf_off + 3);

            auto g_0_xyz = pbuffer.data(sf_off + 4);

            auto g_0_xzz = pbuffer.data(sf_off + 5);

            auto g_0_yyy = pbuffer.data(sf_off + 6);

            auto g_0_yyz = pbuffer.data(sf_off + 7);

            auto g_0_yzz = pbuffer.data(sf_off + 8);

            auto g_0_zzz = pbuffer.data(sf_off + 9);

            /// set up bra offset for contr_buffer_xxsd

            const auto sd_geom_10_off = idx_geom_10_xxsd + (i * bcomps + j) * 6;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xx = cbuffer.data(sd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xy = cbuffer.data(sd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xz = cbuffer.data(sd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_yy = cbuffer.data(sd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_yz = cbuffer.data(sd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_zz = cbuffer.data(sd_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_0_xx, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yz, g_0_zz, g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xx[k] = -g_0_xx[k] * cd_x[k] + g_0_xxx[k];

                g_x_0_0_xy[k] = -g_0_xy[k] * cd_x[k] + g_0_xxy[k];

                g_x_0_0_xz[k] = -g_0_xz[k] * cd_x[k] + g_0_xxz[k];

                g_x_0_0_yy[k] = -g_0_yy[k] * cd_x[k] + g_0_xyy[k];

                g_x_0_0_yz[k] = -g_0_yz[k] * cd_x[k] + g_0_xyz[k];

                g_x_0_0_zz[k] = -g_0_zz[k] * cd_x[k] + g_0_xzz[k];
            }
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xx = cbuffer.data(sd_geom_10_off + 6 * acomps * bcomps + 0);

            auto g_y_0_0_xy = cbuffer.data(sd_geom_10_off + 6 * acomps * bcomps + 1);

            auto g_y_0_0_xz = cbuffer.data(sd_geom_10_off + 6 * acomps * bcomps + 2);

            auto g_y_0_0_yy = cbuffer.data(sd_geom_10_off + 6 * acomps * bcomps + 3);

            auto g_y_0_0_yz = cbuffer.data(sd_geom_10_off + 6 * acomps * bcomps + 4);

            auto g_y_0_0_zz = cbuffer.data(sd_geom_10_off + 6 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_y, g_0_xx, g_0_xxy, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_yy, g_0_yyy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xx[k] = -g_0_xx[k] * cd_y[k] + g_0_xxy[k];

                g_y_0_0_xy[k] = -g_0_xy[k] * cd_y[k] + g_0_xyy[k];

                g_y_0_0_xz[k] = -g_0_xz[k] * cd_y[k] + g_0_xyz[k];

                g_y_0_0_yy[k] = -g_0_yy[k] * cd_y[k] + g_0_yyy[k];

                g_y_0_0_yz[k] = -g_0_yz[k] * cd_y[k] + g_0_yyz[k];

                g_y_0_0_zz[k] = -g_0_zz[k] * cd_y[k] + g_0_yzz[k];
            }
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xx = cbuffer.data(sd_geom_10_off + 12 * acomps * bcomps + 0);

            auto g_z_0_0_xy = cbuffer.data(sd_geom_10_off + 12 * acomps * bcomps + 1);

            auto g_z_0_0_xz = cbuffer.data(sd_geom_10_off + 12 * acomps * bcomps + 2);

            auto g_z_0_0_yy = cbuffer.data(sd_geom_10_off + 12 * acomps * bcomps + 3);

            auto g_z_0_0_yz = cbuffer.data(sd_geom_10_off + 12 * acomps * bcomps + 4);

            auto g_z_0_0_zz = cbuffer.data(sd_geom_10_off + 12 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_0_xx, g_0_xxz, g_0_xy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_0_zzz, g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xx[k] = -g_0_xx[k] * cd_z[k] + g_0_xxz[k];

                g_z_0_0_xy[k] = -g_0_xy[k] * cd_z[k] + g_0_xyz[k];

                g_z_0_0_xz[k] = -g_0_xz[k] * cd_z[k] + g_0_xzz[k];

                g_z_0_0_yy[k] = -g_0_yy[k] * cd_z[k] + g_0_yyz[k];

                g_z_0_0_yz[k] = -g_0_yz[k] * cd_z[k] + g_0_yzz[k];

                g_z_0_0_zz[k] = -g_0_zz[k] * cd_z[k] + g_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

