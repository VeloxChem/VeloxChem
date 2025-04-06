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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xsg(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xsg,
                                        const size_t idx_xsg,
                                        const size_t idx_xsh,
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
        /// Set up components of auxilary buffer : SSG

        const auto sg_off = idx_xsg + i * 15;

        auto g_0_xxxx = cbuffer.data(sg_off + 0);

        auto g_0_xxxy = cbuffer.data(sg_off + 1);

        auto g_0_xxxz = cbuffer.data(sg_off + 2);

        auto g_0_xxyy = cbuffer.data(sg_off + 3);

        auto g_0_xxyz = cbuffer.data(sg_off + 4);

        auto g_0_xxzz = cbuffer.data(sg_off + 5);

        auto g_0_xyyy = cbuffer.data(sg_off + 6);

        auto g_0_xyyz = cbuffer.data(sg_off + 7);

        auto g_0_xyzz = cbuffer.data(sg_off + 8);

        auto g_0_xzzz = cbuffer.data(sg_off + 9);

        auto g_0_yyyy = cbuffer.data(sg_off + 10);

        auto g_0_yyyz = cbuffer.data(sg_off + 11);

        auto g_0_yyzz = cbuffer.data(sg_off + 12);

        auto g_0_yzzz = cbuffer.data(sg_off + 13);

        auto g_0_zzzz = cbuffer.data(sg_off + 14);

        /// Set up components of auxilary buffer : SSH

        const auto sh_off = idx_xsh + i * 21;

        auto g_0_xxxxx = cbuffer.data(sh_off + 0);

        auto g_0_xxxxy = cbuffer.data(sh_off + 1);

        auto g_0_xxxxz = cbuffer.data(sh_off + 2);

        auto g_0_xxxyy = cbuffer.data(sh_off + 3);

        auto g_0_xxxyz = cbuffer.data(sh_off + 4);

        auto g_0_xxxzz = cbuffer.data(sh_off + 5);

        auto g_0_xxyyy = cbuffer.data(sh_off + 6);

        auto g_0_xxyyz = cbuffer.data(sh_off + 7);

        auto g_0_xxyzz = cbuffer.data(sh_off + 8);

        auto g_0_xxzzz = cbuffer.data(sh_off + 9);

        auto g_0_xyyyy = cbuffer.data(sh_off + 10);

        auto g_0_xyyyz = cbuffer.data(sh_off + 11);

        auto g_0_xyyzz = cbuffer.data(sh_off + 12);

        auto g_0_xyzzz = cbuffer.data(sh_off + 13);

        auto g_0_xzzzz = cbuffer.data(sh_off + 14);

        auto g_0_yyyyy = cbuffer.data(sh_off + 15);

        auto g_0_yyyyz = cbuffer.data(sh_off + 16);

        auto g_0_yyyzz = cbuffer.data(sh_off + 17);

        auto g_0_yyzzz = cbuffer.data(sh_off + 18);

        auto g_0_yzzzz = cbuffer.data(sh_off + 19);

        auto g_0_zzzzz = cbuffer.data(sh_off + 20);

        /// set up bra offset for contr_buffer_xxsg

        const auto sg_geom_10_off = idx_geom_10_xsg + i * 15;

        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_0_xxxx = cbuffer.data(sg_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_0_xxxy = cbuffer.data(sg_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_0_xxxz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_0_xxyy = cbuffer.data(sg_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_0_xxyz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_0_xxzz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_0_xyyy = cbuffer.data(sg_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_0_xyyz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_0_xyzz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_0_xzzz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_0_yyyy = cbuffer.data(sg_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_0_yyyz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_0_yyzz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_0_yzzz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_0_zzzz = cbuffer.data(sg_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_0_xxxx, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_0_0_xxxx, g_x_0_0_xxxy, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyz, g_x_0_0_xyzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_0_xxxx[k] = -g_0_xxxx[k] * cd_x[k] + g_0_xxxxx[k];

            g_x_0_0_xxxy[k] = -g_0_xxxy[k] * cd_x[k] + g_0_xxxxy[k];

            g_x_0_0_xxxz[k] = -g_0_xxxz[k] * cd_x[k] + g_0_xxxxz[k];

            g_x_0_0_xxyy[k] = -g_0_xxyy[k] * cd_x[k] + g_0_xxxyy[k];

            g_x_0_0_xxyz[k] = -g_0_xxyz[k] * cd_x[k] + g_0_xxxyz[k];

            g_x_0_0_xxzz[k] = -g_0_xxzz[k] * cd_x[k] + g_0_xxxzz[k];

            g_x_0_0_xyyy[k] = -g_0_xyyy[k] * cd_x[k] + g_0_xxyyy[k];

            g_x_0_0_xyyz[k] = -g_0_xyyz[k] * cd_x[k] + g_0_xxyyz[k];

            g_x_0_0_xyzz[k] = -g_0_xyzz[k] * cd_x[k] + g_0_xxyzz[k];

            g_x_0_0_xzzz[k] = -g_0_xzzz[k] * cd_x[k] + g_0_xxzzz[k];

            g_x_0_0_yyyy[k] = -g_0_yyyy[k] * cd_x[k] + g_0_xyyyy[k];

            g_x_0_0_yyyz[k] = -g_0_yyyz[k] * cd_x[k] + g_0_xyyyz[k];

            g_x_0_0_yyzz[k] = -g_0_yyzz[k] * cd_x[k] + g_0_xyyzz[k];

            g_x_0_0_yzzz[k] = -g_0_yzzz[k] * cd_x[k] + g_0_xyzzz[k];

            g_x_0_0_zzzz[k] = -g_0_zzzz[k] * cd_x[k] + g_0_xzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_0_xxxx = cbuffer.data(sg_geom_10_off + 15 * acomps  + 0);

        auto g_y_0_0_xxxy = cbuffer.data(sg_geom_10_off + 15 * acomps  + 1);

        auto g_y_0_0_xxxz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 2);

        auto g_y_0_0_xxyy = cbuffer.data(sg_geom_10_off + 15 * acomps  + 3);

        auto g_y_0_0_xxyz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 4);

        auto g_y_0_0_xxzz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 5);

        auto g_y_0_0_xyyy = cbuffer.data(sg_geom_10_off + 15 * acomps  + 6);

        auto g_y_0_0_xyyz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 7);

        auto g_y_0_0_xyzz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 8);

        auto g_y_0_0_xzzz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 9);

        auto g_y_0_0_yyyy = cbuffer.data(sg_geom_10_off + 15 * acomps  + 10);

        auto g_y_0_0_yyyz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 11);

        auto g_y_0_0_yyzz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 12);

        auto g_y_0_0_yzzz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 13);

        auto g_y_0_0_zzzz = cbuffer.data(sg_geom_10_off + 15 * acomps  + 14);

        #pragma omp simd aligned(cd_y, g_0_xxxx, g_0_xxxxy, g_0_xxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxxz, g_0_xxyy, g_0_xxyyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xyyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_yyyy, g_0_yyyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_y_0_0_xxxx, g_y_0_0_xxxy, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyz, g_y_0_0_xyzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_0_xxxx[k] = -g_0_xxxx[k] * cd_y[k] + g_0_xxxxy[k];

            g_y_0_0_xxxy[k] = -g_0_xxxy[k] * cd_y[k] + g_0_xxxyy[k];

            g_y_0_0_xxxz[k] = -g_0_xxxz[k] * cd_y[k] + g_0_xxxyz[k];

            g_y_0_0_xxyy[k] = -g_0_xxyy[k] * cd_y[k] + g_0_xxyyy[k];

            g_y_0_0_xxyz[k] = -g_0_xxyz[k] * cd_y[k] + g_0_xxyyz[k];

            g_y_0_0_xxzz[k] = -g_0_xxzz[k] * cd_y[k] + g_0_xxyzz[k];

            g_y_0_0_xyyy[k] = -g_0_xyyy[k] * cd_y[k] + g_0_xyyyy[k];

            g_y_0_0_xyyz[k] = -g_0_xyyz[k] * cd_y[k] + g_0_xyyyz[k];

            g_y_0_0_xyzz[k] = -g_0_xyzz[k] * cd_y[k] + g_0_xyyzz[k];

            g_y_0_0_xzzz[k] = -g_0_xzzz[k] * cd_y[k] + g_0_xyzzz[k];

            g_y_0_0_yyyy[k] = -g_0_yyyy[k] * cd_y[k] + g_0_yyyyy[k];

            g_y_0_0_yyyz[k] = -g_0_yyyz[k] * cd_y[k] + g_0_yyyyz[k];

            g_y_0_0_yyzz[k] = -g_0_yyzz[k] * cd_y[k] + g_0_yyyzz[k];

            g_y_0_0_yzzz[k] = -g_0_yzzz[k] * cd_y[k] + g_0_yyzzz[k];

            g_y_0_0_zzzz[k] = -g_0_zzzz[k] * cd_y[k] + g_0_yzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_0_xxxx = cbuffer.data(sg_geom_10_off + 30 * acomps  + 0);

        auto g_z_0_0_xxxy = cbuffer.data(sg_geom_10_off + 30 * acomps  + 1);

        auto g_z_0_0_xxxz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 2);

        auto g_z_0_0_xxyy = cbuffer.data(sg_geom_10_off + 30 * acomps  + 3);

        auto g_z_0_0_xxyz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 4);

        auto g_z_0_0_xxzz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 5);

        auto g_z_0_0_xyyy = cbuffer.data(sg_geom_10_off + 30 * acomps  + 6);

        auto g_z_0_0_xyyz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 7);

        auto g_z_0_0_xyzz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 8);

        auto g_z_0_0_xzzz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 9);

        auto g_z_0_0_yyyy = cbuffer.data(sg_geom_10_off + 30 * acomps  + 10);

        auto g_z_0_0_yyyz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 11);

        auto g_z_0_0_yyzz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 12);

        auto g_z_0_0_yzzz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 13);

        auto g_z_0_0_zzzz = cbuffer.data(sg_geom_10_off + 30 * acomps  + 14);

        #pragma omp simd aligned(cd_z, g_0_xxxx, g_0_xxxxz, g_0_xxxy, g_0_xxxyz, g_0_xxxz, g_0_xxxzz, g_0_xxyy, g_0_xxyyz, g_0_xxyz, g_0_xxyzz, g_0_xxzz, g_0_xxzzz, g_0_xyyy, g_0_xyyyz, g_0_xyyz, g_0_xyyzz, g_0_xyzz, g_0_xyzzz, g_0_xzzz, g_0_xzzzz, g_0_yyyy, g_0_yyyyz, g_0_yyyz, g_0_yyyzz, g_0_yyzz, g_0_yyzzz, g_0_yzzz, g_0_yzzzz, g_0_zzzz, g_0_zzzzz, g_z_0_0_xxxx, g_z_0_0_xxxy, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyz, g_z_0_0_xyzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_0_xxxx[k] = -g_0_xxxx[k] * cd_z[k] + g_0_xxxxz[k];

            g_z_0_0_xxxy[k] = -g_0_xxxy[k] * cd_z[k] + g_0_xxxyz[k];

            g_z_0_0_xxxz[k] = -g_0_xxxz[k] * cd_z[k] + g_0_xxxzz[k];

            g_z_0_0_xxyy[k] = -g_0_xxyy[k] * cd_z[k] + g_0_xxyyz[k];

            g_z_0_0_xxyz[k] = -g_0_xxyz[k] * cd_z[k] + g_0_xxyzz[k];

            g_z_0_0_xxzz[k] = -g_0_xxzz[k] * cd_z[k] + g_0_xxzzz[k];

            g_z_0_0_xyyy[k] = -g_0_xyyy[k] * cd_z[k] + g_0_xyyyz[k];

            g_z_0_0_xyyz[k] = -g_0_xyyz[k] * cd_z[k] + g_0_xyyzz[k];

            g_z_0_0_xyzz[k] = -g_0_xyzz[k] * cd_z[k] + g_0_xyzzz[k];

            g_z_0_0_xzzz[k] = -g_0_xzzz[k] * cd_z[k] + g_0_xzzzz[k];

            g_z_0_0_yyyy[k] = -g_0_yyyy[k] * cd_z[k] + g_0_yyyyz[k];

            g_z_0_0_yyyz[k] = -g_0_yyyz[k] * cd_z[k] + g_0_yyyzz[k];

            g_z_0_0_yyzz[k] = -g_0_yyzz[k] * cd_z[k] + g_0_yyzzz[k];

            g_z_0_0_yzzz[k] = -g_0_yzzz[k] * cd_z[k] + g_0_yzzzz[k];

            g_z_0_0_zzzz[k] = -g_0_zzzz[k] * cd_z[k] + g_0_zzzzz[k];
        }
    }
}

} // t3ceri namespace

