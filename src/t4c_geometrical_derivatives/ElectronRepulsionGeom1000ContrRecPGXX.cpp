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

#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_pgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_pgxx,
                                            const size_t idx_sgxx,
                                            const size_t idx_geom_10_sgxx,
                                            const size_t idx_geom_10_shxx,
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
            /// Set up components of auxilary buffer : SGSS

            const auto sg_off = idx_sgxx + i * dcomps + j;

            auto g_0_xxxx = cbuffer.data(sg_off + 0 * ccomps * dcomps);

            auto g_0_xxxy = cbuffer.data(sg_off + 1 * ccomps * dcomps);

            auto g_0_xxxz = cbuffer.data(sg_off + 2 * ccomps * dcomps);

            auto g_0_xxyy = cbuffer.data(sg_off + 3 * ccomps * dcomps);

            auto g_0_xxyz = cbuffer.data(sg_off + 4 * ccomps * dcomps);

            auto g_0_xxzz = cbuffer.data(sg_off + 5 * ccomps * dcomps);

            auto g_0_xyyy = cbuffer.data(sg_off + 6 * ccomps * dcomps);

            auto g_0_xyyz = cbuffer.data(sg_off + 7 * ccomps * dcomps);

            auto g_0_xyzz = cbuffer.data(sg_off + 8 * ccomps * dcomps);

            auto g_0_xzzz = cbuffer.data(sg_off + 9 * ccomps * dcomps);

            auto g_0_yyyy = cbuffer.data(sg_off + 10 * ccomps * dcomps);

            auto g_0_yyyz = cbuffer.data(sg_off + 11 * ccomps * dcomps);

            auto g_0_yyzz = cbuffer.data(sg_off + 12 * ccomps * dcomps);

            auto g_0_yzzz = cbuffer.data(sg_off + 13 * ccomps * dcomps);

            auto g_0_zzzz = cbuffer.data(sg_off + 14 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_10_off = idx_geom_10_sgxx + i * dcomps + j;

            auto g_x_0_0_xxxx = cbuffer.data(sg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxy = cbuffer.data(sg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxz = cbuffer.data(sg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxyy = cbuffer.data(sg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxyz = cbuffer.data(sg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxzz = cbuffer.data(sg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xyyy = cbuffer.data(sg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xyyz = cbuffer.data(sg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xyzz = cbuffer.data(sg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xzzz = cbuffer.data(sg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_yyyy = cbuffer.data(sg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_yyyz = cbuffer.data(sg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_yyzz = cbuffer.data(sg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_yzzz = cbuffer.data(sg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_zzzz = cbuffer.data(sg_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_0_xxxx = cbuffer.data(sg_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_0_xxxy = cbuffer.data(sg_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_0_xxxz = cbuffer.data(sg_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_0_xxyy = cbuffer.data(sg_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_0_xxyz = cbuffer.data(sg_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_0_xxzz = cbuffer.data(sg_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xyyy = cbuffer.data(sg_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xyyz = cbuffer.data(sg_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xyzz = cbuffer.data(sg_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xzzz = cbuffer.data(sg_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_0_yyyy = cbuffer.data(sg_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_0_yyyz = cbuffer.data(sg_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_0_yyzz = cbuffer.data(sg_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_yzzz = cbuffer.data(sg_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_zzzz = cbuffer.data(sg_geom_10_off + 29 * ccomps * dcomps);

            auto g_z_0_0_xxxx = cbuffer.data(sg_geom_10_off + 30 * ccomps * dcomps);

            auto g_z_0_0_xxxy = cbuffer.data(sg_geom_10_off + 31 * ccomps * dcomps);

            auto g_z_0_0_xxxz = cbuffer.data(sg_geom_10_off + 32 * ccomps * dcomps);

            auto g_z_0_0_xxyy = cbuffer.data(sg_geom_10_off + 33 * ccomps * dcomps);

            auto g_z_0_0_xxyz = cbuffer.data(sg_geom_10_off + 34 * ccomps * dcomps);

            auto g_z_0_0_xxzz = cbuffer.data(sg_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_0_xyyy = cbuffer.data(sg_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_0_xyyz = cbuffer.data(sg_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_0_xyzz = cbuffer.data(sg_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_0_xzzz = cbuffer.data(sg_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_0_yyyy = cbuffer.data(sg_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_0_yyyz = cbuffer.data(sg_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_0_yyzz = cbuffer.data(sg_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_0_yzzz = cbuffer.data(sg_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_0_zzzz = cbuffer.data(sg_geom_10_off + 44 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_10_off = idx_geom_10_shxx + i * dcomps + j;

            auto g_x_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 53 * ccomps * dcomps);

            auto g_z_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 54 * ccomps * dcomps);

            auto g_z_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 55 * ccomps * dcomps);

            auto g_z_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 56 * ccomps * dcomps);

            auto g_z_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 57 * ccomps * dcomps);

            auto g_z_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 58 * ccomps * dcomps);

            auto g_z_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 62 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pgxx

            const auto pg_geom_10_off = idx_geom_10_pgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_xxxx = cbuffer.data(pg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxy = cbuffer.data(pg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxz = cbuffer.data(pg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxyy = cbuffer.data(pg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxyz = cbuffer.data(pg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxzz = cbuffer.data(pg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xyyy = cbuffer.data(pg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xyyz = cbuffer.data(pg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xyzz = cbuffer.data(pg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xzzz = cbuffer.data(pg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_yyyy = cbuffer.data(pg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_yyyz = cbuffer.data(pg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_yyzz = cbuffer.data(pg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_yzzz = cbuffer.data(pg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_zzzz = cbuffer.data(pg_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_0_0_xxxx, g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxy, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxz, g_x_0_0_xxxzz, g_x_0_0_xxyy, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyz, g_x_0_0_xxyzz, g_x_0_0_xxzz, g_x_0_0_xxzzz, g_x_0_0_xyyy, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyz, g_x_0_0_xyyzz, g_x_0_0_xyzz, g_x_0_0_xyzzz, g_x_0_0_xzzz, g_x_0_0_xzzzz, g_x_0_0_yyyy, g_x_0_0_yyyz, g_x_0_0_yyzz, g_x_0_0_yzzz, g_x_0_0_zzzz, g_x_0_x_xxxx, g_x_0_x_xxxy, g_x_0_x_xxxz, g_x_0_x_xxyy, g_x_0_x_xxyz, g_x_0_x_xxzz, g_x_0_x_xyyy, g_x_0_x_xyyz, g_x_0_x_xyzz, g_x_0_x_xzzz, g_x_0_x_yyyy, g_x_0_x_yyyz, g_x_0_x_yyzz, g_x_0_x_yzzz, g_x_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_xxxx[k] = -g_0_xxxx[k] - g_x_0_0_xxxx[k] * ab_x + g_x_0_0_xxxxx[k];

                g_x_0_x_xxxy[k] = -g_0_xxxy[k] - g_x_0_0_xxxy[k] * ab_x + g_x_0_0_xxxxy[k];

                g_x_0_x_xxxz[k] = -g_0_xxxz[k] - g_x_0_0_xxxz[k] * ab_x + g_x_0_0_xxxxz[k];

                g_x_0_x_xxyy[k] = -g_0_xxyy[k] - g_x_0_0_xxyy[k] * ab_x + g_x_0_0_xxxyy[k];

                g_x_0_x_xxyz[k] = -g_0_xxyz[k] - g_x_0_0_xxyz[k] * ab_x + g_x_0_0_xxxyz[k];

                g_x_0_x_xxzz[k] = -g_0_xxzz[k] - g_x_0_0_xxzz[k] * ab_x + g_x_0_0_xxxzz[k];

                g_x_0_x_xyyy[k] = -g_0_xyyy[k] - g_x_0_0_xyyy[k] * ab_x + g_x_0_0_xxyyy[k];

                g_x_0_x_xyyz[k] = -g_0_xyyz[k] - g_x_0_0_xyyz[k] * ab_x + g_x_0_0_xxyyz[k];

                g_x_0_x_xyzz[k] = -g_0_xyzz[k] - g_x_0_0_xyzz[k] * ab_x + g_x_0_0_xxyzz[k];

                g_x_0_x_xzzz[k] = -g_0_xzzz[k] - g_x_0_0_xzzz[k] * ab_x + g_x_0_0_xxzzz[k];

                g_x_0_x_yyyy[k] = -g_0_yyyy[k] - g_x_0_0_yyyy[k] * ab_x + g_x_0_0_xyyyy[k];

                g_x_0_x_yyyz[k] = -g_0_yyyz[k] - g_x_0_0_yyyz[k] * ab_x + g_x_0_0_xyyyz[k];

                g_x_0_x_yyzz[k] = -g_0_yyzz[k] - g_x_0_0_yyzz[k] * ab_x + g_x_0_0_xyyzz[k];

                g_x_0_x_yzzz[k] = -g_0_yzzz[k] - g_x_0_0_yzzz[k] * ab_x + g_x_0_0_xyzzz[k];

                g_x_0_x_zzzz[k] = -g_0_zzzz[k] - g_x_0_0_zzzz[k] * ab_x + g_x_0_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_xxxx = cbuffer.data(pg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_y_xxxy = cbuffer.data(pg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_y_xxxz = cbuffer.data(pg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_y_xxyy = cbuffer.data(pg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_y_xxyz = cbuffer.data(pg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_y_xxzz = cbuffer.data(pg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_y_xyyy = cbuffer.data(pg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_y_xyyz = cbuffer.data(pg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_y_xyzz = cbuffer.data(pg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_y_xzzz = cbuffer.data(pg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_y_yyyy = cbuffer.data(pg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_y_yyyz = cbuffer.data(pg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_y_yyzz = cbuffer.data(pg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_yzzz = cbuffer.data(pg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_y_zzzz = cbuffer.data(pg_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxx, g_x_0_0_xxxxy, g_x_0_0_xxxy, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxz, g_x_0_0_xxyy, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyz, g_x_0_0_xxyzz, g_x_0_0_xxzz, g_x_0_0_xyyy, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyz, g_x_0_0_xyyzz, g_x_0_0_xyzz, g_x_0_0_xyzzz, g_x_0_0_xzzz, g_x_0_0_yyyy, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyz, g_x_0_0_yyyzz, g_x_0_0_yyzz, g_x_0_0_yyzzz, g_x_0_0_yzzz, g_x_0_0_yzzzz, g_x_0_0_zzzz, g_x_0_y_xxxx, g_x_0_y_xxxy, g_x_0_y_xxxz, g_x_0_y_xxyy, g_x_0_y_xxyz, g_x_0_y_xxzz, g_x_0_y_xyyy, g_x_0_y_xyyz, g_x_0_y_xyzz, g_x_0_y_xzzz, g_x_0_y_yyyy, g_x_0_y_yyyz, g_x_0_y_yyzz, g_x_0_y_yzzz, g_x_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_xxxx[k] = -g_x_0_0_xxxx[k] * ab_y + g_x_0_0_xxxxy[k];

                g_x_0_y_xxxy[k] = -g_x_0_0_xxxy[k] * ab_y + g_x_0_0_xxxyy[k];

                g_x_0_y_xxxz[k] = -g_x_0_0_xxxz[k] * ab_y + g_x_0_0_xxxyz[k];

                g_x_0_y_xxyy[k] = -g_x_0_0_xxyy[k] * ab_y + g_x_0_0_xxyyy[k];

                g_x_0_y_xxyz[k] = -g_x_0_0_xxyz[k] * ab_y + g_x_0_0_xxyyz[k];

                g_x_0_y_xxzz[k] = -g_x_0_0_xxzz[k] * ab_y + g_x_0_0_xxyzz[k];

                g_x_0_y_xyyy[k] = -g_x_0_0_xyyy[k] * ab_y + g_x_0_0_xyyyy[k];

                g_x_0_y_xyyz[k] = -g_x_0_0_xyyz[k] * ab_y + g_x_0_0_xyyyz[k];

                g_x_0_y_xyzz[k] = -g_x_0_0_xyzz[k] * ab_y + g_x_0_0_xyyzz[k];

                g_x_0_y_xzzz[k] = -g_x_0_0_xzzz[k] * ab_y + g_x_0_0_xyzzz[k];

                g_x_0_y_yyyy[k] = -g_x_0_0_yyyy[k] * ab_y + g_x_0_0_yyyyy[k];

                g_x_0_y_yyyz[k] = -g_x_0_0_yyyz[k] * ab_y + g_x_0_0_yyyyz[k];

                g_x_0_y_yyzz[k] = -g_x_0_0_yyzz[k] * ab_y + g_x_0_0_yyyzz[k];

                g_x_0_y_yzzz[k] = -g_x_0_0_yzzz[k] * ab_y + g_x_0_0_yyzzz[k];

                g_x_0_y_zzzz[k] = -g_x_0_0_zzzz[k] * ab_y + g_x_0_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_xxxx = cbuffer.data(pg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_z_xxxy = cbuffer.data(pg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_z_xxxz = cbuffer.data(pg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_z_xxyy = cbuffer.data(pg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_z_xxyz = cbuffer.data(pg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_z_xxzz = cbuffer.data(pg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_z_xyyy = cbuffer.data(pg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_z_xyyz = cbuffer.data(pg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_z_xyzz = cbuffer.data(pg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_z_xzzz = cbuffer.data(pg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_z_yyyy = cbuffer.data(pg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_z_yyyz = cbuffer.data(pg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_z_yyzz = cbuffer.data(pg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_z_yzzz = cbuffer.data(pg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_z_zzzz = cbuffer.data(pg_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxx, g_x_0_0_xxxxz, g_x_0_0_xxxy, g_x_0_0_xxxyz, g_x_0_0_xxxz, g_x_0_0_xxxzz, g_x_0_0_xxyy, g_x_0_0_xxyyz, g_x_0_0_xxyz, g_x_0_0_xxyzz, g_x_0_0_xxzz, g_x_0_0_xxzzz, g_x_0_0_xyyy, g_x_0_0_xyyyz, g_x_0_0_xyyz, g_x_0_0_xyyzz, g_x_0_0_xyzz, g_x_0_0_xyzzz, g_x_0_0_xzzz, g_x_0_0_xzzzz, g_x_0_0_yyyy, g_x_0_0_yyyyz, g_x_0_0_yyyz, g_x_0_0_yyyzz, g_x_0_0_yyzz, g_x_0_0_yyzzz, g_x_0_0_yzzz, g_x_0_0_yzzzz, g_x_0_0_zzzz, g_x_0_0_zzzzz, g_x_0_z_xxxx, g_x_0_z_xxxy, g_x_0_z_xxxz, g_x_0_z_xxyy, g_x_0_z_xxyz, g_x_0_z_xxzz, g_x_0_z_xyyy, g_x_0_z_xyyz, g_x_0_z_xyzz, g_x_0_z_xzzz, g_x_0_z_yyyy, g_x_0_z_yyyz, g_x_0_z_yyzz, g_x_0_z_yzzz, g_x_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_xxxx[k] = -g_x_0_0_xxxx[k] * ab_z + g_x_0_0_xxxxz[k];

                g_x_0_z_xxxy[k] = -g_x_0_0_xxxy[k] * ab_z + g_x_0_0_xxxyz[k];

                g_x_0_z_xxxz[k] = -g_x_0_0_xxxz[k] * ab_z + g_x_0_0_xxxzz[k];

                g_x_0_z_xxyy[k] = -g_x_0_0_xxyy[k] * ab_z + g_x_0_0_xxyyz[k];

                g_x_0_z_xxyz[k] = -g_x_0_0_xxyz[k] * ab_z + g_x_0_0_xxyzz[k];

                g_x_0_z_xxzz[k] = -g_x_0_0_xxzz[k] * ab_z + g_x_0_0_xxzzz[k];

                g_x_0_z_xyyy[k] = -g_x_0_0_xyyy[k] * ab_z + g_x_0_0_xyyyz[k];

                g_x_0_z_xyyz[k] = -g_x_0_0_xyyz[k] * ab_z + g_x_0_0_xyyzz[k];

                g_x_0_z_xyzz[k] = -g_x_0_0_xyzz[k] * ab_z + g_x_0_0_xyzzz[k];

                g_x_0_z_xzzz[k] = -g_x_0_0_xzzz[k] * ab_z + g_x_0_0_xzzzz[k];

                g_x_0_z_yyyy[k] = -g_x_0_0_yyyy[k] * ab_z + g_x_0_0_yyyyz[k];

                g_x_0_z_yyyz[k] = -g_x_0_0_yyyz[k] * ab_z + g_x_0_0_yyyzz[k];

                g_x_0_z_yyzz[k] = -g_x_0_0_yyzz[k] * ab_z + g_x_0_0_yyzzz[k];

                g_x_0_z_yzzz[k] = -g_x_0_0_yzzz[k] * ab_z + g_x_0_0_yzzzz[k];

                g_x_0_z_zzzz[k] = -g_x_0_0_zzzz[k] * ab_z + g_x_0_0_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_xxxx = cbuffer.data(pg_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_x_xxxy = cbuffer.data(pg_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_x_xxxz = cbuffer.data(pg_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_x_xxyy = cbuffer.data(pg_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_x_xxyz = cbuffer.data(pg_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_x_xxzz = cbuffer.data(pg_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_x_xyyy = cbuffer.data(pg_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_x_xyyz = cbuffer.data(pg_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_x_xyzz = cbuffer.data(pg_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_x_xzzz = cbuffer.data(pg_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_x_yyyy = cbuffer.data(pg_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_x_yyyz = cbuffer.data(pg_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_x_yyzz = cbuffer.data(pg_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_x_yzzz = cbuffer.data(pg_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_x_zzzz = cbuffer.data(pg_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxx, g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxy, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxz, g_y_0_0_xxxzz, g_y_0_0_xxyy, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyz, g_y_0_0_xxyzz, g_y_0_0_xxzz, g_y_0_0_xxzzz, g_y_0_0_xyyy, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyz, g_y_0_0_xyyzz, g_y_0_0_xyzz, g_y_0_0_xyzzz, g_y_0_0_xzzz, g_y_0_0_xzzzz, g_y_0_0_yyyy, g_y_0_0_yyyz, g_y_0_0_yyzz, g_y_0_0_yzzz, g_y_0_0_zzzz, g_y_0_x_xxxx, g_y_0_x_xxxy, g_y_0_x_xxxz, g_y_0_x_xxyy, g_y_0_x_xxyz, g_y_0_x_xxzz, g_y_0_x_xyyy, g_y_0_x_xyyz, g_y_0_x_xyzz, g_y_0_x_xzzz, g_y_0_x_yyyy, g_y_0_x_yyyz, g_y_0_x_yyzz, g_y_0_x_yzzz, g_y_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_xxxx[k] = -g_y_0_0_xxxx[k] * ab_x + g_y_0_0_xxxxx[k];

                g_y_0_x_xxxy[k] = -g_y_0_0_xxxy[k] * ab_x + g_y_0_0_xxxxy[k];

                g_y_0_x_xxxz[k] = -g_y_0_0_xxxz[k] * ab_x + g_y_0_0_xxxxz[k];

                g_y_0_x_xxyy[k] = -g_y_0_0_xxyy[k] * ab_x + g_y_0_0_xxxyy[k];

                g_y_0_x_xxyz[k] = -g_y_0_0_xxyz[k] * ab_x + g_y_0_0_xxxyz[k];

                g_y_0_x_xxzz[k] = -g_y_0_0_xxzz[k] * ab_x + g_y_0_0_xxxzz[k];

                g_y_0_x_xyyy[k] = -g_y_0_0_xyyy[k] * ab_x + g_y_0_0_xxyyy[k];

                g_y_0_x_xyyz[k] = -g_y_0_0_xyyz[k] * ab_x + g_y_0_0_xxyyz[k];

                g_y_0_x_xyzz[k] = -g_y_0_0_xyzz[k] * ab_x + g_y_0_0_xxyzz[k];

                g_y_0_x_xzzz[k] = -g_y_0_0_xzzz[k] * ab_x + g_y_0_0_xxzzz[k];

                g_y_0_x_yyyy[k] = -g_y_0_0_yyyy[k] * ab_x + g_y_0_0_xyyyy[k];

                g_y_0_x_yyyz[k] = -g_y_0_0_yyyz[k] * ab_x + g_y_0_0_xyyyz[k];

                g_y_0_x_yyzz[k] = -g_y_0_0_yyzz[k] * ab_x + g_y_0_0_xyyzz[k];

                g_y_0_x_yzzz[k] = -g_y_0_0_yzzz[k] * ab_x + g_y_0_0_xyzzz[k];

                g_y_0_x_zzzz[k] = -g_y_0_0_zzzz[k] * ab_x + g_y_0_0_xzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_xxxx = cbuffer.data(pg_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_y_xxxy = cbuffer.data(pg_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_y_xxxz = cbuffer.data(pg_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_y_xxyy = cbuffer.data(pg_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_y_xxyz = cbuffer.data(pg_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_y_xxzz = cbuffer.data(pg_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_y_xyyy = cbuffer.data(pg_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_y_xyyz = cbuffer.data(pg_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_y_xyzz = cbuffer.data(pg_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_y_xzzz = cbuffer.data(pg_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_y_yyyy = cbuffer.data(pg_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_y_yyyz = cbuffer.data(pg_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_y_yyzz = cbuffer.data(pg_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_y_yzzz = cbuffer.data(pg_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_y_zzzz = cbuffer.data(pg_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_y_0_0_xxxx, g_y_0_0_xxxxy, g_y_0_0_xxxy, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxz, g_y_0_0_xxyy, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyz, g_y_0_0_xxyzz, g_y_0_0_xxzz, g_y_0_0_xyyy, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyz, g_y_0_0_xyyzz, g_y_0_0_xyzz, g_y_0_0_xyzzz, g_y_0_0_xzzz, g_y_0_0_yyyy, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyz, g_y_0_0_yyyzz, g_y_0_0_yyzz, g_y_0_0_yyzzz, g_y_0_0_yzzz, g_y_0_0_yzzzz, g_y_0_0_zzzz, g_y_0_y_xxxx, g_y_0_y_xxxy, g_y_0_y_xxxz, g_y_0_y_xxyy, g_y_0_y_xxyz, g_y_0_y_xxzz, g_y_0_y_xyyy, g_y_0_y_xyyz, g_y_0_y_xyzz, g_y_0_y_xzzz, g_y_0_y_yyyy, g_y_0_y_yyyz, g_y_0_y_yyzz, g_y_0_y_yzzz, g_y_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_xxxx[k] = -g_0_xxxx[k] - g_y_0_0_xxxx[k] * ab_y + g_y_0_0_xxxxy[k];

                g_y_0_y_xxxy[k] = -g_0_xxxy[k] - g_y_0_0_xxxy[k] * ab_y + g_y_0_0_xxxyy[k];

                g_y_0_y_xxxz[k] = -g_0_xxxz[k] - g_y_0_0_xxxz[k] * ab_y + g_y_0_0_xxxyz[k];

                g_y_0_y_xxyy[k] = -g_0_xxyy[k] - g_y_0_0_xxyy[k] * ab_y + g_y_0_0_xxyyy[k];

                g_y_0_y_xxyz[k] = -g_0_xxyz[k] - g_y_0_0_xxyz[k] * ab_y + g_y_0_0_xxyyz[k];

                g_y_0_y_xxzz[k] = -g_0_xxzz[k] - g_y_0_0_xxzz[k] * ab_y + g_y_0_0_xxyzz[k];

                g_y_0_y_xyyy[k] = -g_0_xyyy[k] - g_y_0_0_xyyy[k] * ab_y + g_y_0_0_xyyyy[k];

                g_y_0_y_xyyz[k] = -g_0_xyyz[k] - g_y_0_0_xyyz[k] * ab_y + g_y_0_0_xyyyz[k];

                g_y_0_y_xyzz[k] = -g_0_xyzz[k] - g_y_0_0_xyzz[k] * ab_y + g_y_0_0_xyyzz[k];

                g_y_0_y_xzzz[k] = -g_0_xzzz[k] - g_y_0_0_xzzz[k] * ab_y + g_y_0_0_xyzzz[k];

                g_y_0_y_yyyy[k] = -g_0_yyyy[k] - g_y_0_0_yyyy[k] * ab_y + g_y_0_0_yyyyy[k];

                g_y_0_y_yyyz[k] = -g_0_yyyz[k] - g_y_0_0_yyyz[k] * ab_y + g_y_0_0_yyyyz[k];

                g_y_0_y_yyzz[k] = -g_0_yyzz[k] - g_y_0_0_yyzz[k] * ab_y + g_y_0_0_yyyzz[k];

                g_y_0_y_yzzz[k] = -g_0_yzzz[k] - g_y_0_0_yzzz[k] * ab_y + g_y_0_0_yyzzz[k];

                g_y_0_y_zzzz[k] = -g_0_zzzz[k] - g_y_0_0_zzzz[k] * ab_y + g_y_0_0_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_xxxx = cbuffer.data(pg_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_z_xxxy = cbuffer.data(pg_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_z_xxxz = cbuffer.data(pg_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_z_xxyy = cbuffer.data(pg_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_z_xxyz = cbuffer.data(pg_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_z_xxzz = cbuffer.data(pg_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_z_xyyy = cbuffer.data(pg_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_z_xyyz = cbuffer.data(pg_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_z_xyzz = cbuffer.data(pg_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_z_xzzz = cbuffer.data(pg_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_z_yyyy = cbuffer.data(pg_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_z_yyyz = cbuffer.data(pg_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_z_yyzz = cbuffer.data(pg_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_z_yzzz = cbuffer.data(pg_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_z_zzzz = cbuffer.data(pg_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxx, g_y_0_0_xxxxz, g_y_0_0_xxxy, g_y_0_0_xxxyz, g_y_0_0_xxxz, g_y_0_0_xxxzz, g_y_0_0_xxyy, g_y_0_0_xxyyz, g_y_0_0_xxyz, g_y_0_0_xxyzz, g_y_0_0_xxzz, g_y_0_0_xxzzz, g_y_0_0_xyyy, g_y_0_0_xyyyz, g_y_0_0_xyyz, g_y_0_0_xyyzz, g_y_0_0_xyzz, g_y_0_0_xyzzz, g_y_0_0_xzzz, g_y_0_0_xzzzz, g_y_0_0_yyyy, g_y_0_0_yyyyz, g_y_0_0_yyyz, g_y_0_0_yyyzz, g_y_0_0_yyzz, g_y_0_0_yyzzz, g_y_0_0_yzzz, g_y_0_0_yzzzz, g_y_0_0_zzzz, g_y_0_0_zzzzz, g_y_0_z_xxxx, g_y_0_z_xxxy, g_y_0_z_xxxz, g_y_0_z_xxyy, g_y_0_z_xxyz, g_y_0_z_xxzz, g_y_0_z_xyyy, g_y_0_z_xyyz, g_y_0_z_xyzz, g_y_0_z_xzzz, g_y_0_z_yyyy, g_y_0_z_yyyz, g_y_0_z_yyzz, g_y_0_z_yzzz, g_y_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_xxxx[k] = -g_y_0_0_xxxx[k] * ab_z + g_y_0_0_xxxxz[k];

                g_y_0_z_xxxy[k] = -g_y_0_0_xxxy[k] * ab_z + g_y_0_0_xxxyz[k];

                g_y_0_z_xxxz[k] = -g_y_0_0_xxxz[k] * ab_z + g_y_0_0_xxxzz[k];

                g_y_0_z_xxyy[k] = -g_y_0_0_xxyy[k] * ab_z + g_y_0_0_xxyyz[k];

                g_y_0_z_xxyz[k] = -g_y_0_0_xxyz[k] * ab_z + g_y_0_0_xxyzz[k];

                g_y_0_z_xxzz[k] = -g_y_0_0_xxzz[k] * ab_z + g_y_0_0_xxzzz[k];

                g_y_0_z_xyyy[k] = -g_y_0_0_xyyy[k] * ab_z + g_y_0_0_xyyyz[k];

                g_y_0_z_xyyz[k] = -g_y_0_0_xyyz[k] * ab_z + g_y_0_0_xyyzz[k];

                g_y_0_z_xyzz[k] = -g_y_0_0_xyzz[k] * ab_z + g_y_0_0_xyzzz[k];

                g_y_0_z_xzzz[k] = -g_y_0_0_xzzz[k] * ab_z + g_y_0_0_xzzzz[k];

                g_y_0_z_yyyy[k] = -g_y_0_0_yyyy[k] * ab_z + g_y_0_0_yyyyz[k];

                g_y_0_z_yyyz[k] = -g_y_0_0_yyyz[k] * ab_z + g_y_0_0_yyyzz[k];

                g_y_0_z_yyzz[k] = -g_y_0_0_yyzz[k] * ab_z + g_y_0_0_yyzzz[k];

                g_y_0_z_yzzz[k] = -g_y_0_0_yzzz[k] * ab_z + g_y_0_0_yzzzz[k];

                g_y_0_z_zzzz[k] = -g_y_0_0_zzzz[k] * ab_z + g_y_0_0_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_xxxx = cbuffer.data(pg_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_x_xxxy = cbuffer.data(pg_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_x_xxxz = cbuffer.data(pg_geom_10_off + 92 * ccomps * dcomps);

            auto g_z_0_x_xxyy = cbuffer.data(pg_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_x_xxyz = cbuffer.data(pg_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_x_xxzz = cbuffer.data(pg_geom_10_off + 95 * ccomps * dcomps);

            auto g_z_0_x_xyyy = cbuffer.data(pg_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_x_xyyz = cbuffer.data(pg_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_x_xyzz = cbuffer.data(pg_geom_10_off + 98 * ccomps * dcomps);

            auto g_z_0_x_xzzz = cbuffer.data(pg_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_x_yyyy = cbuffer.data(pg_geom_10_off + 100 * ccomps * dcomps);

            auto g_z_0_x_yyyz = cbuffer.data(pg_geom_10_off + 101 * ccomps * dcomps);

            auto g_z_0_x_yyzz = cbuffer.data(pg_geom_10_off + 102 * ccomps * dcomps);

            auto g_z_0_x_yzzz = cbuffer.data(pg_geom_10_off + 103 * ccomps * dcomps);

            auto g_z_0_x_zzzz = cbuffer.data(pg_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxx, g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxy, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxz, g_z_0_0_xxxzz, g_z_0_0_xxyy, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyz, g_z_0_0_xxyzz, g_z_0_0_xxzz, g_z_0_0_xxzzz, g_z_0_0_xyyy, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyz, g_z_0_0_xyyzz, g_z_0_0_xyzz, g_z_0_0_xyzzz, g_z_0_0_xzzz, g_z_0_0_xzzzz, g_z_0_0_yyyy, g_z_0_0_yyyz, g_z_0_0_yyzz, g_z_0_0_yzzz, g_z_0_0_zzzz, g_z_0_x_xxxx, g_z_0_x_xxxy, g_z_0_x_xxxz, g_z_0_x_xxyy, g_z_0_x_xxyz, g_z_0_x_xxzz, g_z_0_x_xyyy, g_z_0_x_xyyz, g_z_0_x_xyzz, g_z_0_x_xzzz, g_z_0_x_yyyy, g_z_0_x_yyyz, g_z_0_x_yyzz, g_z_0_x_yzzz, g_z_0_x_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_xxxx[k] = -g_z_0_0_xxxx[k] * ab_x + g_z_0_0_xxxxx[k];

                g_z_0_x_xxxy[k] = -g_z_0_0_xxxy[k] * ab_x + g_z_0_0_xxxxy[k];

                g_z_0_x_xxxz[k] = -g_z_0_0_xxxz[k] * ab_x + g_z_0_0_xxxxz[k];

                g_z_0_x_xxyy[k] = -g_z_0_0_xxyy[k] * ab_x + g_z_0_0_xxxyy[k];

                g_z_0_x_xxyz[k] = -g_z_0_0_xxyz[k] * ab_x + g_z_0_0_xxxyz[k];

                g_z_0_x_xxzz[k] = -g_z_0_0_xxzz[k] * ab_x + g_z_0_0_xxxzz[k];

                g_z_0_x_xyyy[k] = -g_z_0_0_xyyy[k] * ab_x + g_z_0_0_xxyyy[k];

                g_z_0_x_xyyz[k] = -g_z_0_0_xyyz[k] * ab_x + g_z_0_0_xxyyz[k];

                g_z_0_x_xyzz[k] = -g_z_0_0_xyzz[k] * ab_x + g_z_0_0_xxyzz[k];

                g_z_0_x_xzzz[k] = -g_z_0_0_xzzz[k] * ab_x + g_z_0_0_xxzzz[k];

                g_z_0_x_yyyy[k] = -g_z_0_0_yyyy[k] * ab_x + g_z_0_0_xyyyy[k];

                g_z_0_x_yyyz[k] = -g_z_0_0_yyyz[k] * ab_x + g_z_0_0_xyyyz[k];

                g_z_0_x_yyzz[k] = -g_z_0_0_yyzz[k] * ab_x + g_z_0_0_xyyzz[k];

                g_z_0_x_yzzz[k] = -g_z_0_0_yzzz[k] * ab_x + g_z_0_0_xyzzz[k];

                g_z_0_x_zzzz[k] = -g_z_0_0_zzzz[k] * ab_x + g_z_0_0_xzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_xxxx = cbuffer.data(pg_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_y_xxxy = cbuffer.data(pg_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_y_xxxz = cbuffer.data(pg_geom_10_off + 107 * ccomps * dcomps);

            auto g_z_0_y_xxyy = cbuffer.data(pg_geom_10_off + 108 * ccomps * dcomps);

            auto g_z_0_y_xxyz = cbuffer.data(pg_geom_10_off + 109 * ccomps * dcomps);

            auto g_z_0_y_xxzz = cbuffer.data(pg_geom_10_off + 110 * ccomps * dcomps);

            auto g_z_0_y_xyyy = cbuffer.data(pg_geom_10_off + 111 * ccomps * dcomps);

            auto g_z_0_y_xyyz = cbuffer.data(pg_geom_10_off + 112 * ccomps * dcomps);

            auto g_z_0_y_xyzz = cbuffer.data(pg_geom_10_off + 113 * ccomps * dcomps);

            auto g_z_0_y_xzzz = cbuffer.data(pg_geom_10_off + 114 * ccomps * dcomps);

            auto g_z_0_y_yyyy = cbuffer.data(pg_geom_10_off + 115 * ccomps * dcomps);

            auto g_z_0_y_yyyz = cbuffer.data(pg_geom_10_off + 116 * ccomps * dcomps);

            auto g_z_0_y_yyzz = cbuffer.data(pg_geom_10_off + 117 * ccomps * dcomps);

            auto g_z_0_y_yzzz = cbuffer.data(pg_geom_10_off + 118 * ccomps * dcomps);

            auto g_z_0_y_zzzz = cbuffer.data(pg_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxx, g_z_0_0_xxxxy, g_z_0_0_xxxy, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxz, g_z_0_0_xxyy, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyz, g_z_0_0_xxyzz, g_z_0_0_xxzz, g_z_0_0_xyyy, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyz, g_z_0_0_xyyzz, g_z_0_0_xyzz, g_z_0_0_xyzzz, g_z_0_0_xzzz, g_z_0_0_yyyy, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyz, g_z_0_0_yyyzz, g_z_0_0_yyzz, g_z_0_0_yyzzz, g_z_0_0_yzzz, g_z_0_0_yzzzz, g_z_0_0_zzzz, g_z_0_y_xxxx, g_z_0_y_xxxy, g_z_0_y_xxxz, g_z_0_y_xxyy, g_z_0_y_xxyz, g_z_0_y_xxzz, g_z_0_y_xyyy, g_z_0_y_xyyz, g_z_0_y_xyzz, g_z_0_y_xzzz, g_z_0_y_yyyy, g_z_0_y_yyyz, g_z_0_y_yyzz, g_z_0_y_yzzz, g_z_0_y_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_xxxx[k] = -g_z_0_0_xxxx[k] * ab_y + g_z_0_0_xxxxy[k];

                g_z_0_y_xxxy[k] = -g_z_0_0_xxxy[k] * ab_y + g_z_0_0_xxxyy[k];

                g_z_0_y_xxxz[k] = -g_z_0_0_xxxz[k] * ab_y + g_z_0_0_xxxyz[k];

                g_z_0_y_xxyy[k] = -g_z_0_0_xxyy[k] * ab_y + g_z_0_0_xxyyy[k];

                g_z_0_y_xxyz[k] = -g_z_0_0_xxyz[k] * ab_y + g_z_0_0_xxyyz[k];

                g_z_0_y_xxzz[k] = -g_z_0_0_xxzz[k] * ab_y + g_z_0_0_xxyzz[k];

                g_z_0_y_xyyy[k] = -g_z_0_0_xyyy[k] * ab_y + g_z_0_0_xyyyy[k];

                g_z_0_y_xyyz[k] = -g_z_0_0_xyyz[k] * ab_y + g_z_0_0_xyyyz[k];

                g_z_0_y_xyzz[k] = -g_z_0_0_xyzz[k] * ab_y + g_z_0_0_xyyzz[k];

                g_z_0_y_xzzz[k] = -g_z_0_0_xzzz[k] * ab_y + g_z_0_0_xyzzz[k];

                g_z_0_y_yyyy[k] = -g_z_0_0_yyyy[k] * ab_y + g_z_0_0_yyyyy[k];

                g_z_0_y_yyyz[k] = -g_z_0_0_yyyz[k] * ab_y + g_z_0_0_yyyyz[k];

                g_z_0_y_yyzz[k] = -g_z_0_0_yyzz[k] * ab_y + g_z_0_0_yyyzz[k];

                g_z_0_y_yzzz[k] = -g_z_0_0_yzzz[k] * ab_y + g_z_0_0_yyzzz[k];

                g_z_0_y_zzzz[k] = -g_z_0_0_zzzz[k] * ab_y + g_z_0_0_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_xxxx = cbuffer.data(pg_geom_10_off + 120 * ccomps * dcomps);

            auto g_z_0_z_xxxy = cbuffer.data(pg_geom_10_off + 121 * ccomps * dcomps);

            auto g_z_0_z_xxxz = cbuffer.data(pg_geom_10_off + 122 * ccomps * dcomps);

            auto g_z_0_z_xxyy = cbuffer.data(pg_geom_10_off + 123 * ccomps * dcomps);

            auto g_z_0_z_xxyz = cbuffer.data(pg_geom_10_off + 124 * ccomps * dcomps);

            auto g_z_0_z_xxzz = cbuffer.data(pg_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_z_xyyy = cbuffer.data(pg_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_z_xyyz = cbuffer.data(pg_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_z_xyzz = cbuffer.data(pg_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_z_xzzz = cbuffer.data(pg_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_z_yyyy = cbuffer.data(pg_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_z_yyyz = cbuffer.data(pg_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_z_yyzz = cbuffer.data(pg_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_z_yzzz = cbuffer.data(pg_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_z_zzzz = cbuffer.data(pg_geom_10_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_z_0_0_xxxx, g_z_0_0_xxxxz, g_z_0_0_xxxy, g_z_0_0_xxxyz, g_z_0_0_xxxz, g_z_0_0_xxxzz, g_z_0_0_xxyy, g_z_0_0_xxyyz, g_z_0_0_xxyz, g_z_0_0_xxyzz, g_z_0_0_xxzz, g_z_0_0_xxzzz, g_z_0_0_xyyy, g_z_0_0_xyyyz, g_z_0_0_xyyz, g_z_0_0_xyyzz, g_z_0_0_xyzz, g_z_0_0_xyzzz, g_z_0_0_xzzz, g_z_0_0_xzzzz, g_z_0_0_yyyy, g_z_0_0_yyyyz, g_z_0_0_yyyz, g_z_0_0_yyyzz, g_z_0_0_yyzz, g_z_0_0_yyzzz, g_z_0_0_yzzz, g_z_0_0_yzzzz, g_z_0_0_zzzz, g_z_0_0_zzzzz, g_z_0_z_xxxx, g_z_0_z_xxxy, g_z_0_z_xxxz, g_z_0_z_xxyy, g_z_0_z_xxyz, g_z_0_z_xxzz, g_z_0_z_xyyy, g_z_0_z_xyyz, g_z_0_z_xyzz, g_z_0_z_xzzz, g_z_0_z_yyyy, g_z_0_z_yyyz, g_z_0_z_yyzz, g_z_0_z_yzzz, g_z_0_z_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_xxxx[k] = -g_0_xxxx[k] - g_z_0_0_xxxx[k] * ab_z + g_z_0_0_xxxxz[k];

                g_z_0_z_xxxy[k] = -g_0_xxxy[k] - g_z_0_0_xxxy[k] * ab_z + g_z_0_0_xxxyz[k];

                g_z_0_z_xxxz[k] = -g_0_xxxz[k] - g_z_0_0_xxxz[k] * ab_z + g_z_0_0_xxxzz[k];

                g_z_0_z_xxyy[k] = -g_0_xxyy[k] - g_z_0_0_xxyy[k] * ab_z + g_z_0_0_xxyyz[k];

                g_z_0_z_xxyz[k] = -g_0_xxyz[k] - g_z_0_0_xxyz[k] * ab_z + g_z_0_0_xxyzz[k];

                g_z_0_z_xxzz[k] = -g_0_xxzz[k] - g_z_0_0_xxzz[k] * ab_z + g_z_0_0_xxzzz[k];

                g_z_0_z_xyyy[k] = -g_0_xyyy[k] - g_z_0_0_xyyy[k] * ab_z + g_z_0_0_xyyyz[k];

                g_z_0_z_xyyz[k] = -g_0_xyyz[k] - g_z_0_0_xyyz[k] * ab_z + g_z_0_0_xyyzz[k];

                g_z_0_z_xyzz[k] = -g_0_xyzz[k] - g_z_0_0_xyzz[k] * ab_z + g_z_0_0_xyzzz[k];

                g_z_0_z_xzzz[k] = -g_0_xzzz[k] - g_z_0_0_xzzz[k] * ab_z + g_z_0_0_xzzzz[k];

                g_z_0_z_yyyy[k] = -g_0_yyyy[k] - g_z_0_0_yyyy[k] * ab_z + g_z_0_0_yyyyz[k];

                g_z_0_z_yyyz[k] = -g_0_yyyz[k] - g_z_0_0_yyyz[k] * ab_z + g_z_0_0_yyyzz[k];

                g_z_0_z_yyzz[k] = -g_0_yyzz[k] - g_z_0_0_yyzz[k] * ab_z + g_z_0_0_yyzzz[k];

                g_z_0_z_yzzz[k] = -g_0_yzzz[k] - g_z_0_0_yzzz[k] * ab_z + g_z_0_0_yzzzz[k];

                g_z_0_z_zzzz[k] = -g_0_zzzz[k] - g_z_0_0_zzzz[k] * ab_z + g_z_0_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

